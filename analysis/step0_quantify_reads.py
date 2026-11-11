"""
Step 0 (raw reads) - Turn a project's raw sequencing reads into a table of reads per gene.

For datasets published only as raw reads (e.g. PRJNA846869). Three parts:

  1. Fetch: for every run of the project (listed by ENA, the European mirror of NCBI's archive),
     stream the first --reads read pairs from NCBI with `fastq-dump` and save them compressed.
     Reads are stored in the order the sequencer read them, not grouped by gene, so the first
     20 million pairs are a fair sample of the ~85 million per run; fewer reads only lowers the
     precision for weakly expressed genes, not the fold changes.
  2. Count: Salmon (Patro et al. 2017) matches each read pair to the transcripts it could come
     from (library type detected automatically; sequence and GC bias corrected), using the index
     from analysis/step0_build_salmon_index.py.
  3. Add up per gene, like tximport's "lengthScaledTPM" (Soneson, Love & Robinson 2015): each gene's
     transcript abundances are summed and converted back to read counts using the gene's average
     length across samples, then rounded. DESeq2 can use these as counts; using the average length
     keeps a change in which isoform is used from looking like a change in expression.

Outputs:
  data/<project>_gene_counts.tsv.gz      gene_id, gene_name, gene_type, chromosome, one column of counts per
                                         sample (named by the sample's alias, e.g. NC_A) -> the dataset's config
  data/<project>_quantification.tsv      per run: sample, reads used, % mapped, library type, software versions
  raw/<project>/fastq/, raw/<project>/salmon/   the reads and Salmon's output (git-ignored; delete when done)

Run:  uv run analysis/step0_quantify_reads.py PRJNA846869
      (needs `fastq-dump` and `salmon`: brew install sratoolkit salmon; --fetch-only downloads without counting).
      Memory: each Salmon job needs ~12 GB, so --jobs 2 (the default) needs ~24 GB. Interrupted runs can be
      restarted: finished downloads and counts are kept, unfinished ones are redone.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from common import REPO_ROOT, show_path
from step0_build_salmon_index import reference_dir

ENA_RUNS = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project}&result=read_run"
            "&fields=run_accession,sample_alias,library_layout,library_strategy,read_count&format=tsv")


def project_runs(project: str) -> pd.DataFrame:
    """The project's sequencing runs (run accession, sample alias, layout, total read pairs) from ENA."""
    with urllib.request.urlopen(ENA_RUNS.format(project=project), timeout=120) as response:
        runs = pd.DataFrame(csv.DictReader(io.StringIO(response.read().decode()), delimiter="\t"))
    if runs.empty:
        raise SystemExit(f"ENA lists no runs for {project}")
    if (runs["library_layout"] != "PAIRED").any():
        raise SystemExit("this script handles paired-end runs only")
    if runs["sample_alias"].duplicated().any():
        raise SystemExit("two runs share a sample alias; decide how to combine them before counting")
    return runs.sort_values("sample_alias").reset_index(drop=True)


def fetch(run: str, n_reads: int, folder: Path) -> int:
    """Save the first n_reads read pairs of a run as <run>_1.fastq.gz and <run>_2.fastq.gz (skipped
    when already done). Returns the number of read pairs saved."""
    done = folder / f"{run}.done"
    if done.exists():
        return int(done.read_text())
    for attempt in range(1, 4):
        result = subprocess.run(["fastq-dump", "-X", str(n_reads), "--split-files", "--gzip", "-O", str(folder), run],
                                capture_output=True, text=True)
        written = re.search(r"Written (\d+) spots", result.stdout + result.stderr)
        if result.returncode == 0 and written:
            done.write_text(written.group(1))
            print(f"  {run}: {int(written.group(1)):,} read pairs")
            return int(written.group(1))
        print(f"  {run}: fastq-dump failed (attempt {attempt}): {result.stderr.strip()[-300:]}")
    raise SystemExit(f"could not fetch {run}")


def quantify(run: str, fastq: Path, index: Path, out: Path, threads: int, expected_pairs: int) -> None:
    """Salmon on one run's read pairs. Skipped when a finished result exists: quant.sf plus run
    information saying Salmon read every expected pair (so output cut short, e.g. by a crash,
    or a damaged read file is redone or reported, never used)."""
    if (out / "quant.sf").exists() and (out / "aux_info" / "meta_info.json").exists():
        if salmon_summary(out)["reads_processed"] == expected_pairs:
            return
        print(f"  {run}: previous result incomplete, counting again")
    command = ["salmon", "quant", "-i", str(index), "-l", "A", "-1", str(fastq / f"{run}_1.fastq.gz"),
               "-2", str(fastq / f"{run}_2.fastq.gz"), "--seqBias", "--gcBias", "-p", str(threads), "-q", "-o", str(out)]
    subprocess.run(command, check=True)
    processed = salmon_summary(out)["reads_processed"]
    if processed != expected_pairs:
        raise SystemExit(f"{run}: Salmon read {processed} pairs but {expected_pairs} were fetched; the read files may be "
                         f"damaged: delete {run}* in {fastq} and {out}, then run again")
    print(f"  {run}: counted ({salmon_summary(out)['percent_mapped']}% of read pairs mapped)")


def salmon_summary(out: Path) -> dict:
    """Reads processed, % mapped and detected library type from Salmon's run information."""
    meta = json.loads((out / "aux_info" / "meta_info.json").read_text())
    library = json.loads((out / "lib_format_counts.json").read_text()) if (out / "lib_format_counts.json").exists() else {}
    return {"reads_processed": meta.get("num_processed"), "percent_mapped": round(float(meta.get("percent_mapped", np.nan)), 2),
            "library_type": library.get("expected_format", meta.get("library_types")), "salmon_version": meta.get("salmon_version")}


def gene_counts(quants: dict[str, pd.DataFrame], tx2gene: pd.DataFrame) -> pd.DataFrame:
    """Transcript estimates of every sample -> gene counts (tximport's lengthScaledTPM, rounded).

    `quants` maps sample -> Salmon's quant.sf table (index: transcript ID; columns TPM,
    EffectiveLength, NumReads). Returns genes x samples whole numbers."""
    tpm = pd.DataFrame({s: q["TPM"] for s, q in quants.items()})
    length = pd.DataFrame({s: q["EffectiveLength"] for s, q in quants.items()})
    reads = pd.DataFrame({s: q["NumReads"] for s, q in quants.items()})
    gene_of = tx2gene.set_index("transcript_id")["gene_id"].reindex(tpm.index)
    if gene_of.isna().any():
        raise SystemExit(f"{gene_of.isna().sum()} transcripts in Salmon's output are not in tx2gene.tsv")

    gene_tpm = tpm.groupby(gene_of).sum()
    gene_reads = reads.groupby(gene_of).sum()
    # each gene's length in each sample: the average of its transcripts' lengths, weighted by their abundance;
    # where a gene is not expressed in a sample, its average over the other samples
    gene_length = (tpm * length).groupby(gene_of).sum() / gene_tpm.replace(0, np.nan)
    average_length = gene_length.mean(axis=1)
    # TPM x average length is proportional to reads; scale each sample to its total number of counted reads
    scaled = gene_tpm.mul(average_length, axis=0).fillna(0.0)
    scaled = scaled * (gene_reads.sum() / scaled.sum())
    return scaled.round().astype("int64")


def main() -> None:
    parser = argparse.ArgumentParser(description="Raw reads of a sequencing project -> gene count table.")
    parser.add_argument("project", help="BioProject or study accession, e.g. PRJNA846869")
    parser.add_argument("--reads", type=int, default=20_000_000, help="read pairs per run (default 20 million)")
    parser.add_argument("--release", type=int, default=50, help="GENCODE release of the Salmon index (default 50)")
    # Each Salmon job loads the whole index: about 12 GB of memory with the GENCODE + genome-decoy index.
    # Keep jobs x 12 GB well under the computer's memory (4 jobs crashed a 48 GB Mac).
    parser.add_argument("--jobs", type=int, default=2, help="runs fetched / counted at the same time (default 2; "
                                                           "each Salmon job needs ~12 GB of memory)")
    parser.add_argument("--threads", type=int, default=6, help="Salmon threads per run (default 6)")
    parser.add_argument("--fetch-only", action="store_true", help="only download the reads")
    args = parser.parse_args()

    raw = REPO_ROOT / "raw" / args.project
    fastq, salmon_out = raw / "fastq", raw / "salmon"
    fastq.mkdir(parents=True, exist_ok=True)
    runs = project_runs(args.project)
    print(f"{args.project}: {len(runs)} paired-end runs; using the first {args.reads:,} read pairs of each")

    # 1. Fetch
    with ThreadPoolExecutor(args.jobs) as pool:
        runs["read_pairs_used"] = list(pool.map(lambda run: fetch(run, args.reads, fastq), runs["run_accession"]))
    if args.fetch_only:
        return

    # 2. Count with Salmon
    references = reference_dir(args.release)
    index = references / "salmon_index"
    if not index.exists():
        raise SystemExit(f"no Salmon index at {show_path(index)}: run analysis/step0_build_salmon_index.py first")
    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(lambda row: quantify(row[0], fastq, index, salmon_out / row[0], args.threads, row[1]),
                      zip(runs["run_accession"], runs["read_pairs_used"])))

    # 3. Add up per gene
    quants = {row.sample_alias: pd.read_csv(salmon_out / row.run_accession / "quant.sf", sep="\t", index_col="Name")
              for row in runs.itertuples()}
    counts = gene_counts(quants, pd.read_csv(references / "tx2gene.tsv", sep="\t"))
    genes = pd.read_csv(references / "genes.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("gene_id")
    counts = counts[~counts.index.str.endswith("_PAR_Y")]  # chrY copies of X/Y-shared genes: their reads go to the X copy
    table = genes.reindex(counts.index).join(counts)
    table.index = table.index.str.replace(r"\.\d+$", "", regex=True)  # ENSG00000141510.18 -> ENSG00000141510
    table.index.name = "gene_id"
    out = REPO_ROOT / "data" / f"{args.project}_gene_counts.tsv.gz"
    table.to_csv(out, sep="\t", compression={"method": "gzip", "mtime": 0})

    summary = pd.DataFrame([{"run": row.run_accession, "sample": row.sample_alias, "read_pairs_in_run": int(row.read_count),
                             "read_pairs_used": row.read_pairs_used, **salmon_summary(salmon_out / row.run_accession),
                             "reference": f"GENCODE v{args.release} transcripts + GRCh38 decoys"} for row in runs.itertuples()])
    summary_path = REPO_ROOT / "data" / f"{args.project}_quantification.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"\n{len(table):,} genes x {len(quants)} samples -> {show_path(out)}; run summary -> {show_path(summary_path)}")
    print(f"The reads and Salmon output in {show_path(raw)} can be deleted.")


if __name__ == "__main__":
    main()
