"""
Step 0 (raw reads) - Download the human reference and build a Salmon index.

Some datasets (e.g. PRJNA846869) are published only as raw sequencing reads, with no table of
reads per gene. To make that table we count the reads ourselves with Salmon (Patro et al. 2017),
which matches each read to the transcripts it could come from. Salmon needs an index of the
reference transcripts; this script builds it once (about an hour; nothing here is committed).

Reference: GENCODE human release (default 50), GRCh38:
  gencode.v<N>.transcripts.fa.gz     every annotated transcript (what reads are counted against)
  GRCh38.primary_assembly.genome.fa.gz  the genome, used as a "decoy": reads from outside annotated
                                     transcripts (e.g. introns, plentiful in total-RNA libraries) map to it
                                     instead of being forced onto a similar transcript (Srivastava et al. 2020)
  gencode.v<N>.annotation.gtf.gz     gene annotation: symbol, type and chromosome of each gene

Outputs (in references/gencode_v<N>/, git-ignored):
  salmon_index/       the index
  tx2gene.tsv         transcript ID -> gene ID (to add transcripts up to genes)
  genes.tsv           gene_id, gene_name, gene_type, chromosome (one row per gene)

Run:  uv run analysis/step0_build_salmon_index.py            (needs `salmon`: brew install salmon)
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
from pathlib import Path

import pandas as pd

from common import REPO_ROOT, show_path

GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}/{name}"
FILES = {"transcripts": "gencode.v{release}.transcripts.fa.gz",
         "genome": "GRCh38.primary_assembly.genome.fa.gz",
         "annotation": "gencode.v{release}.annotation.gtf.gz"}


def reference_dir(release: int) -> Path:
    return REPO_ROOT / "references" / f"gencode_v{release}"


def download(url: str, path: Path) -> None:
    """Download with curl (resumes an interrupted download; retries on network errors)."""
    if path.exists() and not path.with_suffix(path.suffix + ".part").exists():
        print(f"  already downloaded: {show_path(path)}")
        return
    part = path.with_suffix(path.suffix + ".part")
    print(f"Downloading {url}")
    subprocess.run(["curl", "-fL", "--retry", "5", "--retry-delay", "10", "-C", "-", "-o", str(part), url], check=True)
    part.rename(path)


def sequence_names(fasta_gz: Path) -> list[str]:
    """Names of the sequences in a gzipped FASTA file (the first word after '>')."""
    with gzip.open(fasta_gz, "rt") as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def transcript_table(transcripts_fa: Path) -> pd.DataFrame:
    """transcript_id -> gene_id from GENCODE FASTA headers, which look like
    >ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|"""
    rows = []
    with gzip.open(transcripts_fa, "rt") as f:
        for line in f:
            if line.startswith(">"):
                fields = line[1:].strip().split("|")
                rows.append((fields[0], fields[1]))
    return pd.DataFrame(rows, columns=["transcript_id", "gene_id"])


def gene_table(annotation_gtf: Path) -> pd.DataFrame:
    """One row per gene from the GTF's 'gene' lines: gene_id, gene_name, gene_type, chromosome."""
    rows = []
    with gzip.open(annotation_gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 8)
            if fields[2] != "gene":
                continue
            attributes = dict(re.findall(r'(\S+) "([^"]*)"', fields[8]))
            chromosome = "MT" if fields[0] == "chrM" else fields[0].removeprefix("chr")  # chr17 -> 17
            rows.append((attributes["gene_id"], attributes.get("gene_name", ""), attributes.get("gene_type", ""), chromosome))
    return pd.DataFrame(rows, columns=["gene_id", "gene_name", "gene_type", "chromosome"])


def main() -> None:
    parser = argparse.ArgumentParser(description="Download GENCODE and build a decoy-aware Salmon index.")
    parser.add_argument("--release", type=int, default=50, help="GENCODE human release (default 50)")
    parser.add_argument("--threads", type=int, default=0, help="threads for salmon index (0 = all cores)")
    args = parser.parse_args()

    folder = reference_dir(args.release)
    folder.mkdir(parents=True, exist_ok=True)
    paths = {key: folder / name.format(release=args.release) for key, name in FILES.items()}
    for key, path in paths.items():
        download(GENCODE_URL.format(release=args.release, name=path.name), path)

    # Gene annotation and transcript -> gene map
    transcripts = transcript_table(paths["transcripts"])
    genes = gene_table(paths["annotation"])
    transcripts.to_csv(folder / "tx2gene.tsv", sep="\t", index=False)
    genes.to_csv(folder / "genes.tsv", sep="\t", index=False)
    print(f"  {len(transcripts):,} transcripts of {transcripts['gene_id'].nunique():,} genes; "
          f"annotation for {len(genes):,} genes -> {show_path(folder)}/tx2gene.tsv, genes.tsv")

    # Decoy-aware index: transcripts first, then the genome, whose sequence names are listed as decoys
    index = folder / "salmon_index"
    if (index / "versionInfo.json").exists() or any(index.glob("*")):
        print(f"  index already exists: {show_path(index)} (delete it to rebuild)")
        return
    decoys = folder / "decoys.txt"
    decoys.write_text("\n".join(sequence_names(paths["genome"])) + "\n")
    command = ["salmon", "index", "-t", str(paths["transcripts"]), str(paths["genome"]), "-d", str(decoys),
               "--gencode", "-p", str(args.threads), "-i", str(index)]
    print("Building the Salmon index (this takes a while):\n  " + " ".join(command))
    subprocess.run(command, check=True)
    print(f"  index -> {show_path(index)}")


if __name__ == "__main__":
    main()
