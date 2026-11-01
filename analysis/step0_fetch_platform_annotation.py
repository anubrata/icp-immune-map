"""
Step 0 (microarrays) - Download a microarray platform's probe annotation from GEO.

GEO (NCBI's Gene Expression Omnibus, the public archive the datasets come from) describes each
array design as a "platform" with an ID like GPL10332. Its table lists, for every probe on the
array, the gene the probe measures. Microarray data files only name their probes, so step 1 needs
this table to turn probes into genes. Run once per platform; the output is committed.

Output: resources/<platform>_probes.tsv, one row per probe:
          probe_id         the probe name used in the data files (e.g. A_23_P146146)
          gene_symbol      gene measured by the probe ("" if none)
          gene_name        full gene name
          ensembl_id       Ensembl ID if given (for Agilent arrays a transcript ID, ENST...)
          chromosome       e.g. "7" or "Y" (from the probe's chromosomal location)

Run:  uv run analysis/step0_fetch_platform_annotation.py GPL10332
      (GEO's server is slow and sometimes drops the connection; the script retries. If it keeps
      failing, download the URL below with a browser or curl and pass it: --from-file <file>)
"""

from __future__ import annotations

import argparse
import http.client
import io
import urllib.request
from pathlib import Path

import pandas as pd

from common import RESOURCES_DIR, show_path

URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={platform}&targ=self&form=text&view=data"

# Column names in GEO's table for Agilent platforms -> our names
COLUMNS = {"NAME": "probe_id", "GENE_SYMBOL": "gene_symbol", "GENE_NAME": "gene_name",
           "ENSEMBL_ID": "ensembl_id", "CHROMOSOMAL_LOCATION": "chromosome", "CONTROL_TYPE": "control_type"}


def download(platform: str, attempts: int = 3) -> str:
    url = URL.format(platform=platform)
    for attempt in range(1, attempts + 1):
        print(f"Downloading {url} (attempt {attempt}; GEO can be slow, this may take several minutes)")
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "icp-immune-map (python urllib)"})
            with urllib.request.urlopen(request, timeout=900) as response:
                return response.read().decode("utf-8", errors="replace")
        except (OSError, http.client.HTTPException) as error:
            print(f"  failed: {error!r}")
    raise SystemExit(f"Could not download {url}. Download it another way and use --from-file.")


def parse_platform_table(text: str) -> pd.DataFrame:
    """The rows between !platform_table_begin and !platform_table_end, as a table."""
    lines = text.splitlines()
    if "!platform_table_begin" not in lines or "!platform_table_end" not in lines:
        raise SystemExit("The download does not contain a complete platform table (interrupted?). Try again.")
    start, end = lines.index("!platform_table_begin") + 1, lines.index("!platform_table_end")
    return pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", dtype=str, keep_default_na=False)


def probe_table(table: pd.DataFrame) -> pd.DataFrame:
    """Gene probes only (no controls), one row per probe name, our column names."""
    missing = [c for c in COLUMNS if c not in table.columns]
    if missing:
        raise SystemExit(f"Unexpected platform table: missing columns {missing}. Columns: {list(table.columns)}")
    probes = table[list(COLUMNS)].rename(columns=COLUMNS)
    probes = probes[probes["control_type"].str.upper() == "FALSE"].drop(columns="control_type")
    # "chr7:100000-100060" -> "7"; "unmapped" or empty -> ""
    probes["chromosome"] = probes["chromosome"].str.extract(r"^chr([0-9XYM]+|MT)[:_]", expand=False).fillna("")
    return probes.drop_duplicates("probe_id").sort_values("probe_id")


def main() -> None:
    parser = argparse.ArgumentParser(description="Download a GEO platform's probe annotation.")
    parser.add_argument("platform", help="GEO platform ID, e.g. GPL10332")
    parser.add_argument("--from-file", type=Path, help="use an already downloaded copy of the platform table")
    args = parser.parse_args()
    platform = args.platform
    text = args.from_file.read_text(errors="replace") if args.from_file else download(platform)
    probes = probe_table(parse_platform_table(text))
    out = RESOURCES_DIR / f"{platform}_probes.tsv"
    probes.to_csv(out, sep="\t", index=False)
    with_gene = (probes["gene_symbol"] != "").sum()
    print(f"  {len(probes):,} probes, {with_gene:,} with a gene symbol "
          f"({probes.loc[probes['gene_symbol'] != '', 'gene_symbol'].nunique():,} genes) -> {show_path(out)}")


if __name__ == "__main__":
    main()
