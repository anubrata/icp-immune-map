# Reference files

The ImmPort and HGNC files were downloaded on 2026-10-28 by `analysis/step0_fetch_gene_lists.py`;
re-run it to refresh them. The atlas document is kept by hand: after changing it, run
`analysis/step0_extract_atlas_markers.py` to rebuild the marker table.

| File | Source | Contents |
| --- | --- | --- |
| `immport_gene_lists.gmt` | [ImmPort gene lists](https://s3.immport.org/release/genelists/current/all_gene_lists.gmt) ([docs](https://docs.immport.org/apidocumentation/genelist/overview/)) | 153 immune gene lists (3,118 genes) curated from GO "immune system process" (GO:0002376) and Reactome "Immune System" (R-HSA-168256). One list per line: name, link, gene symbols. |
| `hgnc_symbols.tsv` | [HGNC complete set](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) ([download page](https://www.genenames.org/download/archive/)) | Current official symbol, Ensembl gene ID and previous symbols for 45,289 human genes. Used to match genes whose symbol changed. |
| `human_immune_cell_atlas.txt` | Human Immune Cell Atlas (project document, plain-text export) | Candidate marker genes for ~200 human immune cell populations and functional states, in 14 sections. |
| `immune_cell_markers.tsv` | Made from the atlas by `analysis/step0_extract_atlas_markers.py` | The atlas tables, one row per population: section, name, kind, markers, low markers, original text. |
