# Reference files

These files are committed, so building the website never depends on outside servers. Each is made or refreshed by
the script named below (the scripts don't edit this page: update the dates here when you refresh a file).

| File | Made by (date) | Source | Contents |
| --- | --- | --- | --- |
| `immport_gene_lists.gmt` | `analysis/step0_fetch_gene_lists.py` (2026-10-28) | [ImmPort gene lists](https://s3.immport.org/release/genelists/current/all_gene_lists.gmt) ([docs](https://docs.immport.org/apidocumentation/genelist/overview/)) | 153 immune gene lists (3,118 genes) curated from GO "immune system process" (GO:0002376) and Reactome "Immune System" (R-HSA-168256). One list per line: name, link, gene symbols. |
| `hgnc_symbols.tsv` | `analysis/step0_fetch_gene_lists.py` (2026-10-28) | [HGNC complete set](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) ([download page](https://www.genenames.org/download/archive/)) | Current official symbol, Ensembl gene ID and previous symbols for 45,289 human genes. Used to match genes whose symbol changed. |
| `human_immune_cell_atlas.txt` | kept by hand | Human Immune Cell Atlas (project document, plain-text export) | Candidate marker genes for ~200 human immune cell populations and functional states, in 14 sections. |
| `immune_cell_markers.tsv` | `analysis/step0_extract_atlas_markers.py` | the atlas above | The atlas tables, one row per population: section, name, kind, markers, low markers, original text. Rebuild after editing the atlas. |
| `GPL10332_probes.tsv` | `analysis/step0_fetch_platform_annotation.py GPL10332` (2026-10-31) | [GEO platform GPL10332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10332) | Agilent 4x44K v2 probes: probe ID, gene symbol, gene name, Ensembl transcript ID, chromosome (34,127 probes). Used for GSE46157. |
| `mouse_genes.tsv.gz` | `analysis/step0_fetch_mouse_genes.py` (2026-11-10) | MGI reports [MRK_ENSEMBL.rpt](https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt) and [HOM_MouseHumanSequence.rpt](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt) (homology from the Alliance of Genome Resources), both dated 2026-11-09 | 77,645 mouse Ensembl genes: MGI symbol, name, chromosome, type, and the human ortholog when it is one-to-one (17,099 genes; counted over the whole report in both directions). Ambiguous genes keep their human partners in `human_genes` but get no ortholog. Used for GSE314547. |
| `icp_dataset_handoff.md` | written by hand (2026-11-11) | literature and archive searches | Catalog of 18 candidate ICP datasets (human, mouse, rat, pig) with their status and blockers. |
