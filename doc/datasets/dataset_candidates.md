# ICP dataset candidates

Status of every dataset considered for the ICP Immune Map, grouped by what it would take to add it. The full catalog,
with design details and sources for each entry, is in [`resources/icp_dataset_handoff.md`](../../resources/icp_dataset_handoff.md)
(IDs H01–H10 human, M01–M04 mouse, O01–O04 rat and pig).

Last updated: 2026-11-11 (PRJNA846869 added). Entries marked *checked* were verified against the archives on
2026-11-10 or 2026-11-11; the others are as reported in the handoff catalog and have not been re-checked.

## On the dashboard

| ID | Dataset | Species · tissue | Design | Notes |
| --- | --- | --- | --- | --- |
| H01 | [GSE282182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282182) | Human placenta, RNA-seq | 5 ICP vs 5 healthy | Not adjusted for fetal sex (healthy 4 F / 1 M, ICP 2 F / 3 M) |
| H10 | [GSE46157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46157) | Human placenta, Agilent microarray | 4 ICP vs 2 healthy | Each array is a pool of 5 placentas; healthy and ICP arrays were run on separate slides |
| H02 | [PRJNA846869](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA846869) ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9263217/)) | Human placenta, total RNA-seq | 7 ICP vs 5 healthy | Raw reads only, so counted here: first 20M of ~85M read pairs per placenta, Salmon + GENCODE v50 (`analysis/step0_quantify_reads.py`). Groups from the archive's sample names (NC = healthy). Little decidua except ICP_G |
| M03 | [GSE314547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314547) | Mouse liver, RNA-seq | 4 ethinylestradiol vs 4 vehicle | A model, not patients: humanized PXR-CAR-CYP3A4/3A7 mice, pregnancy not stated; the MI-763F groups are left out. *Checked:* integer counts, library names match the GEO sample records |

## Could be added

| ID | Dataset | Species · tissue | Design | What's needed |
| --- | --- | --- | --- | --- |
| H04 | [HRA003387](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003387) ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10570770/)) | Human maternal plasma, cell-free RNA | 74 ICP or pre-ICP vs 171 healthy | GSA-Human (controlled access): likely needs a data-access application. Blood, not placenta, so analysed separately; only 19 of the 74 had raised bile acids at sampling. |
| M01 | [PRJNA1203061](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1203061) ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12214484/)) | Mouse liver, DEHP-induced gestational cholestasis | 3 vs 3 | *Checked:* not public (no runs at NCBI SRA or ENA). The authors would need to release it. |
| M02 | Xue et al. 2026 ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13046547/)) | Mouse placenta, DDC diet | 6 vs 4 placentas | No accession; ask the authors. Placentas may share mothers, so the litter of each is needed. |
| O01 | [PRJCA058946](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA058946) ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13223021/)) | Rat liver and ileum | 18 RNA libraries | Public, but which animals are control, ICP or treated is not stated; needs the paper's supplement or the authors. |

## Single-cell and small-RNA studies (would need new pipelines)

| ID | Dataset | Design | Issue |
| --- | --- | --- | --- |
| H03 | [PRJNA949258](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA949258) ([PubMed](https://pubmed.ncbi.nlm.nih.gov/38835155/)) | Placenta single-cell | 1 ICP vs 1 healthy: exploratory only |
| H05 | Chen et al. 2026, STAT5A macrophages ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13005436/)) | Placenta single-cell, 3 vs 3 | No accession found; the most valuable one if the authors share donor-level data |
| H06 | Yang et al. 2025 ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12351929/)) | Serum small RNA (tRNA fragments), 3 vs 3 | No accession; a different kind of data |
| H07 | Huang et al. 2025 ([PubMed](https://pubmed.ncbi.nlm.nih.gov/40278950/)) | Placenta single-cell | Donor numbers and data source unclear; may reuse H03 or H05 |
| H08 | Tang et al. 2025 ([PubMed](https://pubmed.ncbi.nlm.nih.gov/40796229/)) | Placenta single-cell and bulk | May reuse cohorts already listed |
| H09 | Tang et al. 2026 ([PubMed](https://pubmed.ncbi.nlm.nih.gov/41397582/)) | Placenta single-cell and public bulk | May reuse cohorts already listed |

## Excluded (not an ICP vs healthy comparison of RNA)

| ID | Dataset | Reason |
| --- | --- | --- |
| M04 | [GSE70179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70179) (mouse) | Abcb11 knockout model, liver arm not pregnant, one array per condition |
| O02 | [GSE110414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110414) (pig placenta) | Compares stages of normal pregnancy, no ICP (*checked:* data downloadable) |
| O03 | [GSE214238](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214238) (pig liver) | Chromatin (ChIP-seq), not RNA |
| O04 | [PRJNA1002065](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1002065) (rat) | Gut microbiome sequencing, not host RNA |

## Suggested next steps

1. **Contact the authors of H05 and M02.** They would add the most where the map is thinnest: immune cells in the
   placenta, resolved by cell type (H05) or in an animal model (M02). H05's single-cell data would also show which
   placental compartment (villi or decidua) the NK-cell differences come from, which bulk samples with varying amounts
   of decidua cannot settle (H02 added on 2026-11-11 did not reproduce the lower NK-cell signal of H01).
