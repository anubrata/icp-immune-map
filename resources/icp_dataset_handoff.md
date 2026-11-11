# ICP datasets discussed in this conversation: catalog and Codex handoff

Prepared: 2026-11-11 (UTC)  
Target repository: <https://github.com/anubrata/icp-immune-map>  
Disease: intrahepatic cholestasis of pregnancy (ICP)

## Task for Codex

Extend the existing ICP immune analysis using the eligible datasets below. First add a source-backed dataset catalog and verify accessions, assay types, sample identities, biological replication, and available count matrices. Then integrate accessible, suitable datasets using the repository's conventions. Keep unresolved candidates in the catalog with their blockers. Do not interpret inclusion in this document as confirmation that a dataset is ready for differential expression.

Read the repository's current README, configs, analysis code, and any applicable AGENTS.md before changing it. This handoff reflects the `main` branch inspected on the preparation date; it is not an implementation or a request to duplicate existing datasets.

### Scope and evidence

- This document consolidates the human, mouse, rat, and pig records discussed in the conversation, including exclusions and unresolved study leads.
- There are **18 catalog entries**, not 18 independent, publicly downloadable RNA-seq case-control cohorts. Several papers may reuse the same cohort; several records use another assay or lack a verified accession.
- An accession reported by a paper does not establish that raw reads or count matrices can currently be downloaded. Recheck access and record the verification date.
- Sample counts below distinguish donors/animals, sequenced specimens, arrays, and libraries where possible. Unknown biological replication must remain unknown.
- The animal records are experimental models or related cholestasis studies. They should not be described as naturally diagnosed human ICP or automatically as wild-type-versus-knockout experiments.
- No additional nonhuman ICP single-cell/single-nucleus RNA-seq dataset was verified in the conversation. This is a search result, not proof that none exists.

## Repository starting point

| Item | Current implementation | Consequence for adding data |
| --- | --- | --- |
| GSE282182 | Already configured in `config/GSE282182.yaml`; human placenta bulk RNA-seq; 5 healthy and 5 ICP samples | Audit and retain; do not create a second copy of this cohort |
| GSE46157 | Already configured in `config/GSE46157.yaml`; human placenta Agilent microarrays; 2 healthy and 4 ICP arrays | Retain as a separate assay; do not call it RNA-seq |
| Bulk input | `analysis/step1_load_data.py` accepts a gene-by-sample table of nonnegative integer read counts | TPM, FPKM, fold changes, and DEG lists cannot substitute for raw counts in DESeq2 |
| Microarray input | Raw Agilent Feature Extraction files plus platform annotation | Other assay platforms need their own validated preprocessing |
| Differential expression | `analysis/step2_differential_expression.py`; DESeq2 for counts, moderated t-tests for arrays | Preserve assay-specific models; verify the biological unit and covariates |
| Immune annotation | `analysis/step3_immune_annotation.py` uses the Human Immune Cell Atlas; ImmPort is optional | Nonhuman data require species-aware references or a documented ortholog mapping layer |
| Export and checks | `analysis/step4_export_site_data.py`, `analysis/step5_check_site_data.py`, `run_all.sh` | New datasets must pass the same export consistency checks |

The existing placenta configs check decidual composition using `PRL` and `IGFBP1`. GSE282182 infers fetal sex from chrY expression. GSE46157 supplies fetal-sex labels but disables sex adjustment because sex coincides with severity; chrY genes are excluded from its test. Do not transfer these assumptions to maternal liver, serum, plasma, ileum, or other species.

The current configs use nominal `pvalue < 0.05` and `abs_log2fc >= 1`, while also reporting adjusted p-values. Preserve existing behavior unless explicitly changing and documenting the statistical policy; distinguish nominal from FDR-significant findings in new outputs.

Repository references: [README](https://github.com/anubrata/icp-immune-map/blob/main/README.md), [GSE282182 config](https://github.com/anubrata/icp-immune-map/blob/main/config/GSE282182.yaml), [GSE46157 config](https://github.com/anubrata/icp-immune-map/blob/main/config/GSE46157.yaml), [analysis directory](https://github.com/anubrata/icp-immune-map/tree/main/analysis).

## Catalog overview

`Candidate` means suitable in principle after the stated checks. `Pending` means access, metadata, cohort independence, or assay identity remains unresolved. `Separate` means the existing human bulk pipeline cannot be used unchanged. `Exclude` means excluded from primary ICP RNA-expression case-control analysis, with the reason retained in the catalog.

| ID | Species | Accession or study lead | Assay / tissue | Sequenced comparison or design | Recommended disposition |
| --- | --- | --- | --- | --- | --- |
| H01 | Human | GSE282182; PRJNA1187538 | Bulk RNA-seq, placenta | 5 ICP vs 5 healthy | Existing |
| H02 | Human | PRJNA846869; reported reanalysis identifier SRP378990 | Bulk total RNA-seq and separate small-RNA libraries, placenta | 7 ICP vs 5 healthy donors | Candidate; verify accession linkage and access |
| H03 | Human | PRJNA949258; SRP429602 | scRNA-seq, placenta / maternal-fetal interface | 1 ICP vs 1 healthy donor | Separate; exploratory only for disease comparison |
| H04 | Human | HRA003387 | Bulk cell-free RNA-seq, maternal plasma | 74 ICP/eventual-ICP vs 171 healthy; separate 40-HBV arm | Candidate, separate compartment and sampling-time contrasts |
| H05 | Human | Chen et al., STAT5A macrophage study, 2026 | scRNA-seq, placenta | 3 ICP vs 3 healthy sequenced donors | Pending accession / counts |
| H06 | Human | Yang et al., tRNA-derived fragments, 2025 | Small RNA-seq, serum | 3 ICP vs 3 healthy donors | Pending accession; separate small-RNA workflow |
| H07 | Human | Huang et al., mitochondrial dysfunction, 2025 | Placental scRNA-seq reported | Sequenced donor counts unverified | Pending metadata and provenance |
| H08 | Human | Tang et al., autophagy multi-omics, 2025 | Placental scRNA-seq / bulk RNA-seq analyses | New versus reused cohorts unresolved | Pending cohort deduplication |
| H09 | Human | Tang et al., ACADVL / ferroptosis, 2026 | Placental scRNA-seq and public bulk RNA-seq analysis | New versus reused cohorts unresolved | Pending cohort deduplication |
| H10 | Human | GSE46157; PRJNA197384 | Agilent microarray, placenta | 2 healthy, 2 mild-ICP, 2 severe-ICP arrays | Existing; separate assay |
| M01 | Mouse | PRJNA1203061 | Bulk mRNA-seq, maternal liver | 3 DEHP-model vs 3 control sequenced samples | Candidate animal model; public access unverified |
| M02 | Mouse | Xue et al., gestational cholestasis / fetal growth, 2026 | Bulk RNA-seq, placenta | 6 DDC vs 4 control placentas | Pending accession and dam/litter mapping |
| M03 | Mouse | GSE314547; PRJNA1391489 | Bulk RNA-seq, liver; humanized mice | 4 groups, 4 samples/group | Separate related estrogen-cholestasis model; not an established WT pregnant comparison |
| M04 | Mouse | GSE70179; PRJNA287800 | Affymetrix microarrays, liver and offspring lung | 12 arrays across distinct comparisons | Exclude from RNA-seq; genetic/mechanistic context |
| O01 | Rat | PRJCA058946; CRA039245 / CRA039246; PRJDB44845; DRP022399 | Metadata labels 18 host RNA-seq libraries: 8 liver, 10 ileum; project also contains amplicons | RNA-library group labels unresolved | Pending assay identity and sample-to-group mapping |
| O02 | Pig | GSE110414; PRJNA433661; SRP132585 | Bulk total RNA-seq, placenta | Gestation days 60 and 90, and farrowing; 10 samples | Exclude from primary ICP case-control DE; developmental reference |
| O03 | Pig | GSE214238; PRJNA884413; SRP399730 | H4K12la ChIP-seq, liver / hepatocytes | High vs low/normal bile acids; IP and input libraries | Exclude from RNA-expression pipeline |
| O04 | Rat | PRJNA1002065; associated OMIX004988 deposition | Fecal shotgun metagenomics, with host tissue qPCR | Healthy / estrogen-induced ICP / treatment arms | Exclude from host RNA-expression pipeline |

## Human records

### H01 — GSE282182 / PRJNA1187538

- **Design:** human placenta, 5 ICP and 5 healthy pregnancies; bulk RNA-seq on Illumina NovaSeq 6000.
- **Availability:** public GEO series with SRA linkage and processed gene-count and FPKM files. Relevant count file: `GSE282182_gene_count.txt.gz`; the FPKM file is not a DESeq2 input.
- **Integration:** already present. Preserve the actual matrix column names in the existing YAML; GEO sample titles and matrix column labels differ. Audit fetal-sex and decidual-composition effects.
- **Deduplication:** a later publication/reanalysis of this series does not constitute another cohort.
- **Sources:** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282182), [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1187538).

### H02 — PRJNA846869; SRP378990 reported for the reused bulk cohort

- **Study:** “Immune Dysfunction Mediated by the ceRNA Regulatory Network in Human Placenta Tissue of Intrahepatic Cholestasis Pregnancy” (2022).
- **Design:** placentas from 7 ICP and 5 healthy donors. The study generated ribodepleted total-RNA libraries for coding/noncoding transcript analyses and separate small-RNA libraries. Twelve donors must not be inflated to the number of libraries or runs.
- **Availability:** the original paper deposits data under PRJNA846869. A subsequent reanalysis identifies the bulk cohort as SRP378990. Confirm the project-to-study/run relationship in repository metadata before recording these as verified aliases. Public raw-read download and a reusable raw-count matrix were not fully confirmed in the conversation.
- **Integration:** first priority among new human placental bulk candidates. Use coding-gene raw counts for the existing bulk analysis; handle miRNA/tRF and other noncoding analyses separately.
- **Sources:** [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA846869), [SRA study search](https://www.ncbi.nlm.nih.gov/sra/?term=SRP378990), [original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9263217/), [DOI](https://doi.org/10.3389/fimmu.2022.883971).

### H03 — PRJNA949258 / SRP429602

- **Study:** “Single-cell RNA sequencing reveals heterogeneity and differential expression of the maternal–fetal interface during ICP and normal pregnancy” (Liu et al., 2024).
- **Design:** placental / maternal-fetal interface scRNA-seq from **one ICP donor and one healthy donor**. Multiple experiments or libraries are not additional biological replicates. The SRA platform includes DNBSEQ-T7; confirm the single-cell chemistry rather than assuming a 10x protocol.
- **Availability:** public SRA records link both accessions. One normal experiment/run discussed was `SRX19788456` / `SRR23984534`; obtain the complete run manifest before downloading.
- **Integration:** useful for cell annotation and descriptive marker exploration. Do not report donor-replicated case-control DE or use cells as independent patients. Any pseudobulk disease contrast remains n=1 versus n=1.
- **Sources:** [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA949258), [SRA study](https://www.ncbi.nlm.nih.gov/sra/?term=SRP429602), [example experiment](https://www.ncbi.nlm.nih.gov/sra/SRX19788456), [PubMed](https://pubmed.ncbi.nlm.nih.gov/38835155/), [DOI](https://doi.org/10.1080/14767058.2024.2361278).

### H04 — HRA003387

- **Study:** “Early characterisation and prediction of liver diseases in pregnancy by plasma cell-free RNAs” (Sun et al., 2023).
- **Design:** maternal plasma cfRNA from 74 ICP/eventual-ICP participants and 171 healthy pregnant controls, with a separate 40-person HBV group. At sampling, **19 ICP participants had elevated TBA and 55 were pre-ICP with normal TBA**. Do not label all 74 as active ICP at collection.
- **Availability:** deposited in GSA-Human; current download permissions, sample metadata, and count availability need checking.
- **Integration:** analyze plasma separately from placenta. Distinguish active-disease and prospective prediction contrasts, account for gestational timing, and keep HBV out of the healthy control group. Validate one-to-one participant/sample mapping and any repeated collections.
- **Sources:** [GSA-Human](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003387), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10570770/), [DOI](https://doi.org/10.1002/ctm2.1439).

### H05 — STAT5A / macrophage scRNA-seq study; accession needed

- **Study:** Chen, Huang, Liao and Shao (2026), “Metamorphosing macrophages promote placental fibrosis by upregulating the transcription factor STAT5A.”
- **Design:** placenta scRNA-seq from **3 ICP and 3 healthy donors**. The larger 10-versus-10 experimental tissue cohort is not the sequenced donor count. The methods describe magnetic cell enrichment; verify the selection and captured populations before treating the data as an unbiased whole-placenta cell census.
- **Availability:** no concrete raw-read or count-matrix accession was verified. The article/supplement availability statement alone does not establish a downloadable expression matrix.
- **Integration:** potentially valuable replicated single-cell candidate if donor-resolved counts and metadata can be obtained. Keep pending until then.
- **Sources:** [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13005436/), [DOI](https://doi.org/10.1186/s12967-026-07727-9).

### H06 — Serum tRNA-derived fragment study; accession needed

- **Study:** Yang et al. (2025), “Dysregulated tRNA-derived fragments impair fatty acid metabolism in intrahepatic cholestasis of pregnancy.”
- **Design:** small RNA-seq of serum from 3 ICP and 3 healthy pregnant donors, one sample per donor; focus on tRNA-derived fragments.
- **Availability:** no concrete accession was verified. The paper refers to repository information in the article/supplement; inspect it before claiming public raw data.
- **Integration:** a separate small-RNA workflow and feature annotation are required. These data cannot directly supply the existing protein-coding immune-gene expression analysis.
- **Sources:** [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12351929/), [DOI](https://doi.org/10.3389/fmed.2025.1630677).

### H07 — Mitochondrial dysfunction / trophoblast study; provenance unresolved

- **Study:** Huang et al. (2025), “Mechanism of mitochondrial dysfunction on placental trophoblastic cells in intrahepatic cholestasis of pregnancy.”
- **Assay:** placental single-cell transcriptomic analysis is reported, but sequenced donor numbers and an independent accession were not established.
- **Uncertainty:** the accessible data-availability statement says no datasets were generated or analyzed, which is inconsistent with the reported transcriptomic analysis. This requires clarification; it is not evidence of an additional downloadable cohort.
- **Integration:** retain as a literature lead. Determine whether the analysis used H03, H05, another deposited cohort, or newly generated data before assigning a new dataset identity.
- **Sources:** [PubMed](https://pubmed.ncbi.nlm.nih.gov/40278950/), [DOI](https://doi.org/10.1007/s10735-025-10427-1).

### H08 — Autophagy multi-omics study; cohort independence unresolved

- **Study:** Tang et al. (2025), “Multi-omics characterization of autophagy-related molecular features in intrahepatic cholestasis of pregnancy.”
- **Assays:** placental scRNA-seq / bulk RNA-seq analyses and lipidomics, with tissue, cell, and animal validation.
- **Availability:** no additional independent sequencing accession was verified in the conversation. Sequenced donor counts and new-versus-reused cohorts need full-text/supplement review.
- **Integration:** check overlap with H02/H03 and other human records before adding a cohort. Animal validation does not by itself establish a mouse RNA-seq dataset.
- **Sources:** [PubMed](https://pubmed.ncbi.nlm.nih.gov/40796229/), [DOI](https://doi.org/10.1093/biolre/ioaf184).

### H09 — ACADVL / placental ferroptosis study; cohort independence unresolved

- **Study:** Tang et al., “ACADVL upregulation is associated with placental ferroptosis in intrahepatic cholestasis of pregnancy”; online December 2025, journal issue March 2026.
- **Assays:** placental scRNA-seq, lipidomics, and public bulk RNA-seq analysis.
- **Availability:** no new independent accession was verified. The abstract's use of public RNA-seq is not evidence of a newly recruited bulk cohort.
- **Integration:** recover the original accessions and donor provenance; map reused data to their existing catalog entries. Count a new cohort only if independently established.
- **Sources:** [PubMed](https://pubmed.ncbi.nlm.nih.gov/41397582/), [DOI](https://doi.org/10.1016/j.bbalip.2025.159711).

### H10 — GSE46157 / PRJNA197384; microarray, already included

- **Study:** Du et al. (2014), placental gene-expression profiling in ICP.
- **Assay:** Agilent one-color expression microarray, not RNA-seq.
- **Design:** six arrays: 2 healthy, 2 mild ICP, 2 severe ICP. The reported underlying cohort comprises 30 donors, pooled in groups of five: 10 healthy, 10 mild ICP, 10 severe ICP. Reconfirm pooling against source/sample metadata. The analysis has six pooled expression units, not 30 separately measured transcriptomes.
- **Integration:** already configured. Retain as a separate platform for corroboration, with pooling, severity, and fetal-sex confounding documented.
- **Sources:** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46157), [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA197384), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4105836/), [DOI](https://doi.org/10.1186/1755-8794-7-42).

## Mouse records

### M01 — PRJNA1203061; DEHP-induced gestational cholestasis

- **Study:** “Transcriptomic and network toxicology approaches reveal potential mechanisms of DEHP induced intrahepatic cholestasis during pregnancy” (2025).
- **Design:** pregnant ICR mice; DEHP 200 mg/kg/day by gavage from gestational day 0 to 16 versus controls. Maternal liver mRNA-seq used **3 samples per group**. The overall animal experiment's 6 mice per group is not the sequencing n.
- **Model interpretation:** chemically exposed standard mice, not a knockout-versus-WT experiment. Toxicant-specific effects may differ from human ICP.
- **Availability:** the paper reports PRJNA1203061, but unrestricted public download was not confirmed. Use the normal accession URL and verify release status; do not depend on reviewer-only access links. Obtain raw counts rather than the reported FPKM values for DESeq2.
- **Integration:** eligible in principle as a separate mouse liver comparison after access and species handling are resolved.
- **Sources:** [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1203061), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12214484/), [DOI](https://doi.org/10.1038/s41598-025-05489-w).

### M02 — DDC gestational cholestasis / placental RNA-seq; accession needed

- **Study:** Xue et al. (2026), “Mechanisms of fetal growth restriction in gestational cholestasis: role of gut microbiota and placental redox.”
- **Design:** pregnant C57BL/6 mice, 0.1% DDC diet from E0.5 to E18.5 versus control diet. Placental RNA-seq used **4 control and 6 DDC placentas**. The overall design reports 6 dams per group.
- **Replication:** placentas were selected from fetuses near the litter mean. Recover dam/litter IDs; do not assume that every placenta is an independent maternal exposure replicate.
- **Availability:** no verified raw-read/count accession. Article/supplement availability and author contact are the starting points.
- **Integration:** potentially useful as a separate mouse placental model. Document DDC's toxic/obstructive cholestasis mechanism and distinguish it from clinical ICP and estrogen-based models.
- **Sources:** [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13046547/), [DOI](https://doi.org/10.3389/fvets.2026.1815498).

### M03 — GSE314547 / PRJNA1391489; humanized estrogen-cholestasis model

- **GEO title:** “Estradiol disrupts constitutive androstane receptor (CAR)-mediated bile acid homeostasis in a murine model of intrahepatic cholestasis of pregnancy (ICP), and this effect is reversed by CAR ligands.” Public GEO release: June 2026.
- **Assay/design:** liver bulk RNA-seq, 16 samples across four groups, n=4 each: `CV` vehicle control, `CM_2.5` MI763F alone, `EV` ethinylestradiol plus vehicle, `EM_2.5` ethinylestradiol plus MI763F. Confirm all labels against individual GSM records.
- **Critical eligibility limits:** sequenced animals are humanized PXR-CAR-CYP3A4/3A7 mice. The GEO record does not establish a WT sequenced arm or pregnancy status. Do not present this as verified WT pregnant controls versus pregnant ICP mice.
- **Availability:** GEO lists `GSE314547_raw_data_merged_filtr.txt.gz`; validate that it contains integer gene counts. SRA records are linked.
- **Integration:** related mechanistic cohort. Disease-model contrast: `EV` versus `CV`; intervention contrast: `EM_2.5` versus `EV`. Keep separate from primary gestational comparisons until model eligibility is resolved.
- **Sources:** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314547), [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1391489).

### M04 — GSE70179 / PRJNA287800; Abcb11 genetic context, microarray

- **Study:** “Maternal bile acid transporter deficiency promotes neonatal demise” (2015).
- **Assay:** Affymetrix expression microarrays, not bulk/scRNA-seq; 12 arrays across different tissues and experimental questions.
- **Design:** liver comparison includes 3 WT and 3 Abcb11-knockout young mice; pregnancy is not established for that liver arm. Offspring lung comparisons span E17.5, E18.5, and birth, with offspring of WT versus Abcb11-deficient dams. There is one deposited array per condition/time point in the lung arm; recover pooling/biological-unit details.
- **Integration:** retain as genetic/mechanistic context. Do not merge liver and lung or analyze it as replicated WT-versus-ICP RNA-seq. A separate array-platform workflow would be needed for expression reanalysis.
- **Sources:** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70179), [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA287800), [PubMed](https://pubmed.ncbi.nlm.nih.gov/26416771/).

## Other species

### O01 — Rat PRJCA058946 / PRJDB44845 / DRP022399; provisional RNA-seq labels

- **Associated paper:** “Resveratrol ameliorates intrahepatic cholestasis of pregnancy by modulating the gut-liver axis and FXR-mediated bile acid homeostasis” (2026).
- **Repository metadata:** records associated with the project label **18 libraries as RNA-Seq**, `TRANSCRIPTOMIC` / cDNA / total RNA: 8 liver and 10 ileum, female rats at gestational day 21. The project also contains **24 amplicon libraries**, which must not be counted as host RNA-seq.
- **Identifiers:** the paper lists GSA accessions `CRA039245` and `CRA039246`; verify which assay each contains. DDBJ RNA-labeled run range `DRR1259605`–`DRR1259622` maps to original runs `CRR2750500`–`CRR2750517`. Avoid downloading both mirror copies.
- **Sample labels:** liver `L47`, `L48`, `L51`, `L52`, `L53`, `L82`, `L84`, `L85`; ileum `V4`, `V11`, `V14`, `V16`, `V18`, `V43`, `V44`, `V45`, `V47`, `V50`.
- **Unresolved:** labels do not establish control, ICP, or resveratrol assignment, and the paper describes microbiome sequencing without clearly documenting these host transcriptomes. Confirm actual host-RNA assay identity, read quality/depth, strain, biological replication, and group mapping. Eighteen libraries do not establish 18 independent animals.
- **Disposition:** hold; do not compute ICP DE from guessed labels.
- **Sources:** [CNCB BioProject](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA058946), [GSA CRA039245](https://ngdc.cncb.ac.cn/gsa/browse/CRA039245), [GSA CRA039246](https://ngdc.cncb.ac.cn/gsa/browse/CRA039246), [BioProject mirror](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB44845), [DDBJ study](https://ddbj.nig.ac.jp/resource/sra-study/DRP022399), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13223021/), [DOI](https://doi.org/10.3389/fimmu.2026.1819374).

### O02 — Pig GSE110414 / PRJNA433661 / SRP132585; developmental reference

- **Study:** “Transcriptome Profiling of Placenta through Pregnancy Reveals Dysregulation of Bile Acids Transport and Detoxification Function” (2019).
- **Assay/design:** placental total-RNA sequencing for mRNA/lncRNA, 10 samples: 3 at gestation day 60, 3 at day 90, and 4 at farrowing (`L0`).
- **Availability:** public GEO/SRA with processed expression files; verify whether processed values are counts or normalized expression.
- **Limitation:** the deposited comparison is gestational stage, not documented healthy-versus-ICP case/control status. Developmental bile-acid changes do not justify relabeling late-gestation samples as ICP.
- **Integration:** optional developmental reference, excluded from primary disease DE and disease-cohort counts.
- **Sources:** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110414), [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA433661), [SRA study](https://www.ncbi.nlm.nih.gov/sra/?term=SRP132585), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6747679/), [DOI](https://doi.org/10.3390/ijms20174099).

### O03 — Pig GSE214238 / PRJNA884413 / SRP399730; ChIP-seq

- **Study title in repository:** “H4K12 Lactylation-activated c-JUN signaling aggravates bile acid dysfunction in intrahepatic cholestasis during pregnancy.”
- **Assay:** anti-H4K12la ChIP-seq, DNA/epigenetic profiling of pig liver/hepatocytes in high- versus low/normal-bile-acid late pregnancy.
- **Design:** four libraries: low-BA IP, low-BA input, high-BA IP, high-BA input. Input controls are not additional biological replicates. Establish the true animal and pooling structure if pursuing this study.
- **Integration:** exclude from RNA-expression, DESeq2 gene-count, and immune-cell expression analyses. Keep only as a possible separate epigenetic evidence source.
- **Sources:** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214238), [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA884413), [SRA study](https://www.ncbi.nlm.nih.gov/sra/?term=SRP399730).

### O04 — Rat PRJNA1002065; associated OMIX004988 deposition; metagenomics

- **Study:** “Roseburia intestinalis relieves intrahepatic cholestasis of pregnancy through bile acid/FXR-FGF15 in rats” (2023).
- **Design:** pregnant rat estrogen-induced ICP model with healthy/model and intervention arms.
- **Assay:** fecal shotgun metagenomic DNA sequencing; tissue expression measurements include targeted qPCR. Neither metagenomics nor qPCR is host bulk RNA-seq or scRNA-seq.
- **Availability:** the paper reports PRJNA1002065 and OMIX004988. Verify which files/assays each deposition holds instead of assuming they are exact duplicate repositories.
- **Integration:** exclude from host transcriptomic DE and immune-expression inference; retain for a separate microbiome mechanism analysis if desired.
- **Sources:** [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1002065), [OMIX](https://ngdc.cncb.ac.cn/omix/release/OMIX004988), [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10679810/).

## Reused cohorts and accession handling

- The 2024 reanalysis in [Journal of Inflammation Research, DOI 10.2147/JIR.S471374](https://doi.org/10.2147/JIR.S471374) discusses bulk `SRP378990` and single-cell `SRP429602`. Record it as a secondary analysis linked to H02/H03, not two additional cohorts. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11451404/).
- A paper that reanalyzes GSE282182 is not a new cohort. Retain GSE282182 as the primary dataset identity.
- Resolve H07–H09 against original accessions before counting them as independent evidence. Similar author lists or biological findings do not prove cohort identity; verify it.
- Keep primary accession, repository mirrors, publications, and reanalysis papers in separate metadata fields. Verify H02's proposed alias relationship and O04's associated deposit at the file/run level.
- Do not include incidental liver-cancer, generic cholestasis, normal tissue atlas, or microbiome reference datasets as additional ICP case-control cohorts merely because a cited paper uses them.

## Implementation order and boundaries

1. **Add the catalog and an accession/sample audit.** Suggested locations: `docs/dataset_catalog.md` and a machine-readable manifest under `resources/` or the repository's chosen metadata directory. Preserve the H/M/O IDs above as provenance IDs. Record unresolved and excluded entries too.
2. **Audit existing H01/H10.** Prevent duplicate configs and sample inflation. Preserve their existing outputs and modeling choices unless a documented correction is needed.
3. **Verify and add H02 first** if donor mapping, raw counts, and access are available. If only FASTQs exist, create a reproducible count-generation workflow with a recorded reference genome/annotation before configuring DE.
4. **Assess H04 separately.** Its cfRNA compartment, prospective sampling, access conditions, and disease timing require their own comparisons. Do not pool it with placenta.
5. **Add single-cell support separately** for H03 and any recovered H05 data. Keep donor IDs, cell-level QC, doublet handling, chemistry, and cell annotations. Use donor-level pseudobulk for eligible replicated comparisons. Treat H03 disease findings as exploratory.
6. **Add animal models only with explicit species/model support.** M01/M02 are the closest gestational bulk candidates discussed, but access and dam/litter metadata remain blockers. M03 is a separate related model. Keep O01 pending until its assay and groups are resolved.
7. **Keep alternate assays and reference-only designs separate.** H06, H10, M04, O02, O03, and O04 must not be silently converted into ordinary ICP mRNA-seq case-control inputs.

Proceed with records that can be verified and analyzed appropriately. Document blockers for the others rather than filling missing labels, accessions, or covariates with guesses.

### Minimum manifest fields

| Field group | Required information |
| --- | --- |
| Identity | Catalog ID, primary accession, verified aliases/mirrors, study title, DOI/PMID, primary source URLs, last verification date |
| Eligibility | Existing/new/pending/excluded status, reason, direct ICP case-control versus related model/reference, independent cohort versus reanalysis |
| Biology | Species/taxon, strain/genotype, tissue and compartment, pregnancy/gestational timing, disease criteria or induction, severity, treatment/exposure |
| Replication | Donor/animal ID, dam/litter ID when relevant, specimen ID, pooling, library ID, run ID, technical replicates, group assignment |
| Covariates | Maternal age and gestational age where available, fetal sex for placenta, batch, treatment such as UDCA, delivery/sampling details, relevant comorbidities |
| Assay | Library strategy/source/selection, platform, single-cell chemistry when relevant, raw versus normalized matrix, feature type, reference genome and annotation |
| Access and files | Public/controlled/unverified status, raw/count/metadata file locations, checksums, file-to-sample mapping, restrictions or missing inputs |
| Analysis | Contrast direction, biological unit, covariates, ortholog/reference handling, QC exclusions, limitations |

Use explicit missing values and notes; do not invent data to satisfy this schema. Extend the repository's config/loader only where necessary for verified inputs rather than adding unsupported keys and assuming they work.

### Statistical and species requirements

- Run case-control analyses within each study, species, tissue, and assay. Do not concatenate raw counts across species or incompatible compartments. Use a justified meta-analysis of study-level effects for cross-study synthesis.
- Confirm count matrices are nonnegative integers with stable gene IDs and sample columns matching the audited manifest. Rebuild counts from reads where required; never round TPM/FPKM to manufacture raw counts.
- Collapse technical runs for the same biological library appropriately. Do not treat run count, cell count, pooled donor count, or multiple placentas from one litter as independent donor/dam replication.
- Use species-appropriate immune markers or a versioned ortholog mapping with one-to-many/ambiguous mappings handled explicitly. Do not apply the Human Immune Cell Atlas to nonhuman symbols unchanged or uppercase mouse gene names as a substitute for mapping.
- Report immune marker expression scores as scores. Whole-tissue expression changes and decidual contamination can alter them without establishing changes in immune-cell proportions.
- Preserve study-specific covariates and treatment distinctions. Add a batch or other design term only when metadata and replication make the model identifiable.
- Prefer FDR-aware interpretation for new discovery claims and show nominal versus adjusted results explicitly. Separate descriptive n=1 single-cell findings from replicated inference.

### Acceptance criteria for the follow-on Codex work

- Every newly analyzed dataset has a verified accession/file source, donor/specimen-to-library mapping, count/assay validation, and stated contrast.
- Existing H01/H10 cohorts are not duplicated; pending entries remain visibly pending and exclusions retain their reasons.
- New eligible bulk configs follow `config/GSE282182.yaml`; supported Agilent inputs follow `config/GSE46157.yaml`. Unsupported modalities/species receive an explicit implementation rather than a misleading config-only addition.
- Run `uv run pytest` for relevant code changes and `./run_all.sh config/<NEW_DATASET>.yaml` for each newly supported dataset; investigate failures. The export consistency check must pass.
- Update documentation and the dashboard's dataset labels with species, tissue, assay, model, sample units, and key limitations. Keep processing reproducible and raw data handling consistent with the repository's data policy.
- End the implementation report with datasets added, datasets already present, pending inputs, exclusions, and validation results. Do not claim that all 18 entries were analyzed if only a subset qualifies.

