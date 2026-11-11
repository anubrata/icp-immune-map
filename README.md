# ICP Immune Map

A dashboard of the genes and immune cell types whose expression changes between healthy pregnancies and
pregnancies with **intrahepatic cholestasis of pregnancy (ICP)**, built from public datasets. Immune cells come from the
Human Immune Cell Atlas: the dashboard shows which immune cell populations have changed marker genes, and scores each
population per patient (ImmPort immune gene lists can be switched on per dataset). It is published with
GitHub Pages and rebuilt automatically on every push.

**Dashboard: https://anubrata.github.io/icp-immune-map/**

Datasets (all from [GEO](https://www.ncbi.nlm.nih.gov/geo/), NCBI's public gene-expression archive):

| Dataset | Species | Tissue | Data | Samples | Method |
| --- | --- | --- | --- | --- | --- |
| [GSE282182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282182) | Human | Placenta | RNA-seq read counts | 5 healthy, 5 ICP | DESeq2, not adjusted for fetal sex |
| [GSE46157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46157) | Human | Placenta | Agilent microarray (4x44K v2) | 2 healthy, 4 ICP (2 mild, 2 severe) | moderated t-test, not adjusted for fetal sex |
| [GSE314547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314547) | Mouse (model) | Liver | RNA-seq read counts | 4 vehicle, 4 ethinylestradiol (humanized PXR-CAR-CYP3A4/3A7 mice) | DESeq2; immune markers via one-to-one human orthologs |

`resources/icp_dataset_handoff.md` lists the other candidate ICP datasets (human, mouse, rat, pig) and why each is or
isn't ready for this analysis.

## How it works

```
data/<file>  +  config/<dataset>.yaml
        │
        ├─ step 1  analysis/step1_load_data.py               any file layout → standard tables
        ├─ step 2  analysis/step2_differential_expression.py ICP vs healthy: DESeq2 (counts) or
        │                                                    moderated t-test (microarrays)
        ├─ step 3  analysis/step3_immune_annotation.py       tag immune cell marker genes (atlas; ImmPort optional;
        │                                                    mouse genes through human orthologs),
        │                                                    score immune cell types (atlas markers)
        ├─ step 4  analysis/step4_export_site_data.py        significant genes + cell types → site/data/*.json
        └─ step 5  analysis/step5_check_site_data.py         recompute every number independently and compare
                                                                       │
site/index.html + app.js + style.css  (a fixed page) ◄── reads ────────┘
```

Each step is a short, commented script that can be run on its own and writes readable tables to
`results/<dataset>/` (open them in Excel). The website itself is a fixed HTML page; only its data
files in `site/data/` are generated. Step 5 recomputes everything the dashboard shows with its own code
(significant genes, statistics, per-patient values, arrows, immune tags, bar counts, tiles, cell-type
scores) and fails if anything differs beyond rounding, so a wrong number never gets published.

## Run it on your computer

Install [uv](https://docs.astral.sh/uv/) once (it installs the right Python and packages for you), then:

```bash
./run_all.sh                                        # steps 1-5 for every dataset in config/
uv run python -m http.server --directory site 8000  # then open http://localhost:8000
```

Run one step at a time with, for example, `uv run analysis/step2_differential_expression.py config/GSE282182.yaml`.
Run the tests (synthetic data only) with `uv run pytest`.

## Add a dataset

1. Put the data in `data/`. Two kinds are supported:
   - **RNA-seq raw read counts** (`data_type: counts`): one table (`.txt`, `.tsv`, `.csv` or `.xlsx`, optionally `.gz`),
     one row per gene and one column per sample. Already-normalized values (TPM, FPKM) cannot be used. Human or mouse
     (`species`); a mouse dataset needs `gene_annotation: resources/mouse_genes.tsv.gz` (made with
     `uv run analysis/step0_fetch_mouse_genes.py`).
   - **Agilent one-color microarrays** (`data_type: microarray`): the raw Feature Extraction file of each array
     (as on GEO, e.g. `GSM1125156_Healthy-1.txt.gz`), plus the platform's probe annotation, made once with
     `uv run analysis/step0_fetch_platform_annotation.py <GPL id>`.
2. Copy `config/GSE282182.yaml` (counts) or `config/GSE46157.yaml` (microarray) to `config/<new id>.yaml` and edit it:

   | Setting | Meaning |
   | --- | --- |
   | `id`, `title`, `tissue`, `geo_accession` | Shown on the site |
   | `species` | `human` (default) or `mouse` |
   | `notes` | Optional: a sentence or two on the design and its limits, shown under the dataset's title |
   | `data_type` | `counts` (default) or `microarray` |
   | `file` | (counts) Path of the data file, relative to the repository root. If its header has no name for the first (gene ID) column, as in files written by R, that column is used as `columns.gene_id` |
   | `gene_annotation` | (counts) Optional table adding gene symbols, names, chromosomes and human orthologs by gene ID; required for mouse (`resources/mouse_genes.tsv.gz`) |
   | `files`, `platform_annotation` | (microarray) One file per sample, and the probe annotation from step 0 |
   | `read_options` | Extra pandas options for unusual layouts, e.g. `{sep: ";", skiprows: 2}` |
   | `columns.gene_id` | (counts) Column with a unique ID per row |
   | `columns.gene_symbol`, `description`, `chromosome`, `biotype` | (counts) Optional annotation columns |
   | `samples` | Every sample → its group |
   | `reference_group`, `compare_group` | Fold changes are compare vs reference (e.g. ICP vs Healthy) |
   | `fetal_sex` | `infer` (counts only, from chrY genes), a per-sample map of `M`/`F`, or leave it out |
   | `adjust_for_fetal_sex` | `true` (default): add fetal sex to the model when both sexes are present; `false`: record it only (chrY genes are then left out) |
   | `use_immport` | `true` to also tag genes with ImmPort's immune gene lists (default `false`: atlas only) |
   | `composition_check` | Optional `name` + `markers` of a tissue-composition signature to check cell-type scores against (placenta: `Decidua`, `[PRL, IGFBP1]`) |
   | `min_count` | (counts) Genes need this many reads in at least as many samples as the smallest group |
   | `cutoffs` | `pvalue` (unadjusted) or `padj` (FDR-adjusted), plus `abs_log2fc` |

3. Run `./run_all.sh config/<new id>.yaml` to check it, then push. Step 1 stops with a plain-English
   message if the file and config don't match.

## Publish on GitHub Pages

1. Create an empty **public** repository on GitHub and push this folder to its `main` branch.
2. In the repository, open **Settings → Pages** and set **Source** to **GitHub Actions**.

From then on, every push runs `.github/workflows/pages.yml`: it runs the tests, runs `./run_all.sh` (including the step 5
cross-check), and publishes `site/` only if all of that passes.

## Methods and caveats

- **Differential expression, RNA-seq:** DESeq2 (negative binomial model, Wald test) via
  [pydeseq2](https://github.com/owkin/PyDESeq2). Genes need at least 10 reads in at least 5 samples to be tested.
- **Differential expression, microarrays:** each probe's processed signal is log2-transformed and quantile-normalized
  across arrays; probes clearly above background on at least as many arrays as the smallest group and measuring a named
  gene are kept (the brightest probe per gene); then a limma-style moderated t-test (Smyth 2004, implemented in
  `analysis/moderated_t.py`) compares the groups. R's limma itself is not used.
- **Fetal sex:** the placenta is fetal tissue. For now no dataset is adjusted for fetal sex (`adjust_for_fetal_sex: false`):
  the first goal is immune-cell changes between ICP and healthy, and sex-specific changes come later. Fetal sex is still
  recorded and shown (dot shapes on the site), and chrY genes are left out of the test because they would only reflect
  which placentas are male. In GSE282182 the groups differ in fetal sex (healthy: 4 female, 1 male; ICP: 2 female,
  3 male, inferred from chrY genes); adjusting for it (`~ fetal_sex + condition`, set `adjust_for_fetal_sex: true`)
  changes some genes a lot (e.g. CXCL9 −7.6× with adjustment vs −1.3× without) but not the NK-cell signal. In GSE46157
  the two mild-ICP placentas are male and the other four female (from the DDX3Y signal), so fetal sex coincides with severity.
- **Significance:** a gene is shown when its **unadjusted** p-value is below 0.05 and it changes at least 2-fold
  (|log2 fold change| ≥ 1). With ~17,000 genes tested, some genes pass this by chance, so treat the list as
  candidates; the FDR (adjusted p-value) is shown next to every gene. Cutoffs are set per dataset in its config.
- **Arrow chart:** each arrow runs from the gene's level in healthy to its level in ICP as estimated by the model
  (with any other model terms, such as fetal sex when adjusted for, held at the cohort average), so it always matches
  the fold change in the table. Each patient's actual values are in the gene panel.
- **⚠ genes:** directions follow the model. If a gene's plain group averages point the other way (possible when fetal
  sex is adjusted for and the groups have different fetal-sex mixes, or when one sample is extreme), it is marked ⚠ and
  the gene panel explains why. With the current settings no gene is flagged.
- **Small study:** 5 vs 5 samples. In this dataset the overall difference between the groups is subtle
  (see "Quality checks" on the site).
- **Immune genes:** a gene is tagged immune-related if it is a marker of an immune cell population in the Human Immune
  Cell Atlas (`resources/human_immune_cell_atlas.txt`). The atlas's cross-lineage functional-state signatures
  (proliferation, hypoxia, …) don't count, as the atlas notes they are not immune-specific. The dashboard's grid
  "Immune cells and their changed marker genes" shows each atlas cell population with its significant marker genes.
  Genes are matched by their symbol and by their current HGNC symbol (via Ensembl ID), so renamed genes still match.
  After editing the atlas, run `uv run analysis/step0_extract_atlas_markers.py`.
- **Mouse datasets:** the atlas lists human genes, so a mouse gene counts as a marker only through its one-to-one human
  ortholog (from MGI's mouse-human homology report, built by `analysis/step0_fetch_mouse_genes.py`), never because the
  names look alike. Genes whose ortholog is ambiguous (several genes on either side, e.g. Lyz1/Lyz2 for LYZ, the Ccl3/Ccl4
  cluster) and human genes mice lack (e.g. GNLY) cannot be measured, so some populations have too few markers to score; in
  GSE314547 (mouse liver) the NK-cell panel is down to Nkg7 and is not scored. GSE314547 is a model, not patients:
  ethinylestradiol-treated vs vehicle-treated humanized PXR-CAR-CYP3A4/3A7 mice, pregnancy not stated; the study's two
  MI-763F (CAR ligand) groups are not included. Its strongest cell-type scores are the "proliferating" populations, which
  share proliferation genes (MKI67, TOP2A) that are also high in growing liver cells, so read them with care.
- **ImmPort (optional, off):** with `use_immport: true` in a dataset config, genes are also tagged with ImmPort's ~150
  immune gene lists (GO "immune system process" and Reactome "Immune System" sub-terms), count as immune-related when in
  one, and the dashboard adds an ImmPort card. The lists are in `resources/immport_gene_lists.gmt`; refresh them (and the
  HGNC table) with `uv run analysis/step0_fetch_gene_lists.py`.
- **Immune cell types:** for every atlas population with at least 2 marker genes among the tested genes, each patient
  gets a score (the mean of the markers' variance-stabilized expression, each gene scaled to mean 0 and SD 1). Scores
  are compared with `score ~ group` (`score ~ fetal_sex + group` when fetal sex is adjusted for), giving the difference
  in SD with a 95% CI, a p-value and an FDR across all scored populations. In GSE282182 no population passes FDR < 0.05;
  NK cells show the strongest signal (−1.26 SD in ICP, p = 0.007; p = 0.02 after the decidua check), with no change in
  the pan-immune marker PTPRC, so this points to fewer NK cells rather than fewer immune cells overall. Because placenta
  samples hold varying amounts of decidua (rich in NK cells), each score is also
  checked against a decidua signature (PRL, IGFBP1) and re-tested adjusted for it. Marker panels overlap (e.g. PRF1 and
  GNLY mark both NK and cytotoxic T cells), so related populations move together.

## Project layout

| Path | What it is |
| --- | --- |
| `analysis/` | Steps 1–5, the step 0 scripts that prepare the reference files, and shared helpers |
| `config/` | One YAML file per dataset |
| `data/` | Raw data files |
| `resources/` | ImmPort immune gene lists, HGNC symbol table, the Human Immune Cell Atlas and its marker table, the microarray probe annotation, the mouse gene table with human orthologs, and the dataset catalog handoff (see `resources/README.md`) |
| `results/` | Intermediate tables from each step (generated, not committed) |
| `site/` | The dashboard page; `site/data/` is generated by step 4 |
| `tests/` | Tests on synthetic data |

## Sources

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
  *Genome Biology* 15, 550 (2014).
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential
  expression analysis. *Bioinformatics* 39, btad547 (2023).
- Human Immune Cell Atlas: Gene Markers for Cell Types and Functional States (`resources/human_immune_cell_atlas.txt`).
- [ImmPort gene lists](https://docs.immport.org/apidocumentation/genelist/overview/) · [HGNC](https://www.genenames.org/) ·
  [MGI](https://www.informatics.jax.org/) and [Alliance of Genome Resources](https://www.alliancegenome.org/) (mouse genes and
  mouse-human homology) · [GEO GSE282182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282182) ·
  [GEO GSE46157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46157) · [GEO GSE314547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314547)
