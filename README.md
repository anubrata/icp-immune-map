# ICP Immune Map

A dashboard of the genes and immune cell types whose expression changes between healthy pregnancies and
pregnancies with **intrahepatic cholestasis of pregnancy (ICP)**, built from public datasets. Immune-related genes are
tagged with ImmPort immune gene lists and the Human Immune Cell Atlas markers, and the atlas marker panels are used to
score immune cell types per patient. It is published with GitHub Pages and rebuilt automatically on every push.

**Dashboard: https://anubrata.github.io/icp-immune-map/**

The first dataset is [GEO GSE282182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282182): bulk RNA-seq of
placenta from 5 healthy and 5 ICP pregnancies.

## How it works

```
data/<file>  +  config/<dataset>.yaml
        │
        ├─ step 1  analysis/step1_load_data.py               any file layout → standard tables
        ├─ step 2  analysis/step2_differential_expression.py DESeq2 (pydeseq2): ICP vs healthy
        ├─ step 3  analysis/step3_immune_annotation.py       tag immune genes (ImmPort + atlas),
        │                                                    score immune cell types (atlas markers)
        └─ step 4  analysis/step4_export_site_data.py        significant genes + cell types → site/data/*.json
                                                                       │
site/index.html + app.js + style.css  (a fixed page) ◄── reads ────────┘
```

Each step is a short, commented script that can be run on its own and writes readable tables to
`results/<dataset>/` (open them in Excel). The website itself is a fixed HTML page; only its data
files in `site/data/` are generated.

## Run it on your computer

Install [uv](https://docs.astral.sh/uv/) once (it installs the right Python and packages for you), then:

```bash
./run_all.sh                                        # steps 1-4 for every dataset in config/
uv run python -m http.server --directory site 8000  # then open http://localhost:8000
```

Run one step at a time with, for example, `uv run analysis/step2_differential_expression.py config/GSE282182.yaml`.
Run the tests (synthetic data only) with `uv run pytest`.

## Add a dataset

1. Put the file in `data/` (`.txt`, `.tsv`, `.csv` or `.xlsx`, optionally `.gz`). It must hold **raw read counts**,
   one row per gene and one column per sample.
2. Copy `config/GSE282182.yaml` to `config/<new id>.yaml` and edit it:

   | Setting | Meaning |
   | --- | --- |
   | `id`, `title`, `tissue`, `geo_accession` | Shown on the site |
   | `file` | Path of the data file, relative to the repository root |
   | `read_options` | Extra pandas options for unusual layouts, e.g. `{sep: ";", skiprows: 2}` |
   | `columns.gene_id` | Column with a unique ID per row (required) |
   | `columns.gene_symbol`, `description`, `chromosome`, `biotype` | Optional annotation columns |
   | `samples` | Every sample column → its group |
   | `reference_group`, `compare_group` | Fold changes are compare vs reference (e.g. ICP vs Healthy) |
   | `fetal_sex` | `infer` (from chrY genes), a per-sample map of `M`/`F`, or leave it out for no adjustment |
   | `composition_check` | Optional `name` + `markers` of a tissue-composition signature to check cell-type scores against (placenta: `Decidua`, `[PRL, IGFBP1]`) |
   | `min_count` | Genes need this many reads in at least as many samples as the smallest group |
   | `cutoffs` | `pvalue` (unadjusted) or `padj` (FDR-adjusted), plus `abs_log2fc` |

3. Run `./run_all.sh config/<new id>.yaml` to check it, then push. Step 1 stops with a plain-English
   message if the file and config don't match.

## Publish on GitHub Pages

1. Create an empty **public** repository on GitHub and push this folder to its `main` branch.
2. In the repository, open **Settings → Pages** and set **Source** to **GitHub Actions**.

From then on, every push runs `.github/workflows/pages.yml`: it runs the tests, runs `./run_all.sh`, and publishes `site/`.

## Methods and caveats

- **Differential expression:** DESeq2 (negative binomial model, Wald test) via
  [pydeseq2](https://github.com/owkin/PyDESeq2). Genes need at least 10 reads in at least 5 samples to be tested.
- **Fetal sex:** the placenta is fetal tissue, and in GSE282182 the groups differ in fetal sex
  (healthy: 4 female, 1 male; ICP: 2 female, 3 male, inferred from chrY genes). The model is
  `~ fetal_sex + condition`, so sex-linked genes are not mistaken for ICP effects.
- **Significance:** a gene is shown when its **unadjusted** p-value is below 0.05 and it changes at least 2-fold
  (|log2 fold change| ≥ 1). With ~17,000 genes tested, some genes pass this by chance, so treat the list as
  candidates; the FDR (adjusted p-value) is shown next to every gene. Cutoffs are set per dataset in its config.
- **Arrow chart:** each arrow runs from the gene's level in healthy to its level in ICP as estimated by the DESeq2
  model, with fetal sex held at the cohort average, so it always matches the fold change in the table. Plain averages
  can point the other way when the groups' fetal-sex mix differs; each patient's actual values are in the gene panel.
- **Small study:** 5 vs 5 samples. In this dataset the overall difference between the groups is subtle
  (see "Quality checks" on the site).
- **Immune genes:** a gene is tagged immune-related if it is in any of ImmPort's ~150 immune gene lists (curated from
  GO "immune system process" and Reactome "Immune System") or is a marker of an immune cell population in the Human
  Immune Cell Atlas (`resources/human_immune_cell_atlas.docx`). The atlas's cross-lineage functional-state signatures
  (proliferation, hypoxia, …) don't count, as the atlas notes they are not immune-specific. Genes are matched by their
  symbol and by their current HGNC symbol (via Ensembl ID), so renamed genes still match. Refresh the ImmPort/HGNC
  files with `uv run analysis/step0_fetch_gene_lists.py`; after editing the atlas, run
  `uv run analysis/step0_extract_atlas_markers.py`.
- **Immune cell types:** for every atlas population with at least 2 marker genes among the tested genes, each patient
  gets a score (the mean of the markers' variance-stabilized expression, each gene scaled to mean 0 and SD 1). Scores
  are compared with `score ~ fetal_sex + group`, giving the difference in SD with a 95% CI, a p-value and an FDR across
  all scored populations. In GSE282182 no population passes FDR < 0.05; NK cells show the strongest major-lineage
  signal (lower in ICP). Because placenta samples hold varying amounts of decidua (rich in NK cells), each score is also
  checked against a decidua signature (PRL, IGFBP1) and re-tested adjusted for it. Marker panels overlap (e.g. PRF1 and
  GNLY mark both NK and cytotoxic T cells), so related populations move together.

## Project layout

| Path | What it is |
| --- | --- |
| `analysis/` | The four analysis steps, the step 0 scripts that prepare the reference files, and shared helpers |
| `config/` | One YAML file per dataset |
| `data/` | Raw data files |
| `resources/` | ImmPort immune gene lists, HGNC symbol table, the Human Immune Cell Atlas and its marker table (see `resources/README.md`) |
| `results/` | Intermediate tables from each step (generated, not committed) |
| `site/` | The dashboard page; `site/data/` is generated by step 4 |
| `tests/` | Tests on synthetic data |

## Sources

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
  *Genome Biology* 15, 550 (2014).
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential
  expression analysis. *Bioinformatics* 39, btad547 (2023).
- Human Immune Cell Atlas: Gene Markers for Cell Types and Functional States (`resources/human_immune_cell_atlas.docx`).
- [ImmPort gene lists](https://docs.immport.org/apidocumentation/genelist/overview/) · [HGNC](https://www.genenames.org/) ·
  [GEO GSE282182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282182)
