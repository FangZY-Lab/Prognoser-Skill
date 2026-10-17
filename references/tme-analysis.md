# Tumour-microenvironment analysis

Downstream module 2 of 2. It characterises the microenvironment of the same
cohorts used for `Prognoser()`, and relates it to the signature.

Three views are produced together because each answers a different question and
each fails differently:

| View | Question it answers | Method |
| --- | --- | --- |
| Bulk scores | How much stroma and immune content is there overall? | ESTIMATE |
| Deconvolution | Which cell types, in what proportions? | xCell, optional quanTIseq |
| Marker panels | Is a specific programme (checkpoints, TLS, CAF) switched on? | ssGSEA on curated markers |

Reporting all three is deliberate. A high `estimate:ImmuneScore` with low
`quantiseq:T.cells.CD8` is a contradiction worth resolving before writing a
conclusion, and it is invisible if only one method is run.

## ESTIMATE

`estimate::estimateScore()` yields `StromalScore`, `ImmuneScore` and
`ESTIMATEScore` (stromal + immune). Higher ESTIMATEScore means a
microenvironment-rich, purity-low sample.

Two implementation details matter and are handled in `dx_estimate_scores()`:

1. **Input must pass through `filterCommonGenes()` first.** `estimateScore()`
   expects the GCT file that `filterCommonGenes()` writes. Feeding it a plain
   tab-delimited matrix does not fail loudly; it returns shifted sample names
   and `NA` scores. The helper always stages a temporary GCT file.
2. **`platform` selects a gene list, not a technology guess.** Legitimate values
   are `affymetrix`, `agilent` and `illumina`. For RNA-seq there is no official
   option; `illumina` is the common substitute and is the default here. Say so
   in the methods section, and check that the number of common genes reported
   (logged, and written to `tme_estimate_scores.csv`) is in the expected range.

**Tumour purity is only returned for `platform = "affymetrix"`.** With
`illumina`/`agilent` input the module writes three scores and no purity row. Do
not report a purity column that is not there.

## xCell

`xCell::xCellAnalysis()` returns 64 cell-type enrichment scores plus
`ImmuneScore`, `StromaScore` and `MicroenvironmentScore` (67 rows), with
spillover compensation applied between closely related cell types.

- `rnaseq = TRUE` (default) expects **linear TPM**, not log2. The module
  back-transforms log2 input with `2^x - 1` and documents that this is exact
  only when the matrix is log2(TPM + 1).
- xCell needs a reasonably complete transcriptome. The wrapper refuses to run on
  fewer than 5 000 genes and logs how many signature genes actually matched.
  On a 12 000-gene test matrix about 6 500 genes matched, which is workable; on
  a 5 000-gene panel it is not.
- xCell scores are **relative enrichment**, not fractions, and they are not
  comparable between cohorts run separately.

## quanTIseq (optional)

`quantiseqr::run_quantiseq()` returns ten immune cell fractions (`B.cells`,
`Macrophages.M1`, `Macrophages.M2`, `Monocytes`, `Neutrophils`, `NK.cells`,
`T.cells.CD4`, `T.cells.CD8`, `Tregs`, `Dendritic.cells`) plus an `Other`
residual. Fractions describe the **immune compartment**, not the whole tumour,
so they are not directly comparable with ESTIMATE percentages.

It is off by default (`--quantiseq=true` enables it) for one reason: it needs
linear TPM. Passing log2 values or counts produces a successful return value
with meaningless fractions — there is no error to catch. Enable it only when
the input matrix is known to be TPM, and inspect the proportion of signature
genes matched (the package prints it).

## Marker panels

Nine short curated panels, scored with the same ssGSEA engine as the cell-state
module:

| Panel | Genes |
| --- | --- |
| `cd8_t_cell` | CD8A, CD8B, GZMA, GZMB, GZMK, PRF1, NKG7, GNLY |
| `treg` | FOXP3, IL2RA, CTLA4, IKZF2, TNFRSF18 |
| `m1_macrophage` | NOS2, IL12B, CXCL9, CXCL10, CXCL11, TNF |
| `m2_macrophage` | CD163, MRC1, MSR1, CCL22, IL10, TGFB1, ARG1 |
| `mdsc` | S100A8, S100A9, ITGAM, ARG1, CXCR2 |
| `checkpoint` | PDCD1, CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, TIGIT, IDO1, VSIR, BTLA |
| `tls` | CXCL13, LTB, CCL19, CCL21, CCR7, SELL, CXCR5, CD79B |
| `caf` | FAP, ACTA2, PDGFRB, COL1A1, COL1A2, COL3A1, POSTN, TAGLN, THY1, DCN |
| `endothelial` | PECAM1, VWF, CDH5, CLDN5, EMCN, KDR |

These are **markers, not deconvolution**. They answer "is this programme on",
not "how many of these cells are present". A panel is scored only if at least
five of its genes are measured; otherwise it is listed with `scored = FALSE` in
the coverage output rather than silently dropped.

`checkpoint` deserves a warning: checkpoint transcript level is not a predictor
of immune-checkpoint-inhibitor response. Do not present it as one.

## Immune phenotype

`dx_immune_phenotype()` assigns each sample to `inflamed`, `excluded` or
`desert` from two within-cohort z-scores:

- CD8 T-cell panel score (cytotoxic content)
- IFN-gamma response panel score (T-cell activation)

Rules, with tertiles of CD8 and the median of IFN-gamma:

| CD8 level | IFN-gamma | Call |
| --- | --- | --- |
| bottom tertile | any | `desert` |
| middle / top tertile | below median | `excluded` |
| middle / top tertile | at or above median | `inflamed` |

This is a **heuristic** implementation of the Chen & Mellman framework, not the
TIDE algorithm and not a validated classifier. `tme_immune_phenotype_thresholds.csv`
records the cut points actually used per cohort, so the calls can be audited or
re-cut. Because the thresholds are cohort-relative, the phenotype proportions
are not comparable across cohorts with different immune baselines.

## Confounders worth stating explicitly

- **Tumour purity.** Every immune and stromal score is partly a purity estimate.
  A signature that correlates with `estimate:StromalScore` may simply mark
  low-purity samples. If purity is available (affymetrix platform, or an
  external estimate), test the association with purity as a covariate.
- **Compositional data.** Cell fractions sum to a constant, so a rise in one
  compartment forces a fall elsewhere. A negative correlation between CD8
  fraction and fibroblast fraction can be entirely compositional. The module
  reports fractions as-is; if composition matters, transform to log-ratios
  before modelling.
- **Platform and normalisation.** xCell, quanTIseq and ESTIMATE make different
  assumptions about the input scale. The module logs the detected scale and
  accepts `--expr-scale=log2|linear` to override the guess — check that line in
  the log before interpreting anything.
- **Cross-cohort comparison.** Scores are computed per cohort. Hazard ratios and
  correlations can be compared; absolute scores cannot.

## What the module reports

| File | Content |
| --- | --- |
| `tme_estimate_scores.csv` | Stromal / Immune / ESTIMATE scores per sample |
| `tme_xcell_scores.csv` | 64 cell types + 3 composite scores |
| `tme_quantiseq_fractions.csv` | Immune fractions (when enabled) |
| `tme_marker_scores.csv` | Marker panel scores per sample |
| `tme_immune_phenotype.csv` | Per-sample phenotype call and inputs |
| `tme_immune_phenotype_thresholds.csv` | Cut points used per cohort |
| `tme_survival.csv` | Cox per SD of every continuous feature |
| `tme_signature_correlation.csv` | Spearman of the signature score vs every feature |
| `tme_markers_heatmap_<dataset>.pdf` | Panel scores, z-scored |
| `tme_survival_forest.pdf` | Forest plot for ESTIMATE and marker features |

## Dependency notes

`estimate` and `xCell` ship lazy-loaded data sets that their own functions
resolve from the search path, so `pkg::fun()` alone fails with errors such as
`object 'common_genes' not found`. The module attaches both packages
(`dx_attach()`), which is the documented usage. Missing optional packages
degrade one view rather than failing the run.
