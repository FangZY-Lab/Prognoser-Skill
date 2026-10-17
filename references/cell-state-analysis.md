# Cell-state analysis

Downstream module 1 of 2. It answers a narrow question: **once a prognostic
signature exists, which cell-state programmes does it track, and do those
programmes carry survival information on their own?**

Everything here operates on bulk expression matrices, using the same
`<Dataset>_<Endpoint>` objects that `Prognoser()` consumed. No single-cell data
is required.

## What is scored

Thirteen cell-state panels, all resolved from MSigDB at run time through
`msigdbr`:

| Panel | Definition | What it captures |
| --- | --- | --- |
| `emt` | `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION` | Mesenchymal shift / dedifferentiation |
| `proliferation` | `HALLMARK_E2F_TARGETS`, `HALLMARK_G2M_CHECKPOINT`, `HALLMARK_MYC_TARGETS_V1` | Cycling fraction |
| `hypoxia` | `HALLMARK_HYPOXIA` | Hypoxic adaptation |
| `glycolysis` | `HALLMARK_GLYCOLYSIS` | Warburg-type metabolism |
| `oxphos` | `HALLMARK_OXIDATIVE_PHOSPHORYLATION` | Mitochondrial metabolism |
| `angiogenesis` | `HALLMARK_ANGIOGENESIS` | Vascular / VEGF programme |
| `apoptosis` | `HALLMARK_APOPTOSIS` | Cell-death competence |
| `dna_repair` | `HALLMARK_DNA_REPAIR` | Repair capacity |
| `senescence` | `GOBP_CELLULAR_SENESCENCE` | Senescence programme |
| `tgf_beta` | `HALLMARK_TGF_BETA_SIGNALING` | Immune exclusion / stromal driver |
| `ifn_gamma` | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | T-cell inflamed programme |
| `inflammation` | `HALLMARK_INFLAMMATORY_RESPONSE` | Innate inflammation |
| `stemness` | `MUELLER_PLURINET` | Pluripotency-derived stemness proxy |

Gene lists are **not** hard-coded. They are pulled from the MSigDB release
bundled with the installed `msigdbr`, so the panels follow the database rather
than a stale copy inside this repository. The exact genes used in a run are
written to `panel_definitions.csv`, and the `msigdbr` version is recorded in
`session_info.txt`.

## How scoring works

1. For each cohort, genes are intersected with the panel; panels with fewer
   than `min_size` (default 5) measured genes are **not scored** and are listed
   in `cell_state_coverage.csv` with `scored = FALSE`.
2. Remaining panels are scored with ssGSEA (`GSVA::ssgseaParam`, Gaussian kernel,
   within-sample normalisation). This is rank-based, so it tolerates the
   normalisation differences that survive the preprocessing step.
3. For cross-cohort figures only, panel scores are z-scored **within each
   cohort** before being placed on one heatmap.

`--method=zscore` swaps ssGSEA for a plain mean z-score over the panel genes.
That path needs no Bioconductor dependency and is easier to explain, but it
ignores the rank distribution and is more sensitive to single high-variance
genes.

**ssGSEA scores are relative within a cohort.** A `hypoxia` score of 0.41 in
cohort A and 0.37 in cohort B does not mean A is more hypoxic. Compare
correlations and hazard ratios across cohorts, never raw scores.

## What the module reports

| File | Content |
| --- | --- |
| `cell_state_coverage.csv` | Measured / nominal genes per panel, `scored` flag |
| `cell_state_scores_long.csv` | Long table: dataset, sample, panel, score |
| `cell_state_survival.csv` | Cox per SD of each panel, with BH-adjusted p |
| `cell_state_signature_correlation.csv` | Spearman of the signature score vs each panel |
| `signature_score_survival.csv` | Cox per SD of the HRGS-minus-LRGS score |
| `cell_state_heatmap.pdf` | Panels x samples, z-scored within cohort |
| `cell_state_forest_<dataset>.pdf` | Forest plot of panel hazard ratios |

## Interpreting the result

The signature score is `mean z(HRGS genes) - mean z(LRGS genes)`, so higher
means worse prognosis. A positive `rho` with `emt` or `proliferation` says the
signature is partly a mesenchymal or proliferative readout; a negative `rho`
with `ifn_gamma` says the signature is anti-correlated with T-cell
inflammation.

Read these as **descriptions of what the signature measures**, not as evidence
that the state causes outcome. Two panels are systematically confounded:

- `emt` and `angiogenesis` partly measure stromal content, which moves with
  tumour purity.
- `proliferation` tracks the malignant fraction, which also moves with purity.

Adjusting for purity is discussed in [tme-analysis.md](tme-analysis.md).

## Honest limitations

- **No malignant-cell state deconvolution.** If the biological question is
  "which malignant state does the signature mark", bulk panels are the wrong
  tool; use single-cell references (BayesPrism, MuSiC) on scRNA data and treat
  this module as a first pass.
- **`stemness` is a pluripotency proxy** (`MUELLER_PLURINET`), not a validated
  cancer stem-cell score such as a machine-learning stemness index.
- **`senescence` overlaps inflammation** through the SASP; a joint increase in
  both panels is not two independent findings.
- **Panels are not independent.** Hallmark sets share genes, so a signature
  correlating with `proliferation` usually correlates with `myc` targets too.
  Do not treat thirteen panels as thirteen independent tests without
  acknowledging the correlation structure in the discussion.

## Adding a panel

Add an entry to `DX_STATE_PANELS` in `scripts/downstream_signatures.R`:

```r
my_panel = list(
  label = "My programme",
  sets = c("HALLMARK_X", "GOBP_Y"),   # MSigDB names, unioned
  use  = "One sentence on what it does and does not measure."
)
```

Marker-style panels with an explicit gene vector belong in `DX_MARKER_PANELS`
instead. Run `Rscript downstream_cell_states.R --help` for the CLI flags, or
call `dx_cell_states_run()` directly from R.
