# Downstream statistics

Shared statistical contract for both downstream modules. Reading this first
avoids the most common misinterpretations.

## Data contract

Every downstream function takes an expression matrix with **genes in rows and
samples in columns**, and returns features in rows and samples in columns. No
function modifies `.GlobalEnv`; the modules load objects into a private
environment and pass them around explicitly.

Inputs are the same `<Dataset>_<Endpoint>` and `<Dataset>_<Endpoint>_survival`
objects that `Prognoser()` requires, so a completed Prognoser run already has
everything the downstream step needs.

Sample alignment is enforced when the objects are loaded: samples present in
expression but missing from survival (or the reverse) are dropped, and the count
dropped is logged. Silently keeping unmatched samples is how a survival analysis
ends up with `time` from the wrong patient.

## Signature score

```
signature = mean z(genes in HRGS*) - mean z(genes in LRGS*)
```

- `z` is a per-gene z-score across samples **within one cohort**.
- The two sides are averaged separately, then subtracted, so a signature with 40
  risk and 8 protective genes is not dominated by the larger side.
- Higher values mean worse prognosis, matching the HRGS/LRGS convention.

Because it is built from within-cohort z-scores, the signature score is not
comparable across cohorts. Compare its hazard ratio and its correlations.

## Association statistics

**Correlation** is Spearman, with `cor.test(..., exact = FALSE)` and `n >= 5`
required. Spearman is used because expression-derived scores are routinely
heavy-tailed and a single sample can move a Pearson coefficient. Ties are common
in ssGSEA output, which is why the exact p-value is not requested.

**Survival** is a univariable Cox model per one standard deviation of the
feature:

```r
coxph(Surv(time, status) ~ z)   # z = as.numeric(scale(feature))
```

The hazard ratio is therefore expressed per SD, which makes features measured on
different scales (a fraction between 0 and 1, an ESTIMATE score in the hundreds)
directly comparable in a forest plot.

Guardrails, applied per feature:

| Condition | Behaviour |
| --- | --- |
| `n < 10` or `events < 5` | returns `NA` instead of fitting |
| zero variance | returns `NA` |
| non-finite values | excluded pairwise before fitting |

**Proportional hazards is not tested.** `cox.zph` is not run, and with
`n` in the tens to low hundreds it would have little power anyway. Treat the HR
as a summary association, not a validated prognostic model.

## Multiple testing

Benjamini-Hochberg is applied **within each output table**, i.e. within the
cell-state family and within the TME family separately. `p_adj` in
`cell_state_survival.csv` is therefore not comparable with `p_adj` in
`tme_survival.csv`; they correct different families. If you report both, say so,
or re-adjust the combined set yourself.

The panels are not independent: Hallmark sets share genes and xCell cell types
are spillover-corrected versions of each other. BH controls the false discovery
rate under independence or positive dependence; it is optimistic here. Report
the effect size and the number of tests, not only the adjusted p-value.

## Reporting guidance

Safe to claim:

- "The signature score correlates with EMT programme expression
  (rho = 0.42, BH-adjusted p = 0.009)."
- "Higher stromal score is associated with worse overall survival
  (HR 1.48 per SD, 95% CI 1.13-1.93)."

Not safe to claim without further work:

- Causality of any kind.
- That a checkpoint panel predicts immune-checkpoint-inhibitor response.
- That the immune phenotype calls are a validated classification.
- Clinical utility, or that the signature should guide treatment.
- Cross-cohort comparisons of absolute ssGSEA scores, cell fractions or
  ESTIMATE values.

## Reproducibility

Each run writes:

- `session_info.txt` — R version and every attached package version.
- `panel_definitions.csv` — the exact genes behind every panel, so a panel can
  be reconstructed without the code.
- `cell_state_coverage.csv` — measured genes per panel, including dropped ones.

`msigdbr` pins the MSigDB release; upgrading it changes panel membership. Report
the `msigdbr` version from `session_info.txt` alongside the results.

Scoring is deterministic apart from `xCell`, which uses parallel workers;
`dx_xcell_scores()` defaults to `parallel.sz = 1` so repeated runs agree.
