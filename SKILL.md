---
name: prognoser
description: End-to-end discovery of prognostic gene signatures with the Prognoser R package, including public data download, multi-endpoint survival assembly, expression normalization, and R execution.
---

# Prognoser

Use this skill when the task is to derive high-risk (`HRGS`) and low-risk
(`LRGS`) prognostic gene signatures from multi-cohort, multi-endpoint
survival data with the `Prognoser` R package.

The workflow is:

1. Acquire public expression and survival data.
2. Build one expression matrix and one survival table per dataset-endpoint.
3. Normalize and align gene identifiers across cohorts.
4. Expose data in the exact naming contract `Prognoser` expects.
5. Run `Prognoser()` and export the resulting signatures.

## Critical data contract

`Prognoser()` reads objects from the R global environment with `get()`. This
contract is non-negotiable and is the most common source of failures.

For every entry in `expression_accession_vector`, for example `TCGA_COAD_OS`,
the global environment must contain:

- An expression `data.frame` named exactly `TCGA_COAD_OS`:
  - genes as rows
  - samples as columns
- A survival `data.frame` named exactly `TCGA_COAD_OS_survival`:
  - numeric column `time`
  - numeric column `status` (`0` = alive, `1` = event)
  - row names equal to the expression matrix column names

Entries sharing the same prefix before `_` are grouped and meta-analyzed as
one dataset:

```text
TCGA_COAD_OS
TCGA_COAD_RFS
TCGA_COAD_PFS
TCGA_COAD_DFS
```

The prefix is `TCGA_COAD`; the four endpoints are combined before downstream
selection.

## Workflow

### 1. Acquire data

Prefer public, documented sources:

- GEO: `GEOquery::getGEO` for expression series; supplementary files for
  clinical survival.
- TCGA: `TCGAbiolinks` (GDC) or `UCSCXena` (Xena browser) for expression and
  clinical endpoints.

Read [references/data-acquisition.md](references/data-acquisition.md) before
downloading or assembling survival endpoints.

### 2. Normalize and align

- Genes in rows, samples in columns.
- Map probe IDs or Ensembl IDs to HGNC symbols.
- Collapse duplicate symbols, remove NA genes.
- Log2-transform counts/TPM/RSEM when values are on a linear scale.
- Do not run gene filtering on survival time/status.

Read [references/preprocessing.md](references/preprocessing.md) for the full
procedure and edge cases.

### 3. Format for Prognoser

Write every prepared dataset-endpoint pair into the global environment with
the names above. Use `assign()` or load `.Rdata` files whose objects already
have those names.

### 4. Run Prognoser

Install once:

```r
install.packages(c("survival", "RobustRankAggreg", "meta", "remotes"))
remotes::install_github("FangZY-Lab/Prognoser")
```

Then run, passing a single `linkage` value:

```r
library(Prognoser)

result <- Prognoser(
  expression_accession_vector = datasets,
  knowledge_genesets = knowledge,
  denovo_threshold = 0.05,
  gene_threshold = 0.001,
  min_genes = 5,
  similarity_threshold = 0.5,
  linkage = "ward.D"
)
```

Full API details and the output schema are in
[references/prognoser-api.md](references/prognoser-api.md). Use
[scripts/run_prognoser.R](scripts/run_prognoser.R) as a reproducible template.

### 5. Export results

`result` contains:

- `result$genesets`: final non-redundant `HRGS*` and `LRGS*` signatures.
- `result$combined_matrix_protection`: signed effect matrix for protective genes.
- `result$combined_matrix_risk`: signed effect matrix for risk genes.

Save with `saveRDS()` and write gene sets as tab-delimited GMT-compatible text
when a downstream tool needs them.

## Guardrails

- Do not fabricate survival time or status. If an endpoint is unavailable,
  omit that dataset-endpoint instead of imputing events.
- Verify sample overlap between expression columns and survival row names;
  report and drop unmatched samples.
- Keep `linkage` a single character value; the default is a vector and will
  fail if left unset.
- `meta` may emit non-fatal warnings during meta-analysis; preserve them rather
  than suppressing them unless they flood the log.
