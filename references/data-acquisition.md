# Public data acquisition

`Prognoser` needs paired expression and survival data. This reference covers
reliable ways to obtain both from public repositories.

## GEO

### Expression series

Use `GEOquery` to download a processed series matrix:

```r
library(GEOquery)

gse <- getGEO("GSE39582", GSEMatrix = TRUE)
eset <- gse[[1]]
expr <- exprs(eset)
```

Map probes to gene symbols:

- Use the platform annotation (`fData(eset)` often has `Gene Symbol` or
  `GENE_SYMBOL`).
- For Affymetrix, `affy` or `hgu133plus2.db` provides a probe-to-symbol map.
- Collapse probes per symbol by mean, median, or maximum expression.

### Survival data

Survival is usually in:

- `pData(eset)` columns such as `os.time`, `os.event`, `rfs.time`,
  `rfs.event`, `dfs.time`, `dfs.event`, `pfs.time`, `pfs.event`.
- Or in a supplementary file linked on the GEO page.

Normalize the event coding to `0/1`:

```r
surv <- data.frame(
  time = as.numeric(pData(eset)$os.time),
  status = as.numeric(pData(eset)$os.event),
  row.names = rownames(pData(eset))
)
surv <- surv[!is.na(surv$time) & !is.na(surv$status), ]
surv$status <- as.integer(surv$status != 0)
```

## TCGA

### GDC via TCGAbiolinks

Download expression and clinical data:

```r
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
expr <- GDCprepare(query)
```

For survival, query clinical data:

```r
clinical_query <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
```

Relevant fields include `days_to_death`, `days_to_last_follow_up`,
`vital_status`, and recurrence/progression timing fields.

### Xena via UCSCXena

```r
library(UCSCXena)

expr <- XenaGenerate(subset = XenaDatasets == "TCGA.COAD.sampleMap/HiSeqV2")
expr <- XenaQuery(expr) |> XenaDownload()
expr_mat <- XenaPrepare(expr)
```

TCGA Pan-Cancer clinical data on Xena contains curated OS, PFI, DFI, and DSS
times and event indicators.

## Multi-endpoint assembly

Endpoints are represented as separate entries that share a dataset prefix:

```text
OS  = overall survival
RFS = recurrence-free survival
PFS = progression-free survival
DFS = disease-free survival
```

Build each endpoint independently:

1. Select the relevant time and event columns.
2. Coerce both to numeric.
3. Keep samples with complete values.
4. Set `status` to `0/1`.
5. Match row names to expression sample names exactly.

When an endpoint does not exist for a cohort, leave it out. Do not duplicate
another endpoint to fill the gap.

## Sample alignment

After preparing expression and survival, align samples:

```r
common <- intersect(colnames(expr), rownames(surv))
expr <- expr[, common, drop = FALSE]
surv <- surv[common, , drop = FALSE]
```

Stop and inspect the source when `length(common)` is unexpectedly small; it
usually means identifiers use different sample-name conventions.
