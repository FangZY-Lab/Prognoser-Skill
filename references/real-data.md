# Real data support

The Prognoser skill is validated against two public colorectal cancer
transcriptome cohorts. Both provide expression plus multiple survival
endpoints, matching the multi-cohort, multi-endpoint workflow the package
requires.

## Validation datasets

| Dataset | GEO accession | Samples | Genes | Endpoints |
| --- | --- | --- | --- | --- |
| Marisa cohort | GSE39582 | 566 | 23,348 | OS, RFS |
| Smith cohort | GSE17536 | 177 | 23,348 | OS, DFS, DSS |

Both matrices are Affymetrix HG-U133 Plus 2.0 data collapsed to HGNC symbols
and log2-scaled. The gene dimension is identical (23,348), which allows direct
gene intersection across cohorts.

## Download with GEOquery

Install the downloader once:

```r
BiocManager::install("GEOquery")
```

Download expression series:

```r
library(GEOquery)

gse39582 <- getGEO("GSE39582", GSEMatrix = TRUE)
gse17536 <- getGEO("GSE17536", GSEMatrix = TRUE)
```

Download supplementary clinical files when survival is not embedded in the
series matrix:

```r
getGEOSuppFiles("GSE39582")
getGEOSuppFiles("GSE17536")
```

## Endpoint extraction

GEO stores sample-level annotations in `pData(eset)` under
`characteristics_ch1`. Column names differ by series, so search with regular
expressions rather than hard-coding one convention:

```r
pd <- pData(eset)

surv_time <- function(pd, pattern) {
  col <- grep(pattern, colnames(pd), ignore.case = TRUE, value = TRUE)
  as.numeric(sub("^[^:]*: *", "", pd[[col[1]]]))
}

surv_event <- function(pd, pattern) {
  col <- grep(pattern, colnames(pd), ignore.case = TRUE, value = TRUE)
  as.integer(sub("^[^:]*: *", "", pd[[col[1]]]) != 0)
}
```

Typical field names:

| Endpoint | Time patterns | Event patterns |
| --- | --- | --- |
| OS | `os.time`, `os_time`, `overall survival`, `months to death` | `os.event`, `os_event`, `vital status`, `death` |
| RFS | `rfs.time`, `recurrence free survival` | `rfs.event`, `recurrence` |
| DFS | `dfs.time`, `disease free survival` | `dfs.event` |
| DSS | `dss.time`, `disease specific survival` | `dss.event` |

## Survival coding

Normalize all event columns to `0/1`:

```r
surv <- data.frame(
  time = surv_time(pd, "os.time"),
  status = surv_event(pd, "os.event"),
  row.names = rownames(pd)
)
surv <- surv[!is.na(surv$time) & !is.na(surv$status), ]
surv$status <- as.integer(surv$status != 0)
```

## Build Prognoser objects

Expose each dataset-endpoint pair with the required names:

```r
assign("GSE39582_OS", expr39582, envir = .GlobalEnv)
assign("GSE39582_OS_survival", surv39582_os, envir = .GlobalEnv)
assign("GSE17536_OS", expr17536, envir = .GlobalEnv)
assign("GSE17536_OS_survival", surv17536_os, envir = .GlobalEnv)
```

The bundled [download_real_data.R](../scripts/download_real_data.R) performs
these steps and saves `.Rdata` files ready for `Prognoser()`.

## Reproducibility notes

- Record the GEO accession, platform, sample count, and survival column names
  used for every run.
- Do not drop genes before downloading; run symbol mapping and duplicate
  collapse only after inspecting the source annotation.
- Keep endpoint availability explicit: if a cohort lacks RFS, omit that
  dataset-endpoint rather than substituting DFS.
