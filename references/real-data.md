# Real data support

The Prognoser skill is validated against two public colorectal cancer
transcriptome cohorts. Both provide expression plus multiple survival
endpoints, matching the multi-cohort, multi-endpoint workflow the package
requires.

## Validation datasets

| Dataset | GEO accession | Samples | Genes | Endpoints |
| --- | --- | --- | --- | --- |
| Marisa cohort | GSE39582 | 585 (579 with OS, 577 with RFS) | 23,520 | OS, RFS |
| Smith cohort | GSE17536 | 177 | 23,520 | OS, DFS, DSS |

Both matrices are Affymetrix HG-U133 Plus 2.0 data. Collapsing probes to HGNC
symbols by mean gives 23,520 genes in both cohorts, an identical gene
dimension that allows direct intersection.

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
`characteristics_ch1.N`, and the field *name lives in the values*
(`"os.delay (months): 42"`). Grepping the column names therefore finds
nothing, so the search has to look at both the column name and the label:

```r
pd <- pData(eset)

field_label <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return("")
  first <- trimws(x[1])
  if (grepl("^[^:]+:", first)) tolower(trimws(sub(":.*$", "", first))) else ""
}

field_index <- function(pd, pattern, kind) {
  labels <- vapply(colnames(pd), function(cn) {
    tolower(paste(cn, field_label(pd[[cn]])))
  }, character(1))
  core <- sub("\\(.*$", "", labels)   # drop the parenthetical gloss
  hit <- grepl(pattern, labels, ignore.case = TRUE)
  is_event <- grepl("event|status|death|dead|recurrence|progress", core)
  is_time <- grepl("time|delay|follow|month|day|year|survival", core)
  idx <- if (kind == "time") which(hit & is_time & !is_event) else which(hit & is_event & !is_time)
  if (!length(idx)) NA_integer_ else idx[1]
}
```

Typical field names:

| Endpoint | Time patterns | Event patterns |
| --- | --- | --- |
| OS | `os.time`, `os_time`, `os.delay`, `overall survival` | `os.event`, `os_event`, `death` |
| RFS | `rfs.time`, `rfs.delay`, `recurrence free survival` | `rfs.event`, `recurrence` |
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
