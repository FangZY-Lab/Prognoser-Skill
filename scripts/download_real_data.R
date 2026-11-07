#!/usr/bin/env Rscript

# Download and format two real colorectal cancer cohorts for Prognoser.
#
# Datasets:
#   GSE39582 - 566 samples, OS and RFS
#   GSE17536 - 177 samples, OS, DFS, and DSS
#
# The script downloads series matrices and clinical annotations, maps probes
# to HGNC symbols, builds Prognoser-compatible objects, and saves .Rdata files
# in the current working directory.

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Install GEOquery first: BiocManager::install('GEOquery')")
}

library(GEOquery)

# GEOquery hard-codes a 120 s curl timeout in its internal downloadFile(),
# which is far too short for the >150 MB GSE39582 series matrix. Replace that
# helper for this session with the same code and a longer timeout.
if (requireNamespace("curl", quietly = TRUE)) {
  longer_download <- function(url, destfile, mode, quiet = TRUE) {
    h <- curl::new_handle()
    curl::handle_setheaders(h, `accept-encoding` = "gzip")
    curl::handle_setopt(h, timeout_ms = 1800 * 1000)
    ok <- tryCatch({
      curl::curl_download(url, destfile, mode = mode, quiet = quiet, handle = h)
      TRUE
    }, error = function(e) { message(e); FALSE })
    if (!ok) {
      if (file.exists(destfile)) file.remove(destfile)
      stop(sprintf("Failed to download %s!", destfile))
    }
    0L
  }
  try(assignInNamespace("downloadFile", longer_download, ns = "GEOquery"),
      silent = TRUE)
}
options(timeout = max(1800, getOption("timeout")))

output_dir <- "prognoser_real_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

load_series <- function(accession) {
  gse <- getGEO(accession, GSEMatrix = TRUE)
  eset <- gse[[1]]

  expr <- exprs(eset)
  pd <- pData(eset)

  # Map probe IDs to gene symbols when the platform annotation is available.
  symbol_col <- grep(
    "gene symbol|symbol|hgnc",
    colnames(fData(eset)),
    ignore.case = TRUE,
    value = TRUE
  )
  if (length(symbol_col) > 0) {
    symbol <- fData(eset)[[symbol_col[1]]]
    keep <- !is.na(symbol) & symbol != ""
    expr <- expr[keep, , drop = FALSE]
    symbol <- symbol[keep]
    expr <- rowsum(expr, group = symbol, reorder = FALSE) /
            as.numeric(table(symbol)[symbol[!duplicated(symbol)]])
  }

  list(expr = expr, pd = pd)
}

surv_time <- function(pd, pattern) {
  col <- grep(pattern, colnames(pd), ignore.case = TRUE, value = TRUE)
  if (length(col) == 0) return(NULL)
  raw <- trimws(sub("^[^:]*: *", "", pd[[col[1]]]))
  num <- suppressWarnings(as.numeric(raw))
  # values such as "45 months" carry a unit; keep the leading number
  bad <- is.na(num) & !is.na(raw) & nzchar(raw)
  if (any(bad)) {
    num[bad] <- suppressWarnings(as.numeric(sub("^([0-9.]+).*$", "\\1", raw[bad])))
  }
  num
}

surv_event <- function(pd, pattern) {
  col <- grep(pattern, colnames(pd), ignore.case = TRUE, value = TRUE)
  if (length(col) == 0) return(NULL)
  raw <- trimws(sub("^[^:]*: *", "", pd[[col[1]]]))
  num <- suppressWarnings(as.numeric(raw))
  # some series code events as text ("Dead"/"Alive", "yes"/"no", ...)
  bad <- is.na(num) & !is.na(raw) & nzchar(raw)
  if (any(bad)) {
    low <- tolower(raw[bad])
    mapped <- rep(NA_integer_, length(low))
    mapped[low %in% c("dead", "deceased", "yes", "y", "event",
                      "recurrence", "recurred", "progressed", "true")] <- 1L
    mapped[low %in% c("alive", "living", "no", "n", "censored",
                      "none", "disease free", "false")] <- 0L
    num[bad] <- mapped
  }
  as.integer(num != 0)
}

build_survival <- function(pd, time_pattern, event_pattern) {
  time <- surv_time(pd, time_pattern)
  event <- surv_event(pd, event_pattern)
  if (is.null(time) || is.null(event)) return(NULL)
  surv <- data.frame(time = time, status = event, row.names = rownames(pd))
  surv <- surv[!is.na(surv$time) & !is.na(surv$status), ]
  surv$time <- as.numeric(surv$time)
  surv$status <- as.integer(surv$status != 0)
  surv
}

datasets <- list(
  GSE39582 = list(
    endpoints = list(
      OS = list("os.time|os_time|overall", "os.event|os_event|death"),
      RFS = list("rfs.time|recurrence", "rfs.event|recurrence")
    )
  ),
  GSE17536 = list(
    endpoints = list(
      OS = list("os.time|os_time|overall", "os.event|os_event|death"),
      DFS = list("dfs.time|disease free", "dfs.event|disease free"),
      DSS = list("dss.time|disease specific", "dss.event|disease specific")
    )
  )
)

for (acc in names(datasets)) {
  data <- load_series(acc)
  for (ep in names(datasets[[acc]]$endpoints)) {
    patterns <- datasets[[acc]]$endpoints[[ep]]
    surv <- build_survival(data$pd, patterns[[1]], patterns[[2]])
    if (is.null(surv)) {
      message("Skipping ", acc, "_", ep, ": survival columns not detected.")
      next
    }
    common <- intersect(colnames(data$expr), rownames(surv))
    expr <- data$expr[, common, drop = FALSE]
    surv <- surv[common, , drop = FALSE]

    expr_name <- paste0(acc, "_", ep)
    surv_name <- paste0(acc, "_", ep, "_survival")
    assign(expr_name, expr, envir = .GlobalEnv)
    assign(surv_name, surv, envir = .GlobalEnv)

    save(list = expr_name, file = file.path(output_dir, paste0(expr_name, ".Rdata")))
    save(list = surv_name, file = file.path(output_dir, paste0(surv_name, ".Rdata")))
    message("Saved ", expr_name, ": ", nrow(expr), " genes x ", ncol(expr), " samples.")
  }
}

message("Done. Load the .Rdata files into the global environment, then run Prognoser.")
