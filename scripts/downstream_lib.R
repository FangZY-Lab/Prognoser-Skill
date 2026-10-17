#!/usr/bin/env Rscript

# Shared helpers for the Prognoser downstream modules.
#
# Sourced by downstream_cell_states.R, downstream_tme.R, run_downstream.R and
# example_downstream.R. This file only defines functions; it never touches the
# global environment and never runs analysis on load.
#
# Design notes
#   * Every scoring function takes an expression matrix with genes in rows and
#     samples in columns, and returns a matrix or data frame with samples in
#     columns / rows respectively. No side effects.
#   * Optional dependencies are checked with dx_require() so the caller gets an
#     actionable message instead of a raw "there is no package called" error.
#   * All user-facing progress goes through dx_log() (stderr), so scripts stay
#     safe to call with stdout redirected to a report file.

## --------------------------------------------------------------- paths ---

dx_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) {
    return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  }
  normalizePath(sub("^--file=", "", hit[[1L]]), winslash = "/", mustWork = FALSE)
}

dx_script_dir <- function() dirname(dx_script_path())

# Locate a sibling script regardless of how the entry point was invoked:
# from the scripts directory, from the repository root, or from elsewhere with
# an absolute --file path. PROGNOSER_SCRIPTS_DIR overrides the search.
dx_find_script <- function(file) {
  cand <- c(
    Sys.getenv("PROGNOSER_SCRIPTS_DIR", ""),
    dx_script_dir(),
    file.path(getwd(), "scripts"),
    getwd()
  )
  cand <- cand[nzchar(cand)]
  for (d in cand) {
    p <- file.path(d, file)
    if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  stop(sprintf(paste0("Cannot locate '%s'. Set PROGNOSER_SCRIPTS_DIR to the ",
                      "directory that holds the Prognoser scripts."), file),
       call. = FALSE)
}

dx_source_sibling <- function(file) {
  source(dx_find_script(file), local = FALSE)
}

# TRUE only when the calling script was launched with Rscript, so modules can be
# both sourced (tests, notebooks, run_downstream.R) and executed directly.
dx_is_main <- function(script_file) {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0L) return(FALSE)
  mine <- normalizePath(sub("^--file=", "", hit[[1L]]), winslash = "/", mustWork = FALSE)
  theirs <- normalizePath(file.path(dx_script_dir(), script_file),
                          winslash = "/", mustWork = FALSE)
  identical(mine, theirs)
}

## ------------------------------------------------------------- logging ---

dx_log <- function(fmt, ...) {
  message(format(Sys.time(), "[%H:%M:%S] "), sprintf(fmt, ...))
}

dx_require <- function(pkgs, hint = NULL) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) == 0L) return(invisible(TRUE))
  msg <- sprintf("Required package(s) not installed: %s", paste(missing, collapse = ", "))
  if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
  stop(msg, call. = FALSE)
}

# Some packages (estimate, xCell) ship lazy-loaded data sets that their own
# functions resolve from the search path, so `pkg::fun()` alone fails with
# "object 'common_genes' not found". Attaching the package is the documented
# way to use them; the side effect is contained and reported.
dx_attach <- function(pkg) {
  if (paste0("package:", pkg) %in% search()) return(invisible(TRUE))
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  )
  invisible(TRUE)
}

## ----------------------------------------------------------- arguments ---

# Minimal `--key=value`, `--key value` and `--flag` parser. Dashes in flag names
# are normalised to underscores so `--write-inputs` lands in `args$write_inputs`.
# Values are returned as character; --flag style switches come back as TRUE.
dx_parse_args <- function(defaults = list(), argv = commandArgs(trailingOnly = TRUE)) {
  out <- defaults
  norm_key <- function(x) gsub("-", "_", x, fixed = TRUE)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) {
      stop(sprintf("Unexpected argument '%s'. Use --key=value.", a), call. = FALSE)
    }
    a <- sub("^--", "", a)
    if (grepl("=", a, fixed = TRUE)) {
      out[[norm_key(sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
      i <- i + 1L
    } else if (i < length(argv) && !startsWith(argv[[i + 1L]], "--")) {
      out[[norm_key(a)]] <- argv[[i + 1L]]
      i <- i + 2L
    } else {
      out[[norm_key(a)]] <- TRUE
      i <- i + 1L
    }
  }
  out
}

dx_split <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  if (length(x) > 1L) return(x)
  if (is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

dx_as_flag <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y", "on")
}

## ------------------------------------------------------ expr utilities ---

dx_as_matrix <- function(x, name = "expression") {
  m <- as.matrix(x)
  if (!is.numeric(m)) stop(sprintf("%s must be numeric.", name), call. = FALSE)
  if (is.null(rownames(m)) || is.null(colnames(m))) {
    stop(sprintf("%s needs both rownames (genes) and colnames (samples).", name),
         call. = FALSE)
  }
  if (anyNA(m)) {
    stop(sprintf("%s contains %d NA value(s); impute or drop them first.",
                 name, sum(is.na(m))), call. = FALSE)
  }
  storage.mode(m) <- "double"
  m
}

# Heuristic only: log2-scale matrices have small ranges and are usually < 30.
# Callers can and should override this when they know their input scale.
dx_detect_scale <- function(m) {
  rng <- range(m, na.rm = TRUE)
  if (rng[2L] > 50) "linear" else "log2"
}

dx_to_linear <- function(m) {
  # Exact inverse of log2(x + 1); an approximation for other log bases, which is
  # why the docs require log2(TPM + 1) for TPM-based deconvolution.
  pmax(2^m - 1, 0)
}

dx_zscore <- function(m) {
  z <- t(scale(t(m)))
  z[!is.finite(z)] <- 0
  attr(z, "scaled:center") <- NULL
  attr(z, "scaled:scale") <- NULL
  z
}

# Keep genes present in every matrix, preserving the order of the first one.
dx_common_genes <- function(mats) {
  if (length(mats) == 1L) return(rownames(mats[[1L]]))
  Reduce(intersect, lapply(mats, rownames))
}

dx_align_to <- function(mat, genes) {
  keep <- intersect(genes, rownames(mat))
  if (length(keep) == 0L) {
    stop("No shared genes between the expression matrix and the gene sets.",
         call. = FALSE)
  }
  mat[keep, , drop = FALSE]
}

## ----------------------------------------------------- gene set scores ---

# ssGSEA through GSVA >= 1.50 (parameter-object API). Falls back to the legacy
# positional call on older GSVA, so the skill keeps working on Bioconductor 3.16.
dx_ssgsea <- function(expr, gene_sets, min_size = 5L, max_size = 500L) {
  dx_require("GSVA", "BiocManager::install('GSVA')")
  sizes <- vapply(gene_sets, length, integer(1L))
  keep <- sizes >= min_size & sizes <= max_size
  if (!all(keep)) {
    dx_log("Dropped %d gene set(s) outside [%d, %d] genes: %s",
           sum(!keep), min_size, max_size,
           paste(names(gene_sets)[!keep], collapse = ", "))
    gene_sets <- gene_sets[keep]
  }
  if (length(gene_sets) == 0L) {
    stop("No gene set survived the size filter.", call. = FALSE)
  }
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::ssgseaParam(expr, gene_sets, minSize = min_size,
                               maxSize = max_size, normalize = TRUE)
    # GSVA prints progress with print(), not message(), so capture.output() is
    # what actually keeps script logs readable.
    invisible(utils::capture.output(
      scores <- suppressMessages(GSVA::gsva(param, verbose = FALSE))
    ))
  } else {
    invisible(utils::capture.output(
      scores <- suppressMessages(
        GSVA::gsva(expr, gene_sets, method = "ssgsea",
                   kcdf = "Gaussian", ssgsea.norm = TRUE, verbose = FALSE)
      )
    ))
  }
  as.matrix(scores)
}

# Dependency-free fallback: mean z-score over the genes of each set. Cheaper and
# more robust than ssGSEA on small gene sets, but it ignores the rank
# distribution, so relative panel magnitudes are not comparable.
dx_zscore_sets <- function(expr, gene_sets) {
  z <- dx_zscore(expr)
  out <- vapply(gene_sets, function(gs) {
    hit <- intersect(gs, rownames(z))
    if (length(hit) == 0L) return(rep(NA_real_, ncol(z)))
    colMeans(z[hit, , drop = FALSE])
  }, numeric(ncol(z)))
  out <- t(out)
  rownames(out) <- names(gene_sets)
  colnames(out) <- colnames(expr)
  out
}

dx_score_sets <- function(expr, gene_sets, method = c("ssgsea", "zscore"),
                          min_size = 5L, max_size = 500L) {
  method <- match.arg(method)
  # Filter on *measured* genes, not on the nominal set size. GSVA silently drops
  # sets whose intersection with the matrix is below minSize, which used to make
  # panels vanish from the returned matrix without a trace.
  measured <- lapply(gene_sets, function(gs) intersect(gs, rownames(expr)))
  n_genes <- vapply(gene_sets, length, integer(1L))
  n_measured <- vapply(measured, length, integer(1L))
  ok <- n_measured >= min_size & n_measured <= max_size
  coverage <- data.frame(
    panel = names(gene_sets),
    n_genes = unname(n_genes),
    n_measured = unname(n_measured),
    coverage = round(unname(n_measured) / unname(n_genes), 3),
    scored = unname(ok),
    row.names = NULL
  )
  if (any(!ok)) {
    dx_log("WARNING: %d panel(s) not scored (need %d-%d measured genes): %s",
           sum(!ok), min_size, max_size,
           paste(sprintf("%s [%d]", coverage$panel[!ok], coverage$n_measured[!ok]),
                 collapse = ", "))
  }
  if (!any(ok)) {
    stop(sprintf("No panel has at least %d measured genes.", min_size),
         call. = FALSE)
  }
  score_sets <- measured[ok]
  scores <- if (identical(method, "ssgsea")) {
    dx_ssgsea(expr, score_sets, min_size = 1L, max_size = Inf)
  } else {
    dx_zscore_sets(expr, score_sets)
  }
  list(scores = scores, coverage = coverage)
}

# Safe row access for score matrices: returns NULL instead of erroring when a
# panel was dropped for insufficient coverage.
dx_row <- function(mat, name) {
  if (is.null(mat) || !name %in% rownames(mat)) return(NULL)
  mat[name, ]
}

# Signed signature score: mean z of the risk genes minus mean z of the
# protective genes. Mirrors how HRGS/LRGS are meant to be read (higher = worse).
dx_signature_score <- function(expr, up, down = character(0)) {
  z <- dx_zscore(expr)
  up_hit <- intersect(up, rownames(z))
  down_hit <- intersect(down, rownames(z))
  if (length(up_hit) == 0L && length(down_hit) == 0L) {
    stop("Neither the risk nor the protective genes are measured.", call. = FALSE)
  }
  if (length(up_hit) == 0L) dx_log("WARNING: no measured risk gene.")
  if (length(down) > 0L && length(down_hit) == 0L) dx_log("WARNING: no measured protective gene.")
  score <- rep(0, ncol(z))
  if (length(up_hit)) score <- score + colMeans(z[up_hit, , drop = FALSE])
  if (length(down_hit)) score <- score - colMeans(z[down_hit, , drop = FALSE])
  names(score) <- colnames(z)
  attr(score, "n_risk") <- length(up_hit)
  attr(score, "n_protective") <- length(down_hit)
  score
}

# Split Prognoser gene sets into the risk (HRGS*) and protective (LRGS*) pools.
dx_genesets_to_pools <- function(genesets) {
  if (is.null(names(genesets)) || any(!nzchar(names(genesets)))) {
    stop("genesets must be a named list, e.g. result$genesets.", call. = FALSE)
  }
  risk <- unique(unlist(genesets[grepl("^HRGS", names(genesets))], use.names = FALSE))
  prot <- unique(unlist(genesets[grepl("^LRGS", names(genesets))], use.names = FALSE))
  if (length(risk) == 0L && length(prot) == 0L) {
    stop("No gene set name starts with HRGS or LRGS; check result$genesets.",
         call. = FALSE)
  }
  list(risk = risk, protective = prot)
}

## --------------------------------------------------------------- input ---

# Load <Dataset>_<Endpoint> and <Dataset>_<Endpoint>_survival pairs from a
# directory of .rds / .Rdata files without polluting the global environment.
dx_load_inputs <- function(input_dir, datasets = NULL) {
  if (!dir.exists(input_dir)) {
    stop(sprintf("Input directory not found: %s", input_dir), call. = FALSE)
  }
  files <- list.files(input_dir, pattern = "[.](rds|Rdata|rda)$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0L) {
    stop(sprintf("No .rds/.Rdata files in %s", input_dir), call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  for (f in files) {
    if (grepl("[.]rds$", f, ignore.case = TRUE)) {
      assign(tools::file_path_sans_ext(basename(f)), readRDS(f), envir = env)
    } else {
      load(f, envir = env)
    }
  }
  available <- ls(env)
  if (is.null(datasets) || length(datasets) == 0L) {
    datasets <- sub("_survival$", "", grep("_survival$", available, value = TRUE))
  }
  if (length(datasets) == 0L) {
    stop("No dataset-endpoint with a matching _survival object was found.",
         call. = FALSE)
  }
  out <- vector("list", length(datasets))
  names(out) <- datasets
  for (d in datasets) {
    if (!exists(d, envir = env, inherits = FALSE)) {
      stop(sprintf("Expression object '%s' is missing.", d), call. = FALSE)
    }
    surv_name <- paste0(d, "_survival")
    if (!exists(surv_name, envir = env, inherits = FALSE)) {
      stop(sprintf("Survival object '%s' is missing.", surv_name), call. = FALSE)
    }
    expr <- dx_as_matrix(get(d, envir = env), sprintf("'%s'", d))
    surv <- get(surv_name, envir = env)
    out[[d]] <- dx_prepare_pair(expr, surv, d)
  }
  out
}

# Coerce survival to numeric time/status and intersect samples with expression.
dx_prepare_pair <- function(expr, surv, label = "dataset") {
  surv <- as.data.frame(surv, stringsAsFactors = FALSE)
  need <- c("time", "status")
  if (!all(need %in% names(surv))) {
    stop(sprintf("'%s_survival' needs columns time and status (found: %s).",
                 label, paste(names(surv), collapse = ", ")), call. = FALSE)
  }
  time <- suppressWarnings(as.numeric(surv$time))
  status <- suppressWarnings(as.numeric(surv$status))
  if (anyNA(time) || anyNA(status)) {
    stop(sprintf("'%s_survival' has non-numeric time/status values.", label),
         call. = FALSE)
  }
  rownames(surv) <- rownames(surv) %||% colnames(expr)
  common <- intersect(colnames(expr), rownames(surv))
  if (length(common) == 0L) {
    stop(sprintf("No shared sample IDs between '%s' and '%s_survival'.", label, label),
         call. = FALSE)
  }
  if (length(common) < ncol(expr) || length(common) < nrow(surv)) {
    dx_log("'%s': kept %d samples present in both expression and survival (dropped %d/%d).",
           label, length(common), ncol(expr) - length(common), ncol(expr))
  }
  idx <- match(common, rownames(surv))
  list(
    expr = expr[, common, drop = FALSE],
    time = time[idx],
    status = status[idx],
    samples = common
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

## ----------------------------------------------------------- statistics ---

dx_bh <- function(p) stats::p.adjust(p, method = "BH")

# Spearman correlation of one vector against every row of a matrix.
dx_spearman <- function(x, mat) {
  if (is.null(rownames(mat))) stop("mat needs rownames.", call. = FALSE)
  rows <- lapply(rownames(mat), function(f) {
    y <- as.numeric(mat[f, ])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 5L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
      return(data.frame(feature = f, n = sum(ok), rho = NA_real_, p = NA_real_))
    }
    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman",
                                           exact = FALSE))
    data.frame(feature = f, n = sum(ok), rho = unname(ct$estimate),
               p = ct$p.value)
  })
  out <- do.call(rbind, rows)
  out$p_adj <- dx_bh(out$p)
  out[order(out$p), , drop = FALSE]
}

# Hazard ratio per one standard deviation of `feature`.
dx_cox_per_sd <- function(feature, time, status) {
  ok <- is.finite(feature) & is.finite(time) & is.finite(status) & time > 0
  n <- sum(ok)
  events <- sum(status[ok] == 1)
  empty <- data.frame(n = n, events = events, HR = NA_real_, lower = NA_real_,
                      upper = NA_real_, p = NA_real_)
  if (n < 10L || events < 5L) return(empty)
  z <- as.numeric(scale(feature[ok]))
  if (!is.finite(stats::sd(z)) || stats::sd(z) == 0) return(empty)
  d <- data.frame(time = time[ok], status = status[ok], z = z)
  fit <- tryCatch(survival::coxph(survival::Surv(time, status) ~ z, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(empty)
  s <- summary(fit)
  data.frame(n = n, events = events,
             HR = unname(s$conf.int[1L, 1L]),
             lower = unname(s$conf.int[1L, 3L]),
             upper = unname(s$conf.int[1L, 4L]),
             p = unname(s$coefficients[1L, 5L]))
}

dx_logrank_median <- function(feature, time, status) {
  ok <- is.finite(feature) & is.finite(time) & is.finite(status) & time > 0
  if (sum(ok) < 10L) return(list(p = NA_real_, n = sum(ok)))
  grp <- feature[ok] > stats::median(feature[ok], na.rm = TRUE)
  if (length(unique(grp)) < 2L) return(list(p = NA_real_, n = sum(ok)))
  d <- data.frame(time = time[ok], status = status[ok],
                  grp = factor(grp, levels = c(FALSE, TRUE)))
  fit <- tryCatch(survival::survdiff(survival::Surv(time, status) ~ grp, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, n = sum(ok)))
  list(p = stats::pchisq(fit$chisq, df = 1L, lower.tail = FALSE), n = sum(ok))
}

## ------------------------------------------------- deconvolution tools ---

# ESTIMATE requires a GCT file, so the matrix is staged through disk. The
# filtered file is what estimateScore() actually expects; skipping
# filterCommonGenes() silently produces shifted, all-NA output.
dx_estimate_scores <- function(expr_log2, platform = "illumina",
                               id = "GeneSymbol") {
  dx_require("estimate", "install.packages('estimate') or remotes::install_github('dviraran/estimate')")
  dx_attach("estimate")
  platform <- match.arg(platform, c("affymetrix", "agilent", "illumina"))
  in_f <- tempfile(fileext = ".txt")
  fil_f <- tempfile(fileext = ".txt")
  out_f <- tempfile(fileext = ".txt")
  on.exit(unlink(c(in_f, fil_f, out_f)), add = TRUE)

  df <- data.frame(GeneSymbol = rownames(expr_log2), expr_log2,
                   check.names = FALSE)
  utils::write.table(df, in_f, sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(utils::capture.output(
    estimate::filterCommonGenes(in_f, fil_f, id = id)
  ))
  n_common <- as.integer(strsplit(readLines(fil_f, n = 2L)[2L], "\t")[[1L]][[2L]])
  invisible(utils::capture.output(
    estimate::estimateScore(fil_f, out_f, platform = platform)
  ))

  scores <- utils::read.table(out_f, sep = "\t", skip = 2L, header = TRUE,
                              check.names = FALSE, stringsAsFactors = FALSE)
  rownames(scores) <- scores$NAME
  scores$NAME <- NULL
  scores$Description <- NULL
  scores <- as.matrix(scores)
  attr(scores, "n_common_genes") <- n_common
  attr(scores, "platform") <- platform
  scores
}

dx_xcell_scores <- function(expr, rnaseq = TRUE, min_genes = 5000L,
                            parallel_sz = 1L) {
  dx_require("xCell", "BiocManager::install('xCell')")
  if (nrow(expr) < min_genes) {
    stop(sprintf(paste0("xCell needs a reasonably complete transcriptome; ",
                        "%d genes is too few (min_genes = %d)."),
                 nrow(expr), min_genes), call. = FALSE)
  }
  dx_attach("xCell")
  old <- options(warn = 1)
  on.exit(options(old), add = TRUE)
  invisible(utils::capture.output(
    res <- suppressMessages(
      xCell::xCellAnalysis(expr, rnaseq = rnaseq, parallel.sz = parallel_sz,
                           save.raw = FALSE)
    )
  ))
  as.matrix(res)
}

# quanTIseq expects linear TPM values (not log2) and its own gene symbols.
dx_quantiseq <- function(expr_linear, is_arraydata = FALSE, is_tumordata = FALSE) {
  dx_require("quantiseqr", "BiocManager::install('quantiseqr')")
  res <- quantiseqr::run_quantiseq(
    expression_data = expr_linear,
    is_arraydata = is_arraydata,
    is_tumordata = is_tumordata,
    scale_mRNA = TRUE,
    method = "lsei"
  )
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  rownames(res) <- res$Sample
  res$Sample <- NULL
  as.matrix(res)
}

## --------------------------------------------------------------- output ---

dx_output_dir <- function(path) {
  if (is.null(path) || !nzchar(path)) path <- file.path(getwd(), "output", "downstream")
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

dx_write_csv <- function(df, dir, file) {
  path <- file.path(dir, file)
  utils::write.csv(df, path, row.names = FALSE)
  dx_log("wrote %s (%d rows)", path, nrow(df))
  invisible(path)
}

dx_write_matrix <- function(mat, dir, file, row_name = "feature") {
  df <- data.frame(feature = rownames(mat), mat, check.names = FALSE)
  names(df)[1L] <- row_name
  dx_write_csv(df, dir, file)
}

# pheatmap is optional; a missing plotting stack must never fail an analysis run.
dx_plot_heatmap <- function(mat, dir, file, annotation = NULL, title = NULL) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    dx_log("pheatmap not installed; skipping %s", file)
    return(invisible(FALSE))
  }
  ok <- tryCatch({
    grDevices::pdf(file.path(dir, file), width = 8, height = 7)
    on.exit(grDevices::dev.off(), add = TRUE)
    pheatmap::pheatmap(mat, scale = "none", annotation_col = annotation,
                       show_rownames = nrow(mat) <= 40, fontsize_row = 8,
                       main = title, silent = FALSE)
    TRUE
  }, error = function(e) {
    if (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
    dx_log("plot %s failed: %s", file, conditionMessage(e))
    FALSE
  })
  invisible(ok)
}

dx_plot_forest <- function(df, dir, file, title = "Hazard ratio per SD") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    dx_log("ggplot2 not installed; skipping %s", file)
    return(invisible(FALSE))
  }
  df <- df[is.finite(df$HR), , drop = FALSE]
  if (nrow(df) == 0L) return(invisible(FALSE))
  ok <- tryCatch({
    df$label <- factor(df$label, levels = rev(df$label))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = HR, y = label)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
                             width = 0.2, orientation = "y") +
      ggplot2::geom_point(ggplot2::aes(colour = p < 0.05), size = 2.4) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "grey40"),
                                   guide = "none") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "HR per SD (log scale)", y = NULL, title = title) +
      ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(file.path(dir, file), p, width = 7,
                    height = 0.35 * nrow(df) + 2, units = "in")
    TRUE
  }, error = function(e) {
    dx_log("plot %s failed: %s", file, conditionMessage(e))
    FALSE
  })
  invisible(ok)
}

dx_session_info <- function(dir) {
  path <- file.path(dir, "session_info.txt")
  utils::capture.output(utils::sessionInfo(), file = path)
  invisible(path)
}

## ------------------------------------------------------------- report ---

dx_top <- function(df, by = "p_adj", n = 10L, cols = NULL, group = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  take <- function(d) {
    if (!by %in% names(d) || !is.atomic(d[[by]])) return(utils::head(d, n))
    keep <- is.finite(d[[by]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) == 0L) return(d)
    utils::head(d[order(d[[by]]), , drop = FALSE], n)
  }
  if (!is.null(group) && group %in% names(df)) {
    parts <- lapply(split(df, df[[group]]), take)
    df <- do.call(rbind, parts)
    rownames(df) <- NULL
  } else {
    df <- take(df)
  }
  if (!is.null(cols)) df <- df[, intersect(cols, names(df)), drop = FALSE]
  rownames(df) <- NULL
  df
}

dx_md_table <- function(df, digits = 3L) {
  if (is.null(df) || nrow(df) == 0L) return("_none_")
  df <- as.data.frame(df)
  for (nm in names(df)) {
    x <- df[[nm]]
    if (!is.numeric(x)) next
    # Rounding p-values to three decimals turns 5e-11 into "0"; keep small
    # p-values legible instead of pretending they are exactly zero.
    df[[nm]] <- if (grepl("^p", nm)) {
      ifelse(is.na(x), NA_character_,
             ifelse(x < 0.001, "<0.001", format(signif(x, 3L), scientific = FALSE)))
    } else {
      signif(x, digits)
    }
  }
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1L, function(r) {
    paste0("| ", paste(ifelse(is.na(r), "-", r), collapse = " | "), " |")
  })
  paste(c(header, sep, body), collapse = "\n")
}

dx_write_report <- function(cell, tme, meta, out_dir) {
  lines <- c(
    "# Prognoser downstream report",
    "",
    sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("- R: %s", R.version.string),
    sprintf("- Cohorts: %s", paste(meta$datasets, collapse = ", ")),
    sprintf("- Samples: %s", paste(sprintf("%s=%d", meta$datasets, meta$n_samples),
                                   collapse = ", ")),
    sprintf("- Events: %s", paste(sprintf("%s=%d", meta$datasets, meta$events),
                                  collapse = ", ")),
    ""
  )

  lines <- c(lines, "## Signature score", "")
  if (!is.null(cell$signature_survival)) {
    lines <- c(lines, "Cox model per SD of the HRGS-minus-LRGS score:", "",
               dx_md_table(dx_top(cell$signature_survival,
                                  cols = c("dataset", "HR", "lower", "upper",
                                           "p", "n", "events"))),
               "")  } else {
    lines <- c(lines, "_No Prognoser result supplied, so no signature score was computed._", "")
  }

  lines <- c(lines, "## Cell states", "")
  if (!is.null(cell$correlation)) {
    lines <- c(lines, "Strongest association with the signature score (Spearman):", "",
               dx_md_table(dx_top(cell$correlation, group = "dataset", n = 5L,
                                  cols = c("dataset", "panel", "rho", "p_adj"))),
               "")
  }
  lines <- c(lines, "Cell state versus survival (Cox per SD):", "",
             dx_md_table(dx_top(cell$survival, group = "dataset", n = 5L,
                                cols = c("dataset", "panel", "HR", "lower", "upper",
                                         "p_adj", "n", "events"))), "")

  lines <- c(lines, "## Tumour microenvironment", "")
  if (!is.null(tme$phenotype)) {
    tab <- table(tme$phenotype$dataset, tme$phenotype$phenotype)
    lines <- c(lines, "Immune phenotype calls (heuristic, tertile-based):", "",
               "```", utils::capture.output(print(tab)), "```", "")
  }
  if (!is.null(tme$estimate)) {
    wide <- tapply(tme$estimate$value,
                   list(tme$estimate$dataset, tme$estimate$feature), mean)
    lines <- c(lines, "Mean ESTIMATE scores per cohort:", "",
               dx_md_table(data.frame(dataset = rownames(wide), wide,
                                      check.names = FALSE)), "")
  }
  if (!is.null(tme$correlation)) {
    lines <- c(lines, "Strongest TME association with the signature score (negative",
               "rho means the signature marks low immune/stromal content):", "",
               dx_md_table(dx_top(tme$correlation, group = "dataset", n = 5L,
                                  cols = c("dataset", "feature", "rho", "p_adj"))),
               "")
  }
  lines <- c(lines, "TME feature versus survival (Cox per SD):", "",
             dx_md_table(dx_top(tme$survival, group = "dataset", n = 5L,
                                cols = c("dataset", "feature", "HR", "lower", "upper",
                                         "p_adj", "n", "events"))), "")

  lines <- c(lines, "## Output files", "",
             paste0("- `", sort(list.files(out_dir, pattern = "[.]csv$")), "`"),
             paste0("- `", sort(list.files(out_dir, pattern = "[.]pdf$")), "`"),
             "",
             "Panel definitions, including the exact genes used, are in",
             "`panel_definitions.csv`.", "")

  path <- file.path(out_dir, "downstream_report.md")
  writeLines(lines, path)
  dx_log("wrote %s", path)
  invisible(path)
}
