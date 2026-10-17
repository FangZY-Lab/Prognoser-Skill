#!/usr/bin/env Rscript

# Downstream module 1/2: cell-state analysis.
#
# Scores prognostic signatures and cell-state programmes in the same cohorts
# used for Prognoser, then asks two questions per cohort:
#
#   1. Does the signature score track a cell state (Spearman, BH-adjusted)?
#   2. Does a cell state carry survival information on its own (Cox per SD)?
#
# Usage:
#   Rscript downstream_cell_states.R --input=input --output=output/downstream \
#     --result=output/prognoser_result.rds --datasets=GSE39582_OS,GSE39582_RFS
#
# All arguments are optional except that either --input or --result must point
# at real data.

if (!exists("dx_parse_args", inherits = TRUE)) {
  .dx_boot_args <- commandArgs(trailingOnly = FALSE)
  .dx_boot_hit <- grep("^--file=", .dx_boot_args, value = TRUE)
  .dx_boot_dir <- if (length(.dx_boot_hit) > 0L) {
    dirname(sub("^--file=", "", .dx_boot_hit[[1L]]))
  } else {
    getwd()
  }
  .dx_boot_cand <- c(file.path(.dx_boot_dir, "downstream_lib.R"),
                     file.path(getwd(), "scripts", "downstream_lib.R"),
                     file.path(getwd(), "downstream_lib.R"))
  .dx_boot_lib <- .dx_boot_cand[file.exists(.dx_boot_cand)][1L]
  if (is.na(.dx_boot_lib)) {
    stop("Cannot locate downstream_lib.R; set PROGNOSER_SCRIPTS_DIR.")
  }
  source(.dx_boot_lib)
  rm(.dx_boot_args, .dx_boot_hit, .dx_boot_dir, .dx_boot_cand, .dx_boot_lib)
}
if (!exists("DX_STATE_PANELS")) dx_source_sibling("downstream_signatures.R")

## ------------------------------------------------------------ analysis ---

dx_cell_states_run <- function(cfg) {
  stopifnot(is.list(cfg))
  inputs <- cfg$inputs
  if (is.null(inputs) || length(inputs) == 0L) {
    stop("cfg$inputs is empty; run dx_load_inputs() first.", call. = FALSE)
  }
  method <- match.arg(cfg$method %||% "ssgsea", c("ssgsea", "zscore"))
  species <- cfg$species %||% "Homo sapiens"
  panels <- cfg$panels %||% dx_state_panel_ids()
  out_dir <- dx_output_dir(cfg$output_dir)
  plots <- dx_as_flag(cfg$plots, default = TRUE)
  min_size <- as.integer(cfg$min_size %||% 5L)
  max_size <- as.integer(cfg$max_size %||% 500L)

  dx_log("Cell-state module: %d dataset(s), %d panel(s), method = %s",
         length(inputs), length(panels), method)

  gene_sets <- dx_panel_genes(panels, species = species)
  dx_log("Resolved %d panel(s): %s", length(gene_sets),
         paste(sprintf("%s (%d genes)", names(gene_sets),
                       vapply(gene_sets, length, integer(1L))), collapse = "; "))

  pools <- NULL
  if (!is.null(cfg$genesets) && length(cfg$genesets) > 0L) {
    pools <- dx_genesets_to_pools(cfg$genesets)
    dx_log("Signature pools: %d risk gene(s), %d protective gene(s)",
           length(pools$risk), length(pools$protective))
  }

  coverage <- list()
  score_long <- list()
  surv_rows <- list()
  corr_rows <- list()
  sig_surv_rows <- list()
  combined <- list()

  for (ds in names(inputs)) {
    d <- inputs[[ds]]
    expr <- d$expr
    dx_log("[%s] %d genes x %d samples, %d events", ds, nrow(expr), ncol(expr),
           sum(d$status == 1))

    scored <- dx_score_sets(expr, gene_sets, method = method,
                            min_size = min_size, max_size = max_size)
    cov <- scored$coverage
    cov$dataset <- ds
    coverage[[ds]] <- cov
    scores <- scored$scores

    score_long[[ds]] <- data.frame(
      dataset = ds,
      sample = rep(colnames(scores), each = nrow(scores)),
      panel = rep(rownames(scores), times = ncol(scores)),
      score = as.numeric(scores),
      stringsAsFactors = FALSE
    )

    # Survival of each cell-state score, per SD of the score.
    for (p in rownames(scores)) {
      cox <- dx_cox_per_sd(scores[p, ], d$time, d$status)
      surv_rows[[length(surv_rows) + 1L]] <- data.frame(
        dataset = ds, panel = p, cox, stringsAsFactors = FALSE
      )
    }

    if (!is.null(pools)) {
      sig <- dx_signature_score(expr, pools$risk, pools$protective)
      cox <- dx_cox_per_sd(sig, d$time, d$status)
      sig_surv_rows[[length(sig_surv_rows) + 1L]] <- data.frame(
        dataset = ds, signature = "HRGS_minus_LRGS", cox,
        n_risk_genes = attr(sig, "n_risk"),
        n_protective_genes = attr(sig, "n_protective"),
        stringsAsFactors = FALSE
      )
      corr <- dx_spearman(as.numeric(sig), scores)
      names(corr)[names(corr) == "feature"] <- "panel"
      corr$dataset <- ds
      corr_rows[[length(corr_rows) + 1L]] <- corr
    }

    # Panel scores z-scored within cohort, for cross-cohort display only.
    combined[[ds]] <- dx_zscore(scores)
  }

  coverage_df <- do.call(rbind, coverage)
  rownames(coverage_df) <- NULL
  score_df <- do.call(rbind, score_long)
  rownames(score_df) <- NULL
  surv_df <- do.call(rbind, surv_rows)
  rownames(surv_df) <- NULL
  surv_df$p_adj <- dx_bh(surv_df$p)

  dx_write_csv(coverage_df, out_dir, "cell_state_coverage.csv")
  dx_write_csv(score_df, out_dir, "cell_state_scores_long.csv")
  dx_write_csv(surv_df, out_dir, "cell_state_survival.csv")

  out <- list(coverage = coverage_df, scores = score_df, survival = surv_df,
              plots_dir = out_dir)

  if (!is.null(corr_rows)) {
    corr_df <- do.call(rbind, corr_rows)
    rownames(corr_df) <- NULL
    corr_df$p_adj <- dx_bh(corr_df$p)
    dx_write_csv(corr_df, out_dir, "cell_state_signature_correlation.csv")
    out$correlation <- corr_df
  }

  if (length(sig_surv_rows) > 0L) {
    sig_df <- do.call(rbind, sig_surv_rows)
    rownames(sig_df) <- NULL
    sig_df$p_adj <- dx_bh(sig_df$p)
    dx_write_csv(sig_df, out_dir, "signature_score_survival.csv")
    out$signature_survival <- sig_df
  }

  if (isTRUE(plots) && length(combined) > 0L) {
    panels_all <- rownames(combined[[1L]])
    mat <- do.call(cbind, combined)
    ann <- data.frame(dataset = rep(names(combined),
                                    times = vapply(combined, ncol, integer(1L))),
                      row.names = colnames(mat))
    dx_plot_heatmap(mat, out_dir, "cell_state_heatmap.pdf",
                    annotation = ann,
                    title = "Cell-state scores (z-scored within cohort)")
    # One forest plot per small-multiple group keeps labels readable.
    for (ds in unique(surv_df$dataset)) {
      sub <- surv_df[surv_df$dataset == ds & is.finite(surv_df$HR), , drop = FALSE]
      if (nrow(sub) == 0L) next
      sub$label <- sprintf("%s (n=%d, ev=%d)", sub$panel, sub$n, sub$events)
      dx_plot_forest(sub[, c("label", "HR", "lower", "upper", "p")], out_dir,
                     sprintf("cell_state_forest_%s.pdf", gsub("[^A-Za-z0-9_]", "_", ds)),
                     title = sprintf("Cell-state survival, %s (HR per SD)", ds))
    }
  }

  dx_write_csv(dx_panel_table(), out_dir, "panel_definitions.csv")
  dx_session_info(out_dir)
  invisible(out)
}

## ----------------------------------------------------------------- CLI ---

dx_cell_states_main <- function() {
  defaults <- list(
    input = "input",
    output = "output/downstream",
    result = "",
    datasets = "",
    method = "ssgsea",
    species = "Homo sapiens",
    panels = "",
    plots = "true",
    min_size = "5",
    max_size = "500",
    `help` = "false"
  )
  args <- dx_parse_args(defaults)
  if (dx_as_flag(args$help)) {
    cat("Downstream cell-state analysis\n\n",
        "  --input=DIR        directory with <Dataset>_<Endpoint>.rds/.Rdata\n",
        "  --datasets=A,B     dataset-endpoint names (default: auto-detect)\n",
        "  --result=FILE      Prognoser result .rds (enables signature scoring)\n",
        "  --output=DIR       output directory (default output/downstream)\n",
        "  --method=ssgsea    ssgsea | zscore\n",
        "  --panels=A,B       panel IDs (default: all cell-state panels)\n",
        "  --plots=true       write PDF figures\n", sep = "")
    return(invisible(NULL))
  }
  dx_require("GSVA", "BiocManager::install('GSVA')")
  inputs <- dx_load_inputs(args$input, dx_split(args$datasets))
  genesets <- NULL
  if (nzchar(args$result)) {
    if (!file.exists(args$result)) {
      stop(sprintf("Result file not found: %s", args$result), call. = FALSE)
    }
    res <- readRDS(args$result)
    genesets <- res$genesets %||% res
    dx_log("Loaded %d gene set(s) from %s", length(genesets), args$result)
  } else {
    dx_log("No --result supplied; signature scoring skipped.")
  }
  panels <- dx_split(args$panels)
  if (length(panels) == 0L) panels <- dx_state_panel_ids()
  dx_cell_states_run(list(
    inputs = inputs,
    genesets = genesets,
    output_dir = args$output,
    method = args$method,
    species = args$species,
    panels = panels,
    plots = args$plots,
    min_size = args$min_size,
    max_size = args$max_size
  ))
}

if (dx_is_main("downstream_cell_states.R")) dx_cell_states_main()
