#!/usr/bin/env Rscript

# Downstream module 2/2: tumour-microenvironment (TME) analysis.
#
# Three complementary views of the same cohorts:
#
#   * Bulk scores: ESTIMATE (stromal / immune / overall / purity).
#   * Deconvolution: xCell cell-type enrichment, optional quanTIseq fractions.
#   * Marker panels: CD8, Treg, M1/M2, MDSC, checkpoint, TLS, CAF, endothelial.
#
# Plus an explicit, reproducible immune-phenotype call (inflamed / excluded /
# desert) built from the CD8 and IFN-gamma scores, with its thresholds written
# to disk so the classification can be audited or re-cut.
#
# Usage:
#   Rscript downstream_tme.R --input=input --output=output/downstream \
#     --result=output/prognoser_result.rds --platform=illumina

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

## ----------------------------------------------------------- phenotype ---

# Heuristic immune phenotype. The Chen & Mellman categories are qualitative; the
# cut points here are tertiles of within-cohort z-scores and are reported in
# the output so a reader can reproduce or move them.
dx_immune_phenotype <- function(cd8_z, ifng_z) {
  q <- stats::quantile(cd8_z, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
  ifng_cut <- stats::quantile(ifng_z, probs = 1 / 2, na.rm = TRUE, names = FALSE)
  cd8_level <- ifelse(cd8_z <= q[1L], "low", ifelse(cd8_z >= q[2L], "high", "mid"))
  phenotype <- ifelse(
    cd8_level == "low", "desert",
    ifelse(ifng_z >= ifng_cut, "inflamed", "excluded")
  )
  list(
    table = data.frame(
      sample = names(cd8_z),
      cd8_z = as.numeric(cd8_z),
      ifng_z = as.numeric(ifng_z),
      cd8_level = cd8_level,
      phenotype = phenotype,
      stringsAsFactors = FALSE
    ),
    thresholds = c(cd8_low = q[1L], cd8_high = q[2L], ifng_high = ifng_cut)
  )
}

## ------------------------------------------------------------ analysis ---

dx_tme_run <- function(cfg) {
  stopifnot(is.list(cfg))
  inputs <- cfg$inputs
  if (is.null(inputs) || length(inputs) == 0L) {
    stop("cfg$inputs is empty; run dx_load_inputs() first.", call. = FALSE)
  }
  out_dir <- dx_output_dir(cfg$output_dir)
  plots <- dx_as_flag(cfg$plots, default = TRUE)
  platform <- match.arg(cfg$platform %||% "illumina",
                        c("illumina", "affymetrix", "agilent"))
  expr_scale <- cfg$expr_scale %||% "auto"
  want_estimate <- dx_as_flag(cfg$estimate, default = TRUE)
  want_xcell <- dx_as_flag(cfg$xcell, default = TRUE)
  want_quantiseq <- dx_as_flag(cfg$quantiseq, default = FALSE)
  marker_ids <- dx_split(cfg$markers)
  if (length(marker_ids) == 0L) marker_ids <- dx_marker_panel_ids()

  dx_log("TME module: %d dataset(s), markers = %s",
         length(inputs), paste(marker_ids, collapse = ", "))

  marker_sets <- dx_panel_genes(marker_ids)
  ifn_set <- dx_panel_genes("ifn_gamma")

  pools <- NULL
  if (!is.null(cfg$genesets) && length(cfg$genesets) > 0L) {
    pools <- dx_genesets_to_pools(cfg$genesets)
  }

  estimate_long <- list()
  celltype_long <- list()
  fraction_long <- list()
  marker_long <- list()
  pheno_long <- list()
  thresh_long <- list()
  surv_rows <- list()
  corr_rows <- list()

  for (ds in names(inputs)) {
    d <- inputs[[ds]]
    expr <- d$expr
    scale_used <- if (identical(expr_scale, "auto")) dx_detect_scale(expr) else expr_scale
    if (!scale_used %in% c("log2", "linear")) {
      stop("expr_scale must be auto, log2 or linear.", call. = FALSE)
    }
    if (identical(scale_used, "log2") && max(expr) > 30) {
      dx_log("[%s] WARNING: values look linear but were treated as log2; ",
             ds, "pass --expr-scale=linear if that is wrong.")
    }
    expr_log2 <- if (identical(scale_used, "log2")) expr else log2(pmax(expr, 0) + 1)
    expr_linear <- if (identical(scale_used, "log2")) dx_to_linear(expr) else pmax(expr, 0)

    dx_log("[%s] %d genes x %d samples (input scale: %s)",
           ds, nrow(expr), ncol(expr), scale_used)

    marker_scored <- dx_score_sets(expr_log2, marker_sets)
    ifn_scored <- dx_score_sets(expr_log2, ifn_set)
    markers <- rbind(marker_scored$scores, ifn_scored$scores)

    features <- list()

    if (want_estimate) {
      est <- tryCatch(dx_estimate_scores(expr_log2, platform = platform),
                      error = function(e) {
                        dx_log("[%s] ESTIMATE failed: %s", ds, conditionMessage(e))
                        NULL
                      })
      if (!is.null(est)) {
        rownames(est) <- paste0("estimate:", rownames(est))
        features$estimate <- est
        estimate_long[[ds]] <- data.frame(
          dataset = ds, sample = rep(colnames(est), each = nrow(est)),
          feature = rep(rownames(est), times = ncol(est)),
          value = as.numeric(est),
          n_common_genes = attr(est, "n_common_genes"),
          stringsAsFactors = FALSE
        )
        dx_log("[%s] ESTIMATE: %d feature(s) from %d common genes", ds,
               nrow(est), attr(est, "n_common_genes"))
      }
    }

    if (want_xcell) {
      xc <- tryCatch(
        dx_xcell_scores(expr_linear),
        error = function(e) {
          dx_log("[%s] xCell skipped: %s", ds, conditionMessage(e))
          NULL
        }
      )
      if (!is.null(xc)) {
        rownames(xc) <- paste0("xcell:", rownames(xc))
        features$xcell <- xc
        celltype_long[[ds]] <- data.frame(
          dataset = ds, sample = rep(colnames(xc), each = nrow(xc)),
          feature = rep(rownames(xc), times = ncol(xc)),
          value = as.numeric(xc),
          stringsAsFactors = FALSE
        )
        dx_log("[%s] xCell: %d features", ds, nrow(xc))
      }
    }

    if (want_quantiseq) {
      q <- tryCatch(dx_quantiseq(expr_linear),
                    error = function(e) {
                      dx_log("[%s] quanTIseq skipped: %s", ds, conditionMessage(e))
                      NULL
                    })
      if (!is.null(q)) {
        rownames(q) <- paste0("quantiseq:", rownames(q))
        features$quantiseq <- q
        fraction_long[[ds]] <- data.frame(
          dataset = ds, sample = rep(colnames(q), each = nrow(q)),
          feature = rep(rownames(q), times = ncol(q)),
          value = as.numeric(q),
          stringsAsFactors = FALSE
        )
      }
    }

    features$marker <- markers
    marker_long[[ds]] <- data.frame(
      dataset = ds, sample = rep(colnames(markers), each = nrow(markers)),
      feature = rep(rownames(markers), times = ncol(markers)),
      value = as.numeric(markers),
      stringsAsFactors = FALSE
    )

    featmat <- do.call(rbind, features)
    if (anyDuplicated(rownames(featmat))) {
      featmat <- featmat[!duplicated(rownames(featmat)), , drop = FALSE]
    }

    # Immune phenotype from CD8 and IFN-gamma, both z-scored within cohort.
    cd8 <- dx_row(markers, "cd8_t_cell")
    ifn <- dx_row(markers, "ifn_gamma")
    if (is.null(cd8) || is.null(ifn)) {
      dx_log("[%s] CD8 or IFN-gamma panel was not scored; immune phenotype skipped.", ds)
    } else if (sum(is.finite(cd8)) > 5L) {
      cd8_z <- as.numeric(scale(cd8))
      ifn_z <- as.numeric(scale(ifn))
      names(cd8_z) <- colnames(markers)
      names(ifn_z) <- colnames(markers)
      pheno <- dx_immune_phenotype(cd8_z, ifn_z)
      pheno$table$dataset <- ds
      pheno_long[[ds]] <- pheno$table
      thresh_long[[length(thresh_long) + 1L]] <- data.frame(
        dataset = ds,
        criterion = names(pheno$thresholds),
        z_value = as.numeric(pheno$thresholds),
        stringsAsFactors = FALSE
      )
      dx_log("[%s] immune phenotype: %s", ds,
             paste(names(table(pheno$table$phenotype)), table(pheno$table$phenotype),
                   sep = "=", collapse = ", "))
    }

    # Survival of every continuous TME feature, per SD.
    for (f in rownames(featmat)) {
      cox <- dx_cox_per_sd(featmat[f, ], d$time, d$status)
      surv_rows[[length(surv_rows) + 1L]] <- data.frame(
        dataset = ds, feature = f, cox, stringsAsFactors = FALSE
      )
    }

    if (!is.null(pools)) {
      sig <- dx_signature_score(expr_log2, pools$risk, pools$protective)
      corr <- dx_spearman(as.numeric(sig), featmat)
      corr$dataset <- ds
      corr_rows[[length(corr_rows) + 1L]] <- corr
    }
  }

  out <- list()

  bind <- function(x) {
    if (length(x) == 0L) return(NULL)
    df <- do.call(rbind, x)
    rownames(df) <- NULL
    df
  }

  est_df <- bind(estimate_long)
  if (!is.null(est_df)) {
    dx_write_csv(est_df, out_dir, "tme_estimate_scores.csv")
    out$estimate <- est_df
  }
  cell_df <- bind(celltype_long)
  if (!is.null(cell_df)) {
    dx_write_csv(cell_df, out_dir, "tme_xcell_scores.csv")
    out$xcell <- cell_df
  }
  frac_df <- bind(fraction_long)
  if (!is.null(frac_df)) {
    dx_write_csv(frac_df, out_dir, "tme_quantiseq_fractions.csv")
    out$quantiseq <- frac_df
  }
  marker_df <- bind(marker_long)
  if (!is.null(marker_df)) {
    dx_write_csv(marker_df, out_dir, "tme_marker_scores.csv")
    out$markers <- marker_df
  }
  pheno_df <- bind(pheno_long)
  if (!is.null(pheno_df)) {
    dx_write_csv(pheno_df, out_dir, "tme_immune_phenotype.csv")
    out$phenotype <- pheno_df
  }
  thresh_df <- bind(thresh_long)
  if (!is.null(thresh_df)) {
    dx_write_csv(thresh_df, out_dir, "tme_immune_phenotype_thresholds.csv")
    out$phenotype_thresholds <- thresh_df
  }
  surv_df <- bind(surv_rows)
  if (!is.null(surv_df)) {
    surv_df$p_adj <- dx_bh(surv_df$p)
    dx_write_csv(surv_df, out_dir, "tme_survival.csv")
    out$survival <- surv_df
  }
  corr_df <- bind(corr_rows)
  if (!is.null(corr_df)) {
    corr_df$p_adj <- dx_bh(corr_df$p)
    dx_write_csv(corr_df, out_dir, "tme_signature_correlation.csv")
    out$correlation <- corr_df
  }

  if (isTRUE(plots) && !is.null(marker_df)) {
    for (ds in unique(marker_df$dataset)) {
      sub <- marker_df[marker_df$dataset == ds, , drop = FALSE]
      feats <- sort(unique(sub$feature))
      samps <- unique(sub$sample)
      mat <- matrix(NA_real_, nrow = length(feats), ncol = length(samps),
                    dimnames = list(feats, samps))
      mat[cbind(match(sub$feature, feats), match(sub$sample, samps))] <- sub$value
      dx_plot_heatmap(dx_zscore(mat), out_dir,
                      sprintf("tme_markers_heatmap_%s.pdf",
                              gsub("[^A-Za-z0-9_]", "_", ds)),
                      title = sprintf("TME marker scores, %s", ds))
    }
    sub <- surv_df[grep("^estimate:|^marker:", surv_df$feature), , drop = FALSE]
    if (nrow(sub) > 0L) {
      ds1 <- sub[sub$dataset == sub$dataset[1L], , drop = FALSE]
      ds1$label <- sprintf("%s (n=%d)", sub("^(estimate|marker):", "", ds1$feature), ds1$n)
      dx_plot_forest(ds1[, c("label", "HR", "lower", "upper", "p")], out_dir,
                     "tme_survival_forest.pdf",
                     title = sprintf("TME features, %s (HR per SD)", ds1$dataset[1L]))
    }
  }

  dx_write_csv(dx_panel_table(), out_dir, "panel_definitions.csv")
  dx_session_info(out_dir)
  invisible(out)
}

## ----------------------------------------------------------------- CLI ---

dx_tme_main <- function() {
  defaults <- list(
    input = "input",
    output = "output/downstream",
    result = "",
    datasets = "",
    platform = "illumina",
    expr_scale = "auto",
    estimate = "true",
    xcell = "true",
    quantiseq = "false",
    markers = "",
    plots = "true",
    `help` = "false"
  )
  args <- dx_parse_args(defaults)
  if (dx_as_flag(args$help)) {
    cat("Downstream tumour-microenvironment analysis\n\n",
        "  --input=DIR          directory with <Dataset>_<Endpoint>.rds/.Rdata\n",
        "  --datasets=A,B       dataset-endpoint names (default: auto-detect)\n",
        "  --result=FILE        Prognoser result .rds (enables correlation)\n",
        "  --output=DIR         output directory\n",
        "  --platform=illumina  illumina | affymetrix | agilent (ESTIMATE)\n",
        "  --expr-scale=auto    auto | log2 | linear\n",
        "  --xcell=true         run xCell\n",
        "  --quantiseq=false    run quanTIseq (needs TPM)\n",
        "  --markers=A,B        marker panel IDs\n",
        "  --plots=true         write PDF figures\n", sep = "")
    return(invisible(NULL))
  }
  inputs <- dx_load_inputs(args$input, dx_split(args$datasets))
  genesets <- NULL
  if (nzchar(args$result)) {
    if (!file.exists(args$result)) {
      stop(sprintf("Result file not found: %s", args$result), call. = FALSE)
    }
    res <- readRDS(args$result)
    genesets <- res$genesets %||% res
    dx_log("Loaded %d gene set(s) from %s", length(genesets), args$result)
  }
  dx_tme_run(list(
    inputs = inputs,
    genesets = genesets,
    output_dir = args$output,
    platform = args$platform,
    expr_scale = args$expr_scale,
    estimate = args$estimate,
    xcell = args$xcell,
    quantiseq = args$quantiseq,
    markers = args$markers,
    plots = args$plots
  ))
}

if (dx_is_main("downstream_tme.R")) dx_tme_main()
