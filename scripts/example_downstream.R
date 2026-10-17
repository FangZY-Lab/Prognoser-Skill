#!/usr/bin/env Rscript

# Self-contained demonstration of the downstream modules.
#
# Two simulated cohorts are built from real HGNC symbols so that every method in
# the TME module (ESTIMATE, xCell, marker scoring) can actually run:
#
#   * a latent risk axis drives proliferation genes and shortens survival
#   * a latent immune axis drives CD8/cytotoxic genes and lengthens survival
#   * a latent stromal axis drives EMT and CAF genes
#
# The planted structure means the expected answer is known: the HRGS-minus-LRGS
# signature score should be positively correlated with EMT/proliferation and
# negatively correlated with CD8 content.
#
# Usage:
#   Rscript example_downstream.R
#   Rscript example_downstream.R --output=output/example --write-inputs=rdata

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
dx_source_sibling("downstream_signatures.R")
dx_source_sibling("downstream_cell_states.R")
dx_source_sibling("downstream_tme.R")

## ----------------------------------------------------------- simulation ---

dx_simulate_cohort <- function(n_samples = 60L, n_genes = 12000L, seed = 1L,
                               cohort = "SimA") {
  dx_require("msigdbr", "install.packages('msigdbr')")
  set.seed(seed)

  h <- dx_msigdbr("Homo sapiens", "H")
  c5 <- dx_msigdbr("Homo sapiens", "C5", "GO:BP")
  panels <- dx_panel_genes(c("emt", "proliferation", "cd8_t_cell", "caf"))

  # The simulation universe must contain every planted gene, otherwise the
  # panels this demo is supposed to exercise are incomplete (CD8A is not in the
  # Hallmark collection, for example).
  marker_extra <- unique(unlist(dx_panel_genes(dx_marker_panel_ids())))
  universe <- unique(c(h$gene_symbol, c5$gene_symbol, marker_extra))
  planted <- unique(unlist(panels))
  if (n_genes < length(planted) + 100L) {
    stop(sprintf("n_genes must be at least %d.", length(planted) + 100L),
         call. = FALSE)
  }
  genes <- sort(c(planted, sample(setdiff(universe, planted),
                                 n_genes - length(planted))))

  risk <- stats::rnorm(n_samples)
  immune <- stats::rnorm(n_samples)
  stroma <- 0.4 * risk + stats::rnorm(n_samples)

  mu <- stats::rnorm(length(genes), mean = 7, sd = 1.2)
  load <- matrix(stats::rnorm(length(genes) * 3L, sd = 0.05), ncol = 3L,
                 dimnames = list(genes, c("risk", "immune", "stroma")))

  add <- function(gs, col, weight) {
    hit <- intersect(gs, genes)
    if (length(hit)) load[hit, col] <<- load[hit, col] + weight
    invisible(NULL)
  }
  add(panels$emt, "stroma", 0.9)
  add(panels$emt, "risk", 0.4)
  add(panels$proliferation, "risk", 0.8)
  add(panels$cd8_t_cell, "immune", 0.9)
  add(panels$caf, "stroma", 0.9)

  noise <- matrix(stats::rnorm(length(genes) * n_samples, sd = 0.6),
                  nrow = length(genes))
  expr <- mu + load[, "risk"] %o% risk + load[, "immune"] %o% immune +
    load[, "stroma"] %o% stroma + noise
  rownames(expr) <- genes
  colnames(expr) <- sprintf("%s_S%02d", cohort, seq_len(n_samples))

  lp <- 0.9 * risk - 0.6 * immune
  time <- stats::rweibull(n_samples, shape = 1.2, scale = exp(3.2 - lp))
  cens <- stats::runif(n_samples, 0.5, 2.0) * max(time)
  status <- as.integer(time <= cens)
  time <- pmin(time, cens)
  surv <- data.frame(time = round(time, 4), status = status,
                     row.names = colnames(expr))

  list(expr = expr, surv = surv)
}

## ------------------------------------------------------------------ demo ---

dx_example_main <- function() {
  defaults <- list(output = "output/example_downstream", write_inputs = "",
                   write_result = "", n_genes = "12000", seed = "7",
                   `help` = "false")
  args <- dx_parse_args(defaults)
  if (dx_as_flag(args$help)) {
    cat("Self-contained downstream demo\n\n",
        "  --output=DIR        where tables and figures are written\n",
        "  --write-inputs=DIR  also write <Dataset>_<Endpoint>.rds for the CLI\n",
        "  --write-result=DIR  also write a stand-in prognoser_result.rds\n",
        "  --n_genes=12000     simulated gene count\n", sep = "")
    return(invisible(NULL))
  }

  out_dir <- dx_output_dir(args$output)
  seed <- as.integer(args$seed)
  n_genes <- as.integer(args$n_genes)

  dx_log("Simulating two cohorts (%d genes each)...", n_genes)
  a <- dx_simulate_cohort(60L, n_genes, seed = seed, cohort = "SimA")
  b <- dx_simulate_cohort(60L, n_genes, seed = seed + 1L, cohort = "SimB")

  inputs <- list(
    SimA_OS = list(expr = a$expr, time = a$surv$time, status = a$surv$status,
                   samples = rownames(a$surv)),
    SimB_OS = list(expr = b$expr, time = b$surv$time, status = b$surv$status,
                   samples = rownames(b$surv))
  )

  if (nzchar(args$write_inputs)) {
    dir.create(args$write_inputs, showWarnings = FALSE, recursive = TRUE)
    for (ds in names(inputs)) {
      saveRDS(inputs[[ds]]$expr, file.path(args$write_inputs, paste0(ds, ".rds")))
      saveRDS(data.frame(time = inputs[[ds]]$time, status = inputs[[ds]]$status,
                         row.names = inputs[[ds]]$samples),
              file.path(args$write_inputs, paste0(ds, "_survival.rds")))
    }
    dx_log("Wrote input objects to %s", args$write_inputs)
  }

  # A stand-in for result$genesets with the planted axis genes.
  panels <- dx_panel_genes(c("proliferation", "cd8_t_cell"))
  genesets <- list(
    HRGS1 = utils::head(panels$proliferation, 40L),
    LRGS1 = utils::head(panels$cd8_t_cell, 8L)
  )

  if (nzchar(args$write_result)) {
    dir.create(args$write_result, showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(genesets = genesets),
            file.path(args$write_result, "prognoser_result.rds"))
    dx_log("Wrote a stand-in Prognoser result to %s", args$write_result)
  }

  dx_log("Running cell-state module...")
  cell <- dx_cell_states_run(list(
    inputs = inputs, genesets = genesets, output_dir = out_dir,
    method = "ssgsea", plots = TRUE
  ))

  dx_log("Running TME module (ESTIMATE + xCell)...")
  tme <- dx_tme_run(list(
    inputs = inputs, genesets = genesets, output_dir = out_dir,
    platform = "illumina", estimate = TRUE, xcell = TRUE,
    quantiseq = FALSE, plots = TRUE
  ))

  meta <- list(datasets = names(inputs),
               n_samples = vapply(inputs, function(d) ncol(d$expr), integer(1L)),
               events = vapply(inputs, function(d) sum(d$status == 1), integer(1L)))
  dx_write_report(cell, tme, meta, out_dir)

  cat("\n=== expected vs observed (SimA_OS) ===\n")
  corr <- cell$correlation[cell$correlation$dataset == "SimA_OS", ]
  for (i in seq_len(nrow(corr))) {
    cat(sprintf("  %-14s rho = %+.3f  p_adj = %.3g\n",
                corr$feature[i], corr$rho[i], corr$p_adj[i]))
  }
  cat("  (planted: signature should rise with proliferation, fall with CD8)\n")

  cat("\n=== TME association (SimA_OS) ===\n")
  tcorr <- tme$correlation[tme$correlation$dataset == "SimA_OS", ]
  tcorr <- tcorr[order(tcorr$p), ]
  for (i in seq_len(min(6L, nrow(tcorr)))) {
    cat(sprintf("  %-28s rho = %+.3f  p_adj = %.3g\n",
                tcorr$feature[i], tcorr$rho[i], tcorr$p_adj[i]))
  }

  w <- warnings()
  if (length(w) > 0L) {
    cat("\n", length(w), " warning(s) were raised:\n", sep = "")
    for (nm in unique(names(w))) cat("  - ", nm, "\n", sep = "")
  }
  cat("\nOutputs in: ", out_dir, "\n", sep = "")
  invisible(list(cell = cell, tme = tme))
}

if (dx_is_main("example_downstream.R")) dx_example_main()
