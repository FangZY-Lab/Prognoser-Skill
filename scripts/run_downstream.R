#!/usr/bin/env Rscript

# Downstream driver: run the cell-state and TME modules on the same cohorts and
# write one Markdown report next to the tables.
#
# Usage:
#   Rscript run_downstream.R --input=input --output=output/downstream \
#     --result=output/prognoser_result.rds --datasets=GSE39582_OS,GSE17536_OS
#
# The driver deliberately keeps both modules independent: either can also be
# run on its own, and a failure in one does not silently produce half a report.

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

## ------------------------------------------------------------ analysis ---

dx_downstream_main <- function() {
  defaults <- list(
    input = "input",
    output = "output/downstream",
    result = "output/prognoser_result.rds",
    datasets = "",
    platform = "illumina",
    expr_scale = "auto",
    method = "ssgsea",
    xcell = "true",
    quantiseq = "false",
    plots = "true",
    `help` = "false"
  )
  args <- dx_parse_args(defaults)
  if (dx_as_flag(args$help)) {
    cat("Prognoser downstream driver\n\n",
        "  --input=DIR        directory with <Dataset>_<Endpoint>.rds/.Rdata\n",
        "  --datasets=A,B     dataset-endpoint names (default: auto-detect)\n",
        "  --result=FILE      Prognoser result .rds\n",
        "  --output=DIR       output directory\n",
        "  --platform=illumina  ESTIMATE platform\n",
        "  --quantiseq=false  enable quanTIseq (needs TPM)\n", sep = "")
    return(invisible(NULL))
  }

  inputs <- dx_load_inputs(args$input, dx_split(args$datasets))
  genesets <- NULL
  if (nzchar(args$result) && file.exists(args$result)) {
    res <- readRDS(args$result)
    genesets <- res$genesets %||% res
    dx_log("Loaded %d gene set(s) from %s", length(genesets), args$result)
  } else if (nzchar(args$result)) {
    dx_log("WARNING: result file not found (%s); signature scoring skipped.",
           args$result)
  }

  cfg <- list(
    inputs = inputs,
    genesets = genesets,
    output_dir = args$output,
    platform = args$platform,
    expr_scale = args$expr_scale,
    method = args$method,
    xcell = args$xcell,
    quantiseq = args$quantiseq,
    plots = args$plots
  )

  cell <- tryCatch(dx_cell_states_run(cfg), error = function(e) {
    dx_log("Cell-state module failed: %s", conditionMessage(e))
    NULL
  })
  tme <- tryCatch(dx_tme_run(cfg), error = function(e) {
    dx_log("TME module failed: %s", conditionMessage(e))
    NULL
  })
  if (is.null(cell) && is.null(tme)) {
    stop("Both downstream modules failed; see the messages above.", call. = FALSE)
  }

  meta <- list(
    datasets = names(inputs),
    n_samples = vapply(inputs, function(d) ncol(d$expr), integer(1L)),
    events = vapply(inputs, function(d) sum(d$status == 1), integer(1L))
  )
  dx_write_report(cell, tme, meta, dx_output_dir(args$output))
  invisible(list(cell = cell, tme = tme))
}

if (dx_is_main("run_downstream.R")) dx_downstream_main()
