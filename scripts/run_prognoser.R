#!/usr/bin/env Rscript

# End-to-end Prognoser execution template.
#
# Expected inputs (loaded into the global environment):
#   <Dataset>_<Endpoint>             expression data frame, genes x samples
#   <Dataset>_<Endpoint>_survival    data frame with `time` and `status`

suppressPackageStartupMessages({
  library(Prognoser)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "input"
output_dir <- if (length(args) >= 2) args[[2]] else "output"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

load_prognoser_objects <- function(path) {
  files <- list.files(path, pattern = "[.](rds|RData|Rdata|rda|Rda)$",
                        full.names = TRUE)
  for (file in files) {
    if (grepl("[.]rds$", file)) {
      obj <- readRDS(file)
      assign(tools::file_path_sans_ext(basename(file)), obj, envir = .GlobalEnv)
    } else {
      load(file, envir = .GlobalEnv)
    }
  }
}

load_prognoser_objects(input_dir)

datasets <- c(
  # Add dataset-endpoint names here, e.g. "TCGA_COAD_OS".
)

knowledge_genesets <- list(
  # Add prior gene sets here, e.g. immune = c("CD8A", "GZMB", "PRF1").
)

if (length(datasets) == 0) {
  stop("No dataset names given: fill in `datasets` with the <Dataset>_<Endpoint> ",
       "objects that are loaded in the global environment.", call. = FALSE)
}

result <- Prognoser(
  expression_accession_vector = datasets,
  knowledge_genesets = knowledge_genesets,
  denovo_threshold = 0.05,
  gene_threshold = 0.001,
  min_genes = 5,
  similarity_threshold = 0.5,
  # "agglomeration" (hierarchical clustering) or "igraph_component"
  de_redundant_basis = "agglomeration",
  linkage = "ward.D"
)

saveRDS(result, file.path(output_dir, "prognoser_result.rds"))

for (name in names(result$genesets)) {
  writeLines(
    paste(result$genesets[[name]], collapse = "\t"),
    file.path(output_dir, paste0(name, ".gmt"))
  )
}

message("Saved Prognoser results to: ", output_dir)
