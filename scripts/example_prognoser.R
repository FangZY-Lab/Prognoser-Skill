#!/usr/bin/env Rscript

# Self-contained two-dataset Prognoser example.
#
# Simulates two cohorts (SimDataA, SimDataB), each with OS/RFS/PFS/DFS
# endpoints, planted risk and protective genes, and ten prior gene sets.
# Run with: Rscript example_prognoser.R

suppressPackageStartupMessages({
  library(Prognoser)
})

set.seed(123)

endpoints <- c("OS", "RFS", "PFS", "DFS")

make_dataset <- function(n_samples = 80, n_genes = 100, seed = 1) {
  set.seed(seed)

  # Latent risk score: higher score -> worse survival.
  risk_score <- rnorm(n_samples)

  mat <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)

  # GENE1-10 are risk genes (higher expression -> worse survival).
  mat[1:10, ] <- mat[1:10, ] +
    matrix(rep(risk_score * 2, each = 10), nrow = 10)

  # GENE11-20 are protective genes (higher expression -> better survival).
  mat[11:20, ] <- mat[11:20, ] -
    matrix(rep(risk_score * 2, each = 10), nrow = 10)

  expr <- as.data.frame(mat)
  rownames(expr) <- paste0("GENE", seq_len(n_genes))
  colnames(expr) <- paste0("S", seq_len(n_samples))

  make_surv <- function(risk_score, ep_seed) {
    set.seed(ep_seed)
    eps <- rnorm(length(risk_score), 0, 0.3)
    time <- rweibull(length(risk_score), shape = 1, scale = exp(2 - risk_score + eps))
    status <- rbinom(length(risk_score), 1, 0.7)
    data.frame(
      time = round(time, 3),
      status = status,
      row.names = paste0("S", seq_along(risk_score))
    )
  }

  list(
    expr = expr,
    OS = make_surv(risk_score, seed * 10 + 1),
    RFS = make_surv(risk_score, seed * 10 + 2),
    PFS = make_surv(risk_score, seed * 10 + 3),
    DFS = make_surv(risk_score, seed * 10 + 4)
  )
}

dA <- make_dataset(seed = 11)
dB <- make_dataset(seed = 22)

for (ds in c("A", "B")) {
  d <- get(paste0("d", ds))
  for (ep in endpoints) {
    assign(paste0("SimData", ds, "_", ep), d$expr, envir = .GlobalEnv)
    assign(
      paste0("SimData", ds, "_", ep, "_survival"),
      d[[ep]],
      envir = .GlobalEnv
    )
  }
}

datasets <- c(
  "SimDataA_OS", "SimDataA_RFS", "SimDataA_PFS", "SimDataA_DFS",
  "SimDataB_OS", "SimDataB_RFS", "SimDataB_PFS", "SimDataB_DFS"
)

knowledge <- list(
  immune_response = c("GENE1", "GENE2", "GENE3", "GENE11"),
  cell_cycle = c("GENE4", "GENE5", "GENE6"),
  apoptosis = c("GENE7", "GENE8", "GENE14"),
  dna_repair = c("GENE9", "GENE10", "GENE15"),
  metabolism = c("GENE11", "GENE12", "GENE13"),
  angiogenesis = c("GENE16", "GENE17", "GENE18"),
  inflammation = c("GENE19", "GENE20", "GENE1"),
  proliferation = c("GENE2", "GENE4", "GENE6", "GENE8"),
  signaling = c("GENE3", "GENE5", "GENE7", "GENE9"),
  epithelial = c("GENE10", "GENE12", "GENE14", "GENE16")
)

# Thresholds are relaxed for this toy cohort (100 genes, 3-gene knowledge
# sets). The package defaults are gene_threshold = 0.001 and min_genes = 5,
# which return empty gene sets on data this small. Real analyses should use the
# defaults shown in SKILL.md and references/prognoser-api.md.
result <- Prognoser(
  expression_accession_vector = datasets,
  knowledge_genesets = knowledge,
  denovo_threshold = 0.05,
  gene_threshold = 0.01,
  min_genes = 3,
  similarity_threshold = 0.5,
  linkage = "ward.D"
)

cat("Final gene sets:\n")
print(names(result$genesets))

for (name in names(result$genesets)) {
  cat("\n", name, ":\n", sep = "")
  print(result$genesets[[name]])
}
