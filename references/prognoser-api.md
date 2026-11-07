# Prognoser API

## Installation

```r
install.packages(c("survival", "RobustRankAggreg", "meta", "igraph", "remotes"))
remotes::install_github("FangZY-Lab/Prognoser")
```

Dependencies: `R >= 3.5.0`, `survival`, `RobustRankAggreg`, `meta`, `igraph`,
`stats`.

## Function signature

```r
Prognoser(
  expression_accession_vector,
  denovo_threshold = 0.05,
  knowledge_genesets,
  na_ratio = 0.3,
  gene_threshold = 0.001,
  min_genes = 5,
  similarity_threshold = 0.5,
  linkage = "ward.D",
  de_redundant_basis = c("agglomeration", "igraph_component")
)
```

`linkage` takes exactly one of `"ward.D"`, `"ward.D2"`, `"single"`,
`"complete"`, `"average"`, `"mcquitty"`, `"median"`, `"centroid"`; it is
ignored when `de_redundant_basis = "igraph_component"`.

## Arguments

| Argument | Meaning | Default |
| --- | --- | --- |
| `expression_accession_vector` | Dataset-endpoint names present in the global environment | required |
| `denovo_threshold` | Cox p-value cutoff for protective/risk genes | `0.05` |
| `knowledge_genesets` | Named list of prior gene sets | required |
| `na_ratio` | Maximum NA proportion per gene across datasets | `0.3` |
| `gene_threshold` | RRA p-value cutoff for the final gene pool | `0.001` |
| `min_genes` | Minimum genes in a final clustered set | `5` |
| `similarity_threshold` | Jaccard cutoff for merging redundant sets | `0.5` |
| `linkage` | Hierarchical clustering method (agglomeration only) | `"ward.D"` |
| `de_redundant_basis` | `"agglomeration"` (hierarchical clustering) or `"igraph_component"` (connected components of the Jaccard graph) | `"agglomeration"` |

## Global environment contract

For `TCGA_COAD_OS` in `expression_accession_vector`:

- `TCGA_COAD_OS`: expression data frame, genes x samples.
- `TCGA_COAD_OS_survival`: data frame with `time` and `status`, row names
  matching expression columns.

```r
datasets <- c(
  "TCGA_COAD_OS", "TCGA_COAD_RFS", "TCGA_COAD_PFS", "TCGA_COAD_DFS",
  "GSE39582_OS", "GSE39582_RFS"
)
```

`TCGA_COAD_*` entries are meta-analyzed as one dataset; `GSE39582_*` entries
as another.

## Method flow

1. Load expression and survival by `get()`.
2. Group entries by the prefix before `_`.
3. Remove zero-variance genes.
4. Fit a univariate Cox model per gene.
5. Meta-analyze endpoints within a dataset with random-effects DerSimonian-Laird.
6. Split genes into protective (`HR < 1`) and risk (`HR >= 1`).
7. Build signed effect matrices (`sign(log HR) * -log10(P)`).
8. Rank genes with Robust Rank Aggregation.
9. Intersect de novo genes with prior knowledge sets.
10. Resolve Jaccard similarity between gene sets, by hierarchical clustering
    (`"agglomeration"`) or by the connected components of the similarity graph
    (`"igraph_component"`), and return non-redundant sets.

## Output

```r
result$genesets                 # HRGS* and LRGS* gene sets
result$combined_matrix_protection
result$combined_matrix_risk
```

`combined_matrix_*` contain one column per dataset, plus `p_protection` or
`p_risk` from RRA. With `de_redundant_basis = "igraph_component"` the list
also contains `network_protection` and `network_risk`, data frames with
columns `from`, `to`, and `jaccard`.

## Saving results

```r
saveRDS(result, "prognoser_result.rds")

for (name in names(result$genesets)) {
  writeLines(
    paste(result$genesets[[name]], collapse = "\t"),
    paste0(name, ".gmt")
  )
}
```

Use the same gene-set files as downstream ssGSEA/GSVA input when scoring new
cohorts.

## Common failures

- Missing `_survival` object for an entry.
- Survival `time`/`status` stored as character or factor.
- Duplicate sample names in survival row names.
- Passing several `linkage` values at once.
- `knowledge_genesets` left missing.
