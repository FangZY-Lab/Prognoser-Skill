# Preprocessing and normalization

The goal is one clean, comparable numeric matrix per dataset-endpoint before
passing data to `Prognoser`.

## Required orientation

- Genes are rows.
- Samples are columns.
- Column names are sample IDs.
- Row names are HGNC gene symbols.

## Gene identifier mapping

Map identifiers in this order of preference:

1. HGNC symbols already present.
2. Probe IDs using the platform annotation or an annotation package.
3. Ensembl IDs using `biomaRt` or an `org.*.eg.db` package.

```r
library(org.Hs.eg.db)
library(AnnotationDbi)

map <- select(
  org.Hs.eg.db,
  keys = rownames(expr),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)
```

Keep only rows that map to a non-NA symbol.

## Duplicate symbols

Collapse duplicate symbols before analysis:

```r
expr <- as.matrix(expr)
expr <- collapseRows(expr, rowGroup = symbol, rowID = rownames(expr))
```

If `WGCNA::collapseRows` is unavailable, aggregate by mean:

```r
expr <- rowsum(expr, group = symbol, reorder = FALSE) /
        as.numeric(table(symbol)[symbol[!duplicated(symbol)]])
```

The divisor must have one entry per *unique* symbol, in the order those symbols
first appear, which is the order `rowsum(reorder = FALSE)` returns. Indexing the
table by the full `symbol` vector gives one count per row instead, so the
division is recycled against the wrong rows.

Prefer mean for microarray, mean of log values for RNA-seq.

## Multi-mapped symbols

Some platforms annotate a probe with several symbols, separated by `///`:

```text
CARD16 /// CASP1
```

Keep the first symbol and drop the rest, so the probe is not silently lost and
does not create a compound name that no gene set can match:

```r
symbol <- trimws(sub("\\s*///.*$", "", as.character(symbol)))
```

This happens after the rows are filtered to non-missing symbols and before
the duplicate collapse, so several probes can still land on `CARD16` and be
averaged together.


## Missing values

- Remove genes with `NA` across most samples.
- `Prognoser` itself applies a per-gene `na_ratio` filter during the effect
  matrix stage, but clean input prevents silent drops.

## Log transformation

Apply `log2(x + 1)` when values are:

- raw counts
- TPM
- RSEM
- FPKM

Keep already log-scaled matrices unchanged. Do not log-transform twice.

```r
if (max(expr, na.rm = TRUE) > 100) {
  expr <- log2(expr + 1)
}
```

## Cross-cohort alignment

For multi-cohort analysis, intersect genes across matrices:

```r
common_genes <- Reduce(intersect, lapply(expr_list, rownames))
expr_list <- lapply(expr_list, function(x) x[common_genes, , drop = FALSE])
```

Prognoser meta-analyzes cohorts after per-dataset Cox regression, so exact
gene alignment is required.

## Variance and scale

Prognoser removes zero-variance genes internally during Cox fitting. You do
not need to scale expression globally; keep the same scale convention across
all datasets of a project so hazard ratios remain comparable.

## Survival table

The survival table must have:

- column `time`, numeric and non-negative
- column `status`, numeric and coded `0/1`
- row names equal to expression column names

```r
surv <- data.frame(
  time = as.numeric(time),
  status = as.numeric(status),
  row.names = samples
)
surv$status <- as.integer(surv$status != 0)
```
