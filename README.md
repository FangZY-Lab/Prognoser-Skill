# Prognoser AI Agent Skill

An installable AI-agent skill for end-to-end discovery of prognostic gene
signatures with the [`Prognoser`](https://github.com/FangZY-Lab/Prognoser) R
package.

The skill guides an agent through the complete research workflow:

1. **Data acquisition** from GEO, TCGA (GDC), or the UCSC Xena browser.
2. **Multi-endpoint assembly** for OS, RFS, PFS, and DFS survival data.
3. **Preprocessing and normalization**, including identifier mapping,
   duplicate collapse, and log transformation.
4. **Prognoser execution** with the required global-environment naming
   contract.
5. **Result export** as RDS objects and GMT-compatible gene-set files.

## Repository layout

```text
.
├── SKILL.md                         Skill entrypoint
├── README.md                        This document
├── references/
│   ├── data-acquisition.md          Public data download and endpoints
│   ├── preprocessing.md             Normalization and identifier alignment
│   ├── prognoser-api.md             Full API and output schema
│   └── real-data.md                 Validated real colorectal cohorts
└── scripts/
    ├── run_prognoser.R              Reproducible execution template
    ├── example_prognoser.R          Self-contained two-cohort demo
    └── download_real_data.R         Download two real CRC cohorts
```

## Installation

Copy this repository into your Codex skills directory:

```bash
cp -R Prognoser-Skill ~/.codex/skills/prognoser
```

Or use the Codex Skill Installer with the repository URL:

```text
https://github.com/FangZY-Lab/Prognoser-Skill
```

The skill is selected automatically when a request involves multi-cohort or
multi-endpoint prognostic signature discovery with Prognoser.

## Minimal two-cohort example

Run the bundled example to verify the environment and data contract:

```bash
Rscript scripts/example_prognoser.R
```

It simulates two datasets, `SimDataA` and `SimDataB`, each with `OS`, `RFS`,
`PFS`, and `DFS` endpoints. The script plants `GENE1-10` as risk genes and
`GENE11-20` as protective genes, then returns `HRGS` and `LRGS` signatures:

```r
datasets <- c(
  "SimDataA_OS", "SimDataA_RFS", "SimDataA_PFS", "SimDataA_DFS",
  "SimDataB_OS", "SimDataB_RFS", "SimDataB_PFS", "SimDataB_DFS"
)
```

## Real data support

The skill is validated against two public colorectal cancer cohorts:

| Dataset | Samples | Endpoints |
| --- | --- | --- |
| GSE39582 | 566 | OS, RFS |
| GSE17536 | 177 | OS, DFS, DSS |

Download and format them with:

```bash
Rscript scripts/download_real_data.R
```

This writes Prognoser-compatible `.Rdata` objects such as `GSE39582_OS`,
`GSE39582_OS_survival`, and `GSE17536_OS`. Load them into the global
environment, then run `Prognoser()` on both cohorts.

## What the skill enforces

### Exact data contract

`Prognoser()` uses `get()` on the R global environment. Every dataset-endpoint
must have two objects:

```r
TCGA_COAD_OS            # genes x samples expression data frame
TCGA_COAD_OS_survival   # time + status, row names = expression columns
```

Dataset-endpoint entries sharing the same prefix are grouped and
meta-analyzed:

```r
datasets <- c(
  "TCGA_COAD_OS",
  "TCGA_COAD_RFS",
  "TCGA_COAD_PFS",
  "TCGA_COAD_DFS"
)
```

### Reproducible execution

```r
library(Prognoser)

result <- Prognoser(
  expression_accession_vector = datasets,
  knowledge_genesets = knowledge_genesets,
  denovo_threshold = 0.05,
  gene_threshold = 0.001,
  min_genes = 5,
  similarity_threshold = 0.5,
  de_redundant_basis = "agglomeration",
  linkage = "ward.D"
)
```

### Validated output

- `result$genesets`: non-redundant `HRGS*` (risk) and `LRGS*` (protective)
  signatures.
- `result$combined_matrix_protection`: signed effect matrix for protective
  genes.
- `result$combined_matrix_risk`: signed effect matrix for risk genes.
- `result$network_protection` / `result$network_risk`: Jaccard edge lists, only
  with `de_redundant_basis = "igraph_component"`.

## Guardrails

- Never fabricate survival time or event status.
- Drop unmatched expression/survival samples and report the count.
- `linkage` defaults to `"ward.D"`; pass exactly one value, never the whole
  vector of choices.
- Keep multi-cohort expression matrices gene-aligned before meta-analysis.

## Related projects

- [Prognoser R package](https://github.com/FangZY-Lab/Prognoser)
- [Prognoser-Web](https://github.com/FangZY-Lab/Prognoser-Web)

## License

MIT.
