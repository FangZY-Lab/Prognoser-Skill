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
6. **Downstream cell-state characterisation** of the resulting signature.
7. **Downstream tumour-microenvironment analysis** of the same cohorts.

## Repository layout

```text
.
├── SKILL.md                         Skill entrypoint
├── README.md                        This document
├── references/
│   ├── data-acquisition.md          Public data download and endpoints
│   ├── preprocessing.md             Normalization and identifier alignment
│   ├── prognoser-api.md             Full API and output schema
│   ├── real-data.md                 Validated real colorectal cohorts
│   ├── cell-state-analysis.md       Downstream: cell-state panels and scoring
│   ├── tme-analysis.md              Downstream: ESTIMATE, xCell, markers
│   └── downstream-statistics.md     Shared scoring and statistics contract
└── scripts/
    ├── run_prognoser.R              Reproducible execution template
    ├── example_prognoser.R          Self-contained two-cohort demo
    ├── download_real_data.R         Download two real CRC cohorts
    ├── downstream_lib.R             Shared downstream helpers
    ├── downstream_signatures.R      Panel registry (cell state + TME markers)
    ├── downstream_cell_states.R     Downstream module 1: cell states
    ├── downstream_tme.R             Downstream module 2: microenvironment
    ├── run_downstream.R             Both modules plus a Markdown report
    └── example_downstream.R         Self-contained downstream demo
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

## Downstream analysis

Once `HRGS*` and `LRGS*` signatures exist, two downstream modules characterise
them on the same cohorts. Both accept the result object from `Prognoser()` and
the same `<Dataset>_<Endpoint>` inputs, so nothing is downloaded twice.

```bash
# cell states and the tumour microenvironment, plus a Markdown report
Rscript scripts/run_downstream.R \
  --input=output/objects \
  --result=output/prognoser_result.rds \
  --output=output/downstream
```

**Cell states** ([references/cell-state-analysis.md](references/cell-state-analysis.md)).
Thirteen panels — EMT, proliferation, hypoxia, glycolysis, OXPHOS,
angiogenesis, apoptosis, DNA repair, senescence, TGF-beta, IFN-gamma,
inflammation and a stemness proxy — are resolved from MSigDB at run time via
`msigdbr` (no hard-coded gene lists), scored with ssGSEA, then tested against
the signature score and against survival.

**Tumour microenvironment** ([references/tme-analysis.md](references/tme-analysis.md)).
ESTIMATE for stromal/immune content, xCell for 64 cell types, optional quanTIseq
fractions for ten immune populations, and nine curated marker panels (CD8, Treg,
M1/M2, MDSC, checkpoints, TLS, CAF, endothelial). An
`inflamed`/`excluded`/`desert` phenotype call is emitted together with the
thresholds that produced it.

Both modules follow the scoring and multiple-testing contract in
[references/downstream-statistics.md](references/downstream-statistics.md):
BH within each output family, Cox per standard deviation, and explicit
guardrails for small cohorts.

### Dependencies

Signature discovery needs only CRAN packages (see
[references/prognoser-api.md](references/prognoser-api.md)). Downstream adds:

| Tier | Packages | Needed for |
| --- | --- | --- |
| Required | `msigdbr`, `GSVA`, `survival` | Panel resolution, ssGSEA, survival models |
| Recommended | `estimate`, `xCell` | Bulk TME scores and cell types |
| Optional | `quantiseqr` | Immune fractions (TPM input only) |
| Optional | `pheatmap`, `ggplot2` | Figures; analyses run without them |

```r
install.packages(c("msigdbr", "estimate"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GSVA", "xCell", "quantiseqr"))
```

Missing optional packages degrade one view instead of failing the run.

### Verifying an installation

```bash
Rscript scripts/example_prognoser.R     # signature discovery
Rscript scripts/example_downstream.R    # downstream, on simulated cohorts
```

The downstream demo plants a risk axis, an immune axis and a stromal axis, then
recovers them: the signature score should rise with proliferation and fall with
CD8 content. If it does not, the environment or the input scale is wrong.

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
  na_ratio = 0.3,
  gene_threshold = 0.001,
  min_genes = 5,
  similarity_threshold = 0.5,
  linkage = "ward.D"
)
```

### Validated output

- `result$genesets`: non-redundant `HRGS*` (risk) and `LRGS*` (protective)
  signatures.
- `result$combined_matrix_protection`: signed effect matrix for protective
  genes.
- `result$combined_matrix_risk`: signed effect matrix for risk genes.

## Guardrails

- Never fabricate survival time or event status.
- Drop unmatched expression/survival samples and report the count.
- Use one `linkage` value; the package default is a vector and errors if left
  unchanged.
- Keep multi-cohort expression matrices gene-aligned before meta-analysis.

## Related projects

- [Prognoser R package](https://github.com/FangZY-Lab/Prognoser)
- [Prognoser-Web](https://github.com/FangZY-Lab/Prognoser-Web)

## License

MIT.
