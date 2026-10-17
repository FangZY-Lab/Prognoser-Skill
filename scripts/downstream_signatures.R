#!/usr/bin/env Rscript

# Gene set registry for the Prognoser downstream modules.
#
# Two families are defined here:
#
#   * Cell-state panels, resolved from MSigDB at run time through msigdbr. No
#     gene list is hard-coded, so the panels stay in sync with the MSigDB
#     release that ships with the installed msigdbr version. Provenance is
#     written next to every score in the output tables.
#
#   * Tumour-microenvironment marker panels, which are short curated marker
#     lists. These are marker-based approximations, not deconvolution results;
#     they are kept small and interpretable on purpose.
#
# Both families share one interface: dx_panel_genes() returns named character
# vectors, dx_panel_table() describes them for documentation output.

## ------------------------------------------------------- msigdbr shim ---

dx_require_msigdbr <- function() {
  dx_require("msigdbr", "install.packages('msigdbr')")
}

# msigdbr renamed category/subcategory to collection/subcollection in v10.
# Support both so the skill does not break on older Bioconductor stacks.
dx_msigdbr <- function(species, collection, subcollection = NULL) {
  dx_require_msigdbr()
  fmls <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% fmls) {
    msigdbr::msigdbr(species = species, collection = collection,
                     subcollection = subcollection)
  } else {
    msigdbr::msigdbr(species = species, category = collection,
                     subcategory = subcollection)
  }
}

.dx_msig_cache <- new.env(parent = emptyenv())

dx_msig_cache_key <- function(species, collection, subcollection) {
  paste(species, collection, subcollection %||% "NA", sep = "|")
}

dx_msig_collection_of <- function(set_name) {
  if (startsWith(set_name, "HALLMARK_")) return(c("H", NA_character_))
  if (startsWith(set_name, "GOBP_")) return(c("C5", "GO:BP"))
  if (startsWith(set_name, "GOCC_")) return(c("C5", "GO:CC"))
  if (startsWith(set_name, "GOMF_")) return(c("C5", "GO:MF"))
  c("C2", "CGP")
}

dx_msig_table <- function(collection, subcollection, species = "Homo sapiens") {
  key <- dx_msig_cache_key(species, collection, subcollection)
  if (!is.null(.dx_msig_cache[[key]])) return(.dx_msig_cache[[key]])
  df <- dx_msigdbr(species, collection, subcollection)
  .dx_msig_cache[[key]] <- df
  df
}

# Resolve MSigDB set names to gene vectors. Missing sets raise an error that
# names them, because silently scoring a subset of a panel changes the biology.
dx_msig_sets <- function(set_names, species = "Homo sapiens") {
  set_names <- unique(set_names)
  groups <- split(set_names, vapply(set_names,
                                    function(s) paste(dx_msig_collection_of(s)[1L],
                                                      dx_msig_collection_of(s)[2L],
                                                      sep = "|"),
                                    character(1L)))
  out <- list()
  for (grp in names(groups)) {
    parts <- strsplit(grp, "|", fixed = TRUE)[[1L]]
    sub <- if (identical(parts[2L], "NA")) NULL else parts[2L]
    tab <- dx_msig_table(parts[1L], sub, species)
    for (s in groups[[grp]]) {
      genes <- tab$gene_symbol[tab$gs_name == s]
      if (length(genes) == 0L) {
        stop(sprintf("MSigDB set not found: %s (collection %s).", s, grp),
             call. = FALSE)
      }
      out[[s]] <- unique(genes)
    }
  }
  out[set_names]
}

## --------------------------------------------------------- registry -----

# Cell-state panels. `sets` are MSigDB names; multiple sets are unioned into one
# panel. `use` states what the panel is and is not good for.
DX_STATE_PANELS <- list(
  emt = list(
    label = "Epithelial-mesenchymal transition",
    sets = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    use = "Dedifferentiation / mesenchymal shift; also tracks stromal admixture."
  ),
  proliferation = list(
    label = "Proliferation",
    sets = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
             "HALLMARK_MYC_TARGETS_V1"),
    use = "Cycling fraction; strongly confounded by tumour purity."
  ),
  hypoxia = list(
    label = "Hypoxia",
    sets = "HALLMARK_HYPOXIA",
    use = "Hypoxic adaptation; correlates with necrosis and poor perfusion."
  ),
  glycolysis = list(
    label = "Glycolysis",
    sets = "HALLMARK_GLYCOLYSIS",
    use = "Warburg-type metabolic state."
  ),
  oxphos = list(
    label = "Oxidative phosphorylation",
    sets = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    use = "Mitochondrial metabolism; often inverse to glycolysis."
  ),
  angiogenesis = list(
    label = "Angiogenesis",
    sets = "HALLMARK_ANGIOGENESIS",
    use = "Vascular/VEGF programme; overlaps endothelial content."
  ),
  apoptosis = list(
    label = "Apoptosis",
    sets = "HALLMARK_APOPTOSIS",
    use = "Cell-death competence."
  ),
  dna_repair = list(
    label = "DNA repair",
    sets = "HALLMARK_DNA_REPAIR",
    use = "Repair capacity; relevant to genomic instability and therapy response."
  ),
  senescence = list(
    label = "Cellular senescence",
    sets = "GOBP_CELLULAR_SENESCENCE",
    use = "Senescence programme; note that SASP genes overlap inflammation."
  ),
  tgf_beta = list(
    label = "TGF-beta signalling",
    sets = "HALLMARK_TGF_BETA_SIGNALING",
    use = "Immune exclusion and stromal activation driver."
  ),
  ifn_gamma = list(
    label = "Interferon-gamma response",
    sets = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    use = "T-cell inflamed / antigen-presentation programme."
  ),
  inflammation = list(
    label = "Inflammatory response",
    sets = "HALLMARK_INFLAMMATORY_RESPONSE",
    use = "Innate inflammation; partly overlaps immune infiltration."
  ),
  stemness = list(
    label = "Pluripotency / stemness programme",
    sets = "MUELLER_PLURINET",
    use = "Pluripotency-derived proxy; not a validated cancer stem-cell score."
  )
)

# TME marker panels. Short, interpretable, and explicitly not deconvolution.
DX_MARKER_PANELS <- list(
  cd8_t_cell = list(
    label = "CD8+ T cell / cytotoxic",
    genes = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "GNLY"),
    use = "Cytotoxic T-cell content."
  ),
  treg = list(
    label = "Regulatory T cell",
    genes = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
    use = "Treg content; overlaps activation markers."
  ),
  m1_macrophage = list(
    label = "M1-like macrophage",
    genes = c("NOS2", "IL12B", "CXCL9", "CXCL10", "CXCL11", "TNF"),
    use = "Classically activated macrophage programme."
  ),
  m2_macrophage = list(
    label = "M2-like macrophage",
    genes = c("CD163", "MRC1", "MSR1", "CCL22", "IL10", "TGFB1", "ARG1"),
    use = "Alternatively activated macrophage programme."
  ),
  mdsc = list(
    label = "Myeloid-derived suppressor cell",
    genes = c("S100A8", "S100A9", "ITGAM", "ARG1", "CXCR2"),
    use = "MDSC-associated genes; overlaps neutrophils."
  ),
  checkpoint = list(
    label = "Immune checkpoint expression",
    genes = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
              "TIGIT", "IDO1", "VSIR", "BTLA"),
    use = "Checkpoint transcripts; not a predictor of ICI response on its own."
  ),
  tls = list(
    label = "Tertiary lymphoid structure",
    genes = c("CXCL13", "LTB", "CCL19", "CCL21", "CCR7", "SELL", "CXCR5",
              "CD79B"),
    use = "TLS-associated genes; histology or imaging is the reference standard."
  ),
  caf = list(
    label = "Cancer-associated fibroblast / stroma",
    genes = c("FAP", "ACTA2", "PDGFRB", "COL1A1", "COL1A2", "COL3A1", "POSTN",
              "TAGLN", "THY1", "DCN"),
    use = "Stromal content; overlaps ESTIMATE StromalScore."
  ),
  endothelial = list(
    label = "Endothelial",
    genes = c("PECAM1", "VWF", "CDH5", "CLDN5", "EMCN", "KDR"),
    use = "Vascular content."
  )
)

dx_state_panel_ids <- function() names(DX_STATE_PANELS)
dx_marker_panel_ids <- function() names(DX_MARKER_PANELS)

dx_panel_table <- function() {
  states <- do.call(rbind, lapply(names(DX_STATE_PANELS), function(id) {
    p <- DX_STATE_PANELS[[id]]
    data.frame(family = "cell_state", panel = id, label = p$label,
               definition = paste(p$sets, collapse = "; "),
               source = "MSigDB via msigdbr", note = p$use,
               stringsAsFactors = FALSE)
  }))
  markers <- do.call(rbind, lapply(names(DX_MARKER_PANELS), function(id) {
    p <- DX_MARKER_PANELS[[id]]
    data.frame(family = "tme_marker", panel = id, label = p$label,
               definition = paste(p$genes, collapse = ", "),
               source = "curated marker list", note = p$use,
               stringsAsFactors = FALSE)
  }))
  rbind(states, markers)
}

# Resolve any mix of state and marker panel IDs into gene vectors. Unknown IDs
# are reported together instead of failing on the first one.
dx_panel_genes <- function(ids, species = "Homo sapiens") {
  if (length(ids) == 0L) return(list())
  known <- c(names(DX_STATE_PANELS), names(DX_MARKER_PANELS))
  unknown <- setdiff(ids, known)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown panel(s): %s\nAvailable: %s",
                 paste(unknown, collapse = ", "), paste(known, collapse = ", ")),
         call. = FALSE)
  }
  state_ids <- intersect(ids, names(DX_STATE_PANELS))
  marker_ids <- intersect(ids, names(DX_MARKER_PANELS))

  out <- list()
  if (length(state_ids) > 0L) {
    set_names <- unique(unlist(lapply(state_ids, function(id) DX_STATE_PANELS[[id]]$sets)))
    resolved <- dx_msig_sets(set_names, species = species)
    for (id in state_ids) {
      out[[id]] <- unique(unlist(resolved[DX_STATE_PANELS[[id]]$sets], use.names = FALSE))
    }
  }
  for (id in marker_ids) {
    out[[id]] <- unique(DX_MARKER_PANELS[[id]]$genes)
  }
  out[ids]
}
