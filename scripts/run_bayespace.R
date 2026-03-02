#!/usr/bin/env Rscript
# ============================================================================
# Script to run BayesSpace spatial enhancement
# ============================================================================

suppressPackageStartupMessages({
  library(BayesSpace)
  library(SingleCellExperiment)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Run BayesSpace spatial enhancement")
parser$add_argument("--sample", required = TRUE, help = "Sample ID")
parser$add_argument("--visium_path", required = TRUE, help = "Path to Visium output folder")
parser$add_argument("--output", required = TRUE, help = "Output RDS file")
parser$add_argument("--n_pcs", type = "integer", default = 7, help = "Number of PCs")
parser$add_argument("--n_hvgs", type = "integer", default = 2000, help = "Number of HVGs")
parser$add_argument("--q", type = "integer", default = 4, help = "Number of clusters")
parser$add_argument("--d", type = "integer", default = 7, help = "Dimensions for clustering")
parser$add_argument("--seed", type = "integer", default = 102, help = "Random seed")

args <- parser$parse_args()

# Load Visium data
cat("Loading Visium data for", args$sample, "...\n")
# Use SpaCET function to read Visium data, then convert to SCE
spacet_temp <- create.SpaCET.object.10X(visiumPath = args$visium_path)
sce <- SingleCellExperiment(
  assays = list(counts = spacet_temp@input$counts),
  colData = DataFrame(spacet_temp@input$spotCoordinates)
)
# Add required metadata for BayesSpace
metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "Visium"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
sce <- sce[, colSums(counts(sce)) > 0]

# Preprocess
cat("Preprocessing with BayesSpace...\n")
set.seed(args$seed)
sce <- spatialPreprocess(
  sce,
  platform = "Visium",
  n.PCs = args$n_pcs,
  n.HVGs = args$n_hvgs,
  log.normalize = TRUE
)

# Tune q parameter
cat("Tuning q parameter...\n")
sce <- qTune(sce, qs = seq(2, 10), platform = "Visium", d = args$d)

# Spatial clustering
cat("Running spatial clustering (q =", args$q, ")...\n")
set.seed(args$seed + 1)
sce <- spatialCluster(
  sce,
  q = args$q,
  platform = "Visium",
  d = args$d,
  init.method = "mclust",
  model = "t",
  gamma = 2,
  nrep = 10000,
  burn.in = 100,
  save.chain = TRUE
)

# Spatial enhancement
cat("Running spatial enhancement...\n")
set.seed(args$seed + 2)
sce_enhanced <- spatialEnhance(
  sce,
  q = args$q,
  platform = "Visium",
  d = args$d,
  model = "t",
  gamma = 2,
  jitter_prior = 0.3,
  jitter_scale = 3.5,
  nrep = 100000,
  burn.in = 100,
  save.chain = TRUE
)

# Enhance all features
cat("Enhancing all features...\n")
markers <- rownames(sce_enhanced)
sce_enhanced_all <- enhanceFeatures(
  sce_enhanced,
  sce,
  feature_names = markers,
  nrounds = 0
)

# Convert logcounts to dgCMatrix
logcounts(sce_enhanced_all) <- as(logcounts(sce_enhanced_all, withDimnames = FALSE), "dgCMatrix")

# Save
cat("Saving enhanced object to", args$output, "...\n")
saveRDS(sce_enhanced_all, args$output)
cat("✓ BayesSpace enhancement completed!\n")
