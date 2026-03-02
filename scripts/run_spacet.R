#!/usr/bin/env Rscript
# ============================================================================
# Script to run SpaCET deconvolution on BayesSpace enhanced data
# ============================================================================

suppressPackageStartupMessages({
  library(SpaCET)
  library(SingleCellExperiment)
  library(jsonlite)
  library(argparse)
  source("../bayespacet_integration.R")
})

# Parse command line arguments
parser <- ArgumentParser(description = "Run SpaCET deconvolution")
parser$add_argument("--sample", required = TRUE, help = "Sample ID")
parser$add_argument("--sce_file", required = TRUE, help = "BayesSpace enhanced SCE file")
parser$add_argument("--visium_path", required = TRUE, help = "Path to Visium output folder")
parser$add_argument("--output", required = TRUE, help = "Output RDS file")
parser$add_argument("--cancer_type", required = TRUE, help = "Cancer type (e.g., KIRC, LUAD)")
parser$add_argument("--cores", type = "integer", default = 8, help = "Number of cores")

args <- parser$parse_args()

# Load BayesSpace enhanced object
cat("Loading BayesSpace enhanced object for", args$sample, "...\n")
sce_enhanced <- readRDS(args$sce_file)

# Create SpaCET object from BayesSpace data
cat("Creating SpaCET object from BayesSpace data...\n")
spacet_obj <- create_spacet_from_bayespace(
  sce = sce_enhanced,
  visium_path = args$visium_path
)

# Run deconvolution
cat("Running SpaCET deconvolution (cancer type:", args$cancer_type, ")...\n")
cat("This may take a while...\n")
spacet_obj <- SpaCET.deconvolution(
  spacet_obj,
  cancerType = args$cancer_type,
  coreNo = args$cores
)

# Save
cat("Saving SpaCET object to", args$output, "...\n")
saveRDS(spacet_obj, args$output)
cat("✓ SpaCET deconvolution completed!\n")
