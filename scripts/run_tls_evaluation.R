#!/usr/bin/env Rscript
# ============================================================================
# Script to run complete TLS signature evaluation
# ============================================================================

suppressPackageStartupMessages({
  library(rmarkdown)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Run TLS signature evaluation")
parser$add_argument("--sample", required = TRUE, help = "Sample ID")
parser$add_argument("--visium_path", required = TRUE, help = "Path to Visium output folder")
parser$add_argument("--sce_file", required = TRUE, help = "BayesSpace enhanced SCE file")
parser$add_argument("--spacet_file", required = TRUE, help = "SpaCET object file")
parser$add_argument("--output_dir", required = TRUE, help = "Output directory")

args <- parser$parse_args()

# Set environment variables for R Markdown
Sys.setenv(SAMPLE = args$sample)
Sys.setenv(VISIUM_PATH = args$visium_path)
Sys.setenv(SCE_FILE = args$sce_file)
Sys.setenv(SPACET_FILE = args$spacet_file)
Sys.setenv(OUTPUT_DIR = args$output_dir)

# Render R Markdown
cat("Rendering TLS signature evaluation report for", args$sample, "...\n")
rmarkdown::render(
  input = "../TLS_signature_evaluation.Rmd",
  output_file = file.path(args$output_dir, paste0(args$sample, "_TLS_signature_evaluation.html")),
  output_dir = args$output_dir,
  params = list(
    sample = args$sample,
    visium_path = args$visium_path,
    sce_file = args$sce_file,
    spacet_file = args$spacet_file
  ),
  quiet = FALSE
)

cat("✓ TLS signature evaluation completed!\n")
