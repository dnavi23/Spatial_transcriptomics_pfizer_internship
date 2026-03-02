# ============================================================================
# Snakemake Pipeline for TLS Signature Evaluation
# ============================================================================
# 
# This pipeline orchestrates the complete workflow:
# 1. Data preparation
# 2. BayesSpace spatial enhancement (optional)
# 3. SpaCET deconvolution
# 4. TLS signature evaluation
# 5. Results generation
#
# Usage: snakemake --cores 8
# ============================================================================

import os
import subprocess

# Configuration
configfile: "config.yaml"

# Sample information
SAMPLES = config["samples"]
CANCER_TYPE = config["cancer_type"]  # e.g., "KIRC", "LUAD"
N_CORES = config.get("cores", 8)

# Directories
DATA_DIR = "data"
RESULTS_DIR = "results"
FIGS_DIR = os.path.join(RESULTS_DIR, "figures")

# Create output directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGS_DIR, exist_ok=True)

# ============================================================================
# Rule: All
# ============================================================================

rule all:
    input:
        # Final outputs
        html_report = os.path.join(RESULTS_DIR, "TLS_signature_evaluation.html"),
        seurat_obj = os.path.join(RESULTS_DIR, "{sample}_bayespacet_seurat.rds"),
        performance_csv = os.path.join(RESULTS_DIR, "{sample}_signature_performance.csv")
    params:
        samples = SAMPLES

# ============================================================================
# Rule: Download Visium Data (if needed)
# ============================================================================

rule download_visium_data:
    """
    Download 10x Visium data from GEO or prepare from local source.
    This is a placeholder - users should download data manually or use
    their own download script.
    """
    output:
        visium_dir = directory(os.path.join(DATA_DIR, "{sample}_outs"))
    shell:
        """
        echo "Please download Visium data for {wildcards.sample} to {output.visium_dir}"
        echo "Or update this rule with your download method"
        mkdir -p {output.visium_dir}
        touch {output.visium_dir}/.placeholder
        """

# ============================================================================
# Rule: Run BayesSpace Enhancement
# ============================================================================

rule run_bayespace:
    """
    Run BayesSpace spatial enhancement on Visium data.
    Creates enhanced SingleCellExperiment object with subspot resolution.
    """
    input:
        visium_dir = os.path.join(DATA_DIR, "{sample}_outs")
    output:
        enhanced_sce = os.path.join(RESULTS_DIR, "{sample}_bayespace_enhanced.rds")
    params:
        n_pcs = config.get("bayespace", {}).get("n_pcs", 7),
        n_hvgs = config.get("bayespace", {}).get("n_hvgs", 2000),
        q = config.get("bayespace", {}).get("q", 4),
        d = config.get("bayespace", {}).get("d", 7)
    script:
        "scripts/run_bayespace.R"
    log:
        os.path.join(RESULTS_DIR, "logs", "{sample}_bayespace.log")
    shell:
        """
        mkdir -p {log:h}
        Rscript scripts/run_bayespace.R \
            --sample {wildcards.sample} \
            --visium_path {input.visium_dir} \
            --output {output.enhanced_sce} \
            --n_pcs {params.n_pcs} \
            --n_hvgs {params.n_hvgs} \
            --q {params.q} \
            --d {params.d} \
            > {log} 2>&1
        """

# ============================================================================
# Rule: Run SpaCET Deconvolution
# ============================================================================

rule run_spacet:
    """
    Run SpaCET deconvolution on BayesSpace enhanced data.
    Can use pre-computed BayesSpace object or create from enhanced SCE.
    """
    input:
        enhanced_sce = os.path.join(RESULTS_DIR, "{sample}_bayespace_enhanced.rds"),
        visium_dir = os.path.join(DATA_DIR, "{sample}_outs")
    output:
        spacet_obj = os.path.join(RESULTS_DIR, "{sample}_bayespaCET.rds")
    params:
        cancer_type = CANCER_TYPE,
        cores = N_CORES
    script:
        "scripts/run_spacet.R"
    log:
        os.path.join(RESULTS_DIR, "logs", "{sample}_spacet.log")
    shell:
        """
        mkdir -p {log:h}
        Rscript scripts/run_spacet.R \
            --sample {wildcards.sample} \
            --sce_file {input.enhanced_sce} \
            --visium_path {input.visium_dir} \
            --output {output.spacet_obj} \
            --cancer_type {params.cancer_type} \
            --cores {params.cores} \
            > {log} 2>&1
        """

# ============================================================================
# Rule: Run TLS Signature Evaluation
# ============================================================================

rule run_tls_evaluation:
    """
    Run complete TLS signature evaluation workflow.
    Creates integrated bayespacet-seurat object and evaluates all signatures.
    """
    input:
        enhanced_sce = os.path.join(RESULTS_DIR, "{sample}_bayespace_enhanced.rds"),
        spacet_obj = os.path.join(RESULTS_DIR, "{sample}_bayespaCET.rds"),
        visium_dir = os.path.join(DATA_DIR, "{sample}_outs")
    output:
        html_report = os.path.join(RESULTS_DIR, "{sample}_TLS_signature_evaluation.html"),
        seurat_obj = os.path.join(RESULTS_DIR, "{sample}_bayespacet_seurat.rds"),
        performance_csv = os.path.join(RESULTS_DIR, "{sample}_signature_performance.csv"),
        metadata_csv = os.path.join(RESULTS_DIR, "{sample}_metadata_with_scores.csv")
    params:
        sample = "{sample}",
        cancer_type = CANCER_TYPE
    script:
        "scripts/run_tls_evaluation.R"
    log:
        os.path.join(RESULTS_DIR, "logs", "{sample}_tls_evaluation.log")
    shell:
        """
        mkdir -p {log:h}
        Rscript scripts/run_tls_evaluation.R \
            --sample {params.sample} \
            --visium_path {input.visium_dir} \
            --sce_file {input.enhanced_sce} \
            --spacet_file {input.spacet_obj} \
            --output_dir {RESULTS_DIR} \
            > {log} 2>&1
        """

# ============================================================================
# Rule: Generate Summary Report
# ============================================================================

rule generate_summary:
    """
    Generate summary report comparing all samples.
    """
    input:
        performance_files = expand(
            os.path.join(RESULTS_DIR, "{sample}_signature_performance.csv"),
            sample = SAMPLES
        )
    output:
        summary_report = os.path.join(RESULTS_DIR, "summary_report.html")
    script:
        "scripts/generate_summary.R"
    log:
        os.path.join(RESULTS_DIR, "logs", "summary.log")
