"""
Reference Masker - Phase 4: Production-ready with error handling and validation
Masks virus-like regions in reference genomes using BBMap tools
"""

import os
import glob
import sys
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('workflow.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("reference_masker")

# Load configuration with error handling
try:
    configfile: "config.yaml"
    logger.info("Configuration loaded successfully")
except Exception as e:
    logger.error(f"Failed to load configuration: {e}")
    sys.exit(1)

# Validate configuration and extract reference genomes
try:
    REFERENCES = config["reference_genomes"]
    SAMPLES = [os.path.splitext(os.path.basename(ref))[0] for ref in REFERENCES]
    
    # Validate required files exist
    if not os.path.exists(config["viral_database"]):
        logger.error(f"Viral database not found: {config['viral_database']}")
        sys.exit(1)
    
    for ref in REFERENCES:
        if not os.path.exists(ref):
            logger.error(f"Reference genome not found: {ref}")
            sys.exit(1)
    
    logger.info(f"Processing {len(SAMPLES)} reference genome(s): {', '.join(SAMPLES)}")
    
except KeyError as e:
    logger.error(f"Missing required configuration key: {e}")
    sys.exit(1)
except Exception as e:
    logger.error(f"Configuration validation failed: {e}")
    sys.exit(1)

# Create output directories with error handling
try:
    for directory in [config["output_dir"], config["temp_dir"], config["log_dir"]]:
        os.makedirs(directory, exist_ok=True)
        logger.info(f"Created/verified directory: {directory}")
except Exception as e:
    logger.error(f"Failed to create output directories: {e}")
    sys.exit(1)

# Final output - process all reference genomes and generate reports
rule all:
    input:
        expand(os.path.join(config["output_dir"], "{sample}_masked.fasta"), sample=SAMPLES),
        expand(os.path.join(config["output_dir"], "{sample}_masking_stats.txt"), sample=SAMPLES),
        os.path.join(config["output_dir"], "reports", "masking_report.html"),
        expand(os.path.join(config["output_dir"], "reports", "{sample}_masked_regions.bed"), sample=SAMPLES)

# Rule 1: Shred viral sequences into fragments (shared across all references)
rule shred_viruses:
    input:
        config["viral_database"]
    output:
        os.path.join(config["temp_dir"], "viral_shreds.fasta")
    params:
        overlap = config["shred"]["overlap"],
        median = config["shred"]["median_length"],
        variance = config["shred"]["variance"],
        minlength = config["shred"]["min_length"]
    conda:
        "envs/bbmap.yaml"
    threads: config["threads"]
    log:
        os.path.join(config["log_dir"], "shred.log")
    shell:
        """
        set -euo pipefail
        
        echo "[$(date)] Starting viral sequence shredding" >> {log}
        echo "Input file: {input}" >> {log}
        echo "Output file: {output}" >> {log}
        echo "Parameters: minlength={params.minlength}, overlap={params.overlap}, median={params.median}, variance={params.variance}" >> {log}
        
        # Check input file exists and is not empty
        if [[ ! -f {input} ]]; then
            echo "ERROR: Input file not found: {input}" >> {log}
            exit 1
        fi
        
        if [[ ! -s {input} ]]; then
            echo "ERROR: Input file is empty: {input}" >> {log}
            exit 1
        fi
        
        # Run shredding with error checking
        if shred.sh \
            in={input} \
            out={output} \
            minlength={params.minlength} \
            overlap={params.overlap} \
            median={params.median} \
            variance={params.variance} \
            overwrite=true \
            2>> {log}; then
            echo "[$(date)] Shredding completed successfully" >> {log}
            
            # Verify output was created
            if [[ ! -f {output} ]]; then
                echo "ERROR: Output file was not created: {output}" >> {log}
                exit 1
            fi
            
            if [[ ! -s {output} ]]; then
                echo "ERROR: Output file is empty: {output}" >> {log}
                exit 1
            fi
            
            echo "Output file size: $(wc -l < {output}) lines" >> {log}
        else
            echo "ERROR: Shredding failed" >> {log}
            exit 1
        fi
        """

# Rule 2: Map shredded sequences to each reference genome
rule map_to_reference:
    input:
        shreds = os.path.join(config["temp_dir"], "viral_shreds.fasta"),
        ref = lambda wildcards: [ref for ref in config["reference_genomes"] if os.path.splitext(os.path.basename(ref))[0] == wildcards.sample][0]
    output:
        temp(os.path.join(config["temp_dir"], "{sample}_mapped.sam"))
    params:
        minid = config["mapping"]["min_identity"],
        maxindel = config["mapping"]["max_indel"]
    conda:
        "envs/bbmap.yaml"
    threads: config["threads"]
    log:
        os.path.join(config["log_dir"], "{sample}_map.log")
    shell:
        """
        bbmap.sh \
            ref={input.ref} \
            in={input.shreds} \
            outm={output} \
            minid={params.minid} \
            maxindel={params.maxindel} \
            ambig=all \
            threads={threads} \
            overwrite=true \
            2> {log}
        """

# Rule 3: Mask each reference genome based on alignments
rule mask_reference:
    input:
        sam = os.path.join(config["temp_dir"], "{sample}_mapped.sam"),
        ref = lambda wildcards: [ref for ref in config["reference_genomes"] if os.path.splitext(os.path.basename(ref))[0] == wildcards.sample][0]
    output:
        os.path.join(config["output_dir"], "{sample}_masked.fasta")
    params:
        entropy = config["masking"]["entropy_threshold"]
    conda:
        "envs/bbmap.yaml"
    threads: config["threads"]
    log:
        os.path.join(config["log_dir"], "{sample}_mask.log")
    shell:
        """
        bbmask.sh \
            in={input.ref} \
            out={output} \
            sam={input.sam} \
            entropy={params.entropy} \
            threads={threads} \
            overwrite=true \
            2> {log}
        """

# Rule to generate basic statistics about masking for each reference
rule masking_stats:
    input:
        original = lambda wildcards: [ref for ref in config["reference_genomes"] if os.path.splitext(os.path.basename(ref))[0] == wildcards.sample][0],
        masked = os.path.join(config["output_dir"], "{sample}_masked.fasta")
    output:
        os.path.join(config["output_dir"], "{sample}_masking_stats.txt")
    conda:
        "envs/bbmap.yaml"
    log:
        os.path.join(config["log_dir"], "{sample}_stats.log")
    shell:
        """
        echo "Masking Statistics for {wildcards.sample}" > {output}
        echo "========================================" >> {output}
        echo "" >> {output}
        echo "Original reference: {input.original}" >> {output}
        stats.sh in={input.original} format=3 >> {output} 2> {log}
        echo "" >> {output}
        echo "Masked reference: {input.masked}" >> {output}
        stats.sh in={input.masked} format=3 >> {output} 2>> {log}
        """

# Rule to generate comprehensive analysis report
rule generate_report:
    input:
        masked_files = expand(os.path.join(config["output_dir"], "{sample}_masked.fasta"), sample=SAMPLES),
        mask_logs = expand(os.path.join(config["log_dir"], "{sample}_mask.log"), sample=SAMPLES),
        map_logs = expand(os.path.join(config["log_dir"], "{sample}_map.log"), sample=SAMPLES)
    output:
        html_report = os.path.join(config["output_dir"], "reports", "masking_report.html"),
        summary_plot = os.path.join(config["output_dir"], "reports", "masking_summary.png"),
        bed_files = expand(os.path.join(config["output_dir"], "reports", "{sample}_masked_regions.bed"), sample=SAMPLES)
    params:
        samples = SAMPLES,
        script = "scripts/generate_report.py"
    conda:
        "envs/plotting.yaml"
    log:
        os.path.join(config["log_dir"], "generate_report.log")
    shell:
        """
        mkdir -p {config[output_dir]}/reports
        python {params.script} \
            --log-dir {config[log_dir]} \
            --results-dir {config[output_dir]} \
            --output-dir {config[output_dir]}/reports \
            --samples {params.samples} \
            2> {log}
        """