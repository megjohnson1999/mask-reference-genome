# Reference Masker

A robust Snakemake workflow for masking virus-like regions in reference genomes to prevent false-positive host removal during metagenomic analysis.

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

Reference Masker identifies and masks genomic regions with similarity to viral sequences, preventing the inadvertent removal of genuine viral reads during host decontamination of metagenomic samples. This is particularly important for gut microbiome studies where viral sequences should be retained for analysis.

### Key Features

- **Configuration-driven**: All parameters specified in YAML config file
- **Multi-reference support**: Process multiple reference genomes in parallel
- **Comprehensive reporting**: HTML reports, BED files, and visualizations
- **Production-ready**: Robust error handling, logging, and validation
- **Conservative masking**: Strict parameters to avoid over-masking host sequences

## Quick Start

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/reference-masker.git
cd reference-masker

# Install conda/mamba if not already available
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

# Create and activate snakemake environment
mamba create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

### 2. Configuration

Edit `config.yaml` to specify your input files and parameters:

```yaml
# Input files
viral_database: "path/to/viral_sequences.fasta"
reference_genomes:
  - "path/to/human_genome.fa"
  - "path/to/mouse_genome.fa"

# Output directories
output_dir: "results"
temp_dir: "temp"
log_dir: "logs"

# Parameters (recommended defaults)
shred:
  overlap: 40
  median_length: 80
  variance: 5
  min_length: 60

mapping:
  min_identity: 90    # Strict for confident matches
  max_indel: 2        # Conservative indel allowance

masking:
  entropy_threshold: 0.5

threads: 4
```

### 3. Validation (Recommended)

Validate your configuration and input files before running:

```bash
python scripts/validate_inputs.py --config config.yaml --report validation_report.txt
```

### 4. Run the Workflow

```bash
# Dry run to check the workflow
snakemake --use-conda -j 4 --dry-run

# Run the full workflow
snakemake --use-conda -j 4
```

### 5. View Results

- **Masked genomes**: `results/{sample}_masked.fasta`
- **Statistics**: `results/{sample}_masking_stats.txt`
- **HTML report**: `results/reports/masking_report.html`
- **BED files**: `results/reports/{sample}_masked_regions.bed`
- **Visualizations**: `results/reports/masking_summary.png`

## Workflow Steps

1. **Shred viruses**: Fragment viral sequences into overlapping pieces
2. **Map to reference**: Align viral fragments to reference genomes using BBMap
3. **Mask reference**: Replace aligned regions with 'N' characters using BBMask
4. **Generate reports**: Create comprehensive analysis reports and visualizations

## Input Files

### Viral Database
- **Format**: FASTA file containing viral sequences
- **Source**: NCBI viral genomes, custom viral database
- **Example**: `ncbi_viruses_20250530.fasta`

### Reference Genomes
- **Format**: FASTA files (can be compressed with gzip)
- **Examples**: Human genome (GRCh38), mouse genome (GRCm39)
- **Multiple genomes**: Specify as list in config.yaml

## Output Files

### Masked Genomes
- `{sample}_masked.fasta`: Reference genome with viral-like regions masked as 'N'
- Use these for host decontamination to preserve viral reads

### Statistics
- `{sample}_masking_stats.txt`: Detailed statistics about masking
- Includes total bases, percentage masked, low-complexity regions

### Reports
- `masking_report.html`: Interactive HTML summary with plots
- `masking_summary.png`: Publication-ready visualization plots
- `{sample}_masked_regions.bed`: Genomic coordinates of masked regions

### BED Files
- Compatible with genome browsers (IGV, UCSC Genome Browser)
- Can be used for downstream analysis of masked regions

## Parameters

### Shredding Parameters
- `overlap`: Base pair overlap between fragments (default: 40)
- `median_length`: Target fragment length (default: 80)
- `variance`: Length variation around median (default: 5)
- `min_length`: Minimum fragment length (default: 60)

### Mapping Parameters
- `min_identity`: Minimum sequence identity for alignment (default: 90%)
- `max_indel`: Maximum indels allowed in alignment (default: 2)

### Masking Parameters
- `entropy_threshold`: Low-complexity masking threshold (default: 0.5)

## Advanced Usage

### Multiple Reference Genomes

```yaml
reference_genomes:
  - "genomes/human_GRCh38.fa"
  - "genomes/mouse_GRCm39.fa"
  - "genomes/custom_host.fa"
```

Each genome will be processed independently with sample-specific outputs.

### Custom Viral Database

You can use your own viral database instead of NCBI:

```yaml
viral_database: "custom_viruses.fasta"
```

### Adjusting Stringency

For more sensitive detection (may increase false positives):
```yaml
mapping:
  min_identity: 80
  max_indel: 5
```

For more specific detection (may miss some viral regions):
```yaml
mapping:
  min_identity: 95
  max_indel: 1
```

## Troubleshooting

### Common Issues

1. **Memory errors**: Increase memory allocation in conda environment config
2. **Permission errors**: Ensure write access to output directories
3. **Empty output**: Check viral database format and reference genome paths

### Validation Errors

Run the validation script for detailed error diagnosis:
```bash
python scripts/validate_inputs.py --config config.yaml --strict
```

### Log Files

Check workflow logs for detailed information:
- `workflow.log`: Overall workflow logging
- `logs/{sample}_map.log`: BBMap alignment details
- `logs/{sample}_mask.log`: BBMask masking details

## Performance

### Runtime Estimates
- Human chromosome (~50MB): ~5-10 minutes
- Whole human genome (~3GB): ~2-4 hours
- Runtime scales with genome size and viral database size

### Resource Requirements
- **CPU**: 4+ cores recommended
- **Memory**: 4-8GB RAM for most genomes
- **Storage**: ~2x input genome size for temporary files

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support:
- Open an issue on GitHub
- Check the documentation and troubleshooting guide
- Review existing issues for solutions

## Changelog

### Version 1.0.0 (Phase 4)
- Production-ready release
- Comprehensive error handling and validation
- HTML reporting and visualizations
- Multi-reference genome support
- Configuration-driven workflow

### Version 0.3.0 (Phase 3)
- Added reporting and analysis features
- BED file generation
- Visualization plots

### Version 0.2.0 (Phase 2)
- Configuration-driven workflow
- Multi-reference support
- Improved modularity

### Version 0.1.0 (Phase 1)
- Initial working pipeline
- Basic masking functionality
