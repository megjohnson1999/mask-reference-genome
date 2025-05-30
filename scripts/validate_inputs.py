#!/usr/bin/env python3
"""
Input validation script for Reference Masker
Validates configuration and input files before running the workflow
"""

import os
import sys
import yaml
import argparse
from pathlib import Path
from Bio import SeqIO
import logging

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('validation.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def validate_config_file(config_path):
    """Validate the configuration YAML file"""
    logger = logging.getLogger(__name__)
    errors = []
    warnings = []
    
    logger.info(f"Validating configuration file: {config_path}")
    
    if not os.path.exists(config_path):
        errors.append(f"Configuration file not found: {config_path}")
        return errors, warnings
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        errors.append(f"Invalid YAML syntax in config file: {e}")
        return errors, warnings
    
    # Required sections
    required_sections = ['viral_database', 'reference_genomes', 'output_dir', 'temp_dir', 'log_dir']
    for section in required_sections:
        if section not in config:
            errors.append(f"Missing required configuration section: {section}")
    
    # Required parameters
    required_params = {
        'shred': ['overlap', 'median_length', 'variance', 'min_length'],
        'mapping': ['min_identity', 'max_indel'],
        'masking': ['entropy_threshold']
    }
    
    for section, params in required_params.items():
        if section not in config:
            errors.append(f"Missing required configuration section: {section}")
        else:
            for param in params:
                if param not in config[section]:
                    errors.append(f"Missing required parameter: {section}.{param}")
    
    # Validate parameter ranges
    if 'mapping' in config:
        if 'min_identity' in config['mapping']:
            identity = config['mapping']['min_identity']
            if not (0 <= identity <= 100):
                errors.append(f"min_identity must be between 0-100, got: {identity}")
            elif identity < 80:
                warnings.append(f"min_identity is quite low ({identity}%), may produce many false positives")
        
        if 'max_indel' in config['mapping']:
            indel = config['mapping']['max_indel']
            if not (0 <= indel <= 20):
                errors.append(f"max_indel must be between 0-20, got: {indel}")
    
    if 'shred' in config:
        if 'overlap' in config['shred'] and 'median_length' in config['shred']:
            overlap = config['shred']['overlap']
            median = config['shred']['median_length']
            if overlap >= median:
                errors.append(f"shred overlap ({overlap}) must be less than median_length ({median})")
    
    if 'threads' in config:
        threads = config['threads']
        if not isinstance(threads, int) or threads < 1:
            errors.append(f"threads must be a positive integer, got: {threads}")
        elif threads > 32:
            warnings.append(f"threads is quite high ({threads}), ensure your system can handle this")
    
    logger.info(f"Configuration validation complete. Errors: {len(errors)}, Warnings: {len(warnings)}")
    return errors, warnings

def validate_fasta_file(fasta_path, file_type="FASTA"):
    """Validate a FASTA file"""
    logger = logging.getLogger(__name__)
    errors = []
    warnings = []
    
    logger.info(f"Validating {file_type} file: {fasta_path}")
    
    if not os.path.exists(fasta_path):
        errors.append(f"{file_type} file not found: {fasta_path}")
        return errors, warnings
    
    if os.path.getsize(fasta_path) == 0:
        errors.append(f"{file_type} file is empty: {fasta_path}")
        return errors, warnings
    
    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        
        if len(records) == 0:
            errors.append(f"No sequences found in {file_type} file: {fasta_path}")
            return errors, warnings
        
        total_length = sum(len(record.seq) for record in records)
        
        logger.info(f"{file_type} file validation:")
        logger.info(f"  - Sequences: {len(records):,}")
        logger.info(f"  - Total length: {total_length:,} bp")
        logger.info(f"  - Average length: {total_length//len(records):,} bp")
        
        # Check for common issues
        if file_type == "Reference genome":
            if len(records) > 1000:
                warnings.append(f"Reference has many contigs ({len(records)}), consider using a more contiguous assembly")
            
            if total_length < 1000000:  # < 1Mb
                warnings.append(f"Reference genome is quite small ({total_length:,} bp)")
        
        elif file_type == "Viral database":
            if len(records) < 100:
                warnings.append(f"Viral database has few sequences ({len(records)}), may miss viral regions")
            
            if total_length < 10000000:  # < 10Mb
                warnings.append(f"Viral database is quite small ({total_length:,} bp)")
        
        # Check for unusual characters
        valid_chars = set('ATCGRYSWKMBDHVN-')
        for i, record in enumerate(records[:10]):  # Check first 10 sequences
            seq_chars = set(str(record.seq).upper())
            invalid_chars = seq_chars - valid_chars
            if invalid_chars:
                warnings.append(f"Sequence {record.id} contains unusual characters: {invalid_chars}")
                break
    
    except Exception as e:
        errors.append(f"Error parsing {file_type} file: {e}")
    
    return errors, warnings

def validate_directories(config):
    """Validate output directories can be created"""
    logger = logging.getLogger(__name__)
    errors = []
    warnings = []
    
    directories = [config.get('output_dir'), config.get('temp_dir'), config.get('log_dir')]
    
    for dir_path in directories:
        if dir_path:
            try:
                os.makedirs(dir_path, exist_ok=True)
                # Test write permissions
                test_file = os.path.join(dir_path, '.test_write')
                with open(test_file, 'w') as f:
                    f.write('test')
                os.remove(test_file)
                logger.info(f"Directory validated: {dir_path}")
            except PermissionError:
                errors.append(f"No write permission for directory: {dir_path}")
            except Exception as e:
                errors.append(f"Cannot create/access directory {dir_path}: {e}")
    
    return errors, warnings

def validate_dependencies():
    """Check for required software dependencies"""
    logger = logging.getLogger(__name__)
    errors = []
    warnings = []
    
    # Check for conda/mamba
    import shutil
    
    required_tools = ['snakemake']
    optional_tools = ['mamba', 'conda']
    
    for tool in required_tools:
        if not shutil.which(tool):
            errors.append(f"Required tool not found: {tool}")
        else:
            logger.info(f"Found required tool: {tool}")
    
    for tool in optional_tools:
        if shutil.which(tool):
            logger.info(f"Found optional tool: {tool}")
            break
    else:
        warnings.append("Neither conda nor mamba found - workflow may fail")
    
    return errors, warnings

def generate_validation_report(all_errors, all_warnings, output_file=None):
    """Generate a validation report"""
    logger = logging.getLogger(__name__)
    
    report_lines = []
    report_lines.append("=== REFERENCE MASKER VALIDATION REPORT ===")
    report_lines.append(f"Generated: {pd.Timestamp.now()}")
    report_lines.append("")
    
    if all_errors:
        report_lines.append("ERRORS (must fix before running):")
        for error in all_errors:
            report_lines.append(f"  ❌ {error}")
        report_lines.append("")
    
    if all_warnings:
        report_lines.append("WARNINGS (recommended to address):")
        for warning in all_warnings:
            report_lines.append(f"  ⚠️  {warning}")
        report_lines.append("")
    
    if not all_errors and not all_warnings:
        report_lines.append("✅ All validations passed!")
    elif not all_errors:
        report_lines.append("✅ No critical errors found - workflow can proceed")
    else:
        report_lines.append("❌ Critical errors found - please fix before running workflow")
    
    report = "\n".join(report_lines)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        logger.info(f"Validation report written to: {output_file}")
    
    print("\n" + report)
    
    return len(all_errors) == 0

def main():
    parser = argparse.ArgumentParser(description='Validate Reference Masker inputs')
    parser.add_argument('--config', default='config.yaml', help='Configuration file path')
    parser.add_argument('--report', help='Output validation report file')
    parser.add_argument('--strict', action='store_true', help='Treat warnings as errors')
    
    args = parser.parse_args()
    
    logger = setup_logging()
    logger.info("Starting Reference Masker input validation")
    
    all_errors = []
    all_warnings = []
    
    # Validate configuration file
    errors, warnings = validate_config_file(args.config)
    all_errors.extend(errors)
    all_warnings.extend(warnings)
    
    if not errors:  # Only continue if config is valid
        try:
            with open(args.config, 'r') as f:
                config = yaml.safe_load(f)
            
            # Validate viral database
            if 'viral_database' in config:
                errors, warnings = validate_fasta_file(config['viral_database'], "Viral database")
                all_errors.extend(errors)
                all_warnings.extend(warnings)
            
            # Validate reference genomes
            if 'reference_genomes' in config:
                for ref_genome in config['reference_genomes']:
                    errors, warnings = validate_fasta_file(ref_genome, "Reference genome")
                    all_errors.extend(errors)
                    all_warnings.extend(warnings)
            
            # Validate directories
            errors, warnings = validate_directories(config)
            all_errors.extend(errors)
            all_warnings.extend(warnings)
            
        except Exception as e:
            all_errors.append(f"Error loading configuration: {e}")
    
    # Validate dependencies
    errors, warnings = validate_dependencies()
    all_errors.extend(errors)
    all_warnings.extend(warnings)
    
    # Generate report
    success = generate_validation_report(all_errors, all_warnings, args.report)
    
    if args.strict and all_warnings:
        logger.error("Strict mode: treating warnings as errors")
        success = False
    
    if success:
        logger.info("Validation completed successfully")
        sys.exit(0)
    else:
        logger.error("Validation failed")
        sys.exit(1)

if __name__ == "__main__":
    # Import pandas here to avoid import error in report generation
    import pandas as pd
    main()