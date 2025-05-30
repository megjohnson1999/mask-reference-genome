#!/usr/bin/env python3
"""
Test script to verify the workflow logic without Snakemake
This simulates what the Snakemake pipeline would do
"""

import os
import subprocess
import sys

# Configuration
VIRAL_DB = "ncbi_viruses_20250530.fasta/sequences.fasta"
REFERENCE_GENOME = "../qc-for-hecatomb/test_data/mock_host_ref.fa"
OUTPUT_DIR = "results"
TEMP_DIR = "temp"
LOG_DIR = "logs"

def run_command(cmd, log_file):
    """Run a shell command and log output"""
    print(f"\nRunning: {' '.join(cmd)}")
    try:
        with open(log_file, 'w') as log:
            result = subprocess.run(cmd, stderr=log, stdout=subprocess.PIPE, text=True)
            if result.returncode != 0:
                print(f"Error: Command failed. Check {log_file}")
                return False
            print(f"Success! Log saved to {log_file}")
            return True
    except FileNotFoundError:
        print(f"Error: Command not found. Is BBMap installed?")
        return False

def main():
    print("Reference Masker - Test Run")
    print("===========================\n")
    
    # Check inputs
    print("Checking input files:")
    if not os.path.exists(VIRAL_DB):
        print(f"❌ Viral database not found: {VIRAL_DB}")
        sys.exit(1)
    else:
        print(f"✓ Viral database: {VIRAL_DB}")
        
    if not os.path.exists(REFERENCE_GENOME):
        print(f"❌ Reference genome not found: {REFERENCE_GENOME}")
        sys.exit(1)
    else:
        print(f"✓ Reference genome: {REFERENCE_GENOME}")
    
    # Create directories
    for directory in [OUTPUT_DIR, TEMP_DIR, LOG_DIR]:
        os.makedirs(directory, exist_ok=True)
    
    print("\nChecking for BBMap tools:")
    # Check if BBMap tools are available
    for tool in ['shred.sh', 'bbmap.sh', 'bbmask.sh']:
        result = subprocess.run(['which', tool], capture_output=True)
        if result.returncode != 0:
            print(f"❌ {tool} not found in PATH")
            print("\nBBMap tools are not installed. To install:")
            print("conda install -c bioconda bbmap")
            print("or")
            print("mamba install -c bioconda bbmap")
            sys.exit(1)
        else:
            print(f"✓ {tool} found")
    
    print("\nStarting workflow...")
    
    # Step 1: Shred viruses
    print("\n1. Shredding viral sequences...")
    shred_output = os.path.join(TEMP_DIR, "viral_shreds.fasta")
    shred_cmd = [
        'shred.sh',
        f'in={VIRAL_DB}',
        f'out={shred_output}',
        'minlength=60',
        'overlap=40',
        'median=80',
        'variance=5',
        'overwrite=true'
    ]
    if not run_command(shred_cmd, os.path.join(LOG_DIR, "shred.log")):
        sys.exit(1)
    
    # Step 2: Map to reference
    print("\n2. Mapping shreds to reference genome...")
    map_output = os.path.join(TEMP_DIR, "mapped_to_reference.sam")
    map_cmd = [
        'bbmap.sh',
        f'ref={REFERENCE_GENOME}',
        f'in={shred_output}',
        f'outm={map_output}',
        'minid=90',
        'maxindel=2',
        'ambig=all',
        'threads=4',
        'overwrite=true'
    ]
    if not run_command(map_cmd, os.path.join(LOG_DIR, "map.log")):
        sys.exit(1)
    
    # Step 3: Mask reference
    print("\n3. Masking reference genome...")
    mask_output = os.path.join(OUTPUT_DIR, "reference_masked.fasta")
    mask_cmd = [
        'bbmask.sh',
        f'in={REFERENCE_GENOME}',
        f'out={mask_output}',
        f'sam={map_output}',
        'entropy=0.5',
        'threads=4',
        'overwrite=true'
    ]
    if not run_command(mask_cmd, os.path.join(LOG_DIR, "mask.log")):
        sys.exit(1)
    
    print(f"\n✓ Workflow completed successfully!")
    print(f"✓ Masked reference saved to: {mask_output}")
    
    # Quick stats
    print("\nGenerating quick statistics...")
    with open(REFERENCE_GENOME, 'r') as f:
        original = f.read()
        original_n_count = original.upper().count('N')
        original_length = len([c for c in original if c in 'ATCGNatcgn'])
    
    with open(mask_output, 'r') as f:
        masked = f.read()
        masked_n_count = masked.upper().count('N')
        masked_length = len([c for c in masked if c in 'ATCGNatcgn'])
    
    new_n = masked_n_count - original_n_count
    percent_masked = (new_n / original_length) * 100 if original_length > 0 else 0
    
    print(f"\nMasking Statistics:")
    print(f"  Original N bases: {original_n_count}")
    print(f"  Masked N bases: {masked_n_count}")
    print(f"  New N bases added: {new_n}")
    print(f"  Percent masked: {percent_masked:.2f}%")

if __name__ == "__main__":
    main()