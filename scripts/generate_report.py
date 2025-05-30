#!/usr/bin/env python3
"""
Generate comprehensive masking report from BBMask results
"""

import argparse
import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import SeqIO
import numpy as np

def parse_mask_log(log_file):
    """Parse BBMask log file to extract masking statistics"""
    stats = {}
    
    with open(log_file, 'r') as f:
        content = f.read()
        
    # Extract key statistics
    ref_bases_match = re.search(r'Ref Bases:\s+(\d+)', content)
    if ref_bases_match:
        stats['total_bases'] = int(ref_bases_match.group(1))
    
    low_complexity_match = re.search(r'Low Complexity Bases:\s+(\d+)', content)
    if low_complexity_match:
        stats['low_complexity_bases'] = int(low_complexity_match.group(1))
    
    sam_masked_match = re.search(r'Sam Bases Masked:\s+(\d+)', content)
    if sam_masked_match:
        stats['viral_masked_bases'] = int(sam_masked_match.group(1))
    
    total_masked_match = re.search(r'Total Bases Masked:\s+(\d+)/(\d+)\s+([\d.]+)%', content)
    if total_masked_match:
        stats['total_masked_bases'] = int(total_masked_match.group(1))
        stats['masking_percentage'] = float(total_masked_match.group(3))
    
    return stats

def parse_map_log(log_file):
    """Parse BBMap log file to extract mapping statistics"""
    stats = {}
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Extract mapping results
    reads_used_match = re.search(r'Reads Used:\s+(\d+)', content)
    if reads_used_match:
        stats['total_reads'] = int(reads_used_match.group(1))
    
    mapped_match = re.search(r'mapped:\s+([\d.]+)%\s+(\d+)', content)
    if mapped_match:
        stats['mapped_percentage'] = float(mapped_match.group(1))
        stats['mapped_reads'] = int(mapped_match.group(2))
    
    match_rate_match = re.search(r'Match Rate:\s+.*?\s+([\d.]+)%', content)
    if match_rate_match:
        stats['match_rate'] = float(match_rate_match.group(1))
    
    return stats

def find_masked_regions(fasta_file):
    """Find coordinates of masked (N) regions in FASTA file"""
    regions = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        chrom = record.id
        
        # Find runs of N's
        in_n_region = False
        start = 0
        
        for i, base in enumerate(seq):
            if base.upper() == 'N':
                if not in_n_region:
                    start = i
                    in_n_region = True
            else:
                if in_n_region:
                    regions.append({
                        'chrom': chrom,
                        'start': start,
                        'end': i,
                        'length': i - start
                    })
                    in_n_region = False
        
        # Handle case where sequence ends with N's
        if in_n_region:
            regions.append({
                'chrom': chrom,
                'start': start,
                'end': len(seq),
                'length': len(seq) - start
            })
    
    return regions

def generate_bed_file(regions, output_file):
    """Generate BED file from masked regions"""
    with open(output_file, 'w') as f:
        f.write("# BED file of masked regions\n")
        f.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
        
        for i, region in enumerate(regions):
            f.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t"
                   f"masked_region_{i+1}\t{region['length']}\t.\n")

def create_visualizations(stats_dict, regions_dict, output_dir):
    """Create visualization plots"""
    
    # Set up plotting style
    plt.style.use('seaborn-v0_8')
    sns.set_palette("husl")
    
    # 1. Masking percentage comparison
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Masking percentages
    samples = list(stats_dict.keys())
    mask_pcts = [stats_dict[s]['masking_percentage'] for s in samples]
    
    axes[0,0].bar(samples, mask_pcts, color='lightcoral')
    axes[0,0].set_title('Masking Percentage by Sample')
    axes[0,0].set_ylabel('Percentage Masked (%)')
    axes[0,0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Viral vs Low-complexity masking
    viral_bases = [stats_dict[s].get('viral_masked_bases', 0) for s in samples]
    lowc_bases = [stats_dict[s].get('low_complexity_bases', 0) for s in samples]
    
    x = np.arange(len(samples))
    width = 0.35
    
    axes[0,1].bar(x - width/2, viral_bases, width, label='Viral regions', color='red', alpha=0.7)
    axes[0,1].bar(x + width/2, lowc_bases, width, label='Low complexity', color='blue', alpha=0.7)
    axes[0,1].set_title('Masked Bases by Type')
    axes[0,1].set_ylabel('Bases Masked')
    axes[0,1].set_xticks(x)
    axes[0,1].set_xticklabels(samples, rotation=45)
    axes[0,1].legend()
    
    # Plot 3: Region length distribution
    if regions_dict:
        all_lengths = []
        for sample, regions in regions_dict.items():
            all_lengths.extend([r['length'] for r in regions])
        
        if all_lengths:
            axes[1,0].hist(all_lengths, bins=50, alpha=0.7, color='green')
            axes[1,0].set_xlabel('Region Length (bp)')
            axes[1,0].set_ylabel('Count')
            axes[1,0].set_title('Distribution of Masked Region Lengths')
            axes[1,0].set_yscale('log')
    
    # Plot 4: Mapping statistics
    if any('mapped_percentage' in stats_dict[s] for s in samples):
        map_pcts = [stats_dict[s].get('mapped_percentage', 0) for s in samples]
        axes[1,1].bar(samples, map_pcts, color='orange', alpha=0.7)
        axes[1,1].set_title('Viral Read Mapping Percentage')
        axes[1,1].set_ylabel('Mapped Reads (%)')
        axes[1,1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'masking_summary.png'), dpi=300, bbox_inches='tight')
    plt.close()

def generate_html_report(stats_dict, regions_dict, output_file):
    """Generate HTML summary report"""
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Reference Masker Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
            .sample {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
            .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; }}
            .stat {{ background-color: #f9f9f9; padding: 10px; border-radius: 3px; }}
            .plot {{ text-align: center; margin: 20px 0; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Reference Masker Analysis Report</h1>
            <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>Samples processed: {len(stats_dict)}</p>
        </div>
        
        <div class="plot">
            <h2>Summary Visualizations</h2>
            <img src="masking_summary.png" alt="Masking Summary Plots" style="max-width: 100%; height: auto;">
        </div>
        
        <h2>Sample Details</h2>
    """
    
    for sample, stats in stats_dict.items():
        regions = regions_dict.get(sample, [])
        
        html_content += f"""
        <div class="sample">
            <h3>{sample}</h3>
            <div class="stats">
                <div class="stat">
                    <strong>Total Bases:</strong><br>
                    {stats.get('total_bases', 'N/A'):,}
                </div>
                <div class="stat">
                    <strong>Masking Percentage:</strong><br>
                    {stats.get('masking_percentage', 'N/A'):.3f}%
                </div>
                <div class="stat">
                    <strong>Viral Bases Masked:</strong><br>
                    {stats.get('viral_masked_bases', 'N/A'):,}
                </div>
                <div class="stat">
                    <strong>Low Complexity Bases:</strong><br>
                    {stats.get('low_complexity_bases', 'N/A'):,}
                </div>
                <div class="stat">
                    <strong>Mapped Reads:</strong><br>
                    {stats.get('mapped_reads', 'N/A'):,} ({stats.get('mapped_percentage', 'N/A'):.4f}%)
                </div>
                <div class="stat">
                    <strong>Masked Regions:</strong><br>
                    {len(regions):,}
                </div>
            </div>
        </div>
        """
    
    html_content += """
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Generate masking analysis report')
    parser.add_argument('--log-dir', required=True, help='Directory containing log files')
    parser.add_argument('--results-dir', required=True, help='Directory containing masked FASTA files')
    parser.add_argument('--output-dir', required=True, help='Output directory for reports')
    parser.add_argument('--samples', nargs='+', required=True, help='Sample names to process')
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    stats_dict = {}
    regions_dict = {}
    
    # Process each sample
    for sample in args.samples:
        print(f"Processing sample: {sample}")
        
        # Parse log files
        mask_log = os.path.join(args.log_dir, f"{sample}_mask.log")
        map_log = os.path.join(args.log_dir, f"{sample}_map.log")
        
        stats = {}
        if os.path.exists(mask_log):
            stats.update(parse_mask_log(mask_log))
        if os.path.exists(map_log):
            stats.update(parse_map_log(map_log))
        
        stats_dict[sample] = stats
        
        # Find masked regions
        masked_fasta = os.path.join(args.results_dir, f"{sample}_masked.fasta")
        if os.path.exists(masked_fasta):
            regions = find_masked_regions(masked_fasta)
            regions_dict[sample] = regions
            
            # Generate BED file
            bed_file = os.path.join(args.output_dir, f"{sample}_masked_regions.bed")
            generate_bed_file(regions, bed_file)
            print(f"Generated BED file: {bed_file}")
    
    # Create visualizations
    print("Generating visualizations...")
    create_visualizations(stats_dict, regions_dict, args.output_dir)
    
    # Generate HTML report
    print("Generating HTML report...")
    html_file = os.path.join(args.output_dir, "masking_report.html")
    generate_html_report(stats_dict, regions_dict, html_file)
    
    print(f"Report generated: {html_file}")
    print(f"Plots saved to: {args.output_dir}/masking_summary.png")

if __name__ == "__main__":
    main()