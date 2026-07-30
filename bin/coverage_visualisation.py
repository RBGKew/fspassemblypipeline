#!/usr/bin/env python3
#!/bin/bash
# Author: @LiaOb21
# Date: 30-04-2026
# Description: This script takes the output of samtools coverage and flagstats and produces a summary plus a figure based on these data
# Originally written for LiaOb21/FSP_assembly_benchmarking: https://github.com/LiaOb21/FSP_assembly_benchmarking/blob/main/workflow/scripts/coverage_visualisation.py

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import os

def parse_flagstat(flagstat_file):
    """
    Parse flagstat output to extract mapping statistics.

    Args:
        flagstat_file (str): Path to the flagstat output file

    Returns:
        dict: Dictionary containing mapping statistics
    """
    mapping_stats = {}

    try:
        with open(flagstat_file, 'r') as f:
            lines = f.readlines()

        # Parse each line looking for specific patterns
        for line in lines:
            line = line.strip()

            # Total reads: "13479 + 0 in total (QC-passed reads + QC-failed reads)"
            if 'in total' in line:
                mapping_stats['total_reads'] = int(line.split()[0])

            # Mapped reads: "2495 + 0 mapped (18.51% : N/A)" (but not "primary mapped")
            elif 'mapped (' in line and 'primary' not in line:
                parts = line.split()
                mapping_stats['mapped_reads'] = int(parts[0])
                # Extract percentage from format "mapped (18.51% : N/A)"
                pct_part = line.split('(')[1].split('%')[0]
                mapping_stats['mapping_rate'] = float(pct_part)

            # Properly paired reads: "2340 + 0 properly paired (17.38% : N/A)"
            elif 'properly paired' in line:
                mapping_stats['properly_paired'] = int(line.split()[0])

        # Calculate unmapped reads if we have the required data
        if 'total_reads' in mapping_stats and 'mapped_reads' in mapping_stats:
            mapping_stats['unmapped_reads'] = mapping_stats['total_reads'] - mapping_stats['mapped_reads']

    except Exception as e:
        print(f"Warning: Could not parse flagstat file: {e}")
        return {}

    return mapping_stats

def create_coverage_plots(coverage_file, output_file, sample_name, assembler_name, flagstat_file=None):
    """
    Create comprehensive coverage visualization plots with 6 subplots.

    Args:
        coverage_file (str): Path to samtools coverage output file
        output_file (str): Path for output PNG file
        sample_name (str): Name of the sample
        assembler_name (str): Name of the assembler
        flagstat_file (str, optional): Path to flagstat file for mapping statistics

    Returns:
        dict: Dictionary containing all calculated statistics
    """

    # ==================== DATA LOADING ====================

    # Read coverage data from samtools coverage output
    try:
        coverage_data = pd.read_csv(coverage_file, sep='\t')
        print(f"Loaded {len(coverage_data)} contigs from {coverage_file}")
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Read mapping statistics if flagstat file is provided
    mapping_stats = {}
    if flagstat_file and os.path.exists(flagstat_file):
        mapping_stats = parse_flagstat(flagstat_file)
        if mapping_stats:  # Only print if we successfully parsed the file
            print(f"Loaded mapping statistics from {flagstat_file}")
        else:
            print(f"Warning: Could not parse mapping statistics from {flagstat_file}")

    # ==================== CALCULATE STATISTICS ====================

    # Basic coverage statistics
    avg_coverage = coverage_data['meandepth'].mean()
    median_coverage = coverage_data['meandepth'].median()
    peak_coverage = coverage_data['meandepth'].max()
    min_coverage = coverage_data['meandepth'].min()
    std_coverage = coverage_data['meandepth'].std()
    total_bases = coverage_data['endpos'].sum()

    # Weighted average coverage (larger contigs have more influence)
    weighted_avg_coverage = np.average(coverage_data['meandepth'], weights=coverage_data['endpos'])

    # Overall coverage percentage across all contigs
    covered_bases = (coverage_data['coverage'] * coverage_data['endpos'] / 100).sum()
    overall_coverage_pct = (covered_bases / total_bases) * 100 if total_bases > 0 else 0

    # Longest and shortest contigs
    longest_contig = coverage_data['endpos'].max()
    shortest_contig = coverage_data['endpos'].min()

    # Count contigs smaller than 250bp
    small_contigs_count = len(coverage_data[coverage_data['endpos'] < 250])
    small_contigs_pct = (small_contigs_count / len(coverage_data)) * 100

    # Count contigs longer than 1kb
    long_contigs_count = len(coverage_data[coverage_data['endpos'] >= 1000])
    long_contigs_pct = (long_contigs_count / len(coverage_data)) * 100

    # Count contigs between 250bp and 1kb
    mid_contigs_count = len(coverage_data[(coverage_data['endpos'] >= 250) & (coverage_data['endpos'] < 1000)])
    mid_contigs_pct = (mid_contigs_count / len(coverage_data)) * 100

    # ==================== PRINT SUMMARY STATISTICS ====================

    print(f"\n=== Coverage Summary for {sample_name} ({assembler_name}) ===")
    print(f"Total contigs: {len(coverage_data):,}")
    print(f"Contigs < 250bp: {small_contigs_count:,} ({small_contigs_pct:.1f}%)")
    print(f"Total bases: {total_bases:,}")
    print(f"Average coverage: {avg_coverage:.2f}x")
    print(f"Weighted average coverage: {weighted_avg_coverage:.2f}x")
    print(f"Median coverage: {median_coverage:.2f}x")
    print(f"Peak coverage: {peak_coverage:.2f}x")
    print(f"Min coverage: {min_coverage:.2f}x")
    print(f"Standard deviation: {std_coverage:.2f}x")
    print(f"Overall coverage: {overall_coverage_pct:.2f}%")

    # Print mapping statistics if available
    if mapping_stats:
        print(f"\n=== Mapping Statistics ===")
        print(f"Total reads: {mapping_stats.get('total_reads', 'N/A'):,}")
        print(f"Mapped reads: {mapping_stats.get('mapped_reads', 'N/A'):,}")
        print(f"Unmapped reads: {mapping_stats.get('unmapped_reads', 'N/A'):,}")
        print(f"Mapping rate: {mapping_stats.get('mapping_rate', 'N/A'):.1f}%")
        print(f"Properly paired: {mapping_stats.get('properly_paired', 'N/A'):,}")
    else:
        print(f"\n=== Mapping Statistics ===")
        print("No flagstat file provided - mapping statistics not available")

    # ==================== CREATE VISUALIZATION ====================

    # Create 2x3 subplot grid
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Coverage Analysis - {sample_name} ({assembler_name})',
                fontsize=16, fontweight='bold')

    # Plot 1: Coverage distribution histogram (top-left)
    ax1 = axes[0, 0]

    # Option 1: X-axis log scale
    ax1.set_xscale('log')

    # Create log-spaced bins for better distribution
    min_coverage = max(0.1, coverage_data['meandepth'].min())  # Avoid log(0)
    max_coverage = coverage_data['meandepth'].max()
    log_bins_coverage = np.logspace(np.log10(min_coverage), np.log10(max_coverage), 40)

    n, bins, patches = ax1.hist(coverage_data['meandepth'], bins=log_bins_coverage, alpha=0.7,
                            color='skyblue', edgecolor='black', density=False)

    # Force Y-axis to show integer counts
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax1.set_ylim(0, max(n) + 1)

    # Add vertical lines for key statistics
    ax1.axvline(avg_coverage, color='red', linestyle='--', linewidth=2,
            label=f'Mean: {avg_coverage:.1f}x')
    ax1.axvline(median_coverage, color='green', linestyle='--', linewidth=2,
            label=f'Median: {median_coverage:.1f}x')
    ax1.axvline(peak_coverage, color='orange', linestyle='--', linewidth=2,
            label=f'Peak: {peak_coverage:.1f}x')
    ax1.set_xlabel('Coverage Depth (x) - Log Scale')  # Update label
    ax1.set_ylabel('Number of Contigs')
    ax1.set_title('Coverage Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Coverage vs Contig Length scatter plot (top-middle)
    ax2 = axes[0, 1]

    # ADD LOG SCALE FOR BOTH AXES
    ax2.set_xscale('log')  # Log scale for contig length
    ax2.set_yscale('log')  # Log scale for coverage depth

    # Color points by coverage percentage, size by default
    scatter = ax2.scatter(coverage_data['endpos'], coverage_data['meandepth'],
                         alpha=0.6, c=coverage_data['coverage'], cmap='viridis', s=20)

    # Add horizontal reference lines
    ax2.axhline(avg_coverage, color='red', linestyle='--', linewidth=2,
               label=f'Mean: {avg_coverage:.1f}x')
    ax2.axhline(weighted_avg_coverage, color='purple', linestyle='--', linewidth=2,
               label=f'Weighted Mean: {weighted_avg_coverage:.1f}x')

    ax2.set_xlabel('Contig Length (bp) - Log Scale')  # Update label
    ax2.set_ylabel('Coverage Depth (x) - Log Scale')   # Update label
    ax2.set_title('Coverage vs Contig Length')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax2, label='Coverage %')

    # Plot 3: Cumulative coverage distribution (top-right)
    ax3 = axes[0, 2]

    # ADD LOG SCALE FOR X-AXIS
    ax3.set_xscale('log')

    # Sort coverage values and create cumulative count
    sorted_coverage = np.sort(coverage_data['meandepth'])
    cumulative_count = np.arange(1, len(sorted_coverage) + 1)
    ax3.plot(sorted_coverage, cumulative_count, linewidth=2, color='blue')

    # Add vertical reference lines
    ax3.axvline(avg_coverage, color='red', linestyle='--', linewidth=2,
               label=f'Mean: {avg_coverage:.1f}x')
    ax3.axvline(median_coverage, color='green', linestyle='--', linewidth=2,
               label=f'Median: {median_coverage:.1f}x')

    ax3.set_xlabel('Coverage Depth (x) - Log Scale')  # Update label
    ax3.set_ylabel('Number of Contigs')
    ax3.set_title('Cumulative Coverage Distribution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Plot 4: Four-panel contig length distribution (bottom-left AND bottom-middle) - LINEAR ONLY
    # This will span across both bottom-left and bottom-middle positions
    ax4 = axes[1, 0]
    ax5 = axes[1, 1]  # We'll use this space too

    # Turn off both axes
    ax4.axis('off')
    ax5.axis('off')

    # Create a combined title spanning both quadrants
    fig.text(0.35, 0.45, 'Contig Length Distribution by Size Class',
            fontsize=12, ha='center')

    # Data preparation
    small_contigs = coverage_data[coverage_data['endpos'] < 250]['endpos']
    mid_contigs = coverage_data[(coverage_data['endpos'] >= 250) & (coverage_data['endpos'] < 1000)]['endpos']
    large_contigs = coverage_data[coverage_data['endpos'] >= 1000]['endpos']
    above_250 = coverage_data[coverage_data['endpos'] >= 250]['endpos']

    # Get positions of both quadrants to create a larger plotting area
    pos4 = ax4.get_position()
    pos5 = ax5.get_position()

    # Create a combined plotting area spanning both quadrants
    combined_width = pos5.x1 - pos4.x0
    combined_height = pos4.height

    # Add spacing between subplots
    horizontal_spacing = 0.05  # 3% horizontal spacing between left and right columns
    vertical_spacing = 0.05    # 4% vertical spacing between top and bottom rows
#    downward_shift = 0.02

    # Calculate subplot dimensions accounting for spacing
    subplot_width = (combined_width - horizontal_spacing) / 2
    subplot_height = (combined_height - vertical_spacing - 0.03) / 2

    # Calculate positions for the 4 subplots with spacing
    x1 = pos4.x0
    x2 = pos4.x0 + subplot_width + horizontal_spacing
    y1 = pos4.y0 + subplot_height + vertical_spacing - 0.02 # Top row (with spacing below)
    y2 = pos4.y0 - 0.03 # Bottom row

    # Subplot 1: Small contigs < 250bp (top-left)
    ax4_1 = fig.add_axes([x1, y1, subplot_width, subplot_height])
    if len(small_contigs) > 0:
        # ADD LOG SCALE AND LOG BINS
        min_small = small_contigs.min()
        log_bins_small = np.logspace(np.log10(max(1, min_small)), np.log10(250), 20)

        n1, bins1, patches1 = ax4_1.hist(small_contigs, bins=log_bins_small, alpha=0.7, color='lightcoral', edgecolor='black')
        ax4_1.set_xscale('log')  # Add log scale
        ax4_1.set_xlim(max(1, min_small), 250)
        ax4_1.set_ylim(0, max(n1) * 1.1)
    else:
        ax4_1.text(0.5, 0.5, 'No contigs\n< 250bp', ha='center', va='center', transform=ax4_1.transAxes)
    ax4_1.set_title(f'< 250bp: {len(small_contigs):,} ({small_contigs_pct:.1f}%)', fontsize=10)
    ax4_1.set_ylabel('Count', fontsize=10)
    ax4_1.set_xlabel('Length (bp) - Log Scale', fontsize=10)  # Update label
    ax4_1.tick_params(labelsize=9)
    ax4_1.grid(True, alpha=0.3)

    # Subplot 2: Medium contigs 250-1000bp (top-right)
    ax4_2 = fig.add_axes([x2, y1, subplot_width, subplot_height])
    if len(mid_contigs) > 0:
        # ADD LOG SCALE AND LOG BINS
        log_bins_mid = np.logspace(np.log10(250), np.log10(1000), 20)

        n2, bins2, patches2 = ax4_2.hist(mid_contigs, bins=log_bins_mid, alpha=0.7, color='lightblue', edgecolor='black')
        ax4_2.set_xscale('log')  # Add log scale
        ax4_2.set_xlim(250, 1000)
        ax4_2.set_ylim(0, max(n2) * 1.1)
    else:
        ax4_2.text(0.5, 0.5, 'No contigs\n250-1000bp', ha='center', va='center', transform=ax4_2.transAxes)
    ax4_2.set_title(f'250bp-1kb: {len(mid_contigs):,} ({mid_contigs_pct:.1f}%)', fontsize=10)
    ax4_2.set_xlabel('Length (bp) - Log Scale', fontsize=10)  # Update label
    ax4_2.tick_params(labelsize=9)
    ax4_2.grid(True, alpha=0.3)

    # Subplot 3: Large contigs > 1000bp (bottom-left)
    ax4_3 = fig.add_axes([x1, y2, subplot_width, subplot_height])
    if len(large_contigs) > 0:
        # ADD LOG SCALE AND LOG BINS
        max_large = large_contigs.max()
        log_bins_large = np.logspace(np.log10(1000), np.log10(max_large), 20)

        n3, bins3, patches3 = ax4_3.hist(large_contigs, bins=log_bins_large, alpha=0.7, color='lightgreen', edgecolor='black')
        ax4_3.set_xscale('log')  # Add log scale
        ax4_3.set_xlabel('Length (bp) - Log Scale', fontsize=10)  # Update label
        ax4_3.set_ylim(0, max(n3) * 1.1)
    else:
        ax4_3.text(0.5, 0.5, 'No contigs\n> 1kb', ha='center', va='center', transform=ax4_3.transAxes)
    ax4_3.set_title(f'> 1kb: {len(large_contigs):,} ({long_contigs_pct:.1f}%)', fontsize=10)
    ax4_3.set_ylabel('Count', fontsize=10)
    ax4_3.tick_params(labelsize=9)
    ax4_3.grid(True, alpha=0.3)

    # Subplot 4: All contigs >= 250bp (bottom-right)
    ax4_4 = fig.add_axes([x2, y2, subplot_width, subplot_height])
    if len(above_250) > 0:
        # ADD LOG SCALE AND LOG BINS
        max_above250 = above_250.max()
        log_bins_above250 = np.logspace(np.log10(250), np.log10(max_above250), 25)

        n4, bins4, patches4 = ax4_4.hist(above_250, bins=log_bins_above250, alpha=0.7, color='mediumpurple', edgecolor='black')
        ax4_4.set_xscale('log')  # Add log scale
        ax4_4.axvline(1000, color='red', linestyle='--', linewidth=2, label='1kb threshold')
        ax4_4.set_xlabel('Length (bp) - Log Scale', fontsize=10)  # Update label
        ax4_4.set_ylim(0, max(n4) * 1.1)
        ax4_4.legend(fontsize=9)
    else:
        ax4_4.text(0.5, 0.5, 'No contigs\n≥ 250bp', ha='center', va='center', transform=ax4_4.transAxes)
    ax4_4.set_title(f'≥ 250bp: {len(above_250):,} ({100-small_contigs_pct:.1f}%)', fontsize=10)
    ax4_4.tick_params(labelsize=9)
    ax4_4.grid(True, alpha=0.3)


    # Plot 6: Combined Coverage Statistics and Quality Metrics (bottom-right quadrant)
    ax6 = axes[1, 2]
    ax6.axis('off')

    # Prepare combined statistics text
    coverage_stats_text = f"""Coverage Statistics:

Mean Coverage: {avg_coverage:.2f}x
Weighted Mean: {weighted_avg_coverage:.2f}x
Median Coverage: {median_coverage:.2f}x
Peak Coverage: {peak_coverage:.2f}x
Min Coverage: {min_coverage:.2f}x
Std Deviation: {std_coverage:.2f}x

Assembly Info:
Total Contigs: {len(coverage_data):,}
Contigs < 250bp: {small_contigs_count:,} ({small_contigs_pct:.1f}%)
Contigs 250bp-1kb: {mid_contigs_count:,} ({mid_contigs_pct:.1f}%)
Contigs > 1kb: {long_contigs_count:,} ({long_contigs_pct:.1f}%)
Total Length: {total_bases:,} bp
Longest Contig: {longest_contig:,} bp
Shortest Contig: {shortest_contig:,} bp
Overall Coverage: {overall_coverage_pct:.1f}%
"""

    # Prepare quality metrics text
    if mapping_stats:
        quality_metrics_text = f"""Quality & Mapping Metrics:

Mapping Statistics:
Total Reads: {mapping_stats.get('total_reads', 'N/A'):,}
Mapped Reads: {mapping_stats.get('mapped_reads', 'N/A'):,}
Unmapped Reads: {mapping_stats.get('unmapped_reads', 'N/A'):,}
Mapping Rate: {mapping_stats.get('mapping_rate', 'N/A'):.1f}%
Properly Paired: {mapping_stats.get('properly_paired', 'N/A'):,}

Quality Scores:
Mean Base Quality: {coverage_data['meanbaseq'].mean():.1f}
Median Base Quality: {coverage_data['meanbaseq'].median():.1f}
Mean Mapping Quality: {coverage_data['meanmapq'].mean():.1f}
Median Mapping Quality: {coverage_data['meanmapq'].median():.1f}

Coverage Info:
Total Mapped Reads: {coverage_data['numreads'].sum():,}
Avg Reads per Contig: {coverage_data['numreads'].mean():.1f}
"""
    else:
        quality_metrics_text = f"""Quality Metrics:

Mean Base Quality: {coverage_data['meanbaseq'].mean():.1f}
Median Base Quality: {coverage_data['meanbaseq'].median():.1f}
Mean Mapping Quality: {coverage_data['meanmapq'].mean():.1f}
Median Mapping Quality: {coverage_data['meanmapq'].median():.1f}

Read Mapping:
Total Mapped Reads: {coverage_data['numreads'].sum():,}
Avg Reads per Contig: {coverage_data['numreads'].mean():.1f}
Max Reads per Contig: {coverage_data['numreads'].max():,}

(No flagstat file provided for mapping stats)
"""

    # Get position of ax6 for creating two side-by-side text boxes
    pos6 = ax6.get_position()

    # Calculate centered positioning
    total_width = pos6.width * 0.9  # Use 90% of available width
    box_spacing = 0.05              # 5% spacing between boxes
    box_width = (total_width - box_spacing) / 2
    box_height = pos6.height * 0.85 # Use 85% of available height

    # Center the entire layout in the quadrant
    left_margin = 0.05  # 5% left margin
    start_x = pos6.x0 + left_margin
    start_y = pos6.y0 + (pos6.height - box_height) / 2

    # Calculate positions for both boxes
    left_x = start_x
    right_x = start_x + box_width + box_spacing

    # Coverage Statistics box (left half)
    ax6_left = fig.add_axes([left_x, start_y, box_width, box_height])
    ax6_left.axis('off')
    ax6_left.text(0.5, 0.5, coverage_stats_text, transform=ax6_left.transAxes,
                fontsize=8, verticalalignment='center', horizontalalignment='center',
                fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightblue', alpha=0.9,
                        edgecolor='darkblue', linewidth=1))
#    ax6_left.set_title('Coverage Statistics', fontweight='bold', fontsize=10)

    # Quality Metrics box (right half)
    ax6_right = fig.add_axes([right_x, start_y, box_width, box_height])
    ax6_right.axis('off')
    ax6_right.text(0.5, 0.5, quality_metrics_text, transform=ax6_right.transAxes,
                fontsize=8, verticalalignment='center', horizontalalignment='center',
                fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='lightcoral', alpha=0.9,
                            edgecolor='darkred', linewidth=1))
#    ax6_right.set_title('Quality & Mapping Metrics', fontweight='bold', fontsize=10)

    # ==================== FINALIZE AND SAVE ====================

    # Adjust layout to prevent overlapping
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.03)

    # Save the plot with high resolution
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()  # Close the figure to free memory

    print(f"Plot saved to: {output_file}")

    # ==================== RETURN STATISTICS DICTIONARY ====================

    # Prepare comprehensive statistics dictionary for summary file
    stats_dict = {
        'sample': sample_name,
        'assembler': assembler_name,
        'total_contigs': len(coverage_data),
        'contigs_under_250bp': small_contigs_count,
        'contigs_under_250bp_percent': small_contigs_pct,
        'contigs_250bp_to_1kb': mid_contigs_count,
        'contigs_250bp_to_1kb_percent': mid_contigs_pct,
        'contigs_over_1kb': long_contigs_count,
        'contigs_over_1kb_percent': long_contigs_pct,
        'total_bases': total_bases,
        'mean_coverage': avg_coverage,
        'weighted_mean_coverage': weighted_avg_coverage,
        'median_coverage': median_coverage,
        'peak_coverage': peak_coverage,
        'min_coverage': min_coverage,
        'std_coverage': std_coverage,
        'overall_coverage_pct': overall_coverage_pct,
        'longest_contig': longest_contig,
        'shortest_contig': shortest_contig,
        'mean_base_quality': coverage_data['meanbaseq'].mean(),
        'median_base_quality': coverage_data['meanbaseq'].median(),
        'mean_mapping_quality': coverage_data['meanmapq'].mean(),
        'median_mapping_quality': coverage_data['meanmapq'].median(),
    }

    # Add mapping statistics if available
    if mapping_stats:
        stats_dict.update({
            'total_reads': mapping_stats.get('total_reads', 'N/A'),
            'mapped_reads': mapping_stats.get('mapped_reads', 'N/A'),
            'unmapped_reads': mapping_stats.get('unmapped_reads', 'N/A'),
            'mapping_rate': mapping_stats.get('mapping_rate', 'N/A'),
            'properly_paired': mapping_stats.get('properly_paired', 'N/A')
        })

    return stats_dict

def main():
    parser = argparse.ArgumentParser(description='Create coverage visualization plots')
    parser.add_argument('coverage_file', help='Input coverage file (samtools coverage output)')
    parser.add_argument('output_file', help='Output plot file (PNG/PDF)')
    parser.add_argument('--sample', default='Sample', help='Sample name')
    parser.add_argument('--assembler', default='Assembler', help='Assembler name')
    parser.add_argument('--flagstat', help='Optional flagstat file for mapping statistics')
    parser.add_argument('--summary', help='Output summary file path')  # Add this line

    args = parser.parse_args()

    # Create the visualization
    stats = create_coverage_plots(args.coverage_file, args.output_file,
                                 args.sample, args.assembler, args.flagstat)

    # Use provided summary file path or generate one
    if args.summary:
        summary_file = args.summary
    else:
        summary_file = args.output_file.replace('.png', '_summary.txt').replace('.pdf', '_summary.txt')

    # Save summary file
    with open(summary_file, 'w') as f:
        f.write(f"Coverage Analysis Summary - {args.sample} ({args.assembler})\n")
        f.write("=" * 60 + "\n\n")
        for key, value in stats.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.2f}\n")
            else:
                f.write(f"{key}: {value}\n")

    print(f"Summary statistics saved to: {summary_file}")

# Run the main function if script is executed directly
if __name__ == "__main__":
    main()
