#!/usr/bin/env python3
# Author: @LiaOb21
# Date: 23-10-2026
# Description: This script collects the metrics from the various final QC outputs of the draft assemblies for each sample and puts them in a single table.
# compared to extract_metrics.py it doesn't process coverage results, and reads multiple lines from quast results (one per assembly).


import pandas as pd
import re
import os
import argparse
import sys

def extract_lineage_from_busco_filename(filepath):
    match = re.search(r'short_summary\.specific\.([^.]+)\.', os.path.basename(filepath))
    return match.group(1) if match else None

def parse_busco_file(filepath):
    """
    Parse a BUSCO summary file and extract relevant statistics.
    Excludes assembly statistics (num_scaffolds, num_contigs, total_length, lineage).

    Args:
        filepath (str): Path to BUSCO summary file

    Returns:
        dict: Dictionary containing parsed BUSCO statistics
    """
    results = {}

    with open(filepath, 'r') as f:
        content = f.read()

    # Extract the results line (contains percentages)
    results_pattern = r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)'
    match = re.search(results_pattern, content)

    if match:
        results['complete_percent'] = float(match.group(1))
        results['single_copy_percent'] = float(match.group(2))
        results['duplicated_percent'] = float(match.group(3))
        results['fragmented_percent'] = float(match.group(4))
        results['missing_percent'] = float(match.group(5))
        results['total_buscos'] = int(match.group(6))

    # Extract absolute numbers
    complete_match = re.search(r'(\d+)\s+Complete BUSCOs \(C\)', content)
    single_copy_match = re.search(r'(\d+)\s+Complete and single-copy BUSCOs \(S\)', content)
    duplicated_match = re.search(r'(\d+)\s+Complete and duplicated BUSCOs \(D\)', content)
    fragmented_match = re.search(r'(\d+)\s+Fragmented BUSCOs \(F\)', content)
    missing_match = re.search(r'(\d+)\s+Missing BUSCOs \(M\)', content)

    if complete_match:
        results['complete_count'] = int(complete_match.group(1))
    if single_copy_match:
        results['single_copy_count'] = int(single_copy_match.group(1))
    if duplicated_match:
        results['duplicated_count'] = int(duplicated_match.group(1))
    if fragmented_match:
        results['fragmented_count'] = int(fragmented_match.group(1))
    if missing_match:
        results['missing_count'] = int(missing_match.group(1))

    return results


def parse_quast_file(filepath, assembly_name):
    """
    Parse a QUAST transposed_report.tsv file and extract relevant statistics.
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        df = df[df['Assembly'] == assembly_name]
        if len(df) == 0:
            print(f"Assembly {assembly_name} not found in QUAST file {filepath}")
            sys.exit(1)
            return {}

        row = df.iloc[0]

        results = {
            'total_contigs>250bp': row.get('# contigs', None),
            'contigs>1000bp': row.get('# contigs (>= 1000 bp)', None),
            'contigs>5000bp': row.get('# contigs (>= 5000 bp)', None),
            'contigs>10000bp': row.get('# contigs (>= 10000 bp)', None),
            'contigs>25000bp': row.get('# contigs (>= 25000 bp)', None),
            'total_length>250bp': row.get('Total length', None),
            'largest_contig': row.get('Largest contig', None),
            'gc_percent': row.get('GC (%)', None),
            'aUN': row.get('auN', None),
            'N50': row.get('N50', None),
            'N90': row.get('N90', None),
            'L50': row.get('L50', None),
            'L90': row.get('L90', None),
            'Ns_per_100kbp': row.get('# N\'s per 100 kbp', None)
        }

        # Convert to appropriate types - preserve decimals
        for key, value in results.items():
            if pd.isna(value):
                results[key] = None
            else:
                try:
                    float_val = float(value)
                    if float_val.is_integer():
                        results[key] = int(float_val)
                    else:
                        results[key] = float_val
                except (ValueError, TypeError):
                    results[key] = value

        return results
    except Exception as e:
        print(f"Error parsing QUAST file {filepath}: {e}")
        return {}


def parse_merqury_qv_file(filepath):
    """
    Parse a merquryfk.qv file and extract Error % and QV.
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        if len(df) == 0:
            return {}

        row = df.iloc[0]

        results = {
            'error_percent': row.get('Error %', None),
            'qv': row.get('QV', None)
        }

        for key, value in results.items():
            if pd.isna(value) or value == 'inf':
                results[key] = None if value != 'inf' else float('inf')
            else:
                try:
                    float_val = float(value)
                    results[key] = float_val
                except (ValueError, TypeError):
                    results[key] = value

        return results

    except Exception as e:
        print(f"Error parsing Merqury QV file {filepath}: {e}")
        return {}


def parse_merqury_completeness_file(filepath):
    """
    Parse a merquryfk.completeness.stats file and extract % Covered.
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        if len(df) == 0:
            return {}

        row = df.iloc[0]

        results = {
            'percent_covered': row.get('% Covered', None)
        }

        for key, value in results.items():
            if pd.isna(value):
                results[key] = None
            else:
                try:
                    float_val = float(value)
                    results[key] = float_val
                except (ValueError, TypeError):
                    results[key] = value

        return results

    except Exception as e:
        print(f"Error parsing Merqury completeness file {filepath}: {e}")
        return {}


def main():
    parser = argparse.ArgumentParser(
        description='Extract assembly metrics for a single sample'
    )
    parser.add_argument('--sample-id',            required=True)
    parser.add_argument('--reads-type',           required=True)
    parser.add_argument('--kmer-strategy',        required=True)
    parser.add_argument('--assembler',            required=True)
    parser.add_argument('--busco-general',        required=True)
    parser.add_argument('--busco-specific',       required=True)
    parser.add_argument('--quast-tsv',            required=True)
    parser.add_argument('--merqury-qv',           required=True)
    parser.add_argument('--merqury-completeness', required=True)
    parser.add_argument('--output',               required=True)

    args = parser.parse_args()

    # Lineage still comes from the BUSCO filename (BUSCO sets it, not us)
    lineage_general  = extract_lineage_from_busco_filename(args.busco_general)
    lineage_specific = extract_lineage_from_busco_filename(args.busco_specific)

    # Parse all input files
    busco_general_data  = parse_busco_file(args.busco_general)
    busco_specific_data = parse_busco_file(args.busco_specific)
    assembly_name       = f"{args.sample_id}_{args.reads_type}_{args.kmer_strategy}_{args.assembler}"
    quast_data          = parse_quast_file(args.quast_tsv, assembly_name)
    merqury_qv_data     = parse_merqury_qv_file(args.merqury_qv)
    merqury_comp_data   = parse_merqury_completeness_file(args.merqury_completeness)

    # Build a single flat row — meta fields first, then metrics
    row = {
        'sample_id':        args.sample_id,
        'reads_type':       args.reads_type,
        'kmer_strategy':    args.kmer_strategy,
        'assembler':        args.assembler,
        'lineage_general':  lineage_general,
        'lineage_specific': lineage_specific,
    }

    # BUSCO metrics: add *_general / *_specific suffix to distinguish the two runs
    for key, val in busco_general_data.items():
        row[f'{key}_general'] = val
    for key, val in busco_specific_data.items():
        row[f'{key}_specific'] = val

    # Add all other metrics to the row
    row.update(quast_data)
    row.update(merqury_qv_data)
    row.update(merqury_comp_data)

    pd.DataFrame([row]).to_csv(args.output, index=False, sep='\t')


if __name__ == "__main__":
    main()
