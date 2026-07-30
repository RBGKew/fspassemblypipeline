#!/usr/bin/env python3
# Author: Wu Huang <w.huang@kew.org>
# Date: 2026-03-11
# Description: Extract key columns and classify k-mer peak shapes.

import argparse
import csv
import os


OUTPUT_COLUMNS = [
    "Seq_ID",
    "Total_number_of_kmer",
    "peak_position",
    "Shape",
    "Estimated_genome_size(bp)",
    "Unique_kmer_number",
    "Kmer_automation_note",
]


def to_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def to_int_rounded(value):
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return None


def build_output_row(row):
    peak_count = to_int(row.get("peak_count"))
    peak1_fit_max = to_float(row.get("peak1_model_fit_max"))
    peak2_fit_max = to_float(row.get("peak2_model_fit_max"))

    shape = ""
    est_size = ""
    peak_position = ""
    note = ""

    if peak1_fit_max < 60 and peak2_fit_max < 60:
        shape = "L-shape"
        est_size = "NA"
        peak_position = "NA"
        note = "uncertain"
    elif peak_count == 0:
        shape = "L-shape"
        est_size = "NA"
        peak_position = "NA"
    elif peak1_fit_max > peak2_fit_max:
        shape = "sharp single peak"
        est_size = row.get("peak1_est_genome_size(len)", "")
        peak_position = row.get("peak1_kcov", "")
    elif peak_count == 2:
        shape = "Double diploid peak"
        est_size = row.get("peak2_est_genome_size(len)", "")
        peak_position = row.get("peak2_kcov", "")
        note = "uncertain"
    else:
        shape = "sharp single peak"
        est_size = row.get("peak1_est_genome_size(len)", "")
        peak_position = row.get("peak1_kcov", "")
        note = "uncertain"

    if peak_position not in ("", "NA"):
        pos_int = to_int_rounded(peak_position)
        if pos_int is not None:
            peak_position = str(pos_int)
            if pos_int < 24:
                est_size = "NA"

    return {
        "Seq_ID": row.get("Seq_ID", ""),
        "Total_number_of_kmer": row.get("total_kmers", ""),
        "peak_position": peak_position,
        "Shape": shape,
        "Estimated_genome_size(bp)": est_size,
        "Unique_kmer_number": row.get("unique_kmers", ""),
        "Kmer_automation_note": note,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract key columns and classify k-mer peak shapes."
    )
    parser.add_argument(
        "input_csv",
        help="Path to statistics_all.csv",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output CSV path (default: kmer_profile_statistics_automated.csv next to input)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_path = args.output
    if output_path is None:
        base_dir = os.path.dirname(os.path.abspath(args.input_csv))
        output_path = os.path.join(base_dir, "kmer_profile_statistics_automated.csv")

    with open(args.input_csv, "r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        with open(output_path, "w", newline="", encoding="utf-8") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=OUTPUT_COLUMNS)
            writer.writeheader()
            for row in reader:
                writer.writerow(build_output_row(row))


if __name__ == "__main__":
    main()
