#!/usr/bin/env python3
# Author: Wu Huang <w.huang@kew.org>
# Date: 2026-03-11
# Description: Collect per-sample k-mer statistics into one CSV table.

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect per-sample k-mer statistics into one CSV table."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing per-sample kmer result folders.",
    )
    parser.add_argument(
        "--peak-csv",
        required=True,
        help="CSV produced by kmer_hist_peak_auto_classification.py",
    )
    parser.add_argument(
        "--out-csv",
        required=True,
        help="Output CSV path (statistics_all.csv)",
    )
    return parser.parse_args()


def load_peak_calls(peak_csv: Path) -> dict[str, dict[str, str]]:
    calls: dict[str, dict[str, str]] = {}
    with peak_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                continue
            calls[sample_id] = {
                "peak_count": row.get("peak_count", ""),
                "peak_positions": row.get("peak_positions", ""),
            }
    return calls


def parse_unique_kmers_from_hist(hist_file: Path) -> int:
    unique_kmers = 0

    with hist_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                count = int(float(parts[1]))
            except ValueError:
                continue
            unique_kmers += count

    return unique_kmers


def parse_total_kmers_from_fastk_log(log_file: Path) -> int | None:
    if not log_file.exists():
        return None

    with log_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = re.search(r"\bSum:\s*([0-9][0-9,]*)", line, flags=re.IGNORECASE)
            if match:
                try:
                    return int(match.group(1).replace(",", ""))
                except ValueError:
                    continue

    return None


def _normalize_numeric_token(value: str) -> str:
    token = value.strip().replace(",", "")
    if token.lower() in {"", "inf", "+inf", "-inf", "nan"}:
        return ""
    try:
        float(token)
    except ValueError:
        return ""
    return token


def _prefer_second_numeric(first: str, second: str) -> str:
    second_clean = _normalize_numeric_token(second)
    if second_clean:
        return second_clean
    return _normalize_numeric_token(first)


def parse_genescopefk_summary(summary_file: Path) -> dict[str, str]:
    if not summary_file.exists():
        return {
            "model_fit_max": "",
            "est_genome_size(len)": "",
        }

    summary_text = summary_file.read_text(encoding="utf-8", errors="ignore")
    model_fit = ""
    est_len = ""

    fit_match = re.search(
        r"Model\s*Fit\s+([0-9.eE+\-,]+)%\s+([0-9.eE+\-,]+)%",
        summary_text,
        flags=re.IGNORECASE,
    )
    if fit_match:
        model_fit = _prefer_second_numeric(fit_match.group(1), fit_match.group(2))

    len_match = re.search(
        r"Genome\s+Haploid\s+Length\s+([0-9.eE+\-,]+|Inf)\s*bp\s+([0-9.eE+\-,]+|Inf)\s*bp",
        summary_text,
        flags=re.IGNORECASE,
    )
    if len_match:
        est_len = _prefer_second_numeric(len_match.group(1), len_match.group(2))

    return {
        "model_fit_max": model_fit,
        "est_genome_size(len)": est_len,
    }


def parse_genescopefk_log(log_file: Path) -> dict[int, dict[str, str]]:
    peak_stats: dict[int, dict[str, str]] = {}
    if not log_file.exists():
        return peak_stats

    current_peak: int | None = None
    with log_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.lower().startswith("genomescope analyzing"):
                peak_match = re.search(r"\bp\s*=\s*(\d+)\b", line, flags=re.IGNORECASE)
                if peak_match:
                    current_peak = int(peak_match.group(1))
                    peak_stats.setdefault(
                        current_peak,
                        {
                            "model_fit_max": "",
                            "est_genome_size(len)": "",
                            "kcov": "",
                        },
                    )
                else:
                    current_peak = None
                continue

            if current_peak is None:
                continue

            kcov_match = re.search(r"\bkcov\s*:\s*([0-9.eE+-]+)", line, flags=re.IGNORECASE)
            if kcov_match:
                peak_stats[current_peak]["kcov"] = kcov_match.group(1)

            fit_match = re.search(r"\bmodel\s*fit\s*:\s*([0-9.eE+-]+)", line, flags=re.IGNORECASE)
            if fit_match:
                peak_stats[current_peak]["model_fit_max"] = fit_match.group(1)

            len_match = re.search(r"\blen\s*:\s*([0-9.eE+-]+)", line, flags=re.IGNORECASE)
            if len_match:
                peak_stats[current_peak]["est_genome_size(len)"] = len_match.group(1)

    return peak_stats


def collect_peak_stats(summary_file: Path, log_peak_stats: dict[str, str]) -> dict[str, str]:
    summary_stats = parse_genescopefk_summary(summary_file)
    return {
        "model_fit_max": summary_stats["model_fit_max"] or log_peak_stats.get("model_fit_max", ""),
        "est_genome_size(len)": summary_stats["est_genome_size(len)"] or log_peak_stats.get("est_genome_size(len)", ""),
        "kcov": log_peak_stats.get("kcov", ""),
    }


def collect_rows(input_dir: Path, peak_calls: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for sample_dir in sorted(path for path in input_dir.iterdir() if path.is_dir()):
        sample_id = sample_dir.name
        hist_file = sample_dir / f"{sample_id}.reads.kmer_freq.hist"
        fastk_log_file = sample_dir / f"{sample_id}.fastK.log"
        genescopefk_log_file = sample_dir / f"{sample_id}.genescopefk.log"
        if not hist_file.exists():
            continue

        unique_kmers = parse_unique_kmers_from_hist(hist_file)
        total_kmers = parse_total_kmers_from_fastk_log(fastk_log_file)
        genescopefk_log_stats = parse_genescopefk_log(genescopefk_log_file)

        peak1 = collect_peak_stats(sample_dir / "peak_1" / "summary.txt", genescopefk_log_stats.get(1, {}))
        peak2 = collect_peak_stats(sample_dir / "peak_2" / "summary.txt", genescopefk_log_stats.get(2, {}))
        peak = peak_calls.get(sample_id, {})

        rows.append(
            {
                "Seq_ID": sample_id,
                "total_kmers": "" if total_kmers is None else str(total_kmers),
                "unique_kmers": str(unique_kmers),
                "peak_count": str(peak.get("peak_count", "")),
                "peak_positions": str(peak.get("peak_positions", "")),
                "peak1_model_fit_max": peak1["model_fit_max"],
                "peak1_est_genome_size(len)": peak1["est_genome_size(len)"],
                "peak1_kcov": peak1["kcov"],
                "peak2_model_fit_max": peak2["model_fit_max"],
                "peak2_est_genome_size(len)": peak2["est_genome_size(len)"],
                "peak2_kcov": peak2["kcov"],
            }
        )

    return rows


def write_rows(rows: list[dict[str, str]], out_csv: Path) -> None:
    fieldnames = [
        "Seq_ID",
        "total_kmers",
        "unique_kmers",
        "peak_count",
        "peak_positions",
        "peak1_model_fit_max",
        "peak1_est_genome_size(len)",
        "peak1_kcov",
        "peak2_model_fit_max",
        "peak2_est_genome_size(len)",
        "peak2_kcov",
    ]

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir).expanduser().resolve()
    peak_csv = Path(args.peak_csv).expanduser().resolve()
    out_csv = Path(args.out_csv).expanduser().resolve()

    if not input_dir.exists():
        raise FileNotFoundError(f"Input dir not found: {input_dir}")
    if not peak_csv.exists():
        raise FileNotFoundError(f"Peak CSV not found: {peak_csv}")

    peak_calls = load_peak_calls(peak_csv)
    rows = collect_rows(input_dir, peak_calls)
    write_rows(rows, out_csv)


if __name__ == "__main__":
    main()
