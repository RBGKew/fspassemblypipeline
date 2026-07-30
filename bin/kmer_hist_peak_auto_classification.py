#!/usr/bin/env python3
# Author: Wu Huang <w.huang@kew.org>
# Date: 2026-03-11
# Description: Automatically classify k-mer histogram peaks.

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import argparse


parser = argparse.ArgumentParser(
	description=(
		"Detect peaks in k-mer histograms and output key columns. "
		"Required arguments: -f and -o."
	)
)
parser.add_argument(
	"-f",
	"--folder",
	required=True,
	help="Absolute path to the folder which contains the hist files.",
)
parser.add_argument(
	"-o",
	"--output",
	required=True,
	help="Output CSV file for peak_count and peak_positions.",
)
args = parser.parse_args()

SAMPLES_DIR = Path(args.folder).expanduser().resolve()
OUTPUT_CSV = Path(args.output).expanduser()

K_MIN = 3
K_MAX = 400
SMOOTH_WINDOW = 7
MIN_PROMINENCE_LOW = 0.08
MIN_PROMINENCE_HIGH = 0.25
MIN_DISTANCE = 4
PEAK_RATIO_MIN = 1.6
PEAK_RATIO_MAX = 2.4
MAX_PEAKS = 2
LOW_RANGE_MAX = 100
MIN_PEAK_K = 15


def smooth_series(values: np.ndarray, window: int) -> np.ndarray:
	if window < 3:
		return values
	if window % 2 == 0:
		window += 1
	if len(values) < window:
		return values
	pad = window // 2
	padded = np.pad(values, (pad, pad), mode="edge")
	kernel = np.ones(window, dtype=float) / window
	return np.convolve(padded, kernel, mode="valid")


def find_peak_candidates(x: np.ndarray, y: np.ndarray) -> list[tuple[int, float]]:
	if len(y) < 3:
		return []
	candidates: list[tuple[int, float]] = []
	for i in range(1, len(y) - 1):
		if y[i] > y[i - 1] and y[i] >= y[i + 1]:
			if int(x[i]) < MIN_PEAK_K:
				continue
			left_min = float(np.min(y[:i])) if i > 0 else float(y[i])
			right_min = float(np.min(y[i + 1 :])) if i + 1 < len(y) else float(y[i])
			prominence = float(y[i] - max(left_min, right_min))
			min_prom = MIN_PROMINENCE_LOW if int(x[i]) <= LOW_RANGE_MAX else MIN_PROMINENCE_HIGH
			if prominence >= min_prom:
				candidates.append((i, prominence))
	return candidates


def select_peaks(
	candidates: Iterable[tuple[int, float]],
	x: np.ndarray,
	min_distance: int = MIN_DISTANCE,
	max_peaks: int = MAX_PEAKS,
) -> list[tuple[int, float]]:
	ordered = sorted(candidates, key=lambda item: item[1], reverse=True)
	selected: list[tuple[int, float]] = []
	for idx, prom in ordered:
		if all(abs(int(x[idx]) - int(x[chosen_idx])) >= min_distance for chosen_idx, _ in selected):
			selected.append((idx, prom))
		if len(selected) >= max_peaks:
			break
	selected.sort(key=lambda item: int(x[item[0]]))
	return selected


def classify_histogram(hist_path: Path) -> tuple[int, list[int]]:
	data = pd.read_csv(
		hist_path,
		sep=r"\s+",
		header=None,
		names=["k", "count"],
		usecols=[0, 1],
		comment=">",
	)
	data["k"] = pd.to_numeric(data["k"], errors="coerce")
	data["count"] = pd.to_numeric(data["count"], errors="coerce")
	data = data.dropna(subset=["k", "count"])
	data = data[(data["k"] >= K_MIN) & (data["k"] <= K_MAX)].copy()
	if data.empty:
		return 0, []

	x = data["k"].to_numpy(dtype=int)
	y = data["count"].to_numpy(dtype=float)
	weighted_y = x * y
	log_y = np.log10(weighted_y + 1.0)
	smooth_y = smooth_series(log_y, SMOOTH_WINDOW)

	candidates = find_peak_candidates(x, smooth_y)
	selected = select_peaks(candidates, x)
	if len(selected) >= 2:
		p1_idx, p1_prom = selected[0]
		p2_idx, p2_prom = selected[1]
		p1 = int(x[p1_idx])
		p2 = int(x[p2_idx])
		if p1 > 0:
			ratio = max(p1, p2) / min(p1, p2)
		else:
			ratio = 0.0
		if ratio < PEAK_RATIO_MIN or ratio > PEAK_RATIO_MAX:
			# Keep the stronger peak only when ratio is not diploid-like.
			selected = [(p1_idx, p1_prom)] if p1_prom >= p2_prom else [(p2_idx, p2_prom)]

	peak_positions = [int(x[idx]) for idx, _ in selected]
	peak_count = len(peak_positions)

	return peak_count, peak_positions


def main() -> None:
	if not SAMPLES_DIR.exists():
		raise FileNotFoundError(f"Input folder not found: {SAMPLES_DIR}")
	hist_files = sorted(SAMPLES_DIR.glob("*/*.reads.kmer_freq.hist"))
	if not hist_files:
		raise FileNotFoundError(f"No *.reads.kmer_freq.hist files found under: {SAMPLES_DIR}")

	rows = []
	for hist_path in hist_files:
		sample_id = hist_path.parent.name
		peak_count, peak_positions_list = classify_histogram(hist_path)
		peak_positions = ";".join(str(p) for p in peak_positions_list)

		rows.append(
			{
				"sample_id": sample_id,
				"peak_count": peak_count,
				"peak_positions": peak_positions,
			}
		)

	output_df = pd.DataFrame(rows)
	if OUTPUT_CSV.parent != Path("."):
		OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
	output_df.to_csv(OUTPUT_CSV, index=False)


if __name__ == "__main__":
	main()
