#!/bin/bash
# Author: Wu Huang <w.huang@kew.org>
# Date: 2025-11-05
# Description: Compile FastQC results into a single file.

shopt -s nullglob

input_files=("$@")
fastqc_data_files=()
summary_files=()

for f in "${input_files[@]}"; do
	[[ -e "$f" ]] || continue
	if [[ "$f" == *_fastqc_data.txt || "$f" == *_data.txt ]]; then
		fastqc_data_files+=("$f")
	elif [[ "$f" == *_summary.txt || "$f" == *gz_summary.txt ]]; then
		summary_files+=("$f")
	fi
done

if [[ ${#fastqc_data_files[@]} -eq 0 ]]; then
	for i in */; do
		[[ -d "$i" ]] || continue
		for f in "$i"/*fastqc_data.txt; do
			[[ -e "$f" ]] || continue
			fastqc_data_files+=("$f")
		done
		for f in "$i"/*_data.txt; do
			[[ -e "$f" ]] || continue
			fastqc_data_files+=("$f")
		done
	done
fi

if [[ ${#fastqc_data_files[@]} -eq 0 ]]; then
	echo "ERROR: No *_fastqc_data.txt or *_data.txt files found." >&2
	exit 1
fi

if [[ ${#summary_files[@]} -eq 0 ]]; then
	summary_files=( *_summary.txt )
fi

if [[ ${#summary_files[@]} -eq 0 ]]; then
	for i in */; do
		[[ -d "$i" ]] || continue
		for f in "$i"/*_summary.txt; do
			[[ -e "$f" ]] || continue
			summary_files+=("$f")
		done
	done
fi

if [[ ${#summary_files[@]} -eq 0 ]]; then
	echo "ERROR: No summary files found." >&2
	exit 1
fi

grep -h 'Filename' "${fastqc_data_files[@]}" | awk '{split($2,a,".");print a[1]}' > name.txt
grep -h '%GC' "${fastqc_data_files[@]}" | awk '{print $2}' > GC.txt
grep -h 'Total Sequences' "${fastqc_data_files[@]}" | awk '{print $3}' > total_reads_no.txt
grep -h 'Total Deduplicated Percentage' "${fastqc_data_files[@]}" | awk '{print 100-$4}' > duplication_level.txt
grep -h 'Per sequence GC content' "${summary_files[@]}" | awk '{print $1}' > GC_pass.txt

echo "name	GC	GC_pass	total_reads_no	Duplicated Percentage" > fastQC_result.txt
paste name.txt GC.txt GC_pass.txt total_reads_no.txt duplication_level.txt >> fastQC_result.txt
rm name.txt GC.txt GC_pass.txt total_reads_no.txt duplication_level.txt
