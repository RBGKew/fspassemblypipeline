#!/bin/bash
# Author: @LiaOb21
# Date: 28-04-2026
# Description: This script selects the best assembly based on BUSCO scores with auN as tiebreaker
# Adapted to nextflow from the original version: https://github.com/LiaOb21/FSP_assembly_benchmarking/blob/main/workflow/scripts/select_best_assembly.sh

# Function to show help
show_help() {
    cat << EOF
Usage: $0 -s <sample> -b <busco_files> -q <quast_tsv> -a <assembly_files> -o <output_directory>

This script finds the best assembly for each sample based on:
1. Highest complete and single-copy BUSCOs
2. If multiple assemblies have the highest BUSCOs, the best assembly is chosen among those with the best BUSCOs based on the auN value from QUAST

Options:
  -s <meta.id>  Sample name
  -b <busco_files>   a space-separated list of BUSCO summary .txt files
  -q <quast_tsv>   quast .tsv file
  -a <assembly_files> space-separated list of assembly fasta files
  -h          Show this help message

Example:
  $0 -s sample1 -b "summary1.txt summary2.txt" -q report.tsv -a "asm1.fa.gz asm2.fa.gz" -o results/best_assemblies

EOF
}

# Initialize variables
sample=""
busco_files=""
quast_tsv=""
assembly_files=""
output_dir="."

# Parse command line arguments
while getopts "s:b:q:a:o:h" opt; do
    case $opt in
        s) sample="$OPTARG"
           echo "Sample: $OPTARG"
        ;;
        b) busco_files="$OPTARG"
           echo "BUSCO files: $OPTARG"
        ;;
        q) quast_tsv="$OPTARG"
           echo "QUAST TSV file: $OPTARG"
        ;;
        a) assembly_files="$OPTARG"
           echo "Assembly files: $OPTARG"
        ;;
        o) output_dir="$OPTARG"
           echo "Output directory: $OPTARG"
        ;;
        h) show_help
           exit 0
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
            show_help
            exit 1
        ;;
    esac
done

# Check required arguments
if [[ -z "$sample" ]]; then
    echo "Error: Sample name (-s) is required!" >&2
    show_help
    exit 1
fi

if [[ -z "$busco_files" ]]; then
    echo "Error: BUSCO files (-b) are required!" >&2
    show_help
    exit 1
fi

if [[ -z "$quast_tsv" ]]; then
    echo "Error: QUAST TSV file (-q) is required!" >&2
    show_help
    exit 1
fi

if [[ ! -f "$quast_tsv" ]]; then
    echo "Error: QUAST TSV file '$quast_tsv' does not exist!" >&2
    exit 1
fi

echo "processing sample: $sample"

# Extract BUSCO complete and single-copy for each assembly across all strategies and read types.
# Avoid original directory walk and use the provided list of busco_files instead.
# The busco_files should be in the format: sample_readstype_kmerstrategy_assembler_short_summary.txt

echo "=== Extracting BUSCO scores from provided summary files ==="

# The assembly label is derived from the filename, e.g.:
# short_summary.specific.fungi_odb12.samplename_R1R2_kmergenie_spades.fa.txt
# strip leading fields up to and including the sample prefix, and the trailing .fa.txt
# yields: R1R2_kmergenie_spades

for file in $busco_files; do
    if [[ -f "$file" ]]; then
        # Extract basename, strip extension .fa.txt, then strip the sample prefix
        basename_noext=$(basename "$file" .fa.txt)                         # e.g. short_summary.specific.fungi_odb12.samplename_R1R2_kmergenie_spades
        assembly_full=$(echo "$basename_noext" | awk -F'.' '{print $NF}')  # last dot-field: samplename_R1R2_kmergenie_spades
        assembly_label="${assembly_full#${sample}_}"                       # strip sample prefix: R1R2_kmergenie_spades

        busco_score=$(awk '/Complete and single-copy/ {print $1}' "$file")
        echo "  Found: ${assembly_label} (BUSCO: $busco_score)" >&2
        echo "${assembly_label} ${busco_score}" >> "$output_dir/complete_single_copy_buscos.txt"
    else
        echo "  Warning: BUSCO file not found: $file" >&2
    fi
done

if [[ ! -f "$output_dir/complete_single_copy_buscos.txt" ]]; then
    echo "  No BUSCO data found for sample $sample"
    exit 1
fi

echo ""
echo "=== All assemblies with BUSCO scores ==="
cat "$output_dir/complete_single_copy_buscos.txt"
echo ""

# Sort the results by the number of complete and single-copy BUSCOs in descending order
# extract max value from the first line
# print all lines with that max value to best_buscos.txt
sort -k2,2nr "$output_dir/complete_single_copy_buscos.txt" | awk 'NR==1{max=$2} $2==max' > "$output_dir/best_buscos.txt"

echo "=== Assemblies with highest BUSCO scores ==="
cat "$output_dir/best_buscos.txt"
echo ""

# Check if only one assembly has the highest BUSCO score
best_busco_count=$(wc -l < "$output_dir/best_buscos.txt")

if [[ $best_busco_count -eq 1 ]]; then
    # Single winner - print the assembly name
    best_assembly=$(cut -d' ' -f1 "$output_dir/best_buscos.txt")
    echo "$best_assembly is the best assembly for sample $sample based on BUSCO score."

    # Find the matching assembly file from the provided list
    best_assembly_file=""
    for asm in $assembly_files; do
        if [[ "$(basename "$asm" .fa.gz)" == "${sample}_${best_assembly}" ]]; then
            best_assembly_file="$asm"
            break
        fi
    done
    if [[ -z "$best_assembly_file" ]]; then
        echo "Error: Could not find assembly file for ${sample}_${best_assembly}" >&2
        exit 1
    fi
    echo "copying ${best_assembly_file} to $output_dir/${sample}_best_assembly.fa.gz"
    cp "$best_assembly_file" "$output_dir/${sample}_best_assembly.fa.gz"

    echo "$best_assembly" > "$output_dir/best_assembly.txt"

    # Extract reads_type, kmer_strategy, and assembler from the best assembly name to be added back to meta later
    reads_type=$(echo "$best_assembly" | cut -d'_' -f1)
    assembler=$(echo "$best_assembly" | awk -F'_' '{print $NF}')
    kmer_strategy=$(echo "$best_assembly" | cut -d'_' -f2- | rev | cut -d'_' -f2- | rev)
    printf "reads_type=%s\nkmer_strategy=%s\nassembler=%s\n" \
        "$reads_type" "$kmer_strategy" "$assembler" \
        > "$output_dir/best_assembly_meta.txt"

    exit 0
fi

if [[ -f "$quast_tsv" ]]; then
    echo "  Found QUAST TSV at: $quast_tsv" >&2
    # Parse tab-separated report.tsv - same structure as report.txt but TSV
    # Assembly names format: sample_reads_type_strategy_assembler
    awk -F'\t' '
        /^Assembly/ {
            for(i=2;i<=NF;i++) {
                # Remove .fa.gz extension
                gsub(/\.fa\.gz$/, "", $i);

                # Extract everything starting from R1R2_ or merged_
                # Also handles strategies with underscores like reads_length
                if (match($i, /(R1R2|merged)_.+/)) {
                    extracted = substr($i, RSTART);
                    a[i-1] = extracted;
                }
            }
        }
        /^auN/ {
            for(i=2;i<=NF;i++) {
                if (a[i-1] != "") {
                    print a[i-1], $i
                }
            }
        }
    ' "$quast_tsv" >> "$output_dir/auN_quast.txt"
fi

# Get auN scores only for the tied BUSCO assemblers
# Get the list of assemblers with the highest BUSCOs
tied_assemblers=$(cut -d' ' -f1 "$output_dir/best_buscos.txt")
# set up variables to track best auN and corresponding assembly
best_aun=0
best_assembly=""

# for each assembler in tied_assemblers list, get its auN score and compare
# if higher than current best_aun, update best_aun and best_assembly
while read assembler; do
    aun_score=$(grep "^$assembler " "$output_dir/auN_quast.txt" | cut -d' ' -f2)
    if [[ -n "$aun_score" ]] && awk "BEGIN {exit ($aun_score > $best_aun) ? 0 : 1}"; then
        best_aun=$aun_score
        best_assembly=$assembler
    fi
done <<< "$tied_assemblers"

echo ""
echo "=== auN scores for tied assemblies ==="
grep -F -f <(cut -d' ' -f1 "$output_dir/best_buscos.txt") "$output_dir/auN_quast.txt"
echo ""

if [[ -n "$best_assembly" ]]; then
    echo "$best_assembly is the best assembly for sample $sample based on auN score." >&2
    echo "  Best auN: $best_aun" >&2

    # Find the matching assembly file from the provided list
    best_assembly_file=""
    for asm in $assembly_files; do
        if [[ "$(basename "$asm" .fa.gz)" == "${sample}_${best_assembly}" ]]; then
            best_assembly_file="$asm"
            break
        fi
    done
    if [[ -z "$best_assembly_file" ]]; then
        echo "Error: Could not find assembly file for ${sample}_${best_assembly}" >&2
        exit 1
    fi
    echo "copying ${best_assembly_file} to $output_dir/${sample}_best_assembly.fa.gz" >&2
    cp "$best_assembly_file" "$output_dir/${sample}_best_assembly.fa.gz"

    echo "$best_assembly" > "$output_dir/best_assembly.txt"

    # Extract reads_type, kmer_strategy, and assembler from the best assembly name to be added back to meta later
    reads_type=$(echo "$best_assembly" | cut -d'_' -f1)
    assembler=$(echo "$best_assembly" | awk -F'_' '{print $NF}')
    kmer_strategy=$(echo "$best_assembly" | cut -d'_' -f2- | rev | cut -d'_' -f2- | rev)
    printf "reads_type=%s\nkmer_strategy=%s\nassembler=%s\n" \
        "$reads_type" "$kmer_strategy" "$assembler" \
        > "$output_dir/best_assembly_meta.txt"

    exit 0
fi

echo "Done! Results saved in $output_dir"
