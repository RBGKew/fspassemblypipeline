process KMER_STAT_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/python_pip_numpy_pandas:0f1c9045fe326238' }"

    input:
    path(fastk_hists)
    path(genomescope_summaries)

    output:
    path('statistics/statistics_all.csv'), emit: stats_all
    path('statistics/kmer_profile_statistics_automated.csv'), emit: summary

    script:
    """
    mkdir -p kmer_inputs statistics

    for hist in ${fastk_hists}; do
        sample_id=\$(basename "\$hist" .hist)
        sample_dir="kmer_inputs/\$sample_id"
        mkdir -p "\$sample_dir"
        cp "\$hist" "\$sample_dir/\$sample_id.reads.kmer_freq.hist"
    done

    for summary in ${genomescope_summaries}; do
        sample_id=\$(basename "\$summary" _summary.txt)
        sample_dir="kmer_inputs/\$sample_id/peak_1"
        mkdir -p "\$sample_dir"
        cp "\$summary" "\$sample_dir/summary.txt"
    done

    python ${projectDir}/bin/kmer_hist_peak_auto_classification.py \
        -f kmer_inputs \
        -o statistics/kmer_hist_peak_auto_classification.csv

    python ${projectDir}/bin/kmer_stats_collect.py \
        --input-dir kmer_inputs \
        --peak-csv statistics/kmer_hist_peak_auto_classification.csv \
        --out-csv statistics/statistics_all.csv

    python ${projectDir}/bin/statistics_summary_extract.py \
        statistics/statistics_all.csv \
        -o statistics/kmer_profile_statistics_automated.csv

    rm -f statistics/kmer_hist_peak_auto_classification.csv
    """
}
