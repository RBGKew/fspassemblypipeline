process KMER_STAT_SUMMARY {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81fd170d9d2bd4aa6561ccd1638a79e27895e4884753200237b6695e626794d2/data':
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python:a03b8e12185b7953' }"

    input:
    path(fastk_hists)
    path(fastk_logs)
    path(genescopefk_p1_summaries, stageAs: 'genescopefk_p1_summaries/*')
    path(genescopefk_p2_summaries, stageAs: 'genescopefk_p2_summaries/*')
    path(genescopefk_p1_logs, stageAs: 'genescopefk_p1_logs/*')
    path(genescopefk_p2_logs, stageAs: 'genescopefk_p2_logs/*')

    output:
    path('statistics/statistics_all.csv'), emit: stats_all
    path('statistics/kmer_profile_statistics_automated.csv'), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version"), topic: versions, emit: versions_kmer_stat_summary

    script:
    def fastk_hist_args = (fastk_hists instanceof List ? fastk_hists : [fastk_hists]).collect { "\"${it}\"" }.join(' ')
    def fastk_log_args = (fastk_logs instanceof List ? fastk_logs : [fastk_logs]).collect { "\"${it}\"" }.join(' ')
    def p1_summary_args = (genescopefk_p1_summaries instanceof List ? genescopefk_p1_summaries : [genescopefk_p1_summaries]).collect { "\"${it}\"" }.join(' ')
    def p2_summary_args = (genescopefk_p2_summaries instanceof List ? genescopefk_p2_summaries : [genescopefk_p2_summaries]).collect { "\"${it}\"" }.join(' ')
    def p1_log_args = (genescopefk_p1_logs instanceof List ? genescopefk_p1_logs : [genescopefk_p1_logs]).collect { "\"${it}\"" }.join(' ')
    def p2_log_args = (genescopefk_p2_logs instanceof List ? genescopefk_p2_logs : [genescopefk_p2_logs]).collect { "\"${it}\"" }.join(' ')
    """
    mkdir -p kmer_inputs statistics

    printf "%s\n" ${fastk_hist_args} > fastk_hists.list
    printf "%s\n" ${fastk_log_args} > fastk_logs.list
    printf "%s\n" ${p1_summary_args} > p1_summaries.list
    printf "%s\n" ${p2_summary_args} > p2_summaries.list
    printf "%s\n" ${p1_log_args} > p1_logs.list
    printf "%s\n" ${p2_log_args} > p2_logs.list

    python - <<'PY'
import shutil
from pathlib import Path

def read_items(path):
    p = Path(path)
    if not p.exists():
        return []
    return [Path(line.strip()) for line in p.read_text(encoding='utf-8').splitlines() if line.strip()]

def sample_from_hist(path):
    name = path.name
    if name.endswith('.hist.txt'):
        return name[:-9]
    if name.endswith('.hist'):
        return name[:-5]
    return path.stem

def sample_from_suffix(path, suffix):
    return path.name[:-len(suffix)] if path.name.endswith(suffix) else path.stem

base = Path('kmer_inputs')
base.mkdir(parents=True, exist_ok=True)

for hist in read_items('fastk_hists.list'):
    sample_id = sample_from_hist(hist)
    sample_dir = base / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(hist, sample_dir / f'{sample_id}.reads.kmer_freq.hist')

for log in read_items('fastk_logs.list'):
    sample_id = sample_from_suffix(log, '.fastK.log')
    sample_dir = base / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(log, sample_dir / f'{sample_id}.fastK.log')

for summary in read_items('p1_summaries.list'):
    sample_id = sample_from_suffix(summary, '_summary.txt')
    peak_dir = base / sample_id / 'peak_1'
    peak_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(summary, peak_dir / 'summary.txt')

for summary in read_items('p2_summaries.list'):
    sample_id = sample_from_suffix(summary, '_summary.txt')
    peak_dir = base / sample_id / 'peak_2'
    peak_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(summary, peak_dir / 'summary.txt')

for log in read_items('p1_logs.list'):
    sample_id = sample_from_suffix(log, '_genescopefk.log')
    peak_dir = base / sample_id / 'peak_1'
    peak_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(log, peak_dir / 'genescopefk.log')

for log in read_items('p2_logs.list'):
    sample_id = sample_from_suffix(log, '_genescopefk.log')
    sample_dir = base / sample_id
    peak_dir = sample_dir / 'peak_2'
    peak_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(log, peak_dir / 'genescopefk.log')

for sample_dir in base.iterdir():
    if not sample_dir.is_dir():
        continue
    merged_log = sample_dir / f'{sample_dir.name}.genescopefk.log'
    with merged_log.open('w', encoding='utf-8') as out:
        for candidate in [sample_dir / 'peak_1' / 'genescopefk.log', sample_dir / 'peak_2' / 'genescopefk.log']:
            if candidate.exists():
                out.write(candidate.read_text(encoding='utf-8', errors='ignore'))
PY

    python ${projectDir}/bin/kmer_hist_peak_auto_classification.py \
        -f kmer_inputs \
        -o statistics/kmer_hist_peak_auto_classification.csv

    python ${projectDir}/bin/kmer_stats_collect.py \
        --input-dir kmer_inputs \
        --peak-csv statistics/kmer_hist_peak_auto_classification.csv \
        --out-csv statistics/statistics_all.csv

    python - <<'PY'
import csv
import re
from pathlib import Path

stats_csv = Path('statistics/statistics_all.csv')
rows = []
with stats_csv.open('r', encoding='utf-8', newline='') as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames or []
    for row in reader:
        sample_id = (row.get('Seq_ID') or '').strip()
        if sample_id:
            p1_log = Path('kmer_inputs') / sample_id / 'peak_1' / 'genescopefk.log'
            if p1_log.exists():
                text = p1_log.read_text(encoding='utf-8', errors='ignore')
                match = re.search(r'kcov\s*:\s*([0-9.eE+-]+)', text, flags=re.IGNORECASE)
                if match:
                    token = match.group(1)
                    try:
                        row['peak_positions'] = str(int(round(float(token))))
                    except ValueError:
                        row['peak_positions'] = token
        rows.append(row)

with stats_csv.open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY

    python ${projectDir}/bin/statistics_summary_extract.py \
        statistics/statistics_all.csv \
        -o statistics/kmer_profile_statistics_automated.csv

    rm -f statistics/kmer_hist_peak_auto_classification.csv
    """
}
