// Module: MERGE_PEAKS_TABLE_CONTRAST
// Purpose: Build a merged peak universe for one specific contrast only

process MERGE_PEAKS_TABLE_CONTRAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.8.3 conda-forge::pandas=1.2.3"
    container "quay.io/biocontainers/python:3.8.3"

    publishDir "${params.outdir}/03_peak_calling/11_differential_binding", mode: 'copy', saveAs: { filename -> "${task.tag}/${filename}" }

    input:
    tuple val(meta), path(peak_beds)

    output:
    tuple val(meta), path("merged_peaks_table.txt"), emit: table
    tuple val(meta), path("merged_peaks.bed"), emit: bed
    tuple val(meta), path("${meta.id}.promoter_gc_contrast_qc.tsv"), emit: qc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bed_files = peak_beds instanceof List ? peak_beds.join(' ') : peak_beds
    """
#!/usr/bin/env python3

import os
import sys

peaks_by_sample = {}
all_peaks = []

bed_files = "${bed_files}".split()

for bed_file in bed_files:
    sample = os.path.basename(bed_file).replace('.macs2.peaks.cut.bed', '').replace('.bed', '')
    peaks_by_sample[sample] = []

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\\t')
            if len(fields) < 3:
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            peak = (chrom, start, end)
            peaks_by_sample[sample].append(peak)
            all_peaks.append(peak)

def merge_peaks(peaks):
    if not peaks:
        return []

    sorted_peaks = sorted(peaks, key=lambda x: (x[0], x[1], x[2]))
    merged = [sorted_peaks[0]]

    for current in sorted_peaks[1:]:
        last = merged[-1]
        if current[0] == last[0] and current[1] <= last[2]:
            merged[-1] = (last[0], last[1], max(last[2], current[2]))
        else:
            merged.append(current)

    return merged

merged_peaks = merge_peaks(all_peaks)

def has_overlap(peak, sample_peaks, min_overlap=0.5):
    chrom, start, end = peak
    peak_len = end - start

    for sp in sample_peaks:
        if sp[0] != chrom:
            continue
        overlap_start = max(start, sp[1])
        overlap_end = min(end, sp[2])
        if overlap_start < overlap_end:
            overlap_len = overlap_end - overlap_start
            if overlap_len >= peak_len * min_overlap or overlap_len >= (sp[2] - sp[1]) * min_overlap:
                return True
    return False

with open('merged_peaks.bed', 'w') as out:
    for i, (chrom, start, end) in enumerate(merged_peaks, 1):
        out.write(f"{chrom}\\t{start}\\t{end}\\tpeak_{i}\\t0\\t.\\n")

qc_file = "${meta.id}.promoter_gc_contrast_qc.tsv"
with open(qc_file, 'w') as qc:
    qc.write('contrast_id\\tgroup_a\\tgroup_b\\tn_input_peak_beds\\tn_merged_peaks\\n')
    qc.write(f"${meta.id}\\t${meta.group_a}\\t${meta.group_b}\\t${meta.n_input_peak_beds}\\t{len(merged_peaks)}\\n")

samples = sorted(peaks_by_sample.keys())

with open('merged_peaks_table.txt', 'w') as out:
    header = ['PeakID', 'Chr', 'Start', 'End', 'Length'] + samples + ['Total']
    out.write('\\t'.join(header) + '\\n')

    for i, (chrom, start, end) in enumerate(merged_peaks, 1):
        peak_id = f"peak_{i}"
        length = end - start

        presence = []
        total = 0
        for sample in samples:
            if has_overlap((chrom, start, end), peaks_by_sample[sample]):
                presence.append('1')
                total += 1
            else:
                presence.append('0')

        row = [peak_id, chrom, str(start), str(end), str(length)] + presence + [str(total)]
        out.write('\\t'.join(row) + '\\n')

with open('versions.yml', 'w') as v:
    v.write('"${task.process}":\\n')
    v.write(f'    python: "{sys.version.split()[0]}"\\n')
    """.stripIndent()
}
