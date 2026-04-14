// Module: PROMOTER_GC_CONTRAST_QC
// Purpose: Report per-contrast peak-universe QC stats

process PROMOTER_GC_CONTRAST_QC {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "quay.io/biocontainers/python:3.8.3"

    publishDir "${params.outdir}/03_peak_calling/11_differential_binding", mode: 'copy', saveAs: { filename -> "${task.tag}/${filename}" }

    input:
    tuple val(meta), path(merged_bed), val(n_input_peak_beds)

    output:
    tuple val(meta), path("${meta.id}.promoter_gc_contrast_qc.tsv"), emit: qc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
import csv

meta_id = "${meta.id}"
group_a = "${meta.group_a}"
group_b = "${meta.group_b}"
input_beds = int("${n_input_peak_beds}")
merged_bed = "${merged_bed}"
out_tsv = f"{meta_id}.promoter_gc_contrast_qc.tsv"

merged_count = 0
with open(merged_bed, 'r') as fh:
    for line in fh:
        if line.strip() and not line.startswith('#'):
            merged_count += 1

with open(out_tsv, 'w', newline='') as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow([
        'contrast_id',
        'group_a',
        'group_b',
        'n_input_peak_beds',
        'n_merged_peaks'
    ])
    writer.writerow([
        meta_id,
        group_a,
        group_b,
        input_beds,
        merged_count
    ])
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "$(python --version | sed 's/Python //')"
    END_VERSIONS
    """
}
