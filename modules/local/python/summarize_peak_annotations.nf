// Module: SUMMARIZE_PEAK_ANNOTATIONS
// Author: Katharina Hayer
// Co-created with: GitHub Copilot (GPT-5.3-Codex)
// Purpose: Summarize HOMER annotatePeaks outputs and generate replicate stacked bar plots

process SUMMARIZE_PEAK_ANNOTATIONS {
    tag "summarize_peak_annotations"
    label 'process_single'

    conda "bioconda::multiqc=1.33"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.33--pyhdfd78af_0' :
        'biocontainers/multiqc:1.33--pyhdfd78af_0' }"

    input:
    path annotate_tables

    output:
    path "homer_peak_annotation.raw_annotation_summary.tsv", emit: raw_summary
    path "homer_peak_annotation.feature_summary.tsv", emit: feature_summary
    path "homer_peak_annotation.sample_stats.tsv", emit: sample_stats
    path "homer_peak_annotation.feature_percent_table.tsv", emit: percent_table
    path "homer_peak_annotation.stacked_bar.png", emit: plot_png
    path "homer_peak_annotation.stacked_bar.pdf", emit: plot_pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${workflow.projectDir}/bin/summarize_peak_annotations.py \
        --inputs *.annotatePeaks.txt \
        --out-prefix homer_peak_annotation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -E -o "([0-9]{1,}\\.)+[0-9]{1,}")
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """
}
