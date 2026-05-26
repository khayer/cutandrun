// Module: DESEQ2_PEAK_ANALYSIS
// Purpose: Differential abundance analysis of peaks using DESeq2

process DESEQ2_PEAK_ANALYSIS {
    tag "${meta.treatment}_vs_${meta.control}"
    label 'process_medium'

    conda "bioconda::deseq2=1.44.0 bioconda::bioconductor-genomicranges=1.56.1 bioconda::bioconductor-rtracklayer=1.62.0 conda-forge::r-base=4.4.0 conda-forge::r-tidyverse=2.0.0 conda-forge::r-optparse=1.7.3"
    container "docker://rocker/verse:4.3"

    publishDir "${params.outdir}/03_peak_calling/09_deseq2_analysis/${meta.treatment}_vs_${meta.control}", mode: 'copy'

    input:
    tuple val(meta), path(counts_table), path(narrowpeak_dir)
    path gtf_file
    path samplesheet

    output:
    tuple val(meta), path("*_results.csv"), emit: results
    tuple val(meta), path("*_significant.csv"), emit: significant
    tuple val(meta), path("*_significant_log2fc.bedGraph"), emit: bedgraph
    tuple val(meta), path("*_volcano_plot.pdf"), emit: volcano_plot
    tuple val(meta), path("*_ma_plot.pdf"), emit: ma_plot
    tuple val(meta), path("*_volcano_plot.png"), emit: volcano_plot_png
    tuple val(meta), path("*_ma_plot.png"), emit: ma_plot_png
    path "versions.yml", emit: versions

    when:
    params.run_deseq2_peak_analysis

    script:
    def prefix = "${meta.treatment}_vs_${meta.control}"
    """
    Rscript ${workflow.projectDir}/bin/deseq2_peak_analysis.R \
        --counts-table ${counts_table} \
        --peak-dir ${narrowpeak_dir} \
        --gtf-file ${gtf_file} \
        --samplesheet ${samplesheet} \
        --treatment-group ${meta.treatment} \
        --control-group ${meta.control} \
        --output-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: \$(Rscript -e "library(DESeq2); cat(packageVersion('DESeq2'))")
        R: \$(R --version | head -n 1 | sed 's/R version //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.treatment}_vs_${meta.control}"
    """
    touch ${prefix}_results.csv
    touch ${prefix}_significant.csv
    touch ${prefix}_significant_log2fc.bedGraph
    touch ${prefix}_volcano_plot.pdf
    touch ${prefix}_ma_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deseq2: 1.44.0
        R: 4.4.0
    END_VERSIONS
    """
}
