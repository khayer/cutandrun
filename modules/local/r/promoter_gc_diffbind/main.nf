// Module: PROMOTER_GC_DIFFBIND
// Purpose: Differential binding and promoter GC comparison for a target contrast

process PROMOTER_GC_DIFFBIND {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker://weishwu/chipseeker:1.42.0'

    publishDir "${params.outdir}/03_peak_calling/11_differential_binding", mode: 'copy', saveAs: { filename -> "${task.tag}/${filename}" }

    input:
    tuple val(meta), path(homer_annotation), path(counts_table)

    output:
    tuple val(meta), path("*.merged_peaks.annotatePeaks.txt"), emit: annotation
    tuple val(meta), path("*.deseq2_results.tsv"), emit: results
    tuple val(meta), path("*.promoter_peaks.tsv"), emit: promoter_peaks
    tuple val(meta), path("*.promoter_loss.tsv"), emit: promoter_loss
    tuple val(meta), path("*.promoter_not_affected.tsv"), emit: promoter_not_affected
    tuple val(meta), path("*.promoter_gc_summary.tsv"), emit: promoter_gc_summary
    tuple val(meta), path("*.promoter_gc_test.tsv"), emit: promoter_gc_test
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_loss.tsv"), emit: top_loss
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_gain.tsv"), emit: top_gain
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_promoter_gc_summary.tsv"), emit: top_promoter_gc_summary
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_promoter_gc_test.tsv"), emit: top_promoter_gc_test
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_promoter_gc_boxplot.png"), emit: top_promoter_gc_boxplot_png
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_promoter_gc_boxplot.pdf"), emit: top_promoter_gc_boxplot_pdf
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_promoter_gc_violin.png"), emit: top_promoter_gc_violin_png
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_promoter_gc_violin.pdf"), emit: top_promoter_gc_violin_pdf
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_volcano_raw_pvalue.png"), emit: top_volcano_raw_pvalue_png
    tuple val(meta), path("*.top${params.promoter_gc_top_n}_volcano_raw_pvalue.pdf"), emit: top_volcano_raw_pvalue_pdf
    tuple val(meta), path("*.promoter_gc_boxplot.png"), emit: promoter_gc_boxplot_png
    tuple val(meta), path("*.promoter_gc_boxplot.pdf"), emit: promoter_gc_boxplot_pdf
    tuple val(meta), path("*.promoter_gc_violin.png"), emit: promoter_gc_violin_png
    tuple val(meta), path("*.promoter_gc_violin.pdf"), emit: promoter_gc_violin_pdf
    tuple val(meta), path("*.volcano.png"), emit: volcano_png
    tuple val(meta), path("*.volcano.pdf"), emit: volcano_pdf
    tuple val(meta), path("*.ma_plot.png"), emit: ma_plot_png
    tuple val(meta), path("*.ma_plot.pdf"), emit: ma_plot_pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.group_a}_vs_${meta.group_b}"
    """
    cp ${homer_annotation} ${prefix}.merged_peaks.annotatePeaks.txt

    Rscript ${workflow.projectDir}/bin/promoter_gc_diffbind.R \
        --counts ${counts_table} \
        --annotation ${prefix}.merged_peaks.annotatePeaks.txt \
        --out-prefix ${prefix} \
        --group-a ${meta.group_a} \
        --group-b ${meta.group_b} \
        --fdr ${params.promoter_gc_fdr} \
        --log2fc-cutoff ${params.promoter_gc_log2fc_cutoff} \
        --promoter-window ${params.promoter_gc_tss_dist} \
        --top-n ${params.promoter_gc_top_n}
    """

    stub:
    def prefix = "${meta.group_a}_vs_${meta.group_b}"
    """
    touch ${prefix}.merged_peaks.annotatePeaks.txt
    touch ${prefix}.deseq2_results.tsv
    touch ${prefix}.promoter_peaks.tsv
    touch ${prefix}.promoter_loss.tsv
    touch ${prefix}.promoter_not_affected.tsv
    touch ${prefix}.promoter_gc_summary.tsv
    touch ${prefix}.promoter_gc_test.tsv
    touch ${prefix}.top${params.promoter_gc_top_n}_loss.tsv
    touch ${prefix}.top${params.promoter_gc_top_n}_gain.tsv
    touch ${prefix}.top${params.promoter_gc_top_n}_promoter_gc_summary.tsv
    touch ${prefix}.top${params.promoter_gc_top_n}_promoter_gc_test.tsv
    touch ${prefix}.top${params.promoter_gc_top_n}_promoter_gc_boxplot.png
    touch ${prefix}.top${params.promoter_gc_top_n}_promoter_gc_boxplot.pdf
    touch ${prefix}.top${params.promoter_gc_top_n}_promoter_gc_violin.png
    touch ${prefix}.top${params.promoter_gc_top_n}_promoter_gc_violin.pdf
    touch ${prefix}.top${params.promoter_gc_top_n}_volcano_raw_pvalue.png
    touch ${prefix}.top${params.promoter_gc_top_n}_volcano_raw_pvalue.pdf
    touch ${prefix}.promoter_gc_boxplot.png
    touch ${prefix}.promoter_gc_boxplot.pdf
    touch ${prefix}.promoter_gc_violin.png
    touch ${prefix}.promoter_gc_violin.pdf
    touch ${prefix}.volcano.png
    touch ${prefix}.volcano.pdf
    touch ${prefix}.ma_plot.png
    touch ${prefix}.ma_plot.pdf
    echo "stub" > versions.yml
    """
}