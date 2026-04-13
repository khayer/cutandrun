// Module: PROMOTER_GC_DIFFBIND
// Purpose: Differential binding and promoter GC comparison for a target contrast

process PROMOTER_GC_DIFFBIND {
    tag "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}"
    label 'process_medium'

    container 'docker://weishwu/chipseeker:1.42.0'

    publishDir "${params.outdir}/03_peak_calling/11_differential_binding/${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}", mode: 'copy'

    input:
    path homer_annotation
    path counts_table

    output:
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.merged_peaks.annotatePeaks.txt", emit: annotation
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.deseq2_results.tsv", emit: results
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_peaks.tsv", emit: promoter_peaks
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_loss.tsv", emit: promoter_loss
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_not_affected.tsv", emit: promoter_not_affected
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_gc_summary.tsv", emit: promoter_gc_summary
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_gc_test.tsv", emit: promoter_gc_test
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_loss.tsv", emit: top_loss
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_gain.tsv", emit: top_gain
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_promoter_gc_summary.tsv", emit: top_promoter_gc_summary
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_promoter_gc_test.tsv", emit: top_promoter_gc_test
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_promoter_gc_boxplot.png", emit: top_promoter_gc_boxplot_png
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_promoter_gc_boxplot.pdf", emit: top_promoter_gc_boxplot_pdf
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_promoter_gc_violin.png", emit: top_promoter_gc_violin_png
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_promoter_gc_violin.pdf", emit: top_promoter_gc_violin_pdf
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_volcano_raw_pvalue.png", emit: top_volcano_raw_pvalue_png
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.top${params.promoter_gc_top_n}_volcano_raw_pvalue.pdf", emit: top_volcano_raw_pvalue_pdf
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_gc_boxplot.png", emit: promoter_gc_boxplot_png
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_gc_boxplot.pdf", emit: promoter_gc_boxplot_pdf
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_gc_violin.png", emit: promoter_gc_violin_png
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.promoter_gc_violin.pdf", emit: promoter_gc_violin_pdf
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.volcano.png", emit: volcano_png
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.volcano.pdf", emit: volcano_pdf
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.ma_plot.png", emit: ma_plot_png
    path "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}.ma_plot.pdf", emit: ma_plot_pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}"
    """
    cp ${homer_annotation} ${prefix}.merged_peaks.annotatePeaks.txt

    Rscript ${workflow.projectDir}/bin/promoter_gc_diffbind.R \
        --counts ${counts_table} \
        --annotation ${prefix}.merged_peaks.annotatePeaks.txt \
        --out-prefix ${prefix} \
        --group-a ${params.promoter_gc_group_a} \
        --group-b ${params.promoter_gc_group_b} \
        --fdr ${params.promoter_gc_fdr} \
        --log2fc-cutoff ${params.promoter_gc_log2fc_cutoff} \
        --promoter-window ${params.promoter_gc_tss_dist} \
        --top-n ${params.promoter_gc_top_n}
    """

    stub:
    def prefix = "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}"
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