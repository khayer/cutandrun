/*
 * Module: VIRAL_WINDOW_DB
 * Authors: Katharina E. Hayer, co-created with Claude
 * Version: 1.0.0  |  Date: 2026-11-25
 *
 * Differential coverage on a viral contig with replicate dispersion, plus a
 * permutation null that is always reported. Needs edgeR; the container is a
 * parameter because this repo pins no bioconductor image and guessing one is
 * how a fabricated tag gets into a pipeline.
 */
process VIRAL_WINDOW_DB {
    tag "$contig"
    label 'process_low'

    container "${ params.viral_db_container ?: 'CONTAINER_NOT_SET' }"

    input:
    tuple val(contig), path(counts), path(samples)
    path contrasts
    val  fdr
    val  merge_gap
    val  min_count
    val  norm_mode

    output:
    tuple val(contig), path("${contig}.viral_db*"), emit: results
    path "${contig}.viral_db.summary.txt"         , emit: summary

    script:
    def ctr = contrasts.name != 'NO_FILE' ? "--contrasts ${contrasts}" : ''
    """
    Rscript ${workflow.projectDir}/bin/viral_window_db.R \\
        --counts ${counts} \\
        --samples ${samples} \\
        --out-prefix ${contig}.viral_db \\
        --fdr ${fdr} \\
        --merge-gap ${merge_gap} \\
        --min-count ${min_count} \\
        --norm ${norm_mode} \\
        ${ctr}
    """
}
