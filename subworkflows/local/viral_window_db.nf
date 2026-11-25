/*
 * Subworkflow: VIRAL_WINDOW_DB_WF
 * Authors: Katharina E. Hayer, co-created with Claude
 * Version: 1.0.0  |  Date: 2026-11-25
 *
 * Windowed differential coverage on one or more viral contigs, as an
 * alternative to peak calling on a small genome.
 *
 * Measured on the C14 Ad5 data: MACS2 called 23 peaks between two replicates
 * of the same IgG control, against 15-21 for the real comparisons, and 14 of
 * 18 "reproducible" target peaks were also called by IgG vs IgG. A peak caller
 * has no model of how much coverage varies between replicates by chance, and
 * on a 36 kb contig at ~24,000x that variation is the dominant signal. This
 * counts windows and lets edgeR estimate that variation from the replicates.
 * See claude/weitzman-jack/ad5_peak_calling_limit.md.
 */

include { VIRAL_WINDOW_COUNTS } from "../../modules/local/viral_window_counts"
include { VIRAL_WINDOW_DB     } from "../../modules/local/r/viral_window_db/main"

workflow VIRAL_WINDOW_DB_WF {
    take:
        ch_bam_bai   // [ meta, bam, bai ] every sample, targets and controls
        ch_fai       // PREPARE_GENOME.out.fai
        contigs      // comma-separated contig names, e.g. "Ad5,HSV17_NR"

    main:
        ch_versions = Channel.empty()

        // One ordered list of every sample, so the count columns and the
        // sample table cannot drift out of step.
        ch_bam_bai
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .map { rows ->
                def metas  = rows.collect { it[0] }
                def bams   = rows.collect { it[1] }
                def bais   = rows.collect { it[2] }
                def names  = metas.collect { it.id }
                def lines  = ["sample\tgroup\tis_control\tcontrol_group"]
                metas.each { m ->
                    lines << [ m.id,
                               m.group,
                               m.is_control ? 1 : 0,
                               m.is_control ? '' : (m.control_group ?: '') ].join('\t')
                }
                [ bams, bais, names, lines.join('\n') ]
            }
            .set { ch_all_samples }

        ch_contigs = Channel.from(contigs.split(',').collect { it.trim() }.findAll { it })

        VIRAL_WINDOW_COUNTS (
            ch_contigs,
            ch_fai.collect(),
            ch_all_samples.map { it[0] }.collect(),
            ch_all_samples.map { it[1] }.collect(),
            ch_all_samples.map { it[2] }.first(),
            ch_all_samples.map { it[3] }.first(),
            params.viral_db_window,
            params.viral_db_step,
            params.viral_db_min_mapq
        )
        ch_versions = ch_versions.mix(VIRAL_WINDOW_COUNTS.out.versions)

        ch_contrasts = params.viral_db_contrasts
            ? Channel.fromPath(params.viral_db_contrasts, checkIfExists: true).collect()
            : Channel.value(file("${projectDir}/assets/NO_FILE"))

        VIRAL_WINDOW_DB (
            VIRAL_WINDOW_COUNTS.out.counts,
            ch_contrasts,
            params.viral_db_fdr,
            params.viral_db_merge_gap,
            params.viral_db_min_count,
            params.viral_db_norm
        )

    emit:
        results  = VIRAL_WINDOW_DB.out.results
        summary  = VIRAL_WINDOW_DB.out.summary
        versions = ch_versions
}
