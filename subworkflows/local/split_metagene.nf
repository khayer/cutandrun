/*
 * Subworkflow: SPLIT_METAGENE_WF
 * Authors: Katharina E. Hayer, co-created with Claude
 * Version: 1.0.0  |  Date: 2026-11-23
 *
 * Builds the gene groups once from the run's own GTF, then draws one split
 * metagene per target group using that group's IgG-subtracted replicate
 * bigWigs. Optional pre-built BEDs bypass MAKE_GENE_GROUPS, which is how a
 * non-intron split (DESeq2 up vs unchanged promoters, say) is used instead.
 */

include { MAKE_GENE_GROUPS } from '../../modules/local/python/make_gene_groups'
include { SPLIT_METAGENE   } from '../../modules/local/python/split_metagene'

workflow SPLIT_METAGENE_WF {
    take:
        ch_bigwig          // [ meta, bigwig ] one row per replicate
        ch_gtf             // path to the run's GTF
        bed_a              // optional pre-built group A BED, or null
        bed_b              // optional pre-built group B BED, or null

    main:
        ch_versions = Channel.empty()

        if (bed_a && bed_b) {
            ch_ic = Channel.value(file(bed_a, checkIfExists: true))
            ch_il = Channel.value(file(bed_b, checkIfExists: true))
        }
        else {
            MAKE_GENE_GROUPS(
                ch_gtf,
                params.split_metagene_min_length,
                params.split_metagene_biotype,
                params.split_metagene_exclude_contigs
            )
            ch_ic       = MAKE_GENE_GROUPS.out.intron_containing
            ch_il       = MAKE_GENE_GROUPS.out.intronless
            ch_versions = ch_versions.mix(MAKE_GENE_GROUPS.out.versions)
        }

        // one row per target group, carrying every replicate bigWig it has
        ch_bigwig
            .filter { row -> row[0].is_control == false }
            .map    { row -> [ row[0].group, row[0], row[1] ] }
            .groupTuple(by: 0)
            .map    { group, metas, bws ->
                        [ metas[0] + [ group: group ], bws.flatten().sort { it.getName() } ] }
            .filter { row -> row[1].size() > 0 }
            .set    { ch_split_input }

        SPLIT_METAGENE(
            ch_split_input,
            ch_ic.first(),
            ch_il.first(),
            params.split_metagene_signal,
            params.split_metagene_ci,
            params.split_metagene_bin_size,
            params.split_metagene_tss_up,
            params.split_metagene_tss_down,
            params.split_metagene_tes_up,
            params.split_metagene_tes_down
        )
        ch_versions = ch_versions.mix(SPLIT_METAGENE.out.versions)

    emit:
        png      = SPLIT_METAGENE.out.png
        matrices = SPLIT_METAGENE.out.matrices
        versions = ch_versions
}
