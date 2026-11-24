/*
 * Subworkflow: DEEPTOOLS_HEATMAPS
 * Authors: Katharina E. Hayer, co-created with Claude
 * Version: 1.0.0  |  Date: 2026-11-24
 *
 * deepTools computeMatrix + plotHeatmap over gene bodies and over peak
 * summits. Lifted verbatim out of workflow CUTANDRUN: the body of that
 * workflow is compiled into a single JVM string constant, which is capped at
 * 65535 code units, and it had grown to within a few hundred bytes of the cap.
 * Moving a self-contained block out is the only way to get that budget back.
 *
 * Behaviour is unchanged, including the fact that the *_GENE_ALL processes do
 * not contribute to the versions channel; that is how it was in CUTANDRUN and
 * this move is not the place to change it.
 */

include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE     } from "../../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS    } from "../../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE_ALL } from "../../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE         } from "../../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS        } from "../../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE_ALL     } from "../../modules/nf-core/deeptools/plotheatmap/main"

workflow DEEPTOOLS_HEATMAPS {
    take:
        ch_bigwig         // [ meta, bigwig ] every sample, controls included
        ch_peaks_summits  // [ meta, bed ]
        ch_genome_bed     // PREPARE_GENOME.out.bed
        calc_all_matrix   // boolean, params.dt_calc_all_matrix

    main:
        ch_versions = Channel.empty()

        /*
        * CHANNEL: Remove IgG from bigwig channel
        */
        ch_bigwig.filter { it[0].is_control == false }
        .set { ch_bigwig_no_igg }

        /*
        * MODULE: Compute DeepTools matrix used in heatmap plotting for Genes
        */
        DEEPTOOLS_COMPUTEMATRIX_GENE (
            ch_bigwig_no_igg,
            ch_genome_bed.collect()
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX_GENE.out.versions)

        /*
        * MODULE: Calculate DeepTools heatmap
        */
        DEEPTOOLS_PLOTHEATMAP_GENE (
            DEEPTOOLS_COMPUTEMATRIX_GENE.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP_GENE.out.versions)

        /*
        * CHANNEL: Structure output for join on id
        */
        ch_peaks_summits
        .map { row -> [row[0].id, row ].flatten()}
        .set { ch_peaks_summits_id }

        /*
        * CHANNEL: Join beds and bigwigs on id
        */
        ch_bigwig_no_igg
        .map { row -> [row[0].id, row ].flatten()}
        .join ( ch_peaks_summits_id )
        .filter ( it -> it[-1].size() > 1)
        .set { ch_dt_bigwig_summits }

        ch_dt_bigwig_summits
        .map { row -> row[1,2] }
        .set { ch_ordered_bigwig }

        ch_dt_bigwig_summits
        .map { row -> row[-1] }
        .set { ch_ordered_peaks_max }

        /*
        * MODULE: Compute DeepTools matrix used in heatmap plotting for Peaks
        */
        DEEPTOOLS_COMPUTEMATRIX_PEAKS (
            ch_ordered_bigwig,
            ch_ordered_peaks_max
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.versions)

        /*
        * MODULE: Calculate DeepTools heatmap
        */
        DEEPTOOLS_PLOTHEATMAP_PEAKS (
            DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP_PEAKS.out.versions)

        if (calc_all_matrix) {
            /*
            * MODULE: Run calc gene matrix for all samples
            */
            ch_bigwig_no_igg
            .map { it[1] }
            .toSortedList()
            .filter { it && it.size() > 0 }
            .map { [[id:'all_genes'], it] }
            .set { ch_all_genes_bigwig_list }

            DEEPTOOLS_COMPUTEMATRIX_GENE_ALL (
                ch_all_genes_bigwig_list,
                ch_genome_bed.toSortedList()
            )

            /*
            * MODULE: Calculate DeepTools heatmap for all samples
            */
            DEEPTOOLS_PLOTHEATMAP_GENE_ALL (
                DEEPTOOLS_COMPUTEMATRIX_GENE_ALL.out.matrix
            )
        }

    emit:
        versions = ch_versions
}
