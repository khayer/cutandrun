/*
 * Convert bam files to bedgraph and bigwig with apropriate normalisation
 */

include { BEDTOOLS_GENOMECOV    } from "../../modules/nf-core/bedtools/genomecov/main"
include { DEEPTOOLS_BAMCOVERAGE } from "../../modules/local/for_patch/deeptools/bamcoverage/main"
include { BEDTOOLS_SORT         } from "../../modules/local/for_patch/bedtools/sort/main"
include { UCSC_BEDCLIP          } from "../../modules/nf-core/ucsc/bedclip/main"
include { UCSC_BEDGRAPHTOBIGWIG } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"

workflow PREPARE_PEAKCALLING {
    take:
    ch_bam           // channel: [ val(meta), [ bam ] ]
    ch_bai           // channel: [ val(meta), [ bai ] ]
    ch_chrom_sizes   // channel: [ sizes ]
    ch_dummy_file    // channel: [ dummy ]
    norm_mode        // value:   ["Spikein", "RPKM", "CPM", "BPM", "RPGC", "None" ]
    metadata         // channel  [ csv ] - spike-in metadata
    target_metadata  // channel  [ csv ] - target genome metadata
    mean_target_reads // value: mean target aligned reads (for dual normalization)

    main:
    ch_versions = Channel.empty()
    ch_bedgraph = Channel.empty()

    if (norm_mode == "Spikein") {
        /*
        * CHANNEL: Load up spike-in alignment metadata into channel
        */
        metadata.splitCsv ( header:true, sep:"," )
            .map { row -> [ row[0].id, row[1] ]}
            .set { ch_spikein_metadata }
        //ch_spikein_metadata | view

        /*
        * CHANNEL: Load up target genome alignment metadata into channel
        */
        target_metadata.splitCsv ( header:true, sep:"," )
            .map { row -> [ row[0].id, row[1] ]}
            .set { ch_target_metadata }
        //ch_target_metadata | view

        /*
        * CHANNEL: Calculate scale factor for each sample.
        * If dual normalization is enabled, multiply spike-in factor by target genome abundance factor.
        * Formula: scale_factor = (normalisation_c / spike_in_reads) * (mean_target_reads / sample_target_reads)
        */
        ch_bam.map { row -> [ row[0].id, row[0], row[1] ]}
            .join ( ch_spikein_metadata )
            .join ( ch_target_metadata )
            .combine ( mean_target_reads )
            .map { row ->
                def spikein_reads = row[3].find{ it.key == "bt2_total_aligned" }?.value.toInteger()
                def target_reads = row[4].find{ it.key == "bt2_total_aligned" }?.value.toInteger()
                def mean_target = row[5]
                
                // Calculate spike-in normalization factor
                def spikein_factor = params.normalisation_c / (spikein_reads != 0 ? spikein_reads : params.normalisation_c)
                
                // Calculate final scale factor
                def final_factor = spikein_factor
                if (params.normalisation_mode_dual && target_reads != null && target_reads > 0) {
                    // Apply dual normalization: spike-in × target abundance
                    def target_factor = mean_target / target_reads
                    final_factor = spikein_factor * target_factor
                    log.info "Sample ${row[0]}: spikein_factor=${spikein_factor}, target_factor=${target_factor}, final_factor=${final_factor}"
                }
                
                [ row[1], row[2], final_factor ]
            }
            .set { ch_bam_scale_factor }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, scale_factor]
        //ch_bam_scale_factor | view
    }
    else if (norm_mode == "None") {
        /*
        * CHANNEL: Assign scale factor of 1
        */
        ch_bam.map { row ->
                [ row[0], row[1], 1 ]
            }
            .set { ch_bam_scale_factor }
        //ch_bam_scale_factor | view
    }

    if (norm_mode == "Spikein" || norm_mode == "None") {
        /*
        * MODULE: Convert bam files to bedgraph
        */
        BEDTOOLS_GENOMECOV (
            ch_bam_scale_factor,
            ch_dummy_file,
            "bedGraph"
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
        ch_bedgraph = BEDTOOLS_GENOMECOV.out.genomecov
        //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
        //BEDTOOLS_GENOMECOV.out.genomecov | view

        /*
        * CHANNEL: Dump scale factor values
        */
        if(params.dump_scale_factors) {
            ch_scale_factor = ch_bam_scale_factor
            .map { [it[0].id, it[0].group, it[2]] }
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .map { list ->
                new File('scale-factors.csv').withWriter('UTF-8') { writer ->
                    list.each { item ->
                        str = item[0] + "," + item[1] + "," + item[2]
                        writer.write(str + "\n")
                    }
                }
            }
        }
    } else {
        /*
        * CHANNEL: Combine bam and bai files on id
        */
        ch_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }
        .set { ch_bam_bai }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bam_bai | view

        /*
        * CHANNEL: Split files based on igg or not
        */
        ch_bam_bai.branch { it ->
            target:  it[0].is_control == false
            control: it[0].is_control == true
        }
        .set { ch_bam_bai_split }

        /*
        * CHANNEL: Assign scale factor of 1 to target files
        */
        ch_bam_bai_split.target
            .map { row ->
                [ row[0], row[1], row[2], 1 ]
            }
        .set { ch_bam_bai_split_target }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_split_target | view

        /*
        * CHANNEL: Assign igg scale factor to target files
        */
        ch_bam_bai_split.control
            .map { row ->
                [ row[0], row[1], row[2], params.igg_scale_factor ]
            }
        .set { ch_bam_bai_split_igg }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_split_igg | view

        /*
        * CHANNEL: Mix the split channels back up
        */
        ch_bam_bai_split_target
            .mix(ch_bam_bai_split_igg)
        .set { ch_bam_bai_scale_factor }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
        //ch_bam_bai_scale_factor | view

        /*
        * MODULE: Convert bam files to bedgraph and normalise
        */
        DEEPTOOLS_BAMCOVERAGE (
            ch_bam_bai_scale_factor
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions)
        ch_bedgraph = DEEPTOOLS_BAMCOVERAGE.out.bedgraph
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bedgraph | view

        /*
        * CHANNEL: Dump scale factor values
        */
        if(params.dump_scale_factors) {
            ch_scale_factor = ch_bam_bai_scale_factor
            .map { [it[0].id, it[0].group, it[3]] }
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .map { list ->
                new File('scale-factors.csv').withWriter('UTF-8') { writer ->
                    list.each { item ->
                        str = item[0] + "," + item[1] + "," + item[2]
                        writer.write(str + "\n")
                    }
                }
            }
        }
    }

    /*
    * MODULE: Sort bedgraph
    */
    BEDTOOLS_SORT (
        ch_bedgraph,
        "bedGraph",
        []
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    /*
    * MODULE: Clip off bedgraphs so none overlap beyond chromosome edge
    */
    UCSC_BEDCLIP (
        BEDTOOLS_SORT.out.sorted,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions)
    //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
    //UCSC_BEDCLIP.out.bedgraph | view

    /*
    * MODULE: Convert bedgraph to bigwig
    */
    UCSC_BEDGRAPHTOBIGWIG (
        UCSC_BEDCLIP.out.bedgraph,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)
    //EXAMPLE CHANNEL STRUCT: [[META], BIGWIG]
    //UCSC_BEDGRAPHTOBIGWIG.out.bigwig | view

    emit:
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
