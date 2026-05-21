/*
========================================================================================
    nf-core/cutandrun - Custom modifications for Weitzman Lab
    Author: Katharina Hayer (khayer)
    Co-created with: GitHub Copilot (Claude Sonnet 4.5)
    
    Key additions:
    - Dual normalization (spike-in + target abundance)
    - Homer motif analysis integration
    - Enhanced peak annotation and reporting
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

// Validate input parameters in specialised library
WorkflowCutandrun.initialise(params, log)
def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Check input path parameters to see if the files exist if they have been specified
checkPathParamList = [
    params.blacklist,
    params.bowtie2,
    params.fasta,
    params.gtf,
    params.input
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if(params.normalisation_mode == "Spikein") {
    // Check spike-in only if it is enabled
    checkPathParamList = [
        params.spikein_bowtie2,
        params.spikein_fasta
    ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
}

// Check mandatory parameters that cannot be checked in the groovy lib as we want a channel for them
if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified!" }

ch_blacklist = Channel.empty()
if (params.blacklist) {
    ch_blacklist = Channel.from( file(params.blacklist) )
}
else {
    ch_blacklist = Channel.empty()
    WorkflowCutandrun.blacklistWarn(log)
}

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Stage awk files for parsing log files
ch_bt2_to_csv_awk     = file("$projectDir/bin/bt2_report_to_csv.awk"    , checkIfExists: true)
ch_dt_frag_to_csv_awk = file("$projectDir/bin/dt_frag_report_to_csv.awk", checkIfExists: true)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

// Load up and check multiqc base config and custom configs
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files for MultiQC
ch_frag_len_header_multiqc              = file("$projectDir/assets/multiqc/frag_len_header.txt", checkIfExists: true)
ch_frip_score_header_multiqc            = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_counts_header_multiqc           = file("$projectDir/assets/multiqc/peak_counts_header.txt", checkIfExists: true)
ch_peak_counts_consensus_header_multiqc = file("$projectDir/assets/multiqc/peak_counts_consensus_header.txt", checkIfExists: true)
ch_peak_reprod_header_multiqc           = file("$projectDir/assets/multiqc/peak_reprod_header.txt", checkIfExists: true)
ch_linear_duplication_header_multiqc    = file("$projectDir/assets/multiqc/linear_duplication_header.txt", checkIfExists: true)


/*
========================================================================================
    INIALISE PARAMETERS AND OPTIONS
========================================================================================
*/

// Init aligners
def prepare_tool_indices = ["bowtie2"]

// Check peak caller params
def caller_list = ['seacr', 'macs2']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : ['seacr']
if ((caller_list + callers).unique().size() != caller_list.size()) {
    exit 1, "Invalid variant calller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}"
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { INPUT_CHECK                } from "../subworkflows/local/input_check"
include { CUT as PEAK_TO_BED         } from '../modules/local/linux/cut'
include { AWK as AWK_NAME_PEAK_BED   } from "../modules/local/linux/awk"
include { IGV_SESSION                } from "../modules/local/python/igv_session"
include { IGV_SESSION as IGV_SESSION_DOWNSAMPLED } from "../modules/local/python/igv_session"
include { AWK as AWK_EXTRACT_SUMMITS } from "../modules/local/linux/awk"
include { SAMTOOLS_CUSTOMVIEW        } from "../modules/local/samtools_custom_view"
include { FRAG_LEN_HIST              } from "../modules/local/python/frag_len_hist"
include { MULTIQC                    } from "../modules/local/multiqc"
include { REPORT                     } from "../modules/local/report"
include { MERGE_PEAKS_TABLE          } from "../modules/local/python/merge_peaks_table"
include { MERGE_PEAKS_TABLE_CONTRAST } from "../modules/local/python/merge_peaks_table_contrast"
include { DOWNSAMPLE_BAM             } from "../modules/local/samtools_downsample"
include { HOMER_FINDMOTIFSGENOME as HOMER_FINDMOTIFSGENOME_MERGED     } from "../modules/local/homer/findmotifsgenome/main"
include { HOMER_FINDMOTIFSGENOME as HOMER_FINDMOTIFSGENOME_CONSENSUS  } from "../modules/local/homer/findmotifsgenome/main"
include { HOMER_ANNOTATEPEAKS                                          } from "../modules/local/homer/annotatepeaks/main"
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS         } from "../modules/local/homer/annotatepeaks/main"
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MERGED            } from "../modules/local/homer/annotatepeaks/main"
include { SUMMARIZE_HOMER_MOTIFS     } from "../modules/local/python/summarize_homer_motifs"
include { CREATE_MOTIF_COMPARISON_TABLES } from "../modules/local/python/create_motif_comparison_tables"
include { SUMMARIZE_PEAK_ANNOTATIONS  } from "../modules/local/python/summarize_peak_annotations"
include { PYGENOMETRACKS_TOP10        } from "../modules/local/python/pygenometracks_top10"
include { PROMOTER_GC_CONTRAST_QC     } from "../modules/local/python/promoter_gc_contrast_qc.nf"
include { CHIPSEEKER_ANNOTATE; CHIPSEEKER_ANNOTATE as CHIPSEEKER_ANNOTATE_CONSENSUS; CHIPSEEKER_COMPARE } from "../modules/local/r/chipseeker/main"
include { PROMOTER_GC_DIFFBIND        } from "../modules/local/r/promoter_gc_diffbind/main"

/*
 * SUBWORKFLOWS
 */
include { PREPARE_GENOME                                   } from "../subworkflows/local/prepare_genome"
include { FASTQC_TRIMGALORE                                } from "../subworkflows/local/fastqc_trimgalore"
include { ALIGN_BOWTIE2                                    } from "../subworkflows/local/align_bowtie2"
include { DESEQ2_PEAK_DIFF_ANALYSIS                        } from "../subworkflows/local/deseq2_peak_analysis"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_TARGET_META  } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_SPIKEIN_META } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_PICARD_DUP_META  } from "../subworkflows/local/extract_metadata_awk"
include { MARK_DUPLICATES_PICARD                           } from "../subworkflows/local/mark_duplicates_picard"
include { MARK_DUPLICATES_PICARD as DEDUPLICATE_PICARD     } from "../subworkflows/local/mark_duplicates_picard"
include { CONSENSUS_PEAKS                                  } from "../subworkflows/local/consensus_peaks"
include { CONSENSUS_PEAKS as CONSENSUS_PEAKS_ALL           } from "../subworkflows/local/consensus_peaks"
include { EXTRACT_FRAGMENTS                                } from "../subworkflows/local/extract_fragments"
include { PREPARE_PEAKCALLING                              } from "../subworkflows/local/prepare_peakcalling"
include { PREPARE_PEAKCALLING as PREPARE_PEAKCALLING_VIS   } from "../subworkflows/local/prepare_peakcalling"
include { DEEPTOOLS_QC                                     } from "../subworkflows/local/deeptools_qc"
include { PEAK_QC                                          } from "../subworkflows/local/peak_qc"
include { SAMTOOLS_VIEW_SORT_STATS as FILTER_READS         } from "../subworkflows/local/samtools_view_sort_stats"
include { DEDUPLICATE_LINEAR                               } from "../subworkflows/local/deduplicate_linear"
include { PEAK_SIGNAL_PROFILER                             } from "../subworkflows/local/peak_signal_profiler"

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { CAT_FASTQ                                                    } from "../modules/nf-core/cat/fastq/main"
include { PRESEQ_LCEXTRAP                                              } from "../modules/nf-core/preseq/lcextrap/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_IGG                         } from "../modules/nf-core/seacr/callpeak/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_NOIGG                       } from "../modules/nf-core/seacr/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_IGG                         } from "../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NOIGG                       } from "../modules/nf-core/macs2/callpeak/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE      } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS     } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE          } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS         } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE_ALL  } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS_ALL } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE_ALL      } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS_ALL     } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_BIGWIGCOMPARE                                      } from "../modules/nf-core/deeptools/bigwigcompare/main"
include { DEEPTOOLS_MULTIBAMSUMMARY_BED                                } from "../modules/nf-core/deeptools/multibamsummary_bed/main"
include { DEEPTOOLS_MULTIBAMSUMMARY_BED as DEEPTOOLS_MULTIBAMSUMMARY_BED_CONTRAST } from "../modules/nf-core/deeptools/multibamsummary_bed/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS                                  } from "../modules/local/custom_dumpsoftwareversions"

/*
 * SUBWORKFLOWS
 */


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CUTANDRUN {

    // Init
    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    if(params.run_genome_prep) {
        PREPARE_GENOME (
            prepare_tool_indices,
            ch_blacklist
        )
        ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.versions)
    }

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    if(params.run_input_check) {
        INPUT_CHECK (
            ch_input
        )

        INPUT_CHECK.out.reads
        .map {
            meta, fastq ->
                meta.id = meta.id.split("_")[0..-2].join("_")
                [ meta, fastq ] }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }
    }

    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    if(params.run_cat_fastq) {
        CAT_FASTQ (
            ch_fastq.multiple
        )
        ch_software_versions = ch_software_versions.mix(CAT_FASTQ.out.versions)

        CAT_FASTQ.out.reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [READS]]
    //ch_cat_fastq | view

    /*
     * SUBWORKFLOW: Read QC, trim adapters and perform post-trim read QC
     */
    if(params.run_trim_galore_fastqc) {
        FASTQC_TRIMGALORE (
            ch_cat_fastq,
            params.skip_fastqc,
            params.skip_trimming
        )
        ch_trimmed_reads     = FASTQC_TRIMGALORE.out.reads
        ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [READS]]
    //FASTQC_TRIMGALORE.out.reads | view

    /*
    * SUBWORKFLOW: Alignment to target and spikein genome using botwtie2
    */
    ch_orig_bam                   = Channel.empty()
    ch_orig_spikein_bam           = Channel.empty()
    ch_bowtie2_log                = Channel.empty()
    ch_bowtie2_spikein_log        = Channel.empty()
    ch_samtools_bam               = Channel.empty()
    ch_samtools_bai               = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_samtools_spikein_bam       = Channel.empty()
    ch_samtools_spikein_bai       = Channel.empty()
    ch_samtools_spikein_stats     = Channel.empty()
    ch_samtools_spikein_flagstat  = Channel.empty()
    ch_samtools_spikein_idxstats  = Channel.empty()
    if(params.run_alignment) {
        if (params.aligner == "bowtie2") {
            ALIGN_BOWTIE2 (
                ch_trimmed_reads,
                PREPARE_GENOME.out.bowtie2_index,
                PREPARE_GENOME.out.bowtie2_spikein_index,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.spikein_fasta
            )
            ch_software_versions          = ch_software_versions.mix(ALIGN_BOWTIE2.out.versions)
            ch_orig_bam                   = ALIGN_BOWTIE2.out.orig_bam
            ch_orig_spikein_bam           = ALIGN_BOWTIE2.out.orig_spikein_bam
            ch_bowtie2_log                = ALIGN_BOWTIE2.out.bowtie2_log
            ch_bowtie2_spikein_log        = ALIGN_BOWTIE2.out.bowtie2_spikein_log

            ch_samtools_bam               = ALIGN_BOWTIE2.out.bam
            ch_samtools_bai               = ALIGN_BOWTIE2.out.bai
            ch_samtools_stats             = ALIGN_BOWTIE2.out.stats
            ch_samtools_flagstat          = ALIGN_BOWTIE2.out.flagstat
            ch_samtools_idxstats          = ALIGN_BOWTIE2.out.idxstats

            ch_samtools_spikein_bam       = ALIGN_BOWTIE2.out.spikein_bam
            ch_samtools_spikein_bai       = ALIGN_BOWTIE2.out.spikein_bai
            ch_samtools_spikein_stats     = ALIGN_BOWTIE2.out.spikein_stats
            ch_samtools_spikein_flagstat  = ALIGN_BOWTIE2.out.spikein_flagstat
            ch_samtools_spikein_idxstats  = ALIGN_BOWTIE2.out.spikein_idxstats
        }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    // Preserve flagstats for FRiP before optional target deduplication
    ch_samtools_flagstat_for_frip = ch_samtools_flagstat

    /*
     * SUBWORKFLOW: extract aligner metadata
     */
    ch_metadata_bt2_target  = Channel.empty()
    ch_metadata_bt2_spikein = Channel.empty()
    if (params.aligner == "bowtie2" && params.run_alignment) {
        EXTRACT_BT2_TARGET_META (
            ch_bowtie2_log,
            ch_bt2_to_csv_awk,
            true
        )
        ch_metadata_bt2_target = EXTRACT_BT2_TARGET_META.out.metadata
        ch_software_versions   = ch_software_versions.mix(EXTRACT_BT2_TARGET_META.out.versions)

        EXTRACT_BT2_SPIKEIN_META (
            ch_bowtie2_spikein_log,
            ch_bt2_to_csv_awk,
            true
        )
        ch_metadata_bt2_spikein = EXTRACT_BT2_SPIKEIN_META.out.metadata
    }
    //ch_metadata_bt2_target | view
    //ch_metadata_bt2_spikein | view

    /*
     *  SUBWORKFLOW: Filter reads based some standard measures
     *  - Unmapped reads 0x004
     *  - Mate unmapped 0x0008
     *  - Multi-mapped reads
     *  - Filter out reads aligned to blacklist regions
     *  - Filter out reads below a threshold q score
     *  - Filter out mitochondrial reads (if required)
     */
    if (params.run_read_filter) {
        FILTER_READS (
            ch_samtools_bam,
            PREPARE_GENOME.out.allowed_regions.collect{it[1]}.ifEmpty([]),
            PREPARE_GENOME.out.fasta
        )
        ch_samtools_bam      = FILTER_READS.out.bam
        ch_samtools_bai      = FILTER_READS.out.bai
        ch_samtools_stats    = FILTER_READS.out.stats
        ch_samtools_flagstat = FILTER_READS.out.flagstat
        ch_samtools_idxstats = FILTER_READS.out.idxstats
        ch_software_versions = ch_software_versions.mix(FILTER_READS.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * MODULE: Run preseq on BAM files before de-duplication
    */
    ch_preseq_output = Channel.empty()
    if (params.run_preseq) {
        PRESEQ_LCEXTRAP (
            ch_samtools_bam
        )
        ch_preseq_output = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.versions)
    }

    /*
     * SUBWORKFLOW: Mark duplicates on all samples
     */
    ch_markduplicates_metrics = Channel.empty()
    if (params.run_mark_dups) {
        MARK_DUPLICATES_PICARD (
            ch_samtools_bam,
            ch_samtools_bai,
            true,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect()
        )
        ch_samtools_bam           = MARK_DUPLICATES_PICARD.out.bam
        ch_samtools_bai           = MARK_DUPLICATES_PICARD.out.bai
        ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_metrics = MARK_DUPLICATES_PICARD.out.metrics
        ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: Remove duplicates - default is on IgG controls only
     */
    if (params.run_remove_dups) {
        DEDUPLICATE_PICARD (
            ch_samtools_bam,
            ch_samtools_bai,
            params.dedup_target_reads,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect()
        )
        ch_samtools_bam      = DEDUPLICATE_PICARD.out.bam
        ch_samtools_bai      = DEDUPLICATE_PICARD.out.bai
        ch_samtools_stats    = DEDUPLICATE_PICARD.out.stats
        ch_samtools_flagstat = DEDUPLICATE_PICARD.out.flagstat
        ch_samtools_idxstats = DEDUPLICATE_PICARD.out.idxstats
        ch_software_versions = ch_software_versions.mix(DEDUPLICATE_PICARD.out.versions)

        if (params.dedup_target_reads) {
            ch_samtools_flagstat_for_frip = ch_samtools_flagstat
        }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bai | view

    /*
    * SUBWORKFLOW: extract duplication stats from picard report
    */
    ch_metadata_picard_duplicates = Channel.empty()
    if (params.run_mark_dups) {
        EXTRACT_PICARD_DUP_META (
            ch_markduplicates_metrics,
            ch_dummy_file.collect(),
            false
        )
        ch_metadata_picard_duplicates = EXTRACT_PICARD_DUP_META.out.metadata
        ch_software_versions          = ch_software_versions.mix(EXTRACT_PICARD_DUP_META.out.versions)
    }

    /*
     * SUBWORKFLOW: Remove linear amplification duplicates - default is false
     */
    ch_linear_metrics         = Channel.empty()
    ch_linear_duplication_mqc = Channel.empty()
    if (params.run_remove_linear_dups) {
        DEDUPLICATE_LINEAR (
            ch_samtools_bam,
            ch_samtools_bai,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect(),
            params.dedup_target_reads,
            ch_linear_duplication_header_multiqc
        )
        ch_samtools_bam           = DEDUPLICATE_LINEAR.out.bam
        ch_samtools_bai           = DEDUPLICATE_LINEAR.out.bai
        ch_samtools_stats         = DEDUPLICATE_LINEAR.out.stats
        ch_samtools_flagstat      = DEDUPLICATE_LINEAR.out.flagstat
        ch_samtools_idxstats      = DEDUPLICATE_LINEAR.out.idxstats
        ch_linear_metrics         = DEDUPLICATE_LINEAR.out.metrics
        ch_linear_duplication_mqc = DEDUPLICATE_LINEAR.out.linear_metrics_mqc
        ch_software_versions      = ch_software_versions.mix(DEDUPLICATE_LINEAR.out.versions)
    }

    ch_bedgraph               = Channel.empty()
    ch_bigwig                 = Channel.empty()
    ch_bigwig_original        = Channel.empty()
    ch_bigwig_visual          = Channel.empty()
    ch_bigwig_subtract        = Channel.empty()
    ch_seacr_peaks            = Channel.empty()
    ch_macs2_peaks            = Channel.empty()
    ch_peaks_primary          = Channel.empty()
    ch_peaks_secondary        = Channel.empty()
    ch_peaks_summits          = Channel.empty()
    ch_consensus_peaks        = Channel.empty()
    ch_consensus_peaks_unfilt = Channel.empty()
    def downsample_enabled = params.downsample_target_coverage && params.downsample_target_coverage > 0
    if(params.run_peak_calling) {
        /*
        * CHANNEL: Calculate mean target genome read count across all samples for dual normalization
        */
        ch_mean_target_reads = Channel.value(1.0)
        if (params.normalisation_mode_dual && params.normalisation_mode == "Spikein") {
            ch_metadata_bt2_target
                .splitCsv(header:true, sep:",")
                .map { row -> 
                    def aligned = row.find{ it.key == "bt2_total_aligned" }?.value
                    aligned ? aligned.toDouble() : 0.0
                }
                .toList()
                .map { list -> 
                    def valid_counts = list.findAll { it > 0 }
                    valid_counts.size() > 0 ? valid_counts.sum() / valid_counts.size() : 1.0
                }
                .set { ch_mean_target_reads }
            
            ch_mean_target_reads.view { "Mean target genome reads for dual normalization: ${it}" }
        }

        /*
        * SUBWORKFLOW: Convert BAM files to bedgraph/bigwig and apply configured normalisation strategy
        */
        def ch_samtools_bam_visual = ch_samtools_bam
        def ch_samtools_bai_visual = ch_samtools_bai

        if (downsample_enabled) {
            ch_samtools_bam
                .map { row -> [row[0].id, row ].flatten() }
                .join ( ch_samtools_bai.map { row -> [row[0].id, row ].flatten()} )
                .map { row -> [row[1], row[2], row[4]] }
                .set { ch_bam_bai_downsample_all }

            ch_bam_bai_downsample_all.branch { it ->
                targets:  it[0].is_control == false
                controls: it[0].is_control == true
            }
            .set { ch_bam_bai_downsample_split }

            def downsample_scope = params.downsample_apply ?: 'all'
            def ch_bam_bai_downsample
            def ch_bam_bai_passthrough = Channel.empty()

            if (downsample_scope == 'targets') {
                ch_bam_bai_downsample = ch_bam_bai_downsample_split.targets
                ch_bam_bai_passthrough = ch_bam_bai_downsample_split.controls
            } else if (downsample_scope == 'controls') {
                ch_bam_bai_downsample = ch_bam_bai_downsample_split.controls
                ch_bam_bai_passthrough = ch_bam_bai_downsample_split.targets
            } else {
                ch_bam_bai_downsample = ch_bam_bai_downsample_all
            }

            DOWNSAMPLE_BAM (
                ch_bam_bai_downsample,
                PREPARE_GENOME.out.chrom_sizes.collect(),
                params.downsample_target_coverage,
                params.downsample_seed
            )
            ch_software_versions = ch_software_versions.mix(DOWNSAMPLE_BAM.out.versions)

            def ch_bam_bai_downsampled = DOWNSAMPLE_BAM.out.bam
            if (downsample_scope == 'targets' || downsample_scope == 'controls') {
                ch_bam_bai_downsampled = ch_bam_bai_downsampled.mix(ch_bam_bai_passthrough)
            }

            // Keep original meta.id so spike-in metadata joins work in PREPARE_PEAKCALLING_VIS
            // (adding .viz suffix would break the join on sample ID)
            ch_samtools_bam_visual = ch_bam_bai_downsampled.map { [it[0], it[1]] }
            ch_samtools_bai_visual = ch_bam_bai_downsampled.map { [it[0], it[2]] }
        }

        PREPARE_PEAKCALLING(
            ch_samtools_bam,
            ch_samtools_bai,
            PREPARE_GENOME.out.chrom_sizes.collect(),
            ch_dummy_file,
            params.normalisation_mode,
            ch_metadata_bt2_spikein,
            ch_metadata_bt2_target,
            ch_mean_target_reads,
            false // visualization_mode: false for normal processing
        )
        ch_bedgraph          = PREPARE_PEAKCALLING.out.bedgraph
        PREPARE_PEAKCALLING.out.bigwig
            .toList()
            .set { ch_bigwig_rows }

        ch_bigwig          = ch_bigwig_rows.flatMap { rows -> rows }
        ch_bigwig_original = ch_bigwig_rows.flatMap { rows -> rows }
        ch_software_versions = ch_software_versions.mix(PREPARE_PEAKCALLING.out.versions)

        if (downsample_enabled) {
            // Visualization branch: downsampled BAMs, no normalization or reporting
            PREPARE_PEAKCALLING_VIS(
                ch_samtools_bam_visual,
                ch_samtools_bai_visual,
                PREPARE_GENOME.out.chrom_sizes.collect(),
                ch_dummy_file,
                "None", // Force no normalization for visualization
                Channel.empty(), // No metadata
                Channel.empty(), // No target metadata
                Channel.value(1.0), // Dummy mean_target_reads
                true // visualization_mode: true for downsampled visualization
            )
            ch_bigwig_visual      = PREPARE_PEAKCALLING_VIS.out.bigwig
            ch_software_versions  = ch_software_versions.mix(PREPARE_PEAKCALLING_VIS.out.versions)
            // Only use visualization outputs for IGV session and bigwig
            ch_bigwig_igv         = ch_bigwig_visual
        }

        /*
        * MODULE: Compute log2 ratio of ChIP vs control bigwigs
        */
        log.info "Checking bigwig subtract: use_control=${params.use_control}, run_bigwig_subtract=${params.run_bigwig_subtract}"
        
        if(params.use_control && params.run_bigwig_subtract) {
            // Separate bigwigs into target and control
            ch_bigwig.filter { it -> it[0].is_control == false }
            .set { ch_bigwig_target_compare }
            
            ch_bigwig.filter { it -> it[0].is_control == true }
            .set { ch_bigwig_control_compare }
            
            // Match targets with their controls
            // Target: key = control_group (already includes replicate, e.g., "igg_ctrl_1")
            ch_bigwig_target_compare
            .map { row -> [ row[0].control_group, row[0], row[1] ] }
            .set { ch_bigwig_target_keyed }
            
            // Control: key = group + "_" + replicate (e.g., "igg_ctrl_1")
            ch_bigwig_control_compare
            .map { row -> [ row[0].group + "_" + row[0].replicate, row[0], row[1] ] }
            .set { ch_bigwig_control_keyed }
            
            // Combine all targets with all controls, then filter for matching keys
            ch_bigwig_target_keyed
            .combine( ch_bigwig_control_keyed )
            .filter { it[0] == it[3] }  // target_key == control_key
            .map { row -> [ row[1], row[2], row[5] ] }  // [meta_target, bigwig_target, bigwig_control]
            .set { ch_bigwig_pairs }
            // EXAMPLE CHANNEL STRUCT: [[META_TARGET], BIGWIG_TARGET, BIGWIG_CONTROL]
            
            if (!params.skip_bigwigcompare) {
                DEEPTOOLS_BIGWIGCOMPARE (
                    ch_bigwig_pairs
                )
                ch_bigwig_subtract   = DEEPTOOLS_BIGWIGCOMPARE.out.bigwig
                ch_software_versions = ch_software_versions.mix(DEEPTOOLS_BIGWIGCOMPARE.out.versions)
            }
        }

        /*
         * CHANNEL: Separate bedgraphs into target/control
         */
        ch_bedgraph.filter { it -> it[0].is_control == false }
        .set { ch_bedgraph_target }
        ch_bedgraph.filter { it -> it[0].is_control == true }
        .set { ch_bedgraph_control }
        //ch_bedgraph_target | view
        //ch_bedgraph_control | view

        // Preserve BAMs for promoter-GC diffbind BEFORE consuming them for peak calling
        // Fan out the BAM/BAI channels to avoid consumption
        ch_samtools_bam
            .multiMap { row ->
                for_peakcalling: row
                for_promoter_gc: row
            }
            .set { ch_bam_split }
        
        ch_samtools_bai
            .multiMap { row ->
                for_peakcalling: row
                for_promoter_gc: row
            }
            .set { ch_bai_split }

        /*
        * CHANNEL: Separate bams into target/control
        */
        ch_bam_split.for_peakcalling.filter { it -> it[0].is_control == false }
        .set { ch_bam_target }
        ch_bam_split.for_peakcalling.filter { it -> it[0].is_control == true }
        .set { ch_bam_control }

        ch_bam_target
            .multiMap { row ->
                for_peakcalling: row
                for_promoter_gc: row
            }
            .set { ch_bam_target_split }

        ch_bam_target_split.for_promoter_gc
            .map { meta, bam -> [meta.id, [meta, bam]] }
            .join(ch_bai_split.for_peakcalling.filter { it -> it[0].is_control == false }.map { meta, bai -> [meta.id, [meta, bai]] }, by: 0)
            .map { sid, left_row, right_row -> [left_row[0], left_row[1], right_row[1]] }
            .set { ch_promoter_gc_bam_bai }

        //ch_bam_target | view
        //ch_bam_control | view

        if(params.use_control) {
            /*
            * MODULE: Call peaks using SEACR with IgG control
            */
            if('seacr' in callers) {
                /*
                * CHANNEL: Create target/control pairings
                */
                ch_bedgraph_control.map{ row -> [row[0].control_group + "_" + row[0].replicate, row] }
                .cross( ch_bedgraph_target.map{ row -> [row[0].control_group, row] } )
                .map {
                    row ->
                    [ row[1][1][0], row[1][1][1], row[0][1][1] ]
                }
                .set { ch_bedgraph_paired }
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]

                SEACR_CALLPEAK_IGG (
                    ch_bedgraph_paired,
                    params.seacr_peak_threshold
                )
                ch_seacr_peaks       = SEACR_CALLPEAK_IGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_IGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_CALLPEAK_IGG.out.bed | view
            }

            if('macs2' in callers) {
                /*
                * CHANNEL: Create target/control pairings
                */
                ch_bam_control.map{ row -> [row[0].control_group + "_" + row[0].replicate, row] }
                .cross( ch_bam_target_split.for_peakcalling.map{ row -> [row[0].control_group, row] } )
                .map {
                    row ->
                    [ row[1][1][0], row[1][1][1], row[0][1][1] ]
                }
                .set { ch_bam_paired }
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BAM, CONTROL_BAM]
                //ch_bam_paired | view

                MACS2_CALLPEAK_IGG (
                    ch_bam_paired,
                    params.macs_gsize
                )
                ch_macs2_peaks       = MACS2_CALLPEAK_IGG.out.peak
                ch_peaks_summits     = MACS2_CALLPEAK_IGG.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_IGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //MACS2_CALLPEAK_IGG.out.peak | view
            }
        }
        else {
            /*
            * MODULE: Call peaks without IgG Control
            */
            if('seacr' in callers) {
                /*
                * CHANNEL: Add fake control channel
                */
                ch_bedgraph_target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_bedgraph_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
                // ch_bedgraph_target_fctrl | view

                SEACR_CALLPEAK_NOIGG (
                    ch_bedgraph_target_fctrl,
                    params.seacr_peak_threshold
                )
                ch_seacr_peaks       = SEACR_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_NO_IGG.out.bed | view
            }

            if('macs2' in callers) {
                /*
                * CHANNEL: Add fake control channel
                */
                ch_bam_target_split.for_peakcalling.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_samtools_bam_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]
                //ch_samtools_bam_target_fctrl | view

                MACS2_CALLPEAK_NOIGG (
                    ch_samtools_bam_target_fctrl,
                    params.macs_gsize
                )
                ch_macs2_peaks       = MACS2_CALLPEAK_NOIGG.out.peak
                ch_peaks_summits     = MACS2_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                // MACS2_CALLPEAK_NOIGG.out.peak | view
            }
        }

        if ("macs2" in callers) {
            /*
            * MODULE: Convert narrow or broad peak to bed
            */
            PEAK_TO_BED ( ch_macs2_peaks )
            ch_macs2_peaks       = PEAK_TO_BED.out.file
            ch_software_versions = ch_software_versions.mix(PEAK_TO_BED.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //PEAK_TO_BED.out.file | view
        }

        // Identify the primary peak data stream for downstream analysis
        if(callers[0] == 'seacr') {
            ch_peaks_primary   = ch_seacr_peaks
            ch_peaks_secondary = ch_macs2_peaks
        }
        if(callers[0] == 'macs2') {
            ch_peaks_primary   = ch_macs2_peaks
            ch_peaks_secondary = ch_seacr_peaks
        }

        // Broadcast primary peaks so each downstream consumer gets a full copy.
        ch_peaks_primary
            .multiMap { row ->
                summits: row
                awk_name: row
                merge_table: row
                homer_annot: row
                promoter_gc_match: row
                promoter_gc_debug: row
                peak_qc: row
                igv_a: row
                igv_b: row
            }
            .set { ch_peaks_primary_split }

        if(callers[0] == 'seacr') {
            /*
            * MODULE: Extract summits from seacr peak beds
            */
            AWK_EXTRACT_SUMMITS (
                ch_peaks_primary_split.summits
            )
            ch_peaks_summits     = AWK_EXTRACT_SUMMITS.out.file
            ch_software_versions = ch_software_versions.mix(AWK_EXTRACT_SUMMITS.out.versions)
            //AWK_EXTRACT_SUMMITS.out.file | view
        }

        /*
        * MODULE: Add sample identifier column to peak beds
        */
        AWK_NAME_PEAK_BED (
            ch_peaks_primary_split.awk_name
        )
        ch_software_versions = ch_software_versions.mix(AWK_NAME_PEAK_BED.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //AWK_NAME_PEAK_BED.out.file | view

        /*
        * MODULE: Create merged peaks table showing peak presence across samples
        */
        MERGE_PEAKS_TABLE (
            ch_peaks_primary_split.merge_table.collect{it[1]}.ifEmpty([])
        )
        ch_software_versions = ch_software_versions.mix(MERGE_PEAKS_TABLE.out.versions)
        
        // Save merged peaks to a variable that can be used by multiple downstream processes
        ch_merged_peaks_bed = MERGE_PEAKS_TABLE.out.bed

        /*
        * Shared reference inputs for downstream peak annotation and contrast analyses
        */
        ch_fasta_for_peak_annotation = PREPARE_GENOME.out.fasta.map { it[1] }.first()
        ch_gtf_for_peak_annotation   = PREPARE_GENOME.out.gtf.first()

        /*
        * MODULE: Annotate replicate peaks with HOMER annotatePeaks and plot feature composition
        */
        if(params.run_homer_peak_annotation) {

            HOMER_ANNOTATEPEAKS (
                ch_peaks_primary_split.homer_annot,
                ch_fasta_for_peak_annotation,
                ch_gtf_for_peak_annotation
            )
            ch_software_versions = ch_software_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)

            SUMMARIZE_PEAK_ANNOTATIONS (
                HOMER_ANNOTATEPEAKS.out.annot.map { meta, table -> table }.collect()
            )
            ch_software_versions = ch_software_versions.mix(SUMMARIZE_PEAK_ANNOTATIONS.out.versions)

            /*
            * MODULE: Create enhanced ChIPseeker visualizations
            */
            if(params.run_chipseeker) {
                CHIPSEEKER_ANNOTATE (
                    HOMER_ANNOTATEPEAKS.out.annot,
                    file("${projectDir}/assets/dummy_file.txt"),  // Placeholder for TxDb
                    ch_gtf_for_peak_annotation
                )
                ch_software_versions = ch_software_versions.mix(CHIPSEEKER_ANNOTATE.out.versions)

                // Optional: Aggregate plots across samples
                if(params.run_chipseeker_compare) {
                    CHIPSEEKER_COMPARE (
                        CHIPSEEKER_ANNOTATE.out.annotation.map { meta, anno -> anno }.collect()
                    )
                    ch_software_versions = ch_software_versions.mix(CHIPSEEKER_COMPARE.out.versions)
                }
            }
        }

        /*
        * MODULE: Run Homer motif finding on merged peaks
        */
        if(params.run_homer_motifs) {
            // Create a value channel for the fasta file
            ch_fasta_for_homer_merged = PREPARE_GENOME.out.fasta.map { it[1] }.first()
            
            HOMER_FINDMOTIFSGENOME_MERGED (
                ch_merged_peaks_bed.map { bed -> [ [id: 'merged_peaks'], bed ] },
                ch_fasta_for_homer_merged,
                params.homer_motif_size
            )
            ch_software_versions = ch_software_versions.mix(HOMER_FINDMOTIFSGENOME_MERGED.out.versions)
            
            // Collect all motif directories for summarization
            HOMER_FINDMOTIFSGENOME_MERGED.out.motifs
                .map { meta, dir -> dir }
                .set { ch_motif_dirs_merged }
        }

        /*
        * MODULE: Annotate merged peaks with read counts from all BAMs
        */
        // Prepare BAM channel for multiBamSummary
        ch_samtools_bam
        .join(ch_samtools_bai, by: 0)
        .map { meta, bam, bai -> [ meta, bam, bai, meta.id ] }
        .toSortedList { it[3] }  // Sort by sample ID
        .map { list ->
            def bams = []
            def bais = []
            def labels = []
            list.each { item ->
                bams.add(item[1])
                bais.add(item[2])
                labels.add(item[3])
            }
            [ [id: 'all_samples'], bams, bais, labels ]
        }
        .set { ch_bams_for_summary }
        
        DEEPTOOLS_MULTIBAMSUMMARY_BED (
            ch_bams_for_summary,
            ch_merged_peaks_bed.collect()
        )
        ch_software_versions = ch_software_versions.mix(DEEPTOOLS_MULTIBAMSUMMARY_BED.out.versions)

        if (params.run_promoter_gc_diffbind) {
            def promoter_gc_contrasts = []
            if (params.promoter_gc_contrasts) {
                def contrasts_file = file(params.promoter_gc_contrasts)
                if (!contrasts_file.exists()) {
                    error "promoter_gc_contrasts file not found: ${params.promoter_gc_contrasts}"
                }

                contrasts_file.readLines().eachWithIndex { line, idx ->
                    def trimmed = line.trim()
                    if (!trimmed) {
                        return
                    }

                    // Allow a simple header row such as: groupA groupB
                    if (idx == 0 && trimmed.toLowerCase().replaceAll("\\s+", "").startsWith("groupagroupb")) {
                        return
                    }

                    def fields = trimmed.split(/\s+/)
                    if (fields.size() < 2) {
                        error "Invalid promoter_gc_contrasts row '${line}'. Expected two columns: groupA groupB"
                    }

                    def group_a = fields[0]
                    def group_b = fields[1]
                    promoter_gc_contrasts << [id: "${group_a}_vs_${group_b}", group_a: group_a, group_b: group_b]
                }

                if (promoter_gc_contrasts.isEmpty()) {
                    error "No contrasts found in promoter_gc_contrasts file: ${params.promoter_gc_contrasts}"
                }
            } else {
                promoter_gc_contrasts << [
                    id     : "${params.promoter_gc_group_a}_vs_${params.promoter_gc_group_b}",
                    group_a: params.promoter_gc_group_a,
                    group_b: params.promoter_gc_group_b
                ]
            }

            def canonical_condition = { value ->
                def s = value?.toString()?.trim()
                if (!s) {
                    return null
                }
                // Normalize common replicate and technical-replicate suffixes such as
                // _R1, _R1_T1, _R1_sorted, or _T1.
                def normalized = s
                while (normalized ==~ /^.+?_[RT]\d+(?:[_\.-].*)?$/) {
                    normalized = normalized.replaceFirst(/_[RT]\d+(?:[_\.-].*)?$/, '')
                }
                return normalized
            }

            def unwrap_meta = { value ->
                def current = value
                while (current instanceof List && !current.isEmpty()) {
                    def first = current[0]
                    if (first instanceof Map) {
                        return first
                    }
                    current = first
                }
                return (current instanceof Map) ? current : null
            }

            def peak_group_key = { meta ->
                meta = unwrap_meta(meta)
                def raw = meta?.group ?: meta?.condition ?: meta?.sample ?: meta?.id ?: meta?.control_group
                return canonical_condition(raw)
            }

            Channel
                .fromList(promoter_gc_contrasts)
                .combine(ch_peaks_primary_split.promoter_gc_match)
                .map { row ->
                    def contrast = row[0]
                    def peak_meta = null
                    def peak_bed = null
                    if (row instanceof List && row.size() > 2 && row[1] instanceof Map) {
                        // combine() may flatten tuples as [contrast, meta, bed]
                        peak_meta = row[1]
                        peak_bed = row[2]
                    } else {
                        // Defensive fallback for nested tuple shapes: [contrast, [meta, bed]]
                        def peak_row = (row instanceof List && row.size() > 1) ? row[1] : null
                        peak_meta = unwrap_meta(peak_row)
                        peak_bed = (peak_row instanceof List && peak_row.size() > 1) ? peak_row[1] : null
                    }
                    def peak_key = peak_group_key(peak_meta)
                    def group_a_key = canonical_condition(contrast?.group_a)
                    def group_b_key = canonical_condition(contrast?.group_b)
                    [contrast, peak_key, group_a_key, group_b_key, peak_bed]
                }
                .filter { row -> row[1] && row[4] && (row[1] == row[2] || row[1] == row[3]) }
                .map { row ->
                    def contrast = row[0]
                    def peak_key = row[1]
                    def peak_bed = row[4]
                    def matched_group = (peak_key == row[2]) ? contrast.group_a : contrast.group_b
                    [contrast.id, contrast, matched_group, peak_bed]
                }
                .groupTuple(by: 0)
                .map { contrast_id, contrasts, groups, beds ->
                    def contrast = contrasts[0]
                    if (!groups.contains(contrast.group_a) || !groups.contains(contrast.group_b)) {
                        error "Contrast ${contrast.id} is missing peaks from one or both groups in ch_peaks_primary"
                    }
                    [contrast, beds, beds.size()]
                }
                .ifEmpty {
                    error "No peaks matched promoter GC contrast(s): ${promoter_gc_contrasts}. Check group names and metadata propagation."
                }
                .set { ch_promoter_gc_peak_beds_by_contrast }

            MERGE_PEAKS_TABLE_CONTRAST (
                ch_promoter_gc_peak_beds_by_contrast.map { meta, beds, n_input_beds -> [meta + [n_input_peak_beds: n_input_beds], beds] }
            )
            ch_software_versions = ch_software_versions.mix(MERGE_PEAKS_TABLE_CONTRAST.out.versions)

            PROMOTER_GC_CONTRAST_QC (
                MERGE_PEAKS_TABLE_CONTRAST.out.bed.map { meta, merged_bed -> [meta, merged_bed, meta.n_input_peak_beds] }
            )
            ch_software_versions = ch_software_versions.mix(PROMOTER_GC_CONTRAST_QC.out.versions)

            MERGE_PEAKS_TABLE_CONTRAST.out.bed
                .multiMap { row ->
                    counts: row
                    homer: row
                }
                .set { ch_promoter_gc_merged_bed_split }

            log.info "DEBUG: Promoter GC contrasts: ${promoter_gc_contrasts.size()} items"
            ch_peaks_primary_split.promoter_gc_debug.view { "DEBUG peak row: meta.group=${it[0]?.group}, meta.id=${it[0]?.id}" }

            def bam_group_key = { meta ->
                meta = unwrap_meta(meta)
                def raw = meta?.group ?: meta?.condition ?: meta?.sample ?: meta?.id ?: meta?.control_group
                return canonical_condition(raw)
            }

            Channel
                .fromList(promoter_gc_contrasts)
                .combine(ch_promoter_gc_bam_bai)
                .view { "DEBUG: Combined contrast+BAM pair: ${it}" }
                .map { pair ->
                    def contrast = pair[0]
                    def bam_meta = null
                    def bam_file = null
                    def bai_file = null
                    if (pair instanceof List && pair.size() > 3 && pair[1] instanceof Map) {
                        // combine() may flatten tuples as [contrast, meta, bam, bai]
                        bam_meta = pair[1]
                        bam_file = pair[2]
                        bai_file = pair[3]
                    } else {
                        // Defensive fallback for nested tuple shapes: [contrast, [meta, bam, bai]]
                        def bam_row = (pair instanceof List && pair.size() > 1) ? pair[1] : null
                        bam_meta = unwrap_meta(bam_row)
                        bam_file = (bam_row instanceof List && bam_row.size() > 1) ? bam_row[1] : null
                        bai_file = (bam_row instanceof List && bam_row.size() > 2) ? bam_row[2] : null
                    }
                    def bam_key = bam_group_key(bam_meta)
                    def group_a_key = canonical_condition(contrast?.group_a)
                    def group_b_key = canonical_condition(contrast?.group_b)
                    def grp = null
                    if (bam_key && bam_key == group_a_key) {
                        grp = contrast.group_a
                    } else if (bam_key && bam_key == group_b_key) {
                        grp = contrast.group_b
                    }
                    [contrast.id, contrast, grp, bam_file, bai_file, bam_meta?.id]
                }
                .filter { row -> row[2] != null }
                .groupTuple(by: 0)
                .map { contrast_id, contrasts, groups, bams, bais, labels ->
                    def contrast = contrasts[0]
                    if (!groups.contains(contrast.group_a) || !groups.contains(contrast.group_b)) {
                        error "Contrast ${contrast.id} is missing BAMs from one or both groups"
                    }
                    def group_counts = groups.countBy { it }
                    def group_a_count = group_counts[contrast.group_a] ?: 0
                    def group_b_count = group_counts[contrast.group_b] ?: 0
                    if (group_a_count < 2 || group_b_count < 2) {
                        error "Contrast ${contrast.id} requires at least two BAMs per group for DESeq2, but found ${group_a_count} for ${contrast.group_a} and ${group_b_count} for ${contrast.group_b}. Update --promoter_gc_contrasts or the input samplesheet."
                    }
                    def idx = (0..<labels.size()).toList().sort { labels[it] }
                    def sorted_bams = idx.collect { bams[it] }
                    def sorted_bais = idx.collect { bais[it] }
                    def sorted_labels = idx.collect { labels[it] }
                    [contrast, sorted_bams, sorted_bais, sorted_labels]
                }
                .ifEmpty {
                    error "No BAMs matched promoter GC contrast(s): ${promoter_gc_contrasts}. Check group names and metadata propagation."
                }
                .set { ch_promoter_gc_bams_by_contrast }

            ch_promoter_gc_bams_by_contrast
                .map { meta, bams, bais, labels -> [meta.id, meta, bams, bais, labels] }
                .join(ch_promoter_gc_merged_bed_split.counts.map { meta, bed -> [meta.id, meta, bed] }, by: 0)
                .map { contrast_id, bam_meta, bams, bais, labels, bed_meta, bed -> [bam_meta, bams, bais, labels, bed] }
                .set { ch_promoter_gc_count_inputs }

            DEEPTOOLS_MULTIBAMSUMMARY_BED_CONTRAST (
                ch_promoter_gc_count_inputs.map { meta, bams, bais, labels, bed -> [meta, bams, bais, labels] },
                ch_promoter_gc_count_inputs.map { meta, bams, bais, labels, bed -> bed }
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_MULTIBAMSUMMARY_BED_CONTRAST.out.versions)

            HOMER_ANNOTATEPEAKS_MERGED (
                ch_promoter_gc_merged_bed_split.homer,
                ch_fasta_for_peak_annotation,
                ch_gtf_for_peak_annotation
            )
            ch_software_versions = ch_software_versions.mix(HOMER_ANNOTATEPEAKS_MERGED.out.versions)

            HOMER_ANNOTATEPEAKS_MERGED.out.annot
                .map { meta, table -> [meta.id, meta, table] }
                .join(DEEPTOOLS_MULTIBAMSUMMARY_BED_CONTRAST.out.table.map { meta, table -> [meta.id, meta, table] }, by: 0)
                .map { contrast_id, annot_meta, annot_table, counts_meta, counts_table -> [annot_meta, annot_table, counts_table] }
                .set { ch_promoter_gc_diffbind_inputs }

            PROMOTER_GC_DIFFBIND (
                ch_promoter_gc_diffbind_inputs
            )
            ch_software_versions = ch_software_versions.mix(PROMOTER_GC_DIFFBIND.out.versions)
        }

        if (params.run_deseq2_peak_analysis) {
            DESEQ2_PEAK_DIFF_ANALYSIS (
                DEEPTOOLS_MULTIBAMSUMMARY_BED.out.table,
                ch_gtf_for_peak_annotation
            )
            ch_software_versions = ch_software_versions.mix(DESEQ2_PEAK_DIFF_ANALYSIS.out.versions)
        }

        if(params.run_consensus_all) {
            /*
            * CHANNEL: Group all samples, filter where the number in the group is > 1
            */
            AWK_NAME_PEAK_BED.out.file
            .map { row -> [ 1, row[1] ] }
            .groupTuple(by: [0])
            .map { row ->
                def new_meta = [:]
                new_meta.put( "id", "all_samples" )
                [ new_meta, row[1].flatten() ]
            }
            .map { row ->
                [ row[0], row[1], row[1].size() ]
            }
            .filter { row -> row[2] > 1 }
            .map { row ->
                [ row[0], row[1] ]
            }
            .set { ch_peaks_bed_all }
            // EXAMPLE CHANNEL STRUCT: [[id: all_samples], [BED1, BED2, BEDn...], count]
            //ch_peaks_bed_all | view

            /*
            * SUBWORKFLOW: Construct group consensus peaks
            */
            CONSENSUS_PEAKS_ALL (
                ch_peaks_bed_all
            )
            ch_consensus_peaks        = CONSENSUS_PEAKS_ALL.out.filtered_bed
            ch_consensus_peaks_unfilt = CONSENSUS_PEAKS_ALL.out.merged_bed
            ch_software_versions      = ch_software_versions.mix(CONSENSUS_PEAKS_ALL.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //CONSENSUS_PEAKS_ALL.out.bed | view
        } else {
            /*
            * CHANNEL: Group samples based on group name
            */
            AWK_NAME_PEAK_BED.out.file
            .map { row -> [ row[0].group, row[1] ] }
            .groupTuple(by: [0])
            .map { row -> [ [id: row[0]], row[1].flatten() ] }
            .set { ch_peaks_bed_group }
            // EXAMPLE CHANNEL STRUCT: [[id: <GROUP>], [BED1, BED2, BEDn...], count]
            //ch_peaks_bed_group | view

            /*
            * SUBWORKFLOW: Construct group consensus peaks
            * where there is more than 1 replicate in a group
            */
            CONSENSUS_PEAKS (
                ch_peaks_bed_group
            )
            ch_consensus_peaks        = CONSENSUS_PEAKS.out.filtered_bed
            ch_consensus_peaks_unfilt = CONSENSUS_PEAKS.out.merged_bed
            ch_software_versions      = ch_software_versions.mix(CONSENSUS_PEAKS.out.versions)
            // EXAMPLE CHANNEL STRUCT: [[META], BED]
            //CONSENSUS_PEAKS.out.bed | view

            /*
            * MODULE: Run Homer motif finding on consensus peaks per group
            */
            if(params.run_homer_motifs) {
                // Create a value channel for the fasta file
                ch_fasta_for_homer = PREPARE_GENOME.out.fasta.map { it[1] }.first()
                
                HOMER_FINDMOTIFSGENOME_CONSENSUS (
                    ch_consensus_peaks_unfilt,
                    ch_fasta_for_homer,
                    params.homer_motif_size
                )
                ch_software_versions = ch_software_versions.mix(HOMER_FINDMOTIFSGENOME_CONSENSUS.out.versions)
                
                // Collect all consensus motif directories
                HOMER_FINDMOTIFSGENOME_CONSENSUS.out.motifs
                    .map { meta, dir -> dir }
                    .collect()
                    .set { ch_motif_dirs_consensus }
                
                // Combine with merged peaks motifs and generate summary
                ch_motif_dirs_merged
                    .mix(ch_motif_dirs_consensus.flatten())
                    .collect()
                    .set { ch_all_motif_dirs }
                
                SUMMARIZE_HOMER_MOTIFS (
                    ch_all_motif_dirs
                )
                ch_software_versions = ch_software_versions.mix(SUMMARIZE_HOMER_MOTIFS.out.versions)
                
                // Generate motif comparison tables
                CREATE_MOTIF_COMPARISON_TABLES (
                    ch_all_motif_dirs
                )
                ch_software_versions = ch_software_versions.mix(CREATE_MOTIF_COMPARISON_TABLES.out.versions)
            }
        }
    }

    ch_dt_corrmatrix              = Channel.empty()
    ch_dt_pcadata                 = Channel.empty()
    ch_dt_fpmatrix                = Channel.empty()
    ch_peakqc_frip_mqc            = Channel.empty()
    ch_peakqc_count_mqc           = Channel.empty()
    ch_peakqc_count_consensus_mqc = Channel.empty()
    ch_peakqc_reprod_perc_mqc     = Channel.empty()
    ch_frag_len_hist_mqc          = Channel.empty()
    if(params.run_reporting) {
        if(params.run_igv) {
            ch_igv_fasta = PREPARE_GENOME.out.fasta
                .map { row -> row instanceof List ? row[1] : row }
                .filter { it != null }

            ch_igv_fasta_index = PREPARE_GENOME.out.fasta_index
                .map { row -> row instanceof List ? row[1] : row }
                .filter { it != null }

            /*
            * MODULE: Create igv session (using original, non-downsampled bigwigs)
            */
            // When downsampling is enabled, subtract bigwigs come from downsampled data,
            // so only include them in the main session when NOT downsampling
            if (downsample_enabled) {
                ch_bigwig_for_igv = ch_bigwig_original.collect{it[1]}.ifEmpty([])
                    .collect()
            } else {
                ch_bigwig_for_igv = ch_bigwig_original.collect{it[1]}.ifEmpty([])
                    .mix(ch_bigwig_subtract.collect{it[1]}.flatten().ifEmpty([]))
                    .collect()
            }
            
            IGV_SESSION (
                ch_igv_fasta,
                ch_igv_fasta_index,
                PREPARE_GENOME.out.bed_index,
                //PREPARE_GENOME.out.gtf.collect(),
                ch_peaks_primary_split.igv_a.collect{it[1]}.filter{ it -> it.size() > 1}.ifEmpty([]),
                ch_peaks_secondary.collect{it[1]}.filter{ it -> it.size() > 1}.ifEmpty([]),
                ch_bigwig_for_igv,
                params.igv_sort_by_groups,
                'igv_session'
            )
            //ch_software_versions = ch_software_versions.mix(IGV_SESSION.out.versions)

            /*
            * MODULE: Create downsampled igv session (using downsampled bigwigs + subtract)
            */
            if (downsample_enabled) {
                ch_bigwig_downsampled_for_igv = ch_bigwig_visual.collect{it[1]}.ifEmpty([])
                    .mix(ch_bigwig_subtract.collect{it[1]}.flatten().ifEmpty([]))
                    .collect()

                IGV_SESSION_DOWNSAMPLED (
                    ch_igv_fasta,
                    ch_igv_fasta_index,
                    PREPARE_GENOME.out.bed_index,
                    ch_peaks_primary_split.igv_b.collect{it[1]}.filter{ it -> it.size() > 1}.ifEmpty([]),
                    ch_peaks_secondary.collect{it[1]}.filter{ it -> it.size() > 1}.ifEmpty([]),
                    ch_bigwig_downsampled_for_igv,
                    params.igv_sort_by_groups,
                    'igv_session_downsampled'
                )
            }
        }

        if (params.run_pygenometracks_top10 && params.run_peak_calling) {
            ch_bigwig_rows
                .map { rows -> [bigwig_rows: rows] }
                .set { ch_bigwig_all_for_pygt }

            HOMER_ANNOTATEPEAKS_CONSENSUS (
                ch_consensus_peaks,
                ch_fasta_for_peak_annotation,
                ch_gtf_for_peak_annotation
            )
            ch_software_versions = ch_software_versions.mix(HOMER_ANNOTATEPEAKS_CONSENSUS.out.versions)

            def ch_consensus_peak_chipseeker_annotations
            if (params.run_chipseeker) {
                CHIPSEEKER_ANNOTATE_CONSENSUS (
                    HOMER_ANNOTATEPEAKS_CONSENSUS.out.annot,
                    file("${projectDir}/assets/dummy_file.txt"),
                    ch_gtf_for_peak_annotation
                )
                ch_software_versions = ch_software_versions.mix(CHIPSEEKER_ANNOTATE_CONSENSUS.out.versions)
                ch_consensus_peak_chipseeker_annotations = CHIPSEEKER_ANNOTATE_CONSENSUS.out.annotation
                    .map { meta, anno -> [meta.id, meta, anno] }
            }

            ch_consensus_peak_annotations = HOMER_ANNOTATEPEAKS_CONSENSUS.out.annot
                .map { meta, table -> [meta.id, meta, table] }

            if (!params.run_chipseeker) {
                ch_consensus_peak_chipseeker_annotations = ch_consensus_peak_annotations
                    .map { id, meta, table -> [id, meta, file("${projectDir}/assets/dummy_file.txt")] }
            }

            ch_consensus_peaks
                .combine(ch_bigwig_all_for_pygt)
                .map { row ->
                    def meta = row[0]
                    def bed = row[1]
                    def bigwig_rows = row[2].bigwig_rows ?: []
                    def valid_rows = bigwig_rows.findAll { it instanceof List && it.size() >= 2 && it[0] != null && it[1] != null }

                    def own_rows = valid_rows.findAll { !it[0].is_control && it[0].group == meta.id }
                    def own_bigwigs = own_rows.collect { it[1] }.sort { it.getName() }

                    def control_rows = valid_rows.findAll { it[0].is_control }
                    def control_rep = control_rows ? control_rows.sort { a, b -> a[0].id <=> b[0].id }[0][1] : null

                    def other_rep_rows = valid_rows
                        .findAll { !it[0].is_control && it[0].group != meta.id }
                        .groupBy { it[0].group }
                        .collect { grp, rows -> rows.sort { it[0].id }[0] }
                    def other_bigwigs = other_rep_rows.collect { it[1] }.sort { it.getName() }

                    def tracks = []
                    tracks.addAll(own_bigwigs)
                    tracks.addAll(other_bigwigs)
                    if (control_rep) {
                        tracks << control_rep
                    }

                    [meta.id, meta, bed, tracks, own_bigwigs.size(), control_rep ? 1 : 0]
                }
                .filter { row -> row[2] && row[2].size() > 0 }
                .join(ch_consensus_peak_annotations, by: 0)
                .join(ch_consensus_peak_chipseeker_annotations, by: 0)
                .map { group_id, meta, bed, tracks, own_count, has_control, annot_meta, annot_table, chip_meta, chip_table ->
                    [meta, bed, tracks, own_count, has_control, chip_table, annot_table]
                }
                .set { ch_pygenometracks_top10_input }

            PYGENOMETRACKS_TOP10(
                ch_pygenometracks_top10_input,
                PREPARE_GENOME.out.gtf.first(),
                params.pygt_top_n,
                params.pygt_peak_flank,
                params.pygt_window_bp,
                params.pygt_rank_mode,
                params.pygt_feature_types,
                params.pygt_feature_anchor_window,
                params.pygt_output_format
            )
            ch_software_versions = ch_software_versions.mix(PYGENOMETRACKS_TOP10.out.versions)
        }

        if (params.run_deeptools_heatmaps && params.run_peak_calling) {
            /*
            * CHANNEL: Remove IgG from bigwig channel
            */
            ch_bigwig.filter { it[0].is_control == false }
            .set { ch_bigwig_no_igg }
            // ch_bigwig_no_igg | view

            /*
            * MODULE: Compute DeepTools matrix used in heatmap plotting for Genes
            */
            DEEPTOOLS_COMPUTEMATRIX_GENE (
                ch_bigwig_no_igg,
                PREPARE_GENOME.out.bed.collect()
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX_GENE.out.versions)

            /*
            * MODULE: Calculate DeepTools heatmap
            */
            DEEPTOOLS_PLOTHEATMAP_GENE (
                DEEPTOOLS_COMPUTEMATRIX_GENE.out.matrix
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_PLOTHEATMAP_GENE.out.versions)

            /*
            * CHANNEL: Structure output for join on id
            */
            ch_peaks_summits
            .map { row -> [row[0].id, row ].flatten()}
            .set { ch_peaks_summits_id }
            //ch_peaks_bed_id | view

            /*
            * CHANNEL: Join beds and bigwigs on id
            */
            ch_bigwig_no_igg
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_peaks_summits_id )
            .filter ( it -> it[-1].size() > 1)
            .set { ch_dt_bigwig_summits }
            //ch_dt_peaks | view

            ch_dt_bigwig_summits
            .map { row -> row[1,2] }
            .set { ch_ordered_bigwig }
            //ch_ordered_bigwig | view

            ch_dt_bigwig_summits
            .map { row -> row[-1] }
            .set { ch_ordered_peaks_max }
            //ch_ordered_peaks_max | view

            /*
            * MODULE: Compute DeepTools matrix used in heatmap plotting for Peaks
            */

            DEEPTOOLS_COMPUTEMATRIX_PEAKS (
                ch_ordered_bigwig,
                ch_ordered_peaks_max
            )

            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.versions)
            //EXAMPLE CHANNEL STRUCT: [[META], MATRIX]
            //DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix | view

            /*
            * MODULE: Calculate DeepTools heatmap
            */
            DEEPTOOLS_PLOTHEATMAP_PEAKS (
                DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
            )
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_PLOTHEATMAP_PEAKS.out.versions)

            if(params.dt_calc_all_matrix) {
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
                    PREPARE_GENOME.out.bed.toSortedList()
                )

                /*
                * MODULE: Calculate DeepTools heatmap for all samples
                */
                DEEPTOOLS_PLOTHEATMAP_GENE_ALL (
                    DEEPTOOLS_COMPUTEMATRIX_GENE_ALL.out.matrix
                )
            }
        }

        if(params.run_deeptools_qc) {
            /*
            * SUBWORKFLOW: Run suite of deeptools QC on bam files
            */
            DEEPTOOLS_QC (
                ch_samtools_bam,
                ch_samtools_bai,
                params.dt_qc_corr_method
            )
            ch_dt_corrmatrix     = DEEPTOOLS_QC.out.correlation_matrix
            ch_dt_pcadata        = DEEPTOOLS_QC.out.pca_data
            ch_dt_fpmatrix       = DEEPTOOLS_QC.out.fingerprint_matrix
            ch_software_versions = ch_software_versions.mix(DEEPTOOLS_QC.out.versions)
        }

        /*
        * CHANNEL: Filter bais for target only
        */
        ch_samtools_bai.filter { it -> it[0].is_control == false }
        .set { ch_bai_target }
        //ch_bai_target | view

        if (params.run_peak_qc && params.run_peak_calling) {
            /*
            * CHANNEL: Filter flagstat for target only
            */
            ch_samtools_flagstat_for_frip.filter { it -> it[0].is_control == false }
            .set { ch_flagstat_target }
            //ch_flagstat_target | view

            /*
            * SUBWORKFLOW: Extract fragments from bam files for fragment-based FRiP score
            */
            EXTRACT_FRAGMENTS (
                ch_bam_target
            )

            /*
            * SUBWORKFLOW: Run suite of peak QC on peaks
            */
            PEAK_QC(
                ch_peaks_primary_split.peak_qc,
                AWK_NAME_PEAK_BED.out.file,
                ch_consensus_peaks,
                ch_consensus_peaks_unfilt,
                EXTRACT_FRAGMENTS.out.bed,
                ch_flagstat_target,
                params.min_frip_overlap,
                ch_frip_score_header_multiqc,
                ch_peak_counts_header_multiqc,
                ch_peak_counts_consensus_header_multiqc,
                ch_peak_reprod_header_multiqc
            )
            ch_peakqc_frip_mqc             = PEAK_QC.out.primary_frip_mqc
            ch_peakqc_count_mqc            = PEAK_QC.out.primary_count_mqc
            ch_peakqc_count_consensus_mqc  = PEAK_QC.out.consensus_count_mqc
            ch_peakqc_reprod_perc_mqc      = PEAK_QC.out.reprod_perc_mqc
            ch_software_versions           = ch_software_versions.mix(PEAK_QC.out.versions)
        }

        // Run PeakSignalProfiler (optional)
        if (params.run_peak_signal_profiler && params.run_peak_calling) {
            // If the user supplied a gene_bed on the CLI use that, otherwise
            // pass the BED produced by the PREPARE_GENOME subworkflow. This
            // avoids calling `file(null)` and lets the pipeline use its own
            // generated `genes.bed` when no external BED is provided.
            if (params.gene_bed) {
                PEAK_SIGNAL_PROFILER( file(params.input), file(params.gene_bed), PREPARE_GENOME.out.fasta_index.map{it[1]}.first() )
            } else {
                PEAK_SIGNAL_PROFILER( file(params.input), PREPARE_GENOME.out.bed, PREPARE_GENOME.out.fasta_index.map{it[1]}.first() )
            }
            // NOTE: intentionally do NOT mix PSP `versions` into ch_software_versions
            // to avoid compile-time channel evaluation issues; PSP `versions.yml`
            // will still be published alongside pipeline outputs.
        }

        //ch_peakqc_reprod_perc_mqc | view

        /*
        * CHANNEL: Combine bam and bai files on id
        */

        ch_bam_target.map { row -> [row[0].id, row ].flatten()}
        .join ( ch_bai_target.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4]] }
        .set { ch_bam_bai }
        // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
        //ch_bam_bai | view

        /*
        * MODULE: Calculate fragment lengths
        */
        SAMTOOLS_CUSTOMVIEW (
            ch_bam_bai
        )
        ch_software_versions = ch_software_versions.mix(SAMTOOLS_CUSTOMVIEW.out.versions)
        //SAMTOOLS_CUSTOMVIEW.out.tsv | view

        /*
        * CHANNEL: Prepare data for generate reports
        */
        // Make sure files are always in order for resume
        ch_frag_len = SAMTOOLS_CUSTOMVIEW.out.tsv
        .toSortedList { row -> row[0].id }
        .map {
            list ->
            def output = []
            list.each{ v -> output.add(v[1]) }
            output
        }
        //ch_frag_len | view

        /*
        * MODULE: Calculate fragment length histogram for mqc
        */
        FRAG_LEN_HIST(
            ch_frag_len,
            ch_frag_len_header_multiqc
        )
        ch_frag_len_hist_mqc = FRAG_LEN_HIST.out.frag_len_mqc
        ch_software_versions = ch_software_versions.mix(FRAG_LEN_HIST.out.versions)
    }
    //ch_frag_len_hist_mqc | view


    if (params.run_multiqc) {
        workflow_summary    = WorkflowCutandrun.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        /*
        * MODULE: Collect software versions used in pipeline
        */
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_software_versions.unique().collectFile()
        )

        /*
        * MODULE: Multiqc
        */
        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect(),
            ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yml"),
            FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_spikein_log.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_metrics.collect{it[1]}.ifEmpty([]),
            ch_preseq_output.collect{it[1]}.ifEmpty([]),
            ch_dt_corrmatrix.collect{it[1]}.ifEmpty([]),
            ch_dt_pcadata.collect{it[1]}.ifEmpty([]),
            ch_dt_fpmatrix.collect{it[1]}.ifEmpty([]),
            ch_peakqc_count_mqc.collect{it[1]}.ifEmpty([]),
            ch_peakqc_frip_mqc.collect{it[1]}.ifEmpty([]),
            ch_peakqc_count_consensus_mqc.collect{it[1]}.ifEmpty([]),
            ch_peakqc_reprod_perc_mqc.collect().ifEmpty([]),
            ch_frag_len_hist_mqc.collect().ifEmpty([]),
            ch_linear_duplication_mqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }

    /*
    * MODULE: Generate CUT_and_RUN HTML Report
    */
    if(params.with_report) {
        REPORT (
            Channel.of(
                params.outdir.startsWith('/') ? params.outdir : "${workflow.launchDir}/${params.outdir}"
            )
        )
        ch_software_versions = ch_software_versions.mix(REPORT.out.versions)
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
