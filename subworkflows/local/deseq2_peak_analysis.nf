/*
 * Subworkflow: DESeq2 Peak Differential Analysis
 * Performs DESeq2 differential abundance testing on peak data
 */

include { DESEQ2_PEAK_ANALYSIS } from "../../modules/local/r/deseq2_peak_analysis/main"

workflow DESEQ2_PEAK_DIFF_ANALYSIS {
    take:
    counts_table_channel
    gtf_file

    main:
    ch_versions = Channel.empty()

    // Extract table file from channel
    counts_table_channel
        .map { meta, table -> table }
        .set { counts_table }

    // Parse contrasts file
    if (!params.deseq2_contrasts) {
        error "run_deseq2_peak_analysis=true requires --deseq2_contrasts file"
    }

    def deseq2_contrasts = []
    def contrasts_file = file(params.deseq2_contrasts)
    if (!contrasts_file.exists()) {
        error "deseq2_contrasts file not found: ${params.deseq2_contrasts}"
    }

    contrasts_file.readLines().eachWithIndex { line, idx ->
        def trimmed = line.trim()
        if (!trimmed || trimmed.startsWith("#")) {
            return
        }

        def fields = trimmed.split(/\t/)
        if (fields.size() < 2) {
            error "Invalid deseq2_contrasts row ${idx+1}: '${line}'\nExpected TSV format: treatment_group control_group\nExample: RAD51_B02_26 RAD51_DMSO_26"
        }

        def treatment_group = fields[0].trim()
        def control_group = fields[1].trim()

        deseq2_contrasts << [
            id: "${treatment_group}_vs_${control_group}",
            treatment: treatment_group,
            control: control_group
        ]
    }

    if (!deseq2_contrasts.isEmpty()) {
        log.info "DESeq2: Processing ${deseq2_contrasts.size()} contrast(s)"

        Channel
            .fromList(deseq2_contrasts)
            .combine(counts_table)
            .map { contrast, table ->
                def meta = [
                    id: contrast.id,
                    treatment: contrast.treatment,
                    control: contrast.control
                ]
                [meta, table, file("${params.outdir}/03_peak_calling/04_called_peaks/macs2")]
            }
            .set { ch_deseq2_inputs }

        DESEQ2_PEAK_ANALYSIS (
            ch_deseq2_inputs,
            gtf_file,
            file(params.input)
        )
        ch_versions = ch_versions.mix(DESEQ2_PEAK_ANALYSIS.out.versions)
    }

    emit:
    results = DESEQ2_PEAK_ANALYSIS.out.results
    significant = DESEQ2_PEAK_ANALYSIS.out.significant
    bedgraph = DESEQ2_PEAK_ANALYSIS.out.bedgraph
    volcano_plot = DESEQ2_PEAK_ANALYSIS.out.volcano_plot
    ma_plot = DESEQ2_PEAK_ANALYSIS.out.ma_plot
    versions = ch_versions
}
