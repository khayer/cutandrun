/*
 * Subworkflow wrapper for PeakSignalProfiler multi-sample run
 */

include { PEAKSIGNALPROFILER_RUN } from '../../modules/local/peak_signal_profiler'

workflow PEAK_SIGNAL_PROFILER {
    take:
        samplesheet
        annotation
        genome
        psp_dir

    main:
        PEAKSIGNALPROFILER_RUN( samplesheet, annotation, genome, psp_dir )

    emit:
        PEAKSIGNALPROFILER_RUN.out.psp_out
        PEAKSIGNALPROFILER_RUN.out.versions
}
