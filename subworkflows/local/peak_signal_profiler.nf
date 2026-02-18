/*
 * Subworkflow wrapper for PeakSignalProfiler multi-sample run
 */

include { PEAKSIGNALPROFILER_RUN } from '../../modules/local/peak_signal_profiler'

workflow PEAK_SIGNAL_PROFILER {
    take:
        samplesheet
        annotation
        genome

    main:
        PEAKSIGNALPROFILER_RUN( samplesheet, annotation, genome )

    emit:
        PEAKSIGNALPROFILER_RUN.out.psp_out
        PEAKSIGNALPROFILER_RUN.out.versions
}
