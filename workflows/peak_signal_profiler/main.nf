nextflow.enable.dsl=2

include { PEAK_SIGNAL_PROFILER } from '../../subworkflows/local/peak_signal_profiler.nf'

workflow {
    // CLI params (defaults are placeholders - user should override)
    params.samplesheet = params.samplesheet ?: null
    params.annotation  = params.annotation  ?: null
    params.genome      = params.genome      ?: null
    params.psp_sif     = params.psp_sif     ?: '/path/to/peaksignalprofiler.sif'
    params.psp_dir     = params.psp_dir     ?: '/path/to/peak-signal-profiler'

    if (!params.samplesheet) {
        log.error "--samplesheet is required"
        System.exit(1)
    }

    if (!params.annotation) {
        log.error "--annotation is required"
        System.exit(1)
    }

    if (!params.genome) {
        log.error "--genome is required"
        System.exit(1)
    }

    PEAK_SIGNAL_PROFILER( file(params.samplesheet), file(params.annotation), file(params.genome) )
}
