// Module: PEAKSIGNALPROFILER_RUN
// Author: Katharina Hayer
// Co-created with: GitHub Copilot (Claude Sonnet 4.5)

process PEAKSIGNALPROFILER_RUN {
    tag "${samplesheet.baseName}"
    label 'process_single'

    // Publish PeakSignalProfiler outputs into the canonical results reporting folder
    // (makes `psp_out/` visible under results/<outdir>/04_reporting/peaksignalprofiler)
    publishDir "${params.outdir}/04_reporting/peaksignalprofiler", mode: 'copy'

    // Use the Singularity image supplied by the user (placeholder by default)
    container "${params.psp_sif ?: '/path/to/peaksignalprofiler.sif'}"

    input:
    file samplesheet
    path annotation
    path genome

    output:
    path "psp_out"                , emit: psp_out
    path "*.log"                  , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // Ensure the SIF is set via params when running. The SIF contains the
    // PSP repository under /app so a `psp_dir` parameter is not required.
    def sif     = params.psp_sif ?: '/path/to/peaksignalprofiler.sif'

    script:
    """
    bash bin/psp_preflight.sh "${samplesheet.getName()}" "${params.outdir}" "${annotation}" "${genome}" "${task.cpus}" "${params.psp_sif}" "${params.psp_dir:-}"  
    """
    
    stub:
    """
    mkdir -p psp_out && touch psp_out/sample_multisample_plot.png psp_out/multisample_table.tsv
    echo "stub" > versions.yml
    """
}
