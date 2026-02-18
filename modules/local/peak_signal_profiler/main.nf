process PEAKSIGNALPROFILER_RUN {
    tag "${samplesheet.baseName}"
    label 'process_single'

    // Use the Singularity image supplied by the user (placeholder by default)
    container "${params.psp_sif ?: '/path/to/peaksignalprofiler.sif'}"

    input:
    path samplesheet
    path annotation
    path genome

    output:
    path "psp_out"                , emit: psp_out
    path "*.log"                  , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Ensure the R script repository and sif are set via params when running.
    def psp_dir = params.psp_dir ?: '/path/to/peak-signal-profiler'
    def sif     = params.psp_sif ?: '/path/to/peaksignalprofiler.sif'

    """
    set -euo pipefail || true
    mkdir -p psp_out

    echo "[PEAKSIGNALPROFILER] samplesheet=${samplesheet} annotation=${annotation} genome=${genome} cores=${task.cpus}" > psp_run.info

    # Run the multisample R script inside the provided Singularity image.
    singularity exec ${sif} \
        Rscript ${psp_dir}/scripts/02_run_multisample.R \
        --samplesheet ${samplesheet} \
        --annotation ${annotation} \
        --genome ${genome} \
        --outdir psp_out \
        --cores ${task.cpus} \
        > psp_run.log 2>&1 || true

    # If R profiling was produced, summarize it (optional)
    if [ -f psp_out/multisample_rprof.out ]; then
        singularity exec ${sif} Rscript -e "summaryRprof('psp_out/multisample_rprof.out')" \
            > psp_out/multisample_rprof.summary.txt 2>/dev/null || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peaksignalprofiler_sif: ${sif}
        r_version: \$(singularity exec ${sif} Rscript --version 2>&1 | sed -n '1p')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p psp_out && touch psp_out/sample_multisample_plot.png psp_out/multisample_table.tsv
    echo "stub" > versions.yml
    """
}
