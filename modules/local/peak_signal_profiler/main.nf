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
    path psp_src

    output:
    path "psp_out"                , emit: psp_out
    path "*.log"                  , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Ensure the R script repository and sif are set via params when running.
    // NOTE: prefer `psp_src` if staged by the workflow (handled at shell runtime);
    // here we keep a params fallback for callers that do not stage the dir.
    def psp_dir = params.psp_dir ?: '/path/to/peak-signal-profiler'
    def sif     = params.psp_sif ?: '/path/to/peaksignalprofiler.sif'

    """
    set -euo pipefail || true
    mkdir -p psp_out

    # Prefer staged `psp_src` (workflow-input) or the preserved directory name
    # inside the work dir; otherwise fall back to the parameter value expanded
    # from Groovy above.
    if [ -n "${psp_src}" ] && [ -d "${psp_src}" ]; then
        psp_dir="${psp_src}"
    elif [ -d peak-signal-profiler ]; then
        psp_dir="peak-signal-profiler"
    else
        psp_dir="${psp_dir}"
    fi

    echo "[PEAKSIGNALPROFILER] samplesheet=${samplesheet} annotation=${annotation} genome=${genome} cores=${task.cpus} psp_dir=\${psp_dir}" > psp_run.info

    # Preflight: verify required R/Bioconductor packages are available inside
    # the container. Fail with a clear error if any are missing so the user
    # knows to rebuild the SIF with the missing packages.
    required_pkgs="rtracklayer,GenomicRanges,IRanges,S4Vectors,ggplot2,optparse,cowplot,png,ggnewscale"
    REQUIRED_PKGS="\${required_pkgs}" Rscript -e "pkgs <- strsplit(Sys.getenv('REQUIRED_PKGS'),',')[[1]]; missing <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly=TRUE))]; if(length(missing)){ cat('MISSING_R_PKGS:' , paste(missing, collapse=','), '\n'); quit(status=2) } else {cat('R_PKGS_OK\n') }" > psp_preflight.log 2>&1 || {
        echo "ERROR: Required R packages missing inside container. Rebuild SIF to include: \${required_pkgs}" >&2
        echo "See psp_preflight.log for details." >&2
        cat psp_preflight.log >&2 || true
        exit 2
    }

    # Run the multisample R script from the PSP repository root so relative
    # paths inside the package (e.g. R/config.R) resolve correctly.
    ( cd "\${psp_dir}" && \
        Rscript "scripts/02_run_multisample.R" \
            --samplesheet ../${samplesheet} \
            --annotation ../${annotation} \
            --genome ../${genome} \
            --outdir ../psp_out \
            --cores ${task.cpus} \
    ) > psp_run.log 2>&1

    # If R profiling was produced, summarize it (optional)
    if [ -f psp_out/multisample_rprof.out ]; then
        Rscript -e "summaryRprof('psp_out/multisample_rprof.out')" \
            > psp_out/multisample_rprof.summary.txt 2>/dev/null || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peaksignalprofiler_sif: ${sif}
        r_version: \$(Rscript --version 2>&1 | sed -n '1p')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p psp_out && touch psp_out/sample_multisample_plot.png psp_out/multisample_table.tsv
    echo "stub" > versions.yml
    """
}
