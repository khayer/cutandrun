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

    script:
    // Ensure the SIF is set via params when running. The SIF contains the
    // PSP repository under /app so a `psp_dir` parameter is not required.
    def sif     = params.psp_sif ?: '/path/to/peaksignalprofiler.sif'

    """
    set -euo pipefail || true
    mkdir -p psp_out

    # If a local staged repository exists use that; otherwise run the
    # bundled scripts inside the SIF at /app. The concrete run-mode is
    # chosen later when invoking the Rscript.

    echo "[PEAKSIGNALPROFILER] samplesheet=${samplesheet} annotation=${annotation} genome=${genome} cores=${task.cpus}" > psp_run.info

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

    # Prepare a PSP-formatted samplesheet in the task work dir (do **not**
    # overwrite the original pipeline samplesheet). If a staged repository
    # directory `peak-signal-profiler` exists in the work dir we also copy
    # the prepared sheet into it so repository-root runs still work.
    PSP_SHEET="psp_${samplesheet.getName()}"

    if head -n1 "${samplesheet.getName()}" | grep -q "peaks_bed"; then
        # already PSP-format -> copy the original into PSP_SHEET
        cp "${samplesheet.getName()}" "${PSP_SHEET}"
    else
        results_root="${params.outdir}/03_peak_calling"
        echo "group,replicate,peaks_bed,signal_bw,control_bw" > "${PSP_SHEET}"
        tail -n +2 "${samplesheet.getName()}" | while IFS=, read -r group replicate fastq1 fastq2 control; do
            # Basic normalise: strip double-quotes from fields
            group="${group//\"/}"
            replicate="${replicate//\"/}"
            control="${control//\"/}"
            sample_base="${group}_R${replicate}"

            # 1) peaks: seacr -> macs2 -> consensus (search primary outdir then fallback to repo `results/`)
            peaks=""
            if [ -f "${params.outdir}/03_peak_calling/04_called_peaks/seacr/${sample_base}.seacr.peaks.stringent.bed" ]; then
                peaks="${params.outdir}/03_peak_calling/04_called_peaks/seacr/${sample_base}.seacr.peaks.stringent.bed"
            elif [ -f "${new File('results').absolutePath}/03_peak_calling/04_called_peaks/seacr/${sample_base}.seacr.peaks.stringent.bed" ]; then
                peaks="${new File('results').absolutePath}/03_peak_calling/04_called_peaks/seacr/${sample_base}.seacr.peaks.stringent.bed"
            elif [ -f "${params.outdir}/03_peak_calling/04_called_peaks/macs2/${sample_base}.macs2.peaks.cut.bed" ]; then
                peaks="${params.outdir}/03_peak_calling/04_called_peaks/macs2/${sample_base}.macs2.peaks.cut.bed"
            elif [ -f "${new File('results').absolutePath}/03_peak_calling/04_called_peaks/macs2/${sample_base}.macs2.peaks.cut.bed" ]; then
                peaks="${new File('results').absolutePath}/03_peak_calling/04_called_peaks/macs2/${sample_base}.macs2.peaks.cut.bed"
            elif [ -f "${params.outdir}/03_peak_calling/05_consensus_peaks/${group}.seacr.consensus.peaks.awk.bed" ]; then
                peaks="${params.outdir}/03_peak_calling/05_consensus_peaks/${group}.seacr.consensus.peaks.awk.bed"
            elif [ -f "${new File('results').absolutePath}/03_peak_calling/05_consensus_peaks/${group}.seacr.consensus.peaks.awk.bed" ]; then
                peaks="${new File('results').absolutePath}/03_peak_calling/05_consensus_peaks/${group}.seacr.consensus.peaks.awk.bed"
            fi

            # 2) signal bigWig (search primary outdir then fallback to repo `results/`)
            signal=""
            if [ -f "${params.outdir}/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig" ]; then
                signal="${params.outdir}/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig"
            elif [ -f "${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig" ]; then
                signal="${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig"
            elif [ -f "${params.outdir}/03_peak_calling/03_bed_to_bigwig_visual/${sample_base}.bigWig" ]; then
                signal="${params.outdir}/03_peak_calling/03_bed_to_bigwig_visual/${sample_base}.bigWig"
            elif [ -f "${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig_visual/${sample_base}.bigWig" ]; then
                signal="${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig_visual/${sample_base}.bigWig"
            fi

            # 3) control bigWig (try same replicate, else R1) - search primary outdir then fallback to repo `results/`
            ctrl_bw=""
            if [ -n "${control}" ]; then
                if [ -f "${params.outdir}/03_peak_calling/03_bed_to_bigwig/${control}_R${replicate}.bigWig" ]; then
                    ctrl_bw="${params.outdir}/03_peak_calling/03_bed_to_bigwig/${control}_R${replicate}.bigWig"
                elif [ -f "${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig/${control}_R${replicate}.bigWig" ]; then
                    ctrl_bw="${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig/${control}_R${replicate}.bigWig"
                elif [ -f "${params.outdir}/03_peak_calling/03_bed_to_bigwig/${control}_R1.bigWig" ]; then
                    ctrl_bw="${params.outdir}/03_peak_calling/03_bed_to_bigwig/${control}_R1.bigWig"
                elif [ -f "${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig/${control}_R1.bigWig" ]; then
                    ctrl_bw="${new File('results').absolutePath}/03_peak_calling/03_bed_to_bigwig/${control}_R1.bigWig"
                fi
            fi

            # Fail early if we couldn't find required files
            if [ -z "${peaks}" ]; then
                echo "ERROR: cannot locate peaks BED for sample ${sample_base}. Searched seacr/macs2/consensus under ${params.outdir}/03_peak_calling" >&2
                exit 1
            fi
            if [ -z "${signal}" ]; then
                echo "ERROR: cannot locate signal bigWig for sample ${sample_base} under ${params.outdir}/03_peak_calling/03_bed_to_bigwig" >&2
                exit 1
            fi

            echo "${group},${replicate},${peaks},${signal},${ctrl_bw}" >> "${PSP_SHEET}"
        done
    fi

    # If a staged PSP repo exists, copy the prepared sheet into it so the
    # repository-root run mode still works.
    if [ -d peak-signal-profiler ]; then
        cp "${PSP_SHEET}" peak-signal-profiler/"${samplesheet.getName()}" 2>/dev/null || true
    fi

    # Run the multisample R script. Prefer a staged repository if present,
    # otherwise run the bundled script inside the SIF at /app.
    if [ -d peak-signal-profiler ]; then
        ( cd "peak-signal-profiler" && \
            Rscript "scripts/02_run_multisample.R" \
                --samplesheet "${samplesheet.getName()}" \
                --annotation ../${annotation} \
                --genome ../${genome} \
                --outdir ../psp_out \
                --cores ${task.cpus} \
        ) > psp_run.log 2>&1
    else
        Rscript /app/scripts/02_run_multisample.R \
            --samplesheet "${PSP_SHEET}" \
            --annotation ${annotation} \
            --genome ${genome} \
            --outdir psp_out \
            --cores ${task.cpus} \
        > psp_run.log 2>&1
    fi

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
