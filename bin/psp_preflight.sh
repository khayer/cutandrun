#!/usr/bin/env bash
set -euo pipefail

PROJECTDIR="$1"
OUTDIR="$2"
SAMPLESHEET="$3"
ANNOTATION="$4"
GENOME="$5"
CORES="${6:-1}"
PSP_SIF="${7:-}"
PSP_DIR="${8:-}"
# derive absolute base directory for results; if OUTDIR is relative, prefix projectDir
if [[ "$OUTDIR" = /* ]]; then
    BASEDIR="$OUTDIR"
else
    BASEDIR="$PROJECTDIR/$OUTDIR"
fi

# Preflight: check bigWig files
LOG=psp_preflight_files.log
echo "[PSP PRECHECK] Verifying bigWig files for all samples..." > "$LOG"
missing_files=0
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r group replicate fastq1 fastq2 control; do
    group="${group//\"/}"
    replicate="${replicate//\"/}"
    sample_base="${group}_R${replicate}"
    found=0
    if [ -f "$BASEDIR/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig" ]; then
        found=1
    elif [ -f "$BASEDIR/03_peak_calling/03_bed_to_bigwig/${group}_R1.bigWig" ]; then
        found=1
    fi
    if [ $found -eq 0 ]; then
        echo "MISSING: $OUTDIR/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig or ${group}_R1.bigWig" >> "$LOG"
        missing_files=1
    fi
done
if [ $missing_files -eq 1 ]; then
    echo "ERROR: One or more bigWig files are missing. See $LOG for details." >&2
    exit 2
fi


# convert samplesheet to PSP format if needed
PSP_SHEET="psp_$(basename "$SAMPLESHEET")"
if head -n1 "$SAMPLESHEET" | grep -q "peaks_bed"; then
    cp "$SAMPLESHEET" "$PSP_SHEET"
else
    echo '"group","replicate","peaks_bed","signal_bw","control_bw"' > "$PSP_SHEET"
    tail -n +2 "$SAMPLESHEET" | while IFS=, read -r group replicate fastq1 fastq2 control; do
        group="${group//\"/}"
        replicate="${replicate//\"/}"
        control="${control//\"/}"
        sample_base="${group}_R${replicate}"

        peaks=""
        if [ -f "$BASEDIR/03_peak_calling/04_called_peaks/seacr/${sample_base}.seacr.peaks.stringent.bed" ]; then
            peaks="$BASEDIR/03_peak_calling/04_called_peaks/seacr/${sample_base}.seacr.peaks.stringent.bed"
        elif [ -f "$BASEDIR/03_peak_calling/04_called_peaks/macs2/${sample_base}.macs2.peaks.cut.bed" ]; then
            peaks="$BASEDIR/03_peak_calling/04_called_peaks/macs2/${sample_base}.macs2.peaks.cut.bed"
        elif [ -f "$BASEDIR/03_peak_calling/05_consensus_peaks/${group}.seacr.consensus.peaks.awk.bed" ]; then
            peaks="$BASEDIR/03_peak_calling/05_consensus_peaks/${group}.seacr.consensus.peaks.awk.bed"
        fi

        signal=""
        if [ -f "$BASEDIR/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig" ]; then
            signal="$BASEDIR/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig"
        elif [ -f "$BASEDIR/03_peak_calling/03_bed_to_bigwig/${group}_R1.bigWig" ]; then
            signal="$BASEDIR/03_peak_calling/03_bed_to_bigwig/${group}_R1.bigWig"
        fi

        ctrl_bw=""
        if [ -n "$control" ]; then
            if [ -f "$BASEDIR/03_peak_calling/03_bed_to_bigwig/${control}_R${replicate}.bigWig" ]; then
                ctrl_bw="$BASEDIR/03_peak_calling/03_bed_to_bigwig/${control}_R${replicate}.bigWig"
            elif [ -f "$BASEDIR/03_peak_calling/03_bed_to_bigwig/${control}_R1.bigWig" ]; then
                ctrl_bw="$BASEDIR/03_peak_calling/03_bed_to_bigwig/${control}_R1.bigWig"
            fi
        fi

        if [ -z "$peaks" ] || [ -z "$signal" ]; then
            # some controls (e.g. igg_ctrl) won't produce peaks/signals, skip them
            if [[ "$group" =~ ^igg ]] || [[ "$group" =~ ^GK ]] || [[ "$group" =~ ^GM ]] || [[ "$group" =~ ^H3 ]]; then
                echo "WARNING: skipping control sample ${sample_base} (no peaks or signal)" >&2
                continue
            fi
            echo "ERROR: cannot locate peaks or signal for ${sample_base}" >&2
            exit 1
        fi
        # wrap paths in quotes, make them relative to OUTDIR
        echo "\"${group}\",\"${replicate}\",\"${peaks}\",\"${signal}\",\"${ctrl_bw}\"" >> "$PSP_SHEET"
    done
fi

# run PSP multisample
mkdir -p psp_out
# default to running inside container with staged repo if available

Rscript /app/scripts/02_run_multisample.R \
    --samplesheet "${PSP_SHEET}" \
    --annotation "$ANNOTATION" \
    --genome "$GENOME" \
    --outdir psp_out \
    --cores $CORES \
> psp_run.log 2>&1


# optional profiling summary
if [ -f psp_out/multisample_rprof.out ]; then
    Rscript -e "summaryRprof('psp_out/multisample_rprof.out')" \
        > psp_out/multisample_rprof.summary.txt 2>/dev/null || true
fi

cat <<-END_VERSIONS > versions.yml
"${task.process:-PEAKSIGNALPROFILER_RUN}":
    peaksignalprofiler_sif: ${PSP_SIF:-}
    r_version: $(Rscript --version 2>&1 | sed -n '1p')
END_VERSIONS
