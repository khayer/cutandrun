#!/usr/bin/env bash
set -euo pipefail

SAMPLESHEET="$1"
OUTDIR="$2"
ANNOTATION="$3"
GENOME="$4"
CORES="${5:-1}"
PSP_SIF="${6:-}"
PSP_DIR="${7:-}"

# Preflight: check bigWig files
LOG=psp_preflight_files.log
echo "[PSP PRECHECK] Verifying bigWig files for all samples..." > "$LOG"
missing_files=0
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r group replicate fastq1 fastq2 control; do
    group="${group//\"/}"
    replicate="${replicate//\"/}"
    sample_base="${group}_R${replicate}"
    found=0
    if [ -f "$OUTDIR/03_peak_calling/03_bed_to_bigwig/${sample_base}.bigWig" ]; then
        found=1
    elif [ -f "$OUTDIR/03_peak_calling/03_bed_to_bigwig/${group}_R1.bigWig" ]; then
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

# ...existing code for samplesheet conversion and PSP run can be added here...
