#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_split_metagene.sh
# Authors : Katharina (Nina) Hayer, co-created with Claude Fable 5.1 (Anthropic)
# Version : 1.1.0
# Date    : 2026-10-06
# Purpose : computeMatrix (deepTools) twice - TSS-anchored and TES-anchored -
#           over all replicate bigWigs of one condition and two gene groups,
#           then draw the split figure (replicates averaged per gene).
# Usage   : bash run_split_metagene.sh <bw_dir> <condition> <ic.bed> <il.bed> \
#                <out_prefix> [subtract|log2ratio] [label] [threads]
# Example : bash run_split_metagene.sh \
#             ../results_<condition>_final/03_peak_calling/07_bigwig_minus_igg \
#             COND_A genes.intron_containing.bed genes.intronless.bed \
#             figs/COND_A_host subtract "COND_A (host, minus IgG)" 8
# -----------------------------------------------------------------------------
#SBATCH --job-name=split_metagene
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail

BW_DIR=${1:?dir with *.subtract.bigWig / *.log2ratio.bigWig}
COND=${2:?condition prefix, e.g. COND_A}
BED_IC=${3:?intron-containing BED}
BED_IL=${4:?intronless BED}
OUT=${5:?output prefix}
TYPE=${6:-subtract}                 # subtract | log2ratio
LABEL=${7:-$COND}
THREADS=${8:-${SLURM_CPUS_PER_TASK:-8}}

mapfile -t BWS < <(ls "${BW_DIR}"/"${COND}"_R*."${TYPE}".bigWig 2>/dev/null | sort -V)
if [ "${#BWS[@]}" -eq 0 ]; then
    echo "no files match ${BW_DIR}/${COND}_R*.${TYPE}.bigWig" >&2; exit 1
fi
LABELS=(); for b in "${BWS[@]}"; do LABELS+=("$(basename "$b" ".${TYPE}.bigWig")"); done
echo "replicates (${#BWS[@]}): ${LABELS[*]}" >&2

UP_TSS=1000;  DOWN_TSS=5000     # -1 kb .. +5 kb around TSS
UP_TES=5000;  DOWN_TES=1000     # -5 kb .. +1 kb around TES
BIN=50

common=( -S "${BWS[@]}" -R "$BED_IC" "$BED_IL"
         --binSize "$BIN" --missingDataAsZero
         --samplesLabel "${LABELS[@]}" -p "$THREADS" )

mkdir -p "$(dirname "$OUT")"
echo "[1/3] computeMatrix TSS" >&2
computeMatrix reference-point --referencePoint TSS \
    -b "$UP_TSS" -a "$DOWN_TSS" "${common[@]}" -o "${OUT}.TSS.mat.gz"

echo "[2/3] computeMatrix TES" >&2
computeMatrix reference-point --referencePoint TES \
    -b "$UP_TES" -a "$DOWN_TES" "${common[@]}" -o "${OUT}.TES.mat.gz"

echo "[3/3] plot" >&2
python3 "$(dirname "$0")/plot_split_metagene.py" \
    --tss "${OUT}.TSS.mat.gz" --tes "${OUT}.TES.mat.gz" \
    --signal "$TYPE" --title "$LABEL" --out "${OUT}.split_metagene"

echo "done: ${OUT}.split_metagene.{pdf,png}" >&2
