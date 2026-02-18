#!/usr/bin/env bash
# Local helper: build SIF from published Docker image (requires Singularity/Apptainer)
# Usage: ./build_sif_instructions.sh <ghcr_owner>

OWNER=${1:-"<OWNER>"}
IMAGE=ghcr.io/${OWNER}/peaksignalprofiler:latest
SIF_NAME=psp-1.0.0.sif

echo "Building SIF from ${IMAGE} -> ${SIF_NAME}"
# Example (requires apptainer/singularity installed):
# singularity build ${SIF_NAME} docker://${IMAGE}

echo "Run the above singularity build command on a host with Singularity/Apptainer installed."
