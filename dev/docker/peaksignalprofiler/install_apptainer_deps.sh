#!/usr/bin/env bash
set -euo pipefail

# Install packages required to build Apptainer/Singularity on Ubuntu runners
# (use in CI or on a build host before compiling Apptainer).

sudo apt-get update -y
sudo apt-get install -y \
  build-essential automake autoconf libtool m4 pkg-config \
  libseccomp-dev libgpgme-dev libssl-dev uuid-dev squashfs-tools \
  libarchive-dev libglib2.0-dev libcap-dev

echo "Apptainer build deps installed."