#!/usr/bin/env bash
# ClosBuild.sh
# Usage: bash scripts/ClosBuild.sh 1>results/logs/ClosBuild.log 2>results/logs/ClosBuild.err &

# Creates a reference genome database
gmap_build -D data \
-d ClostridioidesGmapDb \
/home/martinez.jam/GCF_000009205.2_ASM920v2_genomic.fna
