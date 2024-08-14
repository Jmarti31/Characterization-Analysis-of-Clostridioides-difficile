#!/usr/bin/env bash
# analyzeTrinity.sh
# Usage: bash scripts/analyzeTrinity.sh 1>results/trinity_guided_stats.txt 2>results/logs/analyzeTrinity.err

# Generate Transcriptome Assembly Statistics 
TrinityStats.pl results/trinity_guided/Trinity-GG.fasta
