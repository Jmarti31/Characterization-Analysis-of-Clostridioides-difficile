#!/usr/bin/env bash
# analyzeTrinityDeNovo.sh
# Usage: bash scripts/analyzeTrinityDeNovo.sh 1>results/trinity_de_novo_stats.txt 2>results/logs/analyzeTrinityDeNovo.err

# Generate Transcriptome Assembly Statistics
TrinityStats.pl results/trinity_de_novo/Trinity.fasta
