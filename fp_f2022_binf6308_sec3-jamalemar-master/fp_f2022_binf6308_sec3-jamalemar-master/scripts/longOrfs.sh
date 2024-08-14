#!/usr/bin/env bash
# longOrfs.sh
# Usage: bash scripts/longOrfs.sh <transcriptome> <result folder> 1>results/logs/longOrfs.log 2>results/logs/longOrfs.err

# <transcriptome> might be data/trinity_de_novo/Trinity.fasta
# <result folder> might be results/trinity_de_novo.transdecoder_dir

# Using the Trinity Transcriptome file, the program will identify the longest potenital open reading frames
TransDecoder.LongOrfs -t $1 -O $2
