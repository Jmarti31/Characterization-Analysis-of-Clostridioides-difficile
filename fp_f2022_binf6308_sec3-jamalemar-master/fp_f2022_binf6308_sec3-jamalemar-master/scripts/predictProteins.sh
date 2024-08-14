#!/usr/bin/env bash
# predictProteins.sh
# Usage: bash scripts/predictProteins.sh <transcriptome file> <output location> <ORF retention criteria HMM> <ORF retention criteria blast> 1>results/logs/predictProteins.log 2>results/logs/predictProteins.err

# <transcriptome file> may be data/trinity_de_novo/Trinity.fasta
# <output location> may be results/trinity_de_novo.transdecoder_dir
# <ORF retention criteria HMM> may be results/pfam.domtblout
# <ORF retention criteria blast> may be results/blastPep.outfmt6

# Uses Previous hmmscan, blastp, and transdecoder.longorfs results to refine prediction of coding sequences
TransDecoder.Predict \
    -t $1 \
    -O $2 \
    --retain_pfam_hits $3 \
    --retain_blastp_hits $4
