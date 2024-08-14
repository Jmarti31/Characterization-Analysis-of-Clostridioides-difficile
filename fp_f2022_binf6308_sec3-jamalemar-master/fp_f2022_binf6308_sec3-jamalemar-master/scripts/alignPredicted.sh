#!/usr/bin/env bash
# alignPredicted.sh
# Usage: bash scripts/alignPredicted.sh <predicted protein list> <target database> 1>results/alignPredicted.txt 2>results/logs/alignPredicted.err

# <predicted protein list> may be results/predictedProteins/Trinity.fasta.transdecoder.pep
# <target database> may be /work/courses/BINF6308/inputFiles/blastDB/swissprot


blastp -query $1 -db $2 -num_alignments 1 -outfmt "6 qseqid sacc qlen slen length nident pident evalue stitle" -evalue 1e-10 -num_threads 4 
