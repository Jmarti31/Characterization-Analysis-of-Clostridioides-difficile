#!/usr/bin/env bash
# blastPep.sh
# Usage: bash scripts/blastPep.sh <longest orfs> <database> 1>results/blastPep.outfmt6 2>results/logs/blastPep.err

# <longest orfs> may be results/trinity_de_novo.transdecoder_dir/longest_orfs.pep
# <database> may be /work/courses/BINF6308/inputFiles/blastDB/swissprot

# Executes a protein blast alignment on a set of translated open reading frames
blastp -query $1 -db $2 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 
