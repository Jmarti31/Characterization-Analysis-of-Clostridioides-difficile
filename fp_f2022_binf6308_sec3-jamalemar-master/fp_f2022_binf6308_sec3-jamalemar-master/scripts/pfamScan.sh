#!/usr/bin/env bash
# pfamScan.sh
# Usage: bash scripts/pfamScan.sh <tab delimited save file> <target database> <longest orfs> 1>results/logs/pfamScan.log 2>results/logs/pfamScan.err

# <tab delimited save file> may be results/pfam.domtblout
# <target database> may be /work/courses/BINF6308/inputFiles/SampleDataFiles/Pfam-A.hmm
# <longest orfs> may be results/trinity_de_novo.transdecoder_dir/longest_orfs.pep

# Uses Hidden Markov Models to compare sequence to known domains profiles for similarity
hmmscan --cpu 4 --domtblout $1 $2 $3
