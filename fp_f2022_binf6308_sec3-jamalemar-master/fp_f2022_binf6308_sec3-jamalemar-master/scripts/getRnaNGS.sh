#!/usr/bin/env bash
# getNGS.sh


# Retrieve the Clostridioides difficile NGS RNA-seq reads.
fasterq-dump --split-3 SRR21284802 -O data/rawreads/
fasterq-dump --split-3 SRR21284803 -O data/rawreads/
fasterq-dump --split-3 SRR21284804 -O data/rawreads/
fasterq-dump --split-3 SRR21284805 -O data/rawreads/
fasterq-dump --split-3 SRR21284806 -O data/rawreads/
fasterq-dump --split-3 SRR21284807 -O data/rawreads/
