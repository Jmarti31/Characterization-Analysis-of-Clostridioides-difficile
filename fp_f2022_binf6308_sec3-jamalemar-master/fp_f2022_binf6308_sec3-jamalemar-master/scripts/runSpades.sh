#!/usr/bin/env bash
# runSpades.sh

# Uses the trimmed paired reads to create a de-novo genome assembly
function Spades {
    spades.py \
    -1 data/trimmed/paired/Clostridioides.R1.paired.fastq \
    -2 data/trimmed/paired/Clostridioides.R2.paired.fastq \
    -o results/clostridioides
}

Spades 


