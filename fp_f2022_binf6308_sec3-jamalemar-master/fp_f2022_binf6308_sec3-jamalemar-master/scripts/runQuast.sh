#!/usr/bin/env bash
# runQuast.sh

# Quality assessment tool to grade SPAdes genome assembly, provides contig statistics
function Quast {
    quast.py results/clostridioides/contigs.fasta \
    -o results/quast_results
}

Quast
