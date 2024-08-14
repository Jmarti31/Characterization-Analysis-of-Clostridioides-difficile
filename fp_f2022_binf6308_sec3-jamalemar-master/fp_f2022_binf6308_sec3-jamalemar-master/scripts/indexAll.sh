#!/usr/bin/env bash
# indexAll.sh
# Usage: bash scripts/sortAll.sh 1>results/logs/indexAll.log 2>results/logs/indexAll.err &

# Initialize file location and name strings for BAM and BAI files
fastqPath="results/bam/"
bamPrefix=".sorted.bam"
baiPrefix=".sorted.bai"

# Index each available binary alignment file
function indexAll {
    for allBamFiles in $fastqPath*$bamPrefix
    do
        pathRemoved="${allBamFiles/$fastqPath/}"
        sampleName="${pathRemoved/$bamPrefix/}"
        samtools index \
        $fastqPath$sampleName$bamPrefix
        -o $fastqPath$sampleName$baiPrefix 

    done
}
indexAll
