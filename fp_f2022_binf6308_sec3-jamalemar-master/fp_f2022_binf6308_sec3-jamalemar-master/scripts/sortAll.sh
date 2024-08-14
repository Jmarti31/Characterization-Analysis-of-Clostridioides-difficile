#!/usr/bin/env bash
# sortAll.sh
# Usage: bash scripts/sortAll.sh 1>results/logs/sortAll.log 2>results/logs/sortAll.err &

# Initialize file locations and name strings for SAM and BAM files
fastqPath="results/sam/"
samPrefix=".sam"
bamPrefix=".sorted.bam"
bamDirectory="results/bam/"

# Creates the directory to hold all created BAM files
mkdir -p $bamDirectory

function sortAll {
    for allSamFiles in $fastqPath*$samPrefix
    do
        pathRemoved="${allSamFiles/$fastqPath/}"
        sampleName="${pathRemoved/$samPrefix/}"
        samtools sort \
        $fastqPath$sampleName$samPrefix \
	-o $bamDirectory$sampleName$bamPrefix 

    done
}
sortAll
