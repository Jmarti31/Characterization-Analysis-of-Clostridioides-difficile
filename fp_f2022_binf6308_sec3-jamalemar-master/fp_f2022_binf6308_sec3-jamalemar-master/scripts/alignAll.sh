#!/usr/bin/env bash
# alignAll.sh
# Usage: bash scripts/alignAll.sh 1>results/logs/alignAll.log 2>results/logs/alignAll.err &

# Initialize variables to contain name strings and file locations for the left reads
leftSuffix="_1.fastq"
rightSuffix="_2.fastq"
pairedOutPath="data/RNAtrimmed/paired/"
samPrefix=".sam"
resultPath="results/"
samDirectory="results/sam/"

# Creates the directory to hold all SAM files
mkdir -p $samDirectory

function alignAll {
    for leftInFile in $pairedOutPath*$leftSuffix
    do
        # Remove the path from the filename and assign to pathRemoved
        pathRemoved="${leftInFile/$pairedOutPath/}"
        # Remove the left-read suffix from $pathRemoved and assign to suffixRemoved
        sampleName="${pathRemoved/$leftSuffix/}"
        # Print $sampleName to see what it contains after removing the path
        echo $sampleName
        gsnap \
        -A sam \
        -D data \
        -d ClostridioidesGmapDb \
        -N 1 \
        $pairedOutPath$sampleName$leftSuffix \
        $pairedOutPath$sampleName$rightSuffix \
        1>$samDirectory$sampleName$samPrefix
    done
}
alignAll
