#!/usr/bin/env bash
# trim.sh

# Provides quality trimming to downloaded DNA-seq reads
# It is noted that defaults are left due to missing exact adapter sequence used
# Along with that, MINLEN possibly should have been reduced to reflect reads after trimming, possibly 16
# End result without change was good though, so default worked for our purposes

PATH_TO_TRIMMOMATIC="/shared/centos7/anaconda3/2021.11/envs/BINF-12-2021/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2"
function trim {
    trimmomatic PE \
    -threads 1 -phred33 \
    data/SRR22386611_1.fastq \
    data/SRR22386611_2.fastq \
    data/trimmed/paired/Clostridioides.R1.paired.fastq \
    data/trimmed/unpaired/Clostridioides.R1.unpaired.fastq \
    data/trimmed/paired/Clostridioides.R2.paired.fastq \
    data/trimmed/unpaired/Clostridioides.R2.unpaired.fastq \
    HEADCROP:0 \
    ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36
}
trim
