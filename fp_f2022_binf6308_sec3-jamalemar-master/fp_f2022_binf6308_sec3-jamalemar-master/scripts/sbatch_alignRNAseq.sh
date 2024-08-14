#!/usr/bin/bash
#SBATCH --partition=short               # choose from debug, express, or short
#SBATCH --job-name=alignRNAseq
#SBATCH --time=03:00:00                 # the code pieces should run in far less than 1 hour
#SBATCH -N 1                            # nodes requested
#SBATCH -n 1                            # task per node requested
#SBATCH --output="batch-%x-%j.output"   # where to direct standard output; will be batch-jobname-jobID.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=martinez.jam@northeastern.edu # Update to your user name!

# Usage: sbatch sbatch_alignRNAseq.sh 
# Purpose: Generate Binary alignment Files from RNA-Seq Reads for Refrence Based Transcriptome Assembly

echo "Starting our analysis $(date)"  

echo "Loading our BINF6308 Anaconda environment."
module load anaconda3/2021.11
source activate BINF-12-2021
echo "Loading GSNAP and Samtools."
module load gsnap/2021-12-17
module load samtools/1.10

mkdir -p data/rawreads/
echo "Stores raw FASTQ reads for analysis.

echo "Make directory for log files"
mkdir -p results/logs/


echo "Retrieve RNA-Seq Reads $(date)"
bash scripts/getRnaNGS.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-getRnaNGS.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-getRnaNGS.err

echo "Build the reference genome $(date)"
bash scripts/ClosBuild.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-ClosBuild.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-ClosBuild.err

echo "Trim all reads in data/rawreads/ $(date)"
bash scripts/trimRnaAll.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trimRnaAll.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trimRnaAll.err

echo "Align the reads to the reference with GSNAP $(date)"
bash scripts/alignAll.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-alignAll.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-alignAll.err

echo "Sort the resulting SAM files $(date)"
bash scripts/sortAll.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-sortAll.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-sortAll.err

echo "Index the resulting BAM files $(date)"
bash scripts/indexAll.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-indexAll.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-indexAll.err

echo "Alignment complete $(date)"
