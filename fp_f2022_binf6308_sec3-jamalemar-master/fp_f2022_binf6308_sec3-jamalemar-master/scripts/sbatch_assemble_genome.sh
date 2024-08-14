#!/usr/bin/bash
#SBATCH --partition=short               # choose from debug, express, or short
#SBATCH --job-name=assemblegenome
#SBATCH --time=02:00:00                 # the code pieces should run in far less than 1 hour
#SBATCH -N 1                            # nodes requested
#SBATCH -n 1                            # task per node requested
#SBATCH --output="batch-%x-%j.output"   # where to direct standard output; will be batch-jobname-jobID.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=martinez.jam@northeastern.edu # Update to your user name!

# Usage: sbatch sbatch_alignRNAseq.sh 
# Purpose: Using DNA-seq reads create a de-novo genome assembly

echo "Starting our analysis $(date)"  

ORGANISM="Clostridioides difficile"
SRR_ID="SRR22386611"
export SRR_ID="SRR22386611"

echo "$ORGANISM SRR reads to process: $SRR_ID"

mkdir -p data/
mkdir -p results/
mkdir -p results/logs/

echo "Loading our BINF6308 Anaconda environment."
module load anaconda3/2021.11
source activate BINF-12-2021

echo "Downloading $SRR_ID reads $(date)"
bash scripts/getNGS.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-getNGS.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-getNGS.err

mkdir -p data/trimmed/
mkdir -p data/trimmed/paired/
mkdir -p data/trimmed/unpaired/

echo "Trimming $SRR_ID reads $(date)"
bash scripts/trim.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trim.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-trim.err

mkdir -p results/clostridioides/

echo "Assembling genome from trimmed $SRR_ID reads $(date)"
bash scripts/runSpades.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runSpades.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runSpades.err

mkdir -p results/quast_results/

echo "Analyzing genome assembly $(date)"
bash scripts/runQuast.sh 1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runQuast.log 2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-runQuast.err

echo "Assembly and analysis complete $(date)"
