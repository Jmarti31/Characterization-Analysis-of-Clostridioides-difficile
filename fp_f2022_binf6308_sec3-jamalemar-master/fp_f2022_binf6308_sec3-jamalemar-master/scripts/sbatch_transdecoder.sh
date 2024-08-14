#!/usr/bin/bash
#SBATCH --partition=short               # choose from debug, express, or short
#SBATCH --job-name=transdecoder
#SBATCH --time=08:00:00                 
#SBATCH -N 1                            # nodes requested
#SBATCH -n 4                            # task per node requested
#SBATCH --mem=10Gb
#SBATCH --exclusive
#SBATCH --output="batch-%x-%j.output"   # where to direct standard output; will be batch-jobname-jobID.output
#SBATCH --mail-type=ALL
#SBATCH --mail-user="martinez.jam@northeastern.edu" # Update to your user name!

# Usage: sbatch sbatch_transdecoder.sh
# Purpose: Using Transcripts obtained from Trinity.fasta Assembly, Generate a Set of Predicted Proteins

# Assumes follows sequential execution from scratch space, so input data is in /scratch/martinez.jam/fp_f2022_binf6308_sec3-jamalemar/data/ except for de_novo analysis

echo "Starting our analysis $(date)" 
echo 

USER="martinez.jam"

# define key constants
TRANSCRIPTOME=data/trinity_de_novo/Trinity.fasta
SWISSPROT_DB=/work/courses/BINF6308/inputFiles/blastDB/swissprot 
TRANSDECODER_DIR=results/trinity_de_novo.transdecoder_dir
LONGEST_ORFS=$TRANSDECODER_DIR/longest_orfs.pep
OUTFMT=results/blastPep.outfmt6
DOMTBLOUT=results/pfam.domtblout
PFAMA_PATH=/work/courses/BINF6308/inputFiles/SampleDataFiles/Pfam-A.hmm
PREDICTED_PROTEIN_PATH=results/predictedProteins
FINAL_PROTEINS=$PREDICTED_PROTEIN_PATH/*transdecoder.pep

# record these key constants to our batch*.output file by echoing them:
echo "Key parameters"
echo "TRANSCRIPTOME: $TRANSCRIPTOME"
echo "SWISSPROT_DB: $SWISSPROT_DB"
echo "TRANSDECODER_DIR: $TRANSDECODER_DIR"
echo "LONGEST_ORFS: $LONGEST_ORFS"
echo "OUTFMT: $OUTFMT"
echo "DOMTBLOUT: $DOMTBLOUT"
echo "PFAMA_PATH: $PFAMA_PATH"
echo "PREDICTED_PROTEIN_PATH: $PREDICTED_PROTEIN_PATH"
echo "FINAL_PROTEINS: $FINAL_PROTEINS"
echo "USER: $USER"
echo
echo


echo "Loading our BINF6308 Anaconda environment."
module load anaconda3/2021.11
source activate BINF-12-2021

echo "Make directory for data files"
mkdir -p data/

echo "Make directory for results files"
mkdir -p results/


# If done sequentially, the data should have been created in last sbatch script, but else just cp from the home directory.
echo "Look at sbatch documentation if desired to move de novo Trinity transcriptome data from the home directory"
# cp -r /home/$USER/BINF6308/fp_f2022_binf6308_sec3-jamalemar/data/trinity_de_novo data/trinity_de_novo


echo "Move de novo directory from results to data"
mv results/trinity_de_novo/ data/


echo "Make directory for log files"
mkdir -p results/logs/

echo "Starting ORF prediction pipeline $(date)"
echo "Identify longORFs with TransDecoder.LongOrfs on $TRANSCRIPTOME $(date)"
bash scripts/longOrfs.sh $TRANSCRIPTOME $TRANSDECODER_DIR \
  1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-longOrfs.log \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-longOrfs.err

echo "BLASTp of longest_orfs.pep against SwissProt BLAST DB at $SWISSPROT_DB $(date)"
bash scripts/blastPep.sh $LONGEST_ORFS $SWISSPROT_DB \
  1>$OUTFMT \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-blastPep.err

echo "Create pfamScan with hmmscan using the Pfam-A.hmm file found $PFAMA_PATH $(date)"
bash scripts/pfamScan.sh $DOMTBLOUT $PFAMA_PATH $LONGEST_ORFS \
  1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-pfamScan.log \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-pfamScan.err

echo "Predict proteins with TransDecoder.Predict $(date)"
bash scripts/predictProteins.sh $TRANSCRIPTOME $TRANSDECODER_DIR $DOMTBLOUT $OUTFMT \
  1>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-predictProteins.log \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-predictProteins.err
echo "Copy TransDecoder.Predict outputs to $PREDICTED_PROTEIN_PATH"
mkdir -p $PREDICTED_PROTEIN_PATH
mv *transdecoder* $PREDICTED_PROTEIN_PATH

echo "Align predicted proteins to SwissProt DB $(date)"
bash scripts/alignPredicted.sh $FINAL_PROTEINS $SWISSPROT_DB \
  1>results/alignPredicted.txt \
  2>results/logs/$SLURM_JOB_NAME-$SLURM_JOB_ID-alignPredicted.err

echo "ORF prediction pipeline complete $(date)"

echo "Moving key files back to /home"
cp -r results/ /home/$USER/BINF6308/fp_f2022_binf6308_sec3-jamalemar/data/proteinPrediction/

echo "Analysis complete $(date)

