# Forming Bioinformatics Pipelines: Characterization Analysis of Clostridioides difficile Illumina Reads through De-Novo Omics Assemblies and Protein Prediction.


## PLEASE LOOK AT REPRODUCING_BIOINFORMATICS_PIPELINE FOR ANALYSIS PAPER


## Purpose

Next Generation Sequencing has revolutionized the way we process the world of biology.
Using sequencing techonologies such as Illumina and PacBio platforms, researchers the
world over now have the ability to answer questions spanning genomics to proteomics,
all through the use of computational methodology. In-silico methods have greatly 
expanded research capabilities by processing years of data into useful forms for use in
areas such gene annotaton, protein functionality, assembly reconstruction, etc. The
reason for creating a bionformatics pipeline is to provide users, a streamline
set of steps towards analyzing sets of data to generate concise, accurate, and useful
information from them. The current pipeline focuses on analyzing bacterial species,
and generates Genome and Transcriptome Assemblies from DNA & RNA-Seq data for
protein prediction and functional annotation of a given organism. The species being
used to test the following pipeline is known as Clostridioides difficile, and
is known primarily for its association with gastrointestinal diseases. The pipeline
will provide further characterization of the species, and can potentially be used
to analyze future infectious strains that are on the rise. By using the four sequential
sbatch scripts (1) sbatch_assemble_genome.sh , (2) sbatch_alignRNAseq.sh, 
(3) sbatch_trinity.sh, and (4) sbatch_transdecoder.sh, the desired pipeline
functionality can be accomplished 

## Data Dictionary

- DNA-Seq Reads: SRR22386611; https://www.ncbi.nlm.nih.gov/sra/?term=SRR22386611

- RNA-Seq Reads: SRR21284802, SRR21284803, SRR21284804, SRR21284805, SRR21284806, SRR21284807; https://www.ncbi.nlm.nih.gov/bioproject/PRJNA874447

- Genome File: GCF_000009205.2_ASM920v2_genomic.fna; https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000009205.2/

## Program & Tools

- SRA Toolkit 3.0.2 
- Trimmomatic 0.39-2
- SPAdes 3.15.4
- QUAST 5.2.0
- Bandage v0.9.0
- GMAP & GSNAP Version 2012-12-17
- SAMtools 3.0.2
- Trinity v2.15.0
- HMMER v3.3.2
- BLAST 2.9.0
- TransDecoder 5.6.0
- R v4.2.2

## End Result

The current iteration of the pipeline only supports Genome Assembly,
Transcriptome Assembly, and Protein Prediction. Functionality Identification
will be added into the future versions. Each part of the pipeline
ends in the completed compiled file for the next, such that users should
expect a genome gfa file, a Trinity transcriptome fasta file, and a 
Transdecoder.Predict based alignedProteins text file.

# Genome Assembly

## sbatch_assembleGenome.sh

## Overview

The program is meant to obtain the reads from a sequencing run
of the genome Clostridioides difficile NB647, and attempt to create and
assess a genome assembly. The four shell scripts inside of the sbatch script: 
Retrieve the sequence reads, trim the reads to remove leftover adapter 
sequences and determine poor quality areas to only keep good quality pairs, 
assemble the pairs into potential alignments of contigs aligned in this case 
by de Bruijin Graphs, and finally provide quality check statistics on the 
assembly created.

## Running sbatch_assembleGenome.sh

#### Log into Discovery from the terminal
* $ ssh [username]@login.discovery.neu.edu

### First cd to the area of execution in the scratch file 
* $ cd /scratch/[user]/fp_f2022_binf6308_sec3-[username]/

### Next call the sbatch_assembleGenome.sh script
* $ sbatch scripts/sbatch_assembleGenome.sh


## Methods

### getNGS.sh
#### Retrieves the sequence reads from an experiment
* Usage: fasterq-dump [options] SRR[number]
* Uses: fasterq-dump --split-3 SRR22386611

- fasterq-dump is the program that will be used to retrieve an experiments run.
- --split-3 splits the reads in to 1st mate reads, 2nd mate reads, and unmated reads.
- SRR22386611 represents the experiment ID that contains the reads.


### trim.sh
#### Seperates quality paired end reads from excess adapters and low quality reads
* Uses: trimmomatic PE \
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

- trimmomatic PE tells the program that we are using paired end reads.
- -threads 1 tells how many threads.
- -phred33 identifies the quality encoding method, seen in 4th line of fastq file.
- SRR22386611_1.fastq and SRR22386611_2.fastq are the seperate paired end files to be used.
- Clostridioides.R1.paired.fastq, Clostridioides.R1.unpaired.fastq, Clostridioides.R2.paired.fastq, Clostridioides.R2.unpaired.fastq creates the output files seperating quality reads from removable ones.
- HEADCROP:0 tells the program to remove 0 bases from each read.
- ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 is used to determine adapter mismatches.
- LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36 finally are used to indicate qualities and parameters necessary to keep a read, such as bases trailing and leading requiring a score of 20, the size of a window and its neccessary quality score of 36, and the minimum length of a read.


### runSpades.sh
#### Creates an Assembly From Pair End Reads
* Usage: spades.py [options] -o [output_dir]
* Uses: spades.py \
    -1 data/trimmed/paired/Clostridioides.R1.paired.fastq \
    -2 data/trimmed/paired/Clostridioides.R2.paired.fastq \
    -o results/rhodo

- spades.py is the actual file that aligns the reads to a potential assembly.
- -1 data/trimmed/paired/Clostridioides.R1.paired.fastq uses the first read mate of a pair.
- -2 data/trimmed/paired/Clostridioides.R2.paired.fastq uses the second read mate of a pair.
- -o results/clostridioides directs the output toward the results folder.


### runQuast.sh
#### Quality Assessment tool for Assembly

* Usage: python quast.py [options] [contig_file(s)]
* Uses: quast.py clostridioides/contigs.fasta \
    -o results/quast_results

- quast.py is the actual file that works on a given assembly.
- clostridioides/contigs.fasta is the contig assembly generated from SPAdes.
- -o results/quast_results directs the QUAST output to the results folder.


## Expected Output 

Number of contigs, scaffold size, and N50 for our case are the main metrics we will use to determine a good assembly.
There are 43 contigs, with the largest contig being 549,417 bp long, with the assembly having a length of 4,008,790 bp.
The N50 is of 208886 bp and GC content is 28.37%. No mismatches are present. Basing off the data and graph provided,
the genome assembly is great.


# Transcriptome Assembly


## sbatch_alignRNAseq.sh

## Overview

The program is meant to obtain the cDNA reads from a RNA-Seq experiment
on wild type and ddl deficient Clostridioides difficile 630. The sequencing 
was performed on a two by three, treatment/replicate group of 6 Clostridioides 
difficile samples and will be used to create a Reference-Guided Transcriptome 
Assembly. The created sbatch script loads six sequential .sh scripts: (1) getRnaNGS.sh
which retrieves the RNA-Seq Reads (2) ClosBuild.sh which creates the reference genome database to align the reads 
(3) trimAll.sh which trims the reads (4) alignAll.sh uses reads and the reference 
genome database to create sam alignment files. (5) sortAll.sh uses SAMtools to 
convert sam files to bam files. (6) indexAll.sh is finally used to generate a bam 
file index.


## Running sbatch_alignRNAseq.sh

#### Log into Discovery from the terminal
* $ ssh [username]@login.discovery.neu.edu

### First cd to the module in the scratch file(certain directories like data, scripts, results pre-exist)
* $ cd /scratch/[user]/fp_f2022_binf6308_sec3-[username]/

### Call sbatch to generate alignment files
* $ sbatch scripts/sbatch_alignRNAseq.sh


## Methods

### getRnaNGS.sh
#### Retrieves the sequence reads from an experiment
* Usage: fasterq-dump [options] SRR[number]
* Uses: fasterq-dump --split-3 SRR21284802 -O data/rawreads/
	fasterq-dump --split-3 SRR21284803 -O data/rawreads/
	fasterq-dump --split-3 SRR21284804 -O data/rawreads/
	fasterq-dump --split-3 SRR21284805 -O data/rawreads/
	fasterq-dump --split-3 SRR21284806 -O data/rawreads/
	fasterq-dump --split-3 SRR21284807 -O data/rawreads/

- fasterq-dump is the program that will be used to retrieve an experiments run.
- --split-3 splits the reads in to 1st mate reads, 2nd mate reads, and unmated reads.
- SRR2128480[2-7] represents the experiment ID that contains the reads.

Reads were quality trimmed using Trimmomatic(Bolger, 2014):
### trimRnaAll.sh
#### Seperates quality paired end reads from excess adapters and low quality reads. In the case of this project, a function call was utilized to trim retrieved raw reads from the RNA_seq run for the 24 samples. Within the file we utilized bash script notation to create variables to access each raw read by common file name.

* Usage:
function trimAll {
    for leftInFile in $fastqPath*$leftSuffix
    do
        pathRemoved="${leftInFile/$fastqPath/}"
        sampleName="${pathRemoved/$leftSuffix/}"
        echo $sampleName
        trimmomatic PE \
        -threads 1 -phred33 \
        $fastqPath$sampleName$leftSuffix \
        $fastqPath$sampleName$rightSuffix \
        $pairedOutPath$sampleName$leftSuffix \
        $unpairedOutPath$sampleName$leftSuffix \
        $pairedOutPath$sampleName$rightSuffix \
        $unpairedOutPath$sampleName$rightSuffix \
        HEADCROP:0 \
        ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36
    done
}
trimAll

* flags used:
- trimmomatic PE tells the program that we are using paired end reads.
- -threads 1 tells how many threads.
- -phred33 identifies the quality encoding method, seen in 4th line of fastq file.
- $fastqPath$sampleName$leftSuffix and $fastqPath$sampleName$rightSuffix are the seperate location for the paired end files to be used.
-  $pairedOutPath$sampleName$leftSuffix, $unpairedOutPath$sampleName$leftSuffix, $pairedOutPath$sampleName$rightSuffix, $unpairedOutPath$sampleName$rightSuffix, are the path to the output files used for seperating quality reads from removable ones.
- HEADCROP:0 tells the program to remove 0 bases from each read.
- ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 is used to determine adapter mismatches.
- LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36 finally are used to indicate qualities and parameters necessary to keep a read, such as bases trailing and leading requiring a score of 20, the size of a window being 4 and its neccessary quality score of 36, and the minimum length of a read being 36.; 

Reads were aligned with GSNAP using a GMAP database(Wu, 2016):
### ClosBuild.sh & alignAll.sh
#### GMAP creates an index of the reference genome database in order to provide a basis to align given reads. GSNAP was used with a function call in order to use the database to align paired end reads from the 24 studied specimens and return individual SAM alignment files.

* Usage:
gmap_build -D data \
-d ClostridioidesGmapDb \
/home/martinez.jam/GCF_000009205.2_ASM920v2_genomic.fna
    
function alignAll {
    for leftInFile in $pairedOutPath*$leftSuffix
    do
        pathRemoved="${leftInFile/$pairedOutPath//}"
        sampleName="${pathRemoved/$leftSuffix/}"
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

* flags used:
- gmap_build builds an index of a given reference genome
- -D data sets where the database will be compiled
- -d ClostridioidesGmapDb sets the name for the database
- /home/martinez.jam/GCF_000009205.2_ASM920v2_genomic.fna sets a location for the Clostridioides genome fasta file

- gsnap runs alignment of kept paired end reads to reference genome
- -A sam specifies format of the file produced
- -D data sets where the database will be compiled
- -d ClostridioidesGmapDb sets the name for the database
- -N 1 used to specify identification of novel splice sites
- $pairedOutPath$sampleName$leftSuffix & $pairedOutPath$sampleName$rightSuffix are the file paths to the paired end reads of a sample
- 1>$samDirectory$sampleName$samPrefix redirected STDOUT output into sam file
;

Alignments were sorted and indexed with SAMtools(Li, 2009):
### sortAll.sh & indexAll.sh
#### Using SAMtools in two functions, the "sort" command is used to convert sam files into bam files, while the "index" function is used to generate an index of the bam files. We utilize for loops functions to iterate over all sam files correlated with different samples. 

* Usage:
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

* flags used:
- samtools sort specifies "sort" command to convert sam file to bam file
- $fastqPath$sampleName$samPrefix file path toward a given sam file
- -o $bamDirectory$sampleName$bamPrefix output file location to write sorted bam information too

- samtools index specifies "index" command to generate an index of a bam file
- $fastqPath$sampleName$bamPrefix file path toward a given bam file
- -o $fastqPath$sampleName$baiPrefix output file location to write indexed bam information tool

## Expected Output 

The user should have created a directory of raw reads that was then trimmed and 
seperated into fourteen paired and unpaired read files. Following Read Alignment, the user
should have obtained sam files that were converted to binary alignment files,
and then sorted. The resulting bam files would be used in the next step.


## sbatch_trinity.sh

## Overview

The sbatch script executes Trinity software to assemble an RNA-Seq based Transcriptome for Clostridioides difficile
using both reference guided assembly and de-novo assembly methodology. The first reference
genome guided assembly utilizes three scripts: (1) mergeAll.sh for obtaining a single sorted bam file, 
(2) runTrinity.sh for transciptome assembly via the singular compacted bam file, and (3) analyzeTrinity.sh 
for viewing the assembly output statistics. The second transcriptome assembly using de-novo methods uses
two scripts: (1) trinityDeNovo.sh which uses left and right paired reads to generate a potential transcriptome
assembly and (2) analyzeTrinityDeNovo.sh for viewing of the assembly output statistics. The process
was meant to demonstrate the difference of output when using different methodologies on the same data.
Future refinement would see the reference-guided assembly being useful for SNV analysis of future strains.  

#### Log into Discovery from the terminal
* $ ssh [username]@login.discovery.neu.edu

### First cd to the module in the scratch file(certain directories like data, scripts, results pre-exist)
* $ cd /scratch/[user]/fp_f2022_binf6308_sec3-[username]/

### Call sbatch to generate alignment files
* $ sbatch scripts/sbatch_trinity.sh


## Methods

Sorted bam files were merged into a single file by SAMtools(Li, 2009)
### mergeAll.sh
#### Applying SAMtools, the program is able to take a file with all file names
associated to previously sorted bam files and convert those files into a single bam file.

* Usage:

samtools merge -b data/bam/bamIn.txt data/bam/ClosAll.bam

* flags used:
- samtools merge specifies "merge" command to convert bam files into a bam file
- -b data/bam/bamIn.txt provides the file location to retrieve the bam file names to be merged
- data/bam/ClosAll.bam is the output file location for the merged bam file

Transcriptome assembly through reference guided and de novo methods done through Trinity(Haas, 2013)
### runTrinity.sh and TrinityDeNovo.sh
#### Utilizing Trinity, the reference guided assembly uses a singular compiled bam file with parameters
set up to identify max intron size and resources for parallelizaton for the creation of the transcriptome
assembly. The De Novo method utilizes paired RNA reads to generate an transcriptome assembly from scratch.

* Usage for reference guided:
Trinity \
--genome_guided_bam data/bam/ClosAll.bam \
--genome_guided_max_intron 10000 \
--max_memory 10G --CPU 4 \
--output results/trinity_guided

* flags used:
- --genome_guided_bam data/bam/ClosAll.bam Provides sorted bam file to be used for assembly reference
- --genome_guided_max_intron 10000 sets the size limit for introns within the assembly, not necessary in hindsight due to bacteria
- --max_memory 10G --CPU 4 sets limit of memory to be used and the CPU's to be used(for parallelization) 
- --output results/trinity_guided sets the output directory for the final transcriptome assembly


* Usage for de novo:
Trinity \
--seqType fq \
--output results/trinity_de_novo \
--max_memory 10G --CPU 4 \
--left $leftReads \
--right $rightReads

* flags used:
- --seqType fq sets the sequence type of the reads to be in FASTQ format
- --max_memory 10G --CPU 4 sets limit of memory to be used and the CPU's to be used(for parallelization)
- --output results/trinity_guided sets the output directory for the final transcriptome assembly
- --left $leftReads --right $rightReads sets file names for left and right reads.

Transcriptome assembly stats were produced by Trinity
### analyzeTrinityDeNovo.sh and analyzeTrinity.sh
#### Comparison of the transcriptome assembly statistics produced using both assembly methods.


* Usage:
TrinityStats.pl results/trinity_guided/Trinity-GG.fasta
TrinityStats.pl results/trinity_de_novo/Trinity.fasta 

* flags used:
- TrinityStats.pl Trinity program for returning assembly stats results
- results/trinity_de_novo/Trinity.fasta and results/trinity_guided/Trinity-GG.fasta are output file locations

## Expected output

The type of assembly method used can create different transcriptome profiles from the same data. In the case
of the reference guide method the most important statistics produced were the N50 from all transcript
contigs being 1685 and the N50 from the longest isoform per gene being 1690. Meawhile in the de novo method,
the N50 from all transcript contigs was 2753 and the N50 from the longest isoform per gene was 2350. Without,
extra assessment tools, in face value the de novo assessment seems to have done a slightly better job, since
the N50 is slightly higher. However by viewing the total genes and transcripts created, the lower amount 
produced in the reference guided method supports that it may be the more precise assembly between the two. The 
Trinity.fasta and Trinity-GG.fasta files are produced at the end of the script, and both should be viable 
for the proceeding protein prediction step. The Trinity.fasta file that will be used in the next step should
have 2373 transcripts.


# Protein Prediction/Gene Annotation


## sbatch_transdecoder.sh

## Overview

The sbatch script is meant to Trinity associated software, TransDecoder, to identify open reading frames
for protein prediction using a RNA-Seq based Transcriptome Assembly as a template. In the process to create
a refined pipeline to identify significant similar protein, five scripts were executed: (1) longOrfs.sh 
for identifying long ORF's in the transcriptome assembly. (2) blastPep.sh
for using BLAST to find similar proteins. (3) pfamScan.sh for identifying protein domains. 
(4) predictProteins.sh uses ORFS, BLAST output, and discovered protein domain to refine predicted 
protein identification. (5) alignPredicted.sh BLAST's the final predicted proteins to find
similar matches and formats the output in tabular form.

## Running sbatch_trinity.sh

#### Log into Discovery from the terminal
* $ ssh [username]@login.discovery.neu.edu

### First cd to the module in the scratch file(certain directories like data, scripts, results pre-exist)
* $ cd /scratch/[user]/fp_f2022_binf6308_sec3-[username]/

### Call sbatch to generate alignment files
* $ sbatch scripts/sbatch_transdecoder.sh

## Methods

#### Takes the longest open reading fram and translates to protein. (Haas 2013)
### longOrfs.sh
#### Utilizng the TransDecoder.LongOrfs program with the previous weeks transcriptome assembly,
the longest orfs are identified and converted into amino acid sequences to determine
potential proteins in the sample's proteome.

* Usage:

TransDecoder.LongOrfs -t $1 -O $2

* flags used:
- TransDecoder.LongOrfs a Trinity associated program in which identifies long ORF's in an assembly
- -t $1 the user-defined transcriptome assembly Trinity fasta file to be used for identification
- -O $2 the user-defined output location for the result directory

#### Aligns predicted proteins to potentially similar database proteins. (Pertsemlidis 2001) (Altschul 1997)
### blastPep.sh
#### Using the pep file containing predict proteins created from the TransDecoder.LongOrfs program,
the bash script executes a protein blast (blastp) that will attempt to find significant protein matches
to the query.

* Usage:

blastp -query $1 -db $2 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 

* flags used:
- blastp BLAST program using protein query against a protein database
- -query $1 the user-defined pep file containing predicted proteins for BLASTing
- -db $2 the user-defined database to be searched against when BLASTing
- -max_target_seqs 1 returns specified number of significant match results
- -outfmt 6 tab-delimited output type of end BLAST results
- -evalue 1e-5 sets the evalue threshold for filtering out insignificant matches
- -num_threads 4 specifies cores to be used 

#### Idenftifies protein domains using Hidden Markov Models. (Finn 2015)
### pfamScan.sh
#### Using the pfam database, hidden markov models are utilized to identify potential
protein domains within each of the predicted proteins in the originial pep file.

* Usage:

hmmscan --cpu 4 --domtblout $1 $2 $3

* flags used:
- hmmscan Hidden Markov model scanning program for identifying protein domains
- --cpu 4 sets the amount of memory to be allocated to search through the predicted proteins
- --domtblout $1 the user-defined tab delimited save file to track results
- $2 arguement to specify target database used for hmmscan
- $3 arguement to specify the pep file used for hmmscan

#### Takes ORFS, BLAST output, and discovered protein domain to refine predicted protein identification. (Haas 2013)
### predictProteins.sh
#### Trinity associated program that predicts potential proteins from a transcriptome assembly
using previous hmmscan and Blast results as retention criteria to refine the search.

* Usage:

TransDecoder.Predict \
    -t $1 \
    -O $2 \
    --retain_pfam_hits $3 \
    --retain_blastp_hits $4

* flags used:
- TransDecoder.Predict the program used to predict potential proteins within a transcriptome assembly
- -t $1 the user-defined transcriptome assembly to be searched
- -o $2 the user-defined output location
- --retain_pfam_hits $3 the user-defined list of results of found protein domains to be used for significant ORF retention
- --retain_blastp_hits $4 the user-defined list of BLAST results used for significant ORF retention

### alignPredicted.sh
#### aligns the final output of predicted proteins and significant BLAST results.

* Usage: 

blastp -query $1 -db $2 -num_alignments 5 -outfmt "6 qseqid sacc qlen slen length nident pident evalue stitle" -evalue 1e-10 -num_threads 4 

* flags used:
- blastp BLAST program using protein query against a protein database
- -query $1 user-defined predicted protein list pep file to be BLASTed
- -db $2 user-defined database to be used against BLAST search
- -num_alignments 1 specifies the amount of results to return 5 hits
- -outfmt "6 qseid sacc qlen slen length nident pident evalue stitle" tab delimited file output for end BLAST results
- -evalue 1e-10 sets the evalue threshold for filtering out insignificant matches
- -num_threads 4 specifies cores to be used

## Expected Output  
 
Using the longest orf results pep file from Transdecoder.LongOrfs, 
inside th blastp results:(2413 similar proteins, indicative that some orfs didn't match with any proteins),
and the HMMscan(24002 possible domains), the Transdecoder.Predict program was able to predict 3087 potential 
coding sequences. A final blastp search on the predicted CDS's generated an end text file with 2138 significant 
matches, meaning that some of the previous transcripts didn't have coding areas associated to them or weren't 
identified correctly, or possibly that there are currently no known similar proteins.



# Foreword

The completion of the Genome Assembly, Transcriptome Assembly, and Protein Prediction aspects of the
pipeline generated critical information into characterizing Clostridioides difficile. The hope is to 
extend the pipeline into functional annotations of the predicted proteins identified, in order to
further characterize aspects of the sample. A template of the Functional Annotation part of the pipeline
was left in the root directory for potential use by the next project collaborator's.

# References

- Alhakami, H., Mirebrahim, H., & Lonardi, S. (2017). A comparative evaluation of genome assembly reconciliation tools. Genome biology, 18(1), 93. https://doi.org/10.1186/s13059-017-1213-3
- Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology : a journal of computational molecular cell biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021
- Belitsky, B. R. (2022). VanG- and D-Ala-D-Ser-dependent peptidoglycan synthesis and vancomycin resistance in Clostridioides difficile. Molecular Microbiology, 118, 526– 540. https://doi.org/10.1111/mmi.14980
- Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
- Finn, R. D., Clements, J., & Eddy, S. R. (2011). HMMER web server: interactive sequence similarity searching. Nucleic acids research, 39(Web Server issue), W29–W37. https://doi.org/10.1093/nar/gkr367
- Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. Bioinformatics (Oxford, England), 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086
- Haas, B. J., Papanicolaou, A., Yassour, M., Grabherr, M., Blood, P. D., Bowden, J., Couger, M. B., Eccles, D., Li, B., Lieber, M., MacManes, M. D., Ott, M., Orvis, J., Pochet, N., Strozzi, F., Weeks, N., Westerman, R., William, T., Dewey, C. N., … Regev, A. (2013). De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nature Protocols, 8(8), 1494–1512. https://doi.org/10.1038/nprot.2013.084
- Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
- McGinnis, S., & Madden, T. L. (2004). BLAST: at the core of a powerful and diverse set of sequence analysis tools. Nucleic acids research, 32(Web Server issue), W20–W25. https://doi.org/10.1093/nar/gkh435
- Wick, R. R., Schultz, M. B., Zobel, J., & Holt, K. E. (2015). Bandage: interactive visualization of de novo genome assemblies. Bioinformatics (Oxford, England), 31(20), 3350–3352. https://doi.org/10.1093/bioinformatics/btv383
- Wu, T. D., Reeder, J., Lawrence, M., Becker, G., & Brauer, M. J. (2016). GMAP and GSNAP for Genomic Sequence Alignment: Enhancements to Speed, Accuracy, and Functionality. Statistical Genomics, 283–334. https://doi.org/10.1007/978-1-4939-3578-9_15

