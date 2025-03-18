#!/bin/bash

######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
##  - Create a scratch directory for the dog RNAseq project
##  - Define variables for directories and user ID
##  - Download RNAseq data from NCBI SRA for dog liver samples using SRA Toolkit
##  - Use FastQC to evaluate the quality of the downloaded data
## Download from SRA:
##  - Input: SRA Run IDs (SRR/DRR) from the experimental design
##  - Output: FASTQ files (R1 and R2 for paired-end data)
## FastQC:
##  - Input: Downloaded FASTQ files
##  - Output: A folder with .html and .txt files for quality assessment
## For running the script on the Alabama Super Computer (ASC):
##  - Suggested parameters:
##    - Queue: bigmem (dog data may require more resources than Daphnia)
##    - Cores: 20 (faster downloads and FastQC)
##    - Time limit (HH:MM:SS): 04:00:00 (increased for larger files)
##    - Memory: 150gb (increased for safety)
##    - Run on asax
## After making the script executable with "chmod +x script_name.sh",
## submit it using "run_script script_name.sh"
###############################################

########## SLURM Directives
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=download_fastqc_%j.o

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1
########## Define Variables and Make Directories
## Replace [Your_ASC_ID] with your actual ASC ID
MyID=aubclsd0338  ## Example: MyID=aubtss

## Define directories for your project
DD=/scratch/${MyID}/DogRNAseq/RawData  ## Data Directory for raw FASTQ files
WD=/scratch/${MyID}/DogRNAseq         ## Working Directory for the project
RDQ=RawDataQuality                    ## Directory for FastQC results

## Make directories in SCRATCH for holding the raw data
## -p ensures parent directories are created if they don’t exist
mkdir -p ${DD}
## Move to the Data Directory
cd ${DD}

########## Download Data Files from NCBI SRA Using Run IDs
## Using SRA Toolkit to download and convert SRA files to FASTQ
## -F: Defline contains only original sequence name
## --split-files: Splits paired-end data into R1 and R2 files
## Note: All samples in your table are paired-end (confirmed via SRA metadata, e.g., AvgSpotLen 202 for SRR8996966)

## Samples from your experimental design (large and small dog breeds, liver tissue, male, healthy)
## Large breeds
fastq-dump -F --split-files SRR8996966  ## Newfoundland, PRJNA396033
fastq-dump -F --split-files SRR8997009  ## Belgian Malinois, PRJNA396033
fastq-dump -F --split-files SRR5889337  ##Labrador, PRJNA396033
fastq-dump -F --split-files SRR2960319  ## Tibetan Mastiff, PRJNA302374

## Small breeds
fastq-dump -F --split-files SRR8996977  ## Yorkshire Terrier, PRJNA396033
fastq-dump -F --split-files DRR546744  ## Toy Poodle, PRJDB18013
fastq-dump -F --split-files DRR546790  ## Shih Tzu, PRJDB18013
fastq-dump -F --split-files DRR546743  ## West Highland White Terrier, PRJDB18013


########## FASTQC to Assess Quality of the Sequence Data
## Run FastQC on each FASTQ file to check quality
## Output: A folder with .html and .txt files for each sample
mkdir -p ${WD}/${RDQ}
fastqc *.fastq --outdir=${WD}/${RDQ}

########## Tarball the FastQC Results for Download
## Create a tarball to transfer the FastQC results to your local computer
cd ${WD}/${RDQ}
tar cvzf ${RDQ}.tar.gz ${WD}/${RDQ}/*

## After the job finishes, use scp or rsync to transfer the tarball:
## Example: scp ${MyID}@asax.asc.edu:/scratch/${MyID}/DogRNAseq/RawDataQuality/RawDataQuality.tar.gz .
## Then open the .html files locally to evaluate the quality of your raw data.
