#!/bin/bash

######## FunGen Course Instructions ############
## Purpose: Trim sequencing adapters and low-quality regions from dog RNAseq data using Trimmomatic,
##          then use FastQC to evaluate the quality of the trimmed data.
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
## FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Input Data: Raw paired-end FASTQ files (R1 & R2) and adapter sequences.
## Output Data: Trimmed paired and unpaired FASTQ files, FastQC reports.
## For running on the Alabama Super Computer (ASC): https://hpcdocs.asc.edu/content/slurm-queue-system
## After editing, make executable with "chmod +x trim_fastqc_dog.sh" and run with "run_script trim_fastqc_dog.sh".
###############################################

########## SLURM Directives
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=trim_fastqc_%j.o

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.10.1

########## Define Variables
MyID=aubclsd0338              ## Replace with your ASC ID, e.g., aubtss
WD=/scratch/${MyID}/DogRNAseq   ## Working directory
DD=${WD}/RawData                ## Raw data directory
CD=${WD}/CleanData              ## Clean data directory
PCQ=PostCleanQuality            ## FastQC results directory
adapters=TruSeq3-PE.fa          ## Adapter file (update if needed)

########## Create Directories
mkdir -p ${CD}
mkdir -p ${WD}/${PCQ}

########## Generate List of Samples
cd ${DD}
ls | grep ".fastq" | cut -d "_" -f 1 | sort | uniq > list

########## Copy Adapter File
cp /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/share/trimmomatic/adapters/${adapters} .

########## Trim with Trimmomatic and Run FastQC
while read i
do
    java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar \
    PE -threads 6 -phred33 \
    "$i"_1.fastq "$i"_2.fastq \
    ${CD}/"$i"_1_paired.fastq ${CD}/"$i"_1_unpaired.fastq ${CD}/"$i"_2_paired.fastq ${CD}/"$i"_2_unpaired.fastq \
    ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    fastqc ${CD}/"$i"_1_paired.fastq --outdir=${WD}/${PCQ}
    fastqc ${CD}/"$i"_2_paired.fastq --outdir=${WD}/${PCQ}
done < list

########## Tarball the FastQC Results
cd ${WD}/${PCQ}
tar cvzf ${PCQ}.tar.gz ${WD}/${PCQ}/*
