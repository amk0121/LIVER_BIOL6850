# Functional-Genomics
RNAseq codes
Experimental Design: 


Paragraph Description of Sample Selection Rules:
For this study, we aim to explore the role of the IIS molecular network in body size evolution by analyzing gene expression differences in liver tissue between large and small dog breeds using RNA sequencing (RNAseq) data from the NCBI Sequence Read Archive (SRA). To ensure a robust and reproducible analysis, we established the following sample selection criteria: (1) Tissue specificity: Only liver tissue samples were included to maintain a consistent biological context relevant to metabolic regulation by IIS. (2) Health status: We restricted our selection to healthy control samples, excluding those from diseased dogs or experimental treatments to avoid confounding gene expression changes unrelated to size. (3) Size classification: Breeds were classified as large (>50 lbs or >20 inches at the shoulder) or small (<20 lbs or <15 inches at the shoulder) based on standard adult metrics from veterinary and kennel club data. (4) Diversity: To prevent bias from over-representation, we limited samples to no more than two per breed per size group. (5) Sample size: A minimum of four samples per group (large and small) was targeted to provide sufficient statistical power for differential expression analysis using DESeq2. These rules aim to isolate size-related differences in IIS gene expression while minimizing extraneous variability.
Technical Noise
To control technical noise, we prioritized RNAseq samples sequenced on Illumina platforms to ensure consistency in sequencing technology. Where possible, samples were selected from the same BioProject to reduce batch effects; however, when multiple BioProjects were necessary, we verified that library preparation methods (e.g., poly-A enrichment) and read lengths (e.g., paired-end 100 bp) were comparable. All samples will be aligned to the Tasha the Boxer reference genome (GCF_000002285.5) to standardize mapping and further minimize technical variability. Batch effects will be assessed and corrected during DESeq2 analysis if significant.
Biological Noise: Age, Sex, Location/Region
Biological noise was minimized by controlling for key variables:
Age: Samples were restricted to Adult dogs to reduce gene expression variability due to developmental or aging effects.
Sex: We aimed for only males to reduce sex effect.
Location/Region: While geographic origin was not always specified in SRA metadata, we assumed samples from different BioProjects might reflect regional variation. This was mitigated by focusing on breed and size as the primary variables, with any residual variation modeled as a covariate in DESeq2 if needed.




Parameters for Large and Small Dogs
Large Dogs: Breeds with an average adult weight >50 lbs or height >20 inches at the shoulder. Examples: Newfoundland, Belgian Malanois, Labrador Retriever, Tibetan Mastiff.
Small Dogs: Breeds with an average adult weight <20 lbs or height <15 inches at the shoulder. Examples: Yorkshire Terrier, Toy Poodle, Shih Tzu, West highland white terrier. 
Breeds with ambiguous size classifications (e.g., Beagle) were excluded to maintain a clear size distinction.
Summary Table of RNAseq Samples
Below is a table of RNAseq samples selected for this study, sourced from BioProjects in the NCBI SRA database. All are liver tissue samples from healthy control dogs.


# Dog Liver Transcriptome Data

This repository contains metadata for sequencing runs of dog liver transcriptomes from various projects and breeds. The table below summarizes the SRR (Sequence Read Archive) runs, associated BioProjects, breeds, size classes, sexes, ages, groups, and tissues.

| SRR        | Project    | Breed                     | Size Class | Sex  | Age   | Group | Tissue |
|------------|------------|---------------------------|------------|------|-------|-------|--------|
| SRR8996966 | PRJNA396033| Newfoundland              | Large      | Male | 11    |       | Liver  |
| SRR8997009 | PRJNA396033| Belgian Malinois          | Large      | Male | 3     |       | Liver  |
| SRR5889337 | PRJNA396033| Labrador Retriever        | Large      | Male | 16    |       | Liver  |
| SRR2960319 | PRJNA302374| Tibetan Mastiff           | Large      | Male | Adult |       | Liver  |
| SRR8996977 | PRJNA396033| Yorkshire Terrier         | Small      | Male | 11    |       | Liver  |
| DRR546744  | PRJDB18013 | Toy Poodle                | Small      | Male | 7     |       | Liver  |
| DRR546790  | PRJDB18013 | Shih Tzu                  | Small      | Male | 11    |       | Liver  |
| DRR546743  | PRJDB18013 | West Highland White Terrier | Small      | Male | 11    |       | Liver  |

## Notes
- **SRR/DRR**: Sequence Read Archive run accession numbers.
- **Project**: BioProject accession numbers.
- **Breed**: Dog breed for each sample.
- **Size Class**: Size classification of the breed (Large or Small).
- **Sex**: Sex of the dog (all samples are Male in this dataset).
- **Age**: Age of the dog in years, or "Adult" if not specified numerically.
- **Group**: Experimental group (not specified in this dataset).
- **Tissue**: Tissue sampled (all samples are from the liver).


## Project Management Plan
Task Assignments
All group members will contribute to coding, method development, and result interpretation, but primary responsibilities are assigned as follows:


Run the Standard Bioinformatic Pipeline:
Primary: ?
Support: All members will review pipeline outputs (e.g., FASTQC, alignment rates).


Run Differential Gene Expression Using DESeq2:
Primary: ?
Support: All members will validate differentially expressed genes and statistical significance.


Run GSEA for KEGG Pathways:
Primary: Ava
Support: All members will interpret enriched pathways biologically.


Run GSEA for IIS Molecular Pathway:
Primary: Ava
Support: All members will confirm IIS gene set accuracy and relevance.


Evaluate Protein-Coding Sequence Variation:
Primary: ?
Support: All members will analyze variants in top IIS regulators (e.g., using VCF files or Ensembl).


Manuscript Organization:
Lead: ?
Support: All members will draft their sections and review the final manuscript.


Presentation Organization:
Lead: ?
Support: All members will create slides and rehearse the presentation.


GitHub and Readme Management:
Lead: Anne Marie (GitHub: https://github.com/amk0121/LIVER_BIOL6850)
Support: All members will upload annotated code and update the README.


## Timeline for Project Completion
March 10-15: Finalize sample selection, index the Tasha genome, and run the bioinformatic pipeline (e.g., trimming, alignment, counting).
March 16-19: Generate count matrix and perform DESeq2 differential expression analysis.
March 20-25: Run GSEA for KEGG pathways and IIS pathway.
March 26-31: Assess protein-coding sequence variation in top IIS genes.
April 1: Group check-in to review progress and troubleshoot.
April 2-15: Draft manuscript sections (introduction, methods, results, discussion).
April 16-17: Finalize manuscript draft for peer review.
April 18-21: Conduct peer reviews and revise manuscript.
April 22-30: Prepare and rehearse presentation.
May 1: Deliver group presentation.
May 2-5: Finalize manuscript and GitHub repository.
May 5, 5 PM: Submit final manuscript.
May 7: Complete group reflection survey.

## 1_download_fastqc.sh

**Instructions to Run the Script:**

Save the script (e.g., as download_fastqc.sh).

Replace [Your_ASC_ID] with your actual ASC ID.

Make it executable: chmod +x download_fastqc.sh.

Submit it to the ASC: run_script download_fastqc.sh.

Monitor the job with qstat -u aubclsd0338

This script will download the paired-end FASTQ files for all eight samples and generate FastQC reports, which you can then transfer to your local computer for review.

scpaubclsd0338@asax.asc.edu:/scratch/aubclsd0338/DogRNAseq/RawDataQuality/RawDataQuality.tar.gz ~/Desktop

Review Quality: Unzip the tarball and open the HTML files to check the FastQC reports.

## 2_trim_fastqc_dog.sh

**How to Run the Script**

Save It: Save as trim_fastqc_dog.sh.

Update MyID: Replace [Your_ASC_ID] with your ASC username.

Make Executable: Run chmod +x trim_fastqc_dog.sh.

Submit the Job: Run run_script trim_fastqc_dog.sh.

Check Progress: Use qstat -u aubclsd0338.

Get Results: After it finishes, download the tarball:

scpaubclsd0338@asax.asc.edu:/scratch/aubclsd0338/DogRNAseq/PostCleanQuality/PostCleanQuality.tar.gz ~/Desktop

## 3_map_count_dog.sh

**How to Run the Script**

Save the Script:

Copy the script into a file named map_count_dog.sh.

Update Your ID:

Replace [Your_ASC_ID] with your ASC username (e.g., aubtss).

Make Executable:

Run chmod +x map_count_dog.sh in your terminal.

Submit the Job:

Run run_script map_count_dog.sh to submit it to SLURM.

Monitor Progress:

Check the job status with qstat -u aubclsd0338.

**Retrieve Results:**

After completion, find the count matrices .csv and stats files .txt in RESULTSD (e.g., /home/[Your_ASC_ID]/DogRNAseq/Counts_H_S).

Copy these to your local computer using scp or a file transfer tool.

Output Files

Mapping Statistics: ${MAPD}/*_Stats.txt (e.g., SRR8996966_Stats.txt).

Count Files: ${COUNTSD}/*/. Each sample has a subdirectory with .gtf files.

Count Matrices: ${RESULTSD}/*.csv (gene and transcript count matrices).

**Troubleshooting Tips:**

File Extensions: If the genome or annotation files have different extensions (e.g., .fasta or .gff3), update the cp and gffread lines accordingly.

Resource Limits: If the job fails due to memory or time, increase --mem (e.g., to 200G) or --time (e.g., to 24:00:00) in the SLURM directives.

Module Versions: Confirm that the loaded module versions match those available on ASC.

Next Steps

Once the script completes, use the count matrices (.csv) for downstream analysis, such as:

Differential Expression: Analyze with DESeq2 in R to compare gene expression between conditions (e.g., large vs. small dog breeds).

Pathway Analysis: Perform Gene Set Enrichment Analysis (GSEA) to explore pathways like IIS.



## 4_Deseq2_Dog

1. Experimental Design
   
Design Formula: Changed from ~treatment to ~BioProject + size. This accounts for batch effects (via BioProject) and compares gene expression between large and small breeds (size).

Factor Levels: Updated to dds$size with levels "Small" (reference) and "Large", replacing the original "Ad_lib" and "Caloric_restriction".

3. Input Files
   
Count Data: Kept as gene_count_matrix.csv, assuming it’s the output from a previous step (e.g., prepDE.py3). Uncomment preprocessing lines if your file has extra columns (e.g., length).

Metadata: Assumes PHENO_DATA.txt contains columns: sample (row names), size (Large/Small), and BioProject (e.g., PRJNA396033). Adjust the file name or structure if different.

Annotation: Updated to dog_annotation.csv, which should map gene IDs to gene names or symbols for your dog genome.

5. Visualizations
   
MA Plot Title: Updated to "DESeq2: Large vs Small Dog Breeds" for clarity.

Plot Counts: Suggest replacing "YOUR_GENE_ID" with a gene of interest, like INSR (from the IIS pathway), relevant to your project.

Heatmaps and PCA: Updated annotation to use size and BioProject instead of treatment and type.

7. Output Files
   
DGE Results: Changed to DGE_results_dog.csv to reflect your project.

GSEA and Cytoscape Files: Kept the structure but updated to use your dog annotation file.

9. Comments
    
Retained educational questions (e.g., "What does each column mean?") but tailored the context to your project where applicable.

Prerequisites

Metadata File (PHENO_DATA.txt): Create this file with columns sample, size, and BioProject. Example:

| sample     | size  | BioProject  |
|------------|-------|-------------|
| SRR8996966 | Large | PRJNA396033 |
| SRR8996977 | Small | PRJNA396033 |
| DRR546744  | Small | PRJDB18013  |

Count Matrix (gene_count_matrix.csv): Ensure sample names match PHENO_DATA.txt.

Annotation File (dog_annotation.csv): Should have a gene_id column matching your count matrix and a Name column with gene symbols.

*Next Steps* 

Run the Script: Execute in R or RStudio after setting the working directory and preparing input files.

Check Results: Look for differentially expressed genes related to the IIS pathway (e.g., INSR, IGF1) in DGE_results_dog.csv.

Downstream Analysis: Use DGErankName.rnk for GSEA to test pathway enrichment and NormTransExp_Anno_Names.txt for Cytoscape visualization.
