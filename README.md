# CSeQTL

This repositry contains scripts used to process and analyse data for cell-specific (CS), allele-specific (AS) eQTL across in WHI particpants.

### Aims

-   To perform blood cell-type specific eQTL analyses in a multi-ethnic sample of women from WHI collected at the LLS-1 exam. 
-   To investigate cell-specific eQTLs and their potential associations with GWAS loci relevant to stroke and hypertension. 

### Methods

Cell-sepecific eQTL analysis is based on the methods described in Little, P., Liu, S., Zhabotynsky, V. et al. A computational method for cell type-specific expression quantitative trait loci mapping using bulk RNA-seq data. Nat Commun 14, 3030 (2023). https://doi.org/10.1038/s41467-023-38795-w. Below is a breif description of the analysis plan, along with the scripts used.

#### Steps

Pipeline is based on workflow described in https://github.com/pllittle/CSeQTLworkflow#eqtl-mapping. Generally .sh scripts are used to execute SLURM parallel processing steps. To run CSeQTL, phased genotype data are required, and must be formatted to tab-delim format. Phased genome data is then used to estimate allele-specific reads (ASE), in additon to total reads (TREC). Finally, cibersort deconvolution is used to estimate cell fractions, which can be input into CSeQTL function which jointly models cell-specific eQTL's, in conjuction with TREC, ASE, and prinicipal components. 

*wgs*
-   Pre-process phased genome information from TOPMed/WHI whole genome sequence (wgs) samples

*rnaseq*
-   Pre-proccess RNAseq data to obtain, total read counts per gene (TREC), and allele specific read counts
-   Prepare data format for Cibersort cell imputation by obtaining transcripts per million (TPM)
-   Run cibersort deconvolution analyis to obtain cell fractions and plot output
-   Obtain principal components (genetic and RNA seq) and other model covariates

*eqtl*
-   Using allele-specific reads, covariates and cell fractions (per sample) run eQTL mapping using CSeQTL in R
-   Downstream MT correction with eigenMT and BH for FDR
-   Summary plots and functional analysis


### Cohorts and data progress

At the time of starting this project, not all data were readily availbe across both cohorts we plan on using. Therefore this table will be updated accordingly as preprocessing completes:

![Screenshot](datasets.png)

#### File locations

- SCT RNAseq (processed rna-bam files): /fh/scratch/delete90/kooperberg_c/sct_rnaseq/release_files
- SCT WGS (raw wgs-bams): /fh/scratch/delete90/kooperberg_c/wgs_bam  

- LLS RNAseq (raw rna-bam files): /fh/scratch/delete90/kooperberg_c/lls_rna/bam_files
- LLS WGS (processed wgs-bam -> phased bcfs): /fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only  
  
- CSeQTL data: /fh/scratch/delete90/kooperberg_c/mjohnso5/CSeQTL
- CSeQTL test data (local): 〜/Documents/CSeQTL/data
- CSeQTL scripts (local): 〜/Documents/CSeQTL/scripts/CSeQTL
