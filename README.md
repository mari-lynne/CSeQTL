# CSeQTL

This repositry contains scripts used to process and analyse data for cell-specific (CS), allele-specific (AS) eQTL across WHI particpants.

### Aims

-   To perform blood cell-type specific eQTL analyses in a multi-ethnic sample of women from WHI collected at the LLS-1 exam. 
-   To investigate cell-specific eQTLs and their potential associations with GWAS loci relevant to stroke and hypertension. 

### Methods

Cell-sepecific eQTL analysis is based on the methods described in Little, P., Liu, S., Zhabotynsky, V. et al. A computational method for cell type-specific expression quantitative trait loci mapping using bulk RNA-seq data. Nat Commun 14, 3030 (2023). https://doi.org/10.1038/s41467-023-38795-w. Below is a breif description of the analysis plan, along with the scripts used.

#### Steps

This Pipeline is based on workflow described in https://github.com/pllittle/CSeQTLworkflow#eqtl-mapping and includes additional preprocessing/ data formatting scripts. Generally .sh scripts are used to execute SLURM parallel processing steps.

To run CSeQTL, phased genotype data are used to estimate allele-specific reads (ASE), in additon to total reads (TREC) with the [asSeq](https://github.com/Sun-lab/asSeq) package. The [CibersortX](https://cibersortx.stanford.edu/) deconvolution algorithm is used to estimate cell fractions, which can be input into [CSeQTL](https://github.com/pllittle/CSeQTL) which jointly models cell-specific eQTL's, in conjuction with TREC, ASE, prinicial components and other covariates. 

*wgs*
-   Pre-process whole genome sequence (WGS) data from TOPMed/WHI samples (Freeze 12)
-   Phase genomes using Eagle 2 and TOPMed freeze 10 data as the reference genome

*rnaseq*
-   Pre-proccess RNAseq data to obtain, total read counts per gene (TREC), and allele specific read counts
-   Prepare data format for Cibersort cell imputation by obtaining transcripts per million (TPM)
-   Run cibersort deconvolution analyis to obtain cell fractions and plot output
-   Obtain principal components (genetic and RNA seq) and other model covariates

*eqtl*
-   Using allele-specific reads, covariates and cell fractions (per sample) run eQTL mapping using CSeQTL in R
-   Downstream MT correction with eigenMT and BH for FDR
-   Summary plots and downstream functional analysis


### Cohorts and data progress

2063 post-menopausal women aged between 65-95 years who were originally recruited as part of the Women’s Health Initiative study. 
WGS and RNAseq was collected from whole blood in two batches (LLS-1 and LLS-2 a.k.a SCT study)
WGS harmonized as part of TOPMed Freeze 12, jointly variant called.
RNAseq bam files processed identically as part of the ASeq/CSeQTL pipeline, with RNA-batch included as a covariate in the final model.

| Variable           | LLS-1          | LLS-2          | Total           |
|--------------------|----------------|----------------|-----------------|
| **Race/Ethnicity** |                |                |                 |
| African            | 326            | 718            | 1044            |
| European           | 870            | 0              | 870             |
| Hispanic           | 92             | 57             | 149             |
| **SCT Status**     |                |                |                 |
| 1                  | 1288           | 647            | 1935            |
| 0                  | 0              | 128            | 128             |
| **Age** Mean (SD)  | 81.0 (6.4)     | 75.6 (5.8)     | 79.0 (6.7)      |
| **BMI** Mean (SD)  | 28.2 (5.8)     | 29.5 (6.2)     | 28.7 (6.0)      |


#### File locations

**WGS**
- All TOPMed bam files: (raw-bam files): /fh/scratch/delete90/kooperberg_c/topmed_freeze12/minDP10
- LLS-1/2 bam files: /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study
- LLS-1/2 phased bcf files: /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased
- WGS PC's: /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/phased/pca

**RNA-seq**  
- LLS-1 raw bam files: /fh/scratch/delete90/kooperberg_c/lls_rna/bam_files/bam_files/
- LLS-2 raw bam files: /fh/scratch/delete90/kooperberg_c/sct_rnaseq/release_files/
- LLS-1/2 allele specific read counts: /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/rnaseq/ASE

- LLS-1/2 pre-procesed/QC rnaseq data: /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/sct_lls_merged.rds

**eQTL**
- CSeQTL results and summaries: /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/merged_study/eQTL/combined
- EigenMT: /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/eqtl/merged_study/eigenMT

**Scripts** 
- ~/CSeQTL/scripts_new
