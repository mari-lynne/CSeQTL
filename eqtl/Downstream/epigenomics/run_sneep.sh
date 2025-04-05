#!/bin/bash

#SBATCH --job-name=mono
#SBATCH --time=1:00:00      
#SBATCH --mem-per-cpu=8G            
#SBATCH --cpus-per-task=20
#SBATCH --output=mono.out
#SBATCH --error=mono.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Set up files and directories
module load BEDTools
module load OpenMPI

cd /fh/working/hsu_l/Mari/SNEEP
export PATH=$PATH:/fh/working/hsu_l/Mari/SNEEP/src 

SNP_DIR=/fh/working/hsu_l/Mari/SNEEP/snp_bed_data
REF_DIR=/fh/working/hsu_l/Mari/SNEEP/resources
OUT_DIR="$PWD/results/mono/"

SNP_FILE="${SNP_DIR}/mono_snps_eqtl_nc.bed"
ATAC_FILE="${REF_DIR}/ATC.Bld.05.AllAg.Monocytes.bed"
# BKG_SNP_FILE="${REF_DIR}/dbSNPs_sorted.txt"
TRANSFAC_FILE="examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt"
REF_FILE="/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

# Run SNEEP
echo "Starting SNEEP"

differentialBindingAffinity_multipleSNPs \
-o ${OUT_DIR} \
-p 0.5 -c 0.001 -n 20 \
-b necessaryInputFiles/frequency.txt \
-x necessaryInputFiles/transition_matrix.txt \
-f ${ATAC_FILE} \
${TRANSFAC_FILE} \
${SNP_FILE} \
${REF_FILE} \
necessaryInputFiles/estimatedScalesPerMotif_1.9.txt

echo "SNEEP complete"

# -n 20 -j 500 -k ${BKG_SNP_FILE} -l 2 -q 0 \

# # SNP_position	var1	var2	rsID	MAF	peakPosition	TF	TF-binding_position	strand	effectedPositionInMotif	pvalue_BindAff_var1	pvalue_BindAff_var2	log_pvalueBindAffVar1_pvalueBindAffVar2	pvalue_DiffBindAff
# chr19:49534188-49534189	C	G	rs73582463	chr19	49533871:49534267-ID=SRX6777372	ZNF335	chr19:49534170-49534192	(f)	19	0.000879	0.025145	-3.353895	2.120003e-04	Name=ATAC-Seq%20(@%20T%20cells)
# chr6:42979274-42979275	G	C	rs9462860	chr6	42979026:42979401-ID=SRX10339025	KLF6	chr6:42979270-42979280	(r)	6	0.005700	0.000122	3.842087	8.175470e-04	Name=ATAC-Seq%20(@%20T%20cells)
# chr6:42979274-42979275	G	C	rs9462860	chr6	42979026:42979401-ID=SRX10339025	VEZF1	chr6:42979268-42979276	(r)	2	0.000147	0.004490	-3.422333	7.457274e-04	Name=ATAC-Seq%20(@%20T%20cells)

# Notes:
# Run per cell type, i.e filter T-cell GWAS SNPs then run with ATAC Seq data from T-cells
# Please use varying random seeds for runs with different input SNPs

# T-cell eQTLs are filtered for only eQTLs which are in chromatin open regions in T-cells (could also filter by sc-RNA expression of TFs)
# Then SNEEP calculates if the eQTL SNP is likely to modify transcription factor binding sites (loss or gain)
# Based on this score it defines regulatory or rSNPs and TFs with associated gain/loss of function from said SNPs

# Enrichment statistics
# SNPs randomly sampled from dbSNP database (this isn't providing a stat if cs-eqtl rSNPs are enriched for their corresponding cell type,
# would have to run analysis with other cell type chromatin data, at this point only descriptive)

# SNEEP can also link the disrupted TFs to target genes by using enhancer/gene interaction datatbases
# But we can do this ourselves by rejoining with eQTL gene data

