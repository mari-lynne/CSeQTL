#!/bin/bash
#SBATCH --job-name=annovar
#SBATCH --time=1:00:00      
#SBATCH --mem-per-cpu=19G            
#SBATCH --cpus-per-task=4
#SBATCH --output=mono.out
#SBATCH --error=mono.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org



# Downloaded using annovar wget email url into annovar dir

cd /fh/working/hsu_l/Mari/annovar

chmod+ x *.pl

ml annovar


curl -sI http://www.openbioinformatics.org/annovar/download/avsnp151.txt.gz | grep -i Content-Length


annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene humandb/ \
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c humandb/ \
&& annotate_variation.pl -buildver hg38 -downdb -webfrom annovar 1000g2015aug  humandb/ \
&& annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp151 humandb/ \
&& annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/ \
&& annotate_variation.pl -buildver hg38 -downdb -webfrom annovar regsnpintron humandb/ \
&& annotate_variation.pl -buildver hg38 -downdb -webfrom annovar GTEx_v8_eQTL humandb/




# Run ANNOVAR ------------------------------------------------------------------

apptainer run --bind /fh/working/hsu_l/Mari/annovar:/data annovar_hg38.sif \
  table_annovar.pl /data/annovar_input.txt humandb/ \
  -buildver hg38 \
  -out /data/annotated \
  -remove \
  -protocol refGene \
  -operation g \
  -nastring .

  
apptainer run --bind /fh/working/hsu_l/Mari/annovar:/data annovar_hg38.sif \
table_annovar.pl /data/annovar_input.txt humandb/ \
-buildver hg38 \
-out /data/annotated_snp \
-remove \
-protocol esp6500siv2_all,clinvar_20180603,gnomad211_exome \
-operation f,f,f \
-nastring .

# esp6500siv2_al,

# table_annovar.pl mywork/VCF_files/proband.vcf\
#   humandb/ \
#   -buildver hg19 \
#   -out mywork/proband.annovar \
#   -remove \
#   -protocol refGeneWithVer,clinvar_20240917,gnomad211_exome,dbnsfp47a \
#   -operation g,f,f,f \
#   -arg '-hgvs',,, \
#   -polish -nastring . \
#   -vcfinput \
#   -intronhgvs 100