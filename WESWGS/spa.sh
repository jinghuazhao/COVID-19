#!/usr/bin/bash

# --- step1 ---
cd burden_testing
export TMPDIR=${HPC_WORK}/work
export COHORT=INTERVAL
module load singularity
for panel in cvd2 cvd3 inf neu
do
  export panel=${panel}
  burden.1.6.5 step1 ${COHORT}-${panel}-wes ../work/${panel}-wes.vcf.gz
  for chr in chr{1..22}
  do
    export chr=${chr}
    sbatch --job-name=_${panel}_step1 --account CARDIO-SL0-CPU --partition cardio --qos=cardio --mem=40800 --time=5-00:00:00 --export ALL \
           --output=${TMPDIR}/_${panel}_${chr}.out --error=${TMPDIR}/_${panel}_${chr}.err \
           --wrap ". burden.1.6.5 step1 ${COHORT}-${panel}-wgs ../work/${panel}-wgs-${chr}.vcf.gz"
  done
  for chr in chrX chrY
  do
    export chr=${chr}
    burden.1.6.5 step1 ${COHORT}-${panel}-wgs ../work/${panel}-wgs-${chr}.vcf.gz
  done
  for chr in chr{1..22} chrX chrY
  do zcat ${COHORT}-${panel}.${chr}.variantlist.gz; done | sed 's/^chr//' | gzip -f > ${COHORT}-${panel}-wgs.variantlist.gz
done
single_cohort_munge_variantlist 1 1

# --- geneset (annotation) data ---
# install axel, moreutils
# http://www.tucows.com/preview/231886/Axel
# https://joeyh.name/code/moreutils/

ln -sf geneset_data/ensembl-vep/INSTALL.pl
prepare-regions -o $(pwd)/geneset_data

# --- make (annotation) group files ---
#1. exon severe
#   variants with a "high" predicted consequence according to Ensembl (roughly equivalent to more severe than missense)
#2. exon CADD
#   all exonic variants (+50 bp outside of exons) weighted by CADD scores
#3. exon regulatory 
#   same as above, but weighted by Eigen scores (phred-scaled). Variants in regulatory regions that overlap with eQTL for that gene are also included.
#4. regulatory only
#   same as above but excluding exonic variants

export OPTS1="-g exon"
export OPTS2="-g exon -x 50 -s CADD"
export OPTS3="-g exon -x 50 -e promoter,enhancer,TF_bind -l promoter,enhancer,TF_bind -s EigenPhred"
export OPTS4="-e promoter,enhancer,TF_bind -l promoter,enhancer,TF_bind -s EigenPhred"

for group in OPTS{1..4}
do
  for i in {1..10}; do make-group-file -L chr21.genes -C config.txt -i ${COHORT_NAME}.chr21.variantlist.gz ${!group} -o -w $(pwd) -d 10 -c $i & done
done

find -name "group_file*.txt" -exec cat \{} \+ > group_file.txt
cut -f1 concat.group.file.txt | sort | uniq -c | awk '$1==1{print $2}'> singlesnp.genes.txt
fgrep -wvf singlesnp.genes.txt concat.group.file.txt > concat.group.file.filtered.txt

# --- GRM generation ---


# --- step2 ---

# <olink_protein>_<cohort>_<date_of_analysis>_<analyst_initials>.txt.bgz
# ACE2_INTERVAL_02112020_JHZ.txt.bgz

#Column no	Column name	Description
#1 	SNP	rsID, or NA if unavailable
#2 	CHR	Chromosome number
#3 	POS	Physical base pair position on the chromosome (in b38 coordinates)
#4 	N 	Number of non-missing observations
#5 	EFF_ALLELE 	Allele whose effect is reported (beta estimates)
#6	OTHER_ALLELE 	The other allele at the SNP
#7 	EFF_ALLELE_FREQ	Allele frequency of EFF_ALLELE
#8 	BETA	Effect size estimate, to at least 5 d.p.
#9	SE	Standard error of the beta estimates, to at least 5 d.p.
#10	P_LRT	P-value LRT
#11	P_SCORE	P-value score
