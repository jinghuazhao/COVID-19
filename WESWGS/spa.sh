#!/usr/bin/bash

cd burden_testing
export COHORT_NAME=INTERVAL
for panel in cv2 cvd3
do
  singularity exec -B ${pwd} burden.1.6.5 step1 ${COHORT_NAME} work/${panel}-wes.vcf.gz
  for chr in chr{1.22} chrX chrY
  do
    export chr=chr22
    export WGS=~/COVID-19/WESWGS/work/${panel}-wgs-${chr}.vcf.gz
    singularity exec -B $(pwd) burden.1.6.5 step1 $COHORT_NAME ${WGS}
  done
done

# SLURM
# sbatch --mem 2G -e step1.e -o step1.o --wrap "singularity exec -B $(pwd) burden.1.6.5 step1 $COHORT_NAME $VCF_PATH"

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
