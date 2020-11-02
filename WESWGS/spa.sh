#!/usr/bin/bash

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

# <olink_protein>_<cohort>_<date_of_analysis>_<analyst_initials>.txt.bgz
# ACE2_pooled_MANOLIS_28102019_GP.txt.bgz

cd burden_testing
export COHORT_NAME=INTERVAL
export chr=22
export VCF_PATH=~/COVID-19/WESWGS/wgs/chr${chr}/chr${chr}.intervalwgs_v2_GT_only.vcf.bgz
singularity exec -B $(pwd) burden.1.6.5 step1 $COHORT_NAME $VCF_PATH

# slurm
# sbatch --mem 2G -e step1.e -o step1.o --wrap "singularity exec -B $(pwd) burden.1.6.5 step1 $COHORT_NAME $VCF_PATH"
