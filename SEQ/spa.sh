#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export SEQ=${SCALLOP}/SEQ

function bgen()
{
  for chr in chr{{1..22},X}
  do
    export SLURM_ARRAY_TASK_ID=${chr}
    sbatch --job-name=_${weswgs}_${chr} --account CARDIO-SL0-CPU --partition cardio --qos=cardio \
           --mem=40800 --time=5-00:00:00 --export ALL \
           --output=${TMPDIR}/_${weswgs}_bgen-${chr}_%A_%a.out --error=${TMPDIR}/_${weswgs}_bgen-${chr}_%A_%a.err \
           --wrap ". ${SEQ}/bgen.sb"
  done
}

for weswgs in wes wgs
do
  export weswgs=${weswgs}
# bgen
  cut -f1,2 --complement ${SEQ}/work/${weswgs}.pheno | head -1 | tr '\t' '\n' > ${SEQ}/work/${weswgs}.varlist
  sbatch --export=ALL ${SEQ}/spa.sb
done

# <olink_protein>_<cohort>_<date_of_analysis>_<analyst_initials>.txt.bgz
# ACE2_INTERVAL_02112020_JHZ.txt.bgz

#Column no      Column name     Description
#1      SNP     rsID, or NA if unavailable
#2      CHR     Chromosome number
#3      POS     Physical base pair position on the chromosome (in b38 coordinates)
#4      N       Number of non-missing observations
#5      EFF_ALLELE      Allele whose effect is reported (beta estimates)
#6      OTHER_ALLELE    The other allele at the SNP
#7      EFF_ALLELE_FREQ Allele frequency of EFF_ALLELE
#8      BETA    Effect size estimate, to at least 5 d.p.
#9      SE      Standard error of the beta estimates, to at least 5 d.p.
#10     P_LRT   P-value LRT
#11     P_SCORE P-value score
