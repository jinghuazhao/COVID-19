#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export SEQ=${SCALLOP}/SEQ

bgen <- function()
{
  export weswgs=${weswgs}
  sbatch ${SEQ}/bgen.sb
  for i in X Y
  do
    export SLURM_ARRAY_TASK_ID=${i}
    sbatch --job-name=_${weswgs}_chr${i} --account CARDIO-SL0-CPU --partition cardio --qos=cardio \
           --mem=40800 --time=5-00:00:00 --export ALL \
           --output=${TMPDIR}/_${weswgs}_bgen-chr${i}_%A_%a.out --error=${TMPDIR}/_${weswgs}_bgen-chr${i}_%A_%a.err \
           --wrap ". ${SCALLOP}/SEQ/bgen.sb"
  done
}

for weswgs in wes
do
  export weswgs=${weswgs}
# bgen()
  cut -f1 --complement ${SEQ}/work/${weswgs}.pheno | head -1 | tr '\t' '\n' > ${SEQ}/work/${weswgs}.varlist
  if [ ! -f ${SEQ}/work/${weswgs}.fam2 ]; then
    awk '{$1=$2};1' ${SEQ}/work/${weswgs}.fam > ${SEQ}/work/${weswgs}.fam2
  fi
  sbatch --export=ALL,weswgs ${SEQ}/spa.sb
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
