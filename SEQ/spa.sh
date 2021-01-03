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

function fastGWAsetup()
{
# a sparse GRM from SNP data
  gcta-1.9 --bfile ${SEQ}/work/${weswgs} --make-grm --out ${SEQ}/work/spa/${weswgs}
  gcta-1.9 --grm ${SEQ}/work/spa/${weswgs} --make-bK-sparse 0.05 --out ${SEQ}/work/spa/${weswgs}-spgrm
  echo ${SEQ}/work/${weswgs}-chr{1..22}.bgen | tr ' ' '\n' > ${SEQ}/work/${weswgs}.bgenlist
  (
    head -2 ${SEQ}/work/${weswgs}-chrX.sample
    sed '1,2d' ${SEQ}/work/${weswgs}-chrX.sample | \
    cut -d' ' -f1-3 | \
    join -a1 -e NA -o1.1,1.2,1.3,2.2 - ${SEQ}/work/${weswgs}.sex
  ) > ${SEQ}/work/${weswgs}.sample

  bcftools query -l ${SEQ}/work/${weswgs}-chrX.vcf.gz | awk '{print $1,$1}' > ${SEQ}/work/spa/${weswgs}-chrX.idlist
  bcftools query -f "%ID\n" ${SEQ}/work/${weswgs}-chrX.vcf.gz > ${SEQ}/work/spa/${weswgs}-chrX.snplist

# fastGWA mixed model
  sed '1d' ${SEQ}/work/${weswgs}.pheno > ${SEQ}/work/spa/${weswgs}.pheno
}

for weswgs in wes wgs
do
  export weswgs=${weswgs}
  bgen
  cut -f1,2 --complement ${SEQ}/work/${weswgs}.pheno | head -1 | tr '\t' '\n' > ${SEQ}/work/${weswgs}.varlist
  fastGWAsetup
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
