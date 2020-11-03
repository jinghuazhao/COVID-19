#!/usr/bin/bash

module load ceuadmin/stata
stata -b do weswgs.do
(
  head -1 work/weswgs.txt
  bcftools query -l wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz | grep -f - work/weswgs.txt
) > work/wes.txt
(
  head -1 work/weswgs.txt
  bcftools query -l wgs/chr22/chr22.intervalwgs_v1_all_info.vcf.bgz | grep -f - work/weswgs.txt
) > work/wgs.txt

export TMPDIR=${HPC_WORK}/work
for panel in cvd2 cvd3 inf neu
do
  paste -d, work/${panel}-covariates.txt work/${panel}-protein.txt > work/${panel}.txt
  cut -d, -f1 work/${panel}.txt | awk 'NR>1 && !/-999/' > work/${panel}-wes.samples
  cut -d, -f2 work/${panel}.txt | awk 'NR>1 && !/-999/' > work/${panel}-wgs.samples
  export WES=wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz
  bcftools query -l ${WES} | \
  grep -f work/${panel}-wes.samples | \
  bcftools view -S - ${WES} -O z -o work/${panel}-wes.vcf.gz
  export panel=${panel}
  sbatch --job-name=_${panel} --account CARDIO-SL0-CPU --partition cardio --qos=cardio --array=1-22 --mem=40800 --time=5-00:00:00 --export ALL \
         --output=${TMPDIR}/_${panel}_%A_%a.out --error=${TMPDIR}/_${panel}_%A_%a.err --wrap ". ${HOME}/COVID-19/WESWGS/weswgs.wrap"
  for chr in chrX chrY
  do
    export chr=${chr}
    export WGS=~/COVID-19/WESWGS/wgs/${chr}/${chr}.intervalwgs_v2_GT_only.vcf.bgz
    bcftools query -l ${WGS} | \
    grep -f work/${panel}-wgs.samples | \
    bcftools view -S - ${WGS} -O z -o work/${panel}-wgs-${chr}.vcf.gz
  done
done

cd work
ls *vcf.gz | parallel  -j5 -C' ' 'tabix -f {}' &

function chopped()
# This version avoids Stata but has lines chopped
{
panels="cvd2 cvd3 inf neurology"
export c=1
for i in $(for panel in ${panels}; do export dir=high_dimensional_data/Olink_proteomics_${panel}/qc;export f=$(ls ${dir});echo ${dir}/$f;done)
do
  export panel=$(basename ${i} | awk '{gsub(/olink|_|qc|.csv/,"");print}')
  echo ${i} ${panel}
  (
    awk -vFS="," 'NR==1{$1="";$2="";$(NF)=$(NF) "WES WGS sex age"; print}' ${i} | awk '{$1=$1};1'
    awk -vFS="," 'NR>1{$2="";print}' ${i} | sort -k1,1 | join - <(sed '1d' work/weswgs.txt | cut -f${c},5-8 | tr '\t' ' ' | sort -k1,1) | \
    awk '{$1=""};1' | \
    awk '$1=$1'
  ) > work/${panel}.txt
  export c=$(($c+1))
done
}
