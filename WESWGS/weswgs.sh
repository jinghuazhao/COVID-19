#!/usr/bin/bash

module load ceuadmin/stata

stata -b do weswgs.do

# This version replaces those from Stata by putting covariates before proteins
for panel in cvd2 cvd3 inf neu
do
  paste -d, work/${panel}-covariates.txt work/${panel}-protein.txt > work/${panel}.txt
done

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

bcftools query -l wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz | wc -l
bcftools query -l wgs/chr22/chr22.intervalwgs_v1_all_info.vcf.bgz | wc -l

(
  head -1 work/weswgs.txt
  bcftools query -l wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz | \
  grep -f - work/weswgs.txt
) > work/wes.txt

(
  head -1 work/weswgs.txt
  bcftools query -l wgs/chr22/chr22.intervalwgs_v1_all_info.vcf.bgz | \
  grep -f - work/weswgs.txt
) > work/wgs.txt