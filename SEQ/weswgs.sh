#!/usr/bin/bash

module load ceuadmin/stata
stata -b do idmap.do
sed '1d' work/wes.txt | cut -f2 > work/wes.samples
sed '1d' work/wgs.txt | cut -f2 > work/wgs.samples

# dropouts
# WES
# EGAN00001217281 from omicsMap.csv does not appear in INTERVALdata_28FEB2020.csv
head -1 omicsMap.csv
grep EGAN00001217281 omicsMap.csv
grep EGAN00001217281 INTERVALdata_28FEB2020.csv
# WGS
# EGAN00001240484
# EGAN00001240488
# EGAN00001585278
# EGAN00001586195
# EGAN00001586196
# EGAN00001586198
head -1 omicsMap.csv
grep -e EGAN00001240484 -e EGAN00001240488 -e EGAN00001585278 -e EGAN00001586195 -e EGAN00001586196 -e EGAN00001586198 omicsMap.csv
grep -e EGAN00001240484 -e EGAN00001240488 -e EGAN00001585278 -e EGAN00001586195 -e EGAN00001586196 -e EGAN00001586198 INTERVALdata_28FEB2020.csv

# 50 overlaps between WES and WGS but only 5 with proteomics
awk '$2==$4' WGS-WES-Olink_ID_map_INTERVAL_release_28FEB2020.txt | wc -l
awk '$2==$4 && ($5!="NA"||$6!="NA"||$7!="NA"||$8!="NA")' WGS-WES-Olink_ID_map_INTERVAL_release_28FEB2020.txt | cut -f2 > work/weswgs.overlap
wc -l work/weswgs.overlap

echo WES -- only trimmed by availability of genotypic data
bcftools query -l wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz | grep -f work/weswgs.overlap | grep -f - work/wes.txt | wc -l
bcftools query -l wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz | grep -f work/wes.samples > work/wes.idmap
wc -l work/wes.idmap
grep -f work/wes.idmap work/weswgs.txt | cut -f2,14 | sort -k1,1 > work/wes.sex

echo WGS
bcftools query -l wgs/chr22/chr22.intervalwgs_v2_GT_only.vcf.bgz | grep -f - work/weswgs.overlap | grep -f - work/wgs.txt | wc -l
bcftools query -l wgs/chr22/chr22.intervalwgs_v2_GT_only.vcf.bgz | grep -f work/wgs.samples > work/wgs.idmap
wc -l work/wgs.idmap
grep -f work/wgs.idmap work/wgs.txt | wc -l
grep -f work/wgs.idmap work/weswgs.txt| cut -f3,14 | sort -k1,1 > work/wgs.sex

awk '$2!="NA" && ($5!="NA"||$6!="NA"||$7!="NA"||$8!="NA")' WGS-WES-Olink_ID_map_INTERVAL_release_28FEB2020.txt | wc -l
awk '$2!="NA" && ($5!="NA"||$6!="NA"||$7!="NA"||$8!="NA")' WGS-WES-Olink_ID_map_INTERVAL_release_28FEB2020.txt | cut -f2 | grep -f - work/wgs.idmap | wc -l
grep -f work/wgs.idmap WGS-WES-Olink_ID_map_INTERVAL_release_28FEB2020.txt | cut -f2| grep -v -f - work/wgs.idmap > work/wgs.exclude
grep -f work/wgs.exclude work/weswgs.txt
grep -f work/wgs.exclude work/wgs.txt

export TMPDIR=${HPC_WORK}/work
export WES=wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz
bcftools query -l ${WES} | \
grep -f work/wes.samples | \
bcftools view -S - ${WES} -O z -o work/wes.vcf.gz
tabix -f work/wes.vcf.gz
for chr in chr{1..22} chrX chrY
do
  bcftools view --regions ${chr} work/wes.vcf.gz -O z -o - | \
  bcftools sort -O z -o - | \
  bcftools annotate --set-id +'%CHROM:%POS\_%REF\_%FIRST_ALT' -O z -o work/wes-${chr}.vcf.gz
  bcftools index -f -t ${SCALLOP}/SEQ/work/wes-${chr}.vcf.gz
done
sbatch --job-name=_wgs --account CARDIO-SL0-CPU --partition cardio --qos=cardio --array=1-22 --mem=40800 --time=5-00:00:00 --export ALL \
       --output=${TMPDIR}/_wgs_%A_%a.out --error=${TMPDIR}/_wgs_%A_%a.err --wrap ". ${HOME}/COVID-19/SEQ/wgs.wrap"
export SLURM_ARRAY_TASK_ID=X
${HOME}/COVID-19/SEQ/weswgs.wrap
export SLURM_ARRAY_TASK_ID=Y
${HOME}/COVID-19/SEQ/weswgs.wrap

cd work
ls *chr*vcf.gz | parallel  -j5 -C' ' 'tabix -f {}' &

# https://sylabs.io/singularity/
module load singularity
singularity pull shub://hmgu-itg/burden_testing

sftp INTERVAL_VariantListUpload@146.107.169.40

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
