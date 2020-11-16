#!/usr/bin/bash

# module load ceuadmin/stata
# stata -b do weswgs.do
(
  head -1 work/weswgs.txt
  bcftools query -l wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz | grep -f work/weswgs.overlap -v | grep -w -f - work/weswgs.txt
) > work/wes.txt
sed '1d' work/wes.txt | cut -f2 > work/wes.samples
(
  head -1 work/weswgs.txt
  bcftools query -l wgs/chr22/chr22.intervalwgs_v2_GT_only.vcf.bgz | grep -w -f - work/weswgs.txt
) > work/wgs.txt
sed '1d' work/wgs.txt | cut -f1 > work/wgs.samples

export TMPDIR=${HPC_WORK}/work
export WES=wes/WES_QCed_Info_updated_4006_FINAL.vcf.gz
bcftools query -l ${WES} | \
grep -f work/wes.samples | \
bcftools view -S - ${WES} -O z -o work/wes.vcf.gz
tabix -f work/wes.vcf.gz
for chr in chr{1..22} chrX chrY; do bcftools view --regions ${chr} work/wes.vcf.gz -O z -o work/wes-${chr}.vcf.gz; done
sbatch --job-name=_wgs --account CARDIO-SL0-CPU --partition cardio --qos=cardio --array=1-22 --mem=40800 --time=5-00:00:00 --export ALL \
       --output=${TMPDIR}/_wgs_%A_%a.out --error=${TMPDIR}/_wgs_%A_%a.err --wrap ". ${HOME}/COVID-19/WESWGS/weswgs.wrap"
export SLURM_ARRAY_TASK_ID=X
${HOME}/COVID-19/WESWGS/weswgs.wrap
export SLURM_ARRAY_TASK_ID=Y
${HOME}/COVID-19/WESWGS/weswgs.wrap

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
