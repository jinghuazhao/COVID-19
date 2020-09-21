#!/usr/bin/bash

module load gcc/6
module load ceuadmin/stata

export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
export HGI=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI

function phenotype()
# step 1. SARS_CoV, age, sex, PCs
{
  export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
  export ev=$ref/annot_INT_50PCs_pcs.txt
  sed '1,2d' $autosomes/interval.samples | \
  cut -d' ' -f1 > work/INTERVAL.samples
  awk '{print $1 "_" $1}' work/INTERVAL.samples > work/INTERVAL-X.samples
  stata -b do INTERVAL.do
  export d=20200731
  for dir in ${d}-ANA_C1_V2 ${d}-ANA_C2_V2 \
             ${d}-male-ANA_C1_V2 ${d}-male-ANA_C2_V2 ${d}-female-ANA_C1_V2 ${d}-female-ANA_C2_V2 \
             ${d}-male-60-ANA_C1_V2 ${d}-male-60-ANA_C2_V2 ${d}-female-60-ANA_C1_V2 ${d}-female-60-ANA_C2_V2
  do
    cd ${dir}
    sed '1d' work/INTERVAL-covid.txt | \
    cut -d' ' -f1 > work/INTERVAL-covid.samples
    awk '{print $1, $1}' work/INTERVAL-covid.samples > work/INTERVAL-covid.samples2
    bcftools query -l ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz | \
    join - <(awk '{print $1"_"$1}' work/INTERVAL-covid.samples) -v2 > work/INTERVAL-covid-X.excl-samples
    awk '{if (NR>1) $1=$1 "_" $1};1' work/INTERVAL-covid.txt | \
    grep -v -f work/INTERVAL-covid-X.excl-samples > work/INTERVAL-covid-X.txt
    awk 'NR>1{print $1}' work/INTERVAL-covid-X.txt > work/INTERVAL-covid-X.samples
    awk '{print $1,$1}' work/INTERVAL-covid-X.samples > work/INTERVAL-covid-X.samples2
    grep -v -f work/INTERVAL-covid-X.excl-samples work/INTERVAL-X.FM > work/INTERVAL-covid-X.FM
    awk '$2=="M"' work/INTERVAL-covid-X.FM | \
    cut -f1 > work/INTERVAL-covid-X.samples-male
    cd -
  done
}

function make_bed()
# step 2. inference from directly genotyped variants 
{
# GRM
  module load plink/2.00-alpha
  export merged_imputation=/home/jhz22/rds/post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/merged_imputation
  plink2 --bfile ${merged_imputation} --indep-pairwise 1000kb 1 0.1 --out ${SCALLOP}/HGI/work/INTERVAL-covid
  export d=20200731
  for dir in ${d}-ANA_C1_V2 ${d}-ANA_C2_V2 \
             ${d}-male-ANA_C1_V2 ${d}-male-ANA_C2_V2 ${d}-female-ANA_C1_V2 ${d}-female-ANA_C2_V2 \
             ${d}-male-60-ANA_C1_V2 ${d}-male-60-ANA_C2_V2 ${d}-female-60-ANA_C1_V2 ${d}-female-60-ANA_C2_V2
  do
    cd ${dir}
    plink2 --bfile ${merged_imputation} --make-bed --extract ${SCALLOP}/HGI/work/INTERVAL-covid.prune.in \
           --keep work/INTERVAL-covid.samples2 --out work/INTERVAL-covid
    awk '{$1=$1 "_"_ $1; $2=$2 "_" $2;};1' work/INTERVAL-covid.fam > work/INTERVAL-covid.fam2
    plink --bed work/INTERVAL-covid.bed --bim work/INTERVAL-covid.bim --fam work/INTERVAL-covid.fam2 \
          --make-bed -keep work/INTERVAL-covid-X.samples2 --out work/INTERVAL-covid-X
    cd -
  done
}

# step 3. sample-specific bgen files
function autosomes_sbatch()
{
#!/usr/bin/bash

#SBATCH --job-name=_bgen
#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio
#SBATCH --array=1-22
#SBATCH --mem=28800
#SBATCH --time=5-00:00:00
#SBATCH --output=/rds/user/jhz22/hpc-work/work/_bgen_%A_%a.out
#SBATCH --error=/rds/user/jhz22/hpc-work/work/_bgen_%A_%a.err
#SBATCH --export ALL

export job=$SLURM_ARRAY_TASK_ID
export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37
export TMPDIR=/rds/user/jhz22/hpc-work/work

qctool -g ${autosomes}/imputed/impute_${job}_interval.bgen -s ${autosomes}/imputed/interval.samples \
       -incl-samples work/INTERVAL-covid.samples \
       -bgen-bits 8 -og output/INTERVAL-${job}.bgen -os output/INTERVAL-${job}.samples

bgenix -g output/INTERVAL-${job}.bgen -index -clobber
}

function X()
{
# Chromosome X, 58 samples do not exist in INTERVAL_X_imp_ann_filt_v2.vcf.gz
# Change of sample IDs
# bcftools query -l ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz | \
# tr '_' ' ' | \
# awk '{print $1"_"$1, $2}' | \
# bcftools reheader -s - ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -o work/INTERVAL-X-src.vcf.gz --threads 12
# bcftools annotate --set-id '%CHROM:%POS\_%REF\/%FIRST_ALT' work/INTERVAL-X-src.vcf.gz -O z -o output/INTERVAL-X-src.vcf.gz
# bcftools index -tf output/INTERVAL-X-src.vcf.gz
  bcftools annotate --set-id '%CHROM:%POS\_%REF\/%FIRST_ALT' ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -O z -o work/INTERVAL-X-src.vcf.gz
  bcftools index -tf work/INTERVAL-X-src.vcf.gz
  export d=20200731
  for dir in ${d}-ANA_C1_V2 ${d}-ANA_C2_V2 \
             ${d}-male-ANA_C1_V2 ${d}-male-ANA_C2_V2 ${d}-female-ANA_C1_V2 ${d}-female-ANA_C2_V2 \
             ${d}-male-60-ANA_C1_V2 ${d}-male-60-ANA_C2_V2 ${d}-female-60-ANA_C1_V2 ${d}-female-60-ANA_C2_V2
  do
    cd ${dir}
    grep -v -f work/INTERVAL-covid-X.excl-samples work/INTERVAL-covid-X.samples | \
    bcftools view -S - ${SCALLOP}/HGI/work/INTERVAL-X-src.vcf.gz -O v | \
    bgzip -cf > output/INTERVAL-X-ploidy.vcf.gz
    bcftools index -tf output/INTERVAL-X-ploidy.vcf.gz
    cd -
  done
}
#   tr '_' '\t' | \
#   cut -f1 | \
#   vcf-fix-ploidy -s work/INTERVAL-covid-X.FM | \

# step 4. single variant association tests + gene annodation/association analysis
# INTERVAL.sh + glist_hg19.sh

function Cx_V2_step1()
{
  step1_fitNULLGLMM.R \
     --plinkFile=${dir}/work/INTERVAL-covid \
     --phenoFile=${dir}/work/INTERVAL-covid.txt \
     --phenoCol=SARS_CoV \
     --covarColList=${covlist} \
     --sampleIDColinphenoFile=ID \
     --traitType=binary \
     --outputPrefix=${dir}/output/INTERVAL-covid \
     --nThreads=4 \
     --IsOverwriteVarianceRatioFile=TRUE
}

function Cx_V2_step1_X()
{
  step1_fitNULLGLMM.R \
     --plinkFile=${dir}/work/INTERVAL-covid-X \
     --phenoFile=${dir}/work/INTERVAL-covid-X.txt \
     --phenoCol=SARS_CoV \
     --covarColList=${covlist} \
     --sampleIDColinphenoFile=ID \
     --traitType=binary \
     --outputPrefix=${dir}/output/INTERVAL-covid-X \
     --nThreads=8 \
     --IsOverwriteVarianceRatioFile=TRUE
}

export d=20200731
export covlist=sex,age,age2,sexage,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20
for dir in ${d}-ANA_C1_V2 ${d}-ANA_C2_V2
do
  export dir=${dir}
  Cx_V2_step1
  Cx_V2_step1_X
done

export covlist=age,age2,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20
for dir in ${d}-male-ANA_C1_V2 ${d}-male-ANA_C2_V2 ${d}-female-ANA_C1_V2 ${d}-female-ANA_C2_V2
do
  export dir=${dir}
  Cx_V2_step1
  Cx_V2_step1_X
done

export covlist=PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20
for dir in ${d}-male-60-ANA_C1_V2 ${d}-male-60-ANA_C2_V2 ${d}-female-60-ANA_C1_V2 ${d}-female-60-ANA_C2_V2
do
  export dir=${dir}
  Cx_V2_step1
  Cx_V2_step1_X
done

for dir in ${d}-ANA_C1_V2 ${d}-ANA_C2_V2 \
           ${d}-male-ANA_C1_V2 ${d}-male-ANA_C2_V2 ${d}-female-ANA_C1_V2 ${d}-female-ANA_C2_V2 \
           ${d}-male-60-ANA_C1_V2 ${d}-male-60-ANA_C2_V2 ${d}-female-60-ANA_C1_V2 ${d}-female-60-ANA_C2_V2
do
  export dir=${dir}
  sbatch --wait ${SCALLOP}/HGI/autosomes.sb
  export check=$(echo ${dir} | grep female)
  if [ "${check}" == "" ]; then
     sbatch ${SCALLOP}/HGI/X.sb
  else
     sbatch ${SCALLOP}/HGI/X-female.sb
  fi
done

function aggregate()
# step 5. summary of results
{
# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[gwas software].[YYYYMMDD].txt.gz
# sex=M/MALE, F/FEMALE, ALL
# [dataset].[last name].[analysis_name].[freeze_number].[age].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# export snvResults=output/INTERVAL.Zhao.ANA_C2_V2.5.ALL.EUR.144.612.SAIGE.20200617.txt.gz
  for dir in ${d}-ANA_C1_V2 ${d}-ANA_C2_V2 \
             ${d}-male-ANA_C1_V2 ${d}-male-ANA_C2_V2 ${d}-female-ANA_C1_V2 ${d}-female-ANA_C2_V2 \
             ${d}-male-60-ANA_C1_V2 ${d}-male-60-ANA_C2_V2 ${d}-female-60-ANA_C1_V2 ${d}-female-60-ANA_C2_V2
  do
  export dir=${dir}
  cd ${dir}
  export snvResults=${SCALLOP}/HGI/${dir}.txt.gz
  (
    gunzip -c output/INTERVAL-*.txt.gz | \
    cut -d' ' -f1-3,5-9,11-14,17,21-26 | \
    head -1
    seq 22 | \
    tr ' ' '\n' | \
    parallel -j1 -C' ' '
      gunzip -c output/INTERVAL-{}.txt.gz | \
      sed "1d" | \
      cut -d" " -f1-3,5-9,11-14,17,21-28
    '
    echo X | \
    parallel -j1 -C' ' '
      gunzip -c output/INTERVAL-{}.txt.gz | \
      sed "1d" | \
      cut -d" " -f1-3,4-8,10-13,16,20-27 | \
      sed "s/X/23/"
    '
  ) | \
  awk '{if(NR>1){$1=$1+0;$3=$1+0 ":" $2 "_" $4 "/" $5}};1' | \
  gzip -f > ${snvResults}
  cd -
  done
  export geneResults=output/INTERVAL.Zhao.ANA_C2_V2.6.ALL.EUR.144.612.SAIGE.gene.20200617.txt.gz
  (
    cut -d' ' -f1,2,11-14 output/INTERVAL-*.gene.txt | head -1
    echo $(seq 22) X | \
    tr ' ' '\n' | \
    parallel -j1 -C' ' '
      sed '1d' output/INTERVAL-{}.gene.txt | \
      cut -d" " -f1,2,11-14
    '
  ) | gzip -f > ${geneResults}
# qqman
  export results=${snvResults}
  R --no-save -q <<\ \ END
    analysis <- "ANA_C1_V2"
    results <- Sys.getenv("results")
    require(qqman);
    tbl <- read.table(results,as.is=TRUE,header=TRUE);
    tbl <- within(tbl,{
       SNP <- rsid
       BP <- POS
       P <- p.value
    })
    tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
    qq <- paste0("output/",analysis,"_qq.png");
    png(qq,width=12,height=10,units="in",pointsize=4,res=300)
    qq(with(tbl,P))
    dev.off()
    manhattan <- paste0("output/",analysis,"_manhattan.png");
    png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
    manhattan(tbl,genomewideline=-log10(5e-8),suggestiveline=FALSE,ylim=c(0,25));
    dev.off();
  END
}

function upload()
{
# Request an account
# https://docs.google.com/forms/d/1eAaf-4XNYkplBo5Appbf8LHl2KHJyks9R4t0E3h0jII/viewform?edit_requested=true
# https://console.cloud.google.com/storage/browser/covid19-hg-upload-bugbank
  cd ${HGI}
  module load python/3.7
  virtualenv py37
  source py37/bin/activate
  pip install gsutil==4.50
  gsutil ls gs://covid19-hg-upload-uk--blood-donors-cohort
  gsutil cp $1 gs://covid19-hg-upload-uk--blood-donors-cohort 
  gsutil cp 20200731*/output/INTERVAL.Zhao* gs://covid19-hg-upload-uk--blood-donors-cohort
# HGI spreadsheet
  ls 20200731*/output/INTERVAL.Zhao* | xargs -l basename | xsel -i
# web: https://console.cloud.google.com/storage/browser/covid19-hg-upload-uk--blood-donors-cohort?project=covid-19-hg
# Fill the form (now uses tab in the spreadsheet),
# https://airtable.com/shrdJDwSHRZKXv45H
# Download data
# https://console.cloud.google.com/storage/browser/covid19-hg-analysis
# HGI results
# gs://covid19-hg-analysis/20200619/results/full
# gs://covid19-hg-analysis/20200619/results
# gsutil cp $1 gs://covid19-hg-upload-bugbank
}

# srun -A CARDIO-SL0-CPU -p cardio_intr --qos=cardio_intr -N1 -n1 -c4 --mem=50G -t 12:0:0 --pty bash -i
