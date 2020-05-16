# 15-5-2020 JHZ

export dir=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
export merged_imputation=/home/jhz22/rds/post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/merged_imputation
export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export ev=$ref/annot_INT_50PCs_pcs.txt

function phenofile()
{
  module load ceuadmin/stata
  stata -b do INTERVAL.do
}
rm INTERVAL.log

function genofile()
{
  sed '1d' work/INTERVAL-covid.txt | \
  cut -d' ' -f1 > work/INTERVAL.samples
  paste -d' ' work/INTERVAL.samples work/INTERVAL.samples > work/INTERAL.samples2
# GRM
  module load plink/2.00-alpha
  plink2 --bfile ${merged_imputation} --indep-pairwise 1000kb 1 0.1 --out work/INTERVAL
  plink2 --bfile ${merged_imputation} --make-bed --extract work/INTERVAL.prune.in -keep work/INTERVAL.samples2 --out work/INTERVAL
# Autosomes
  seq 22 | \
  parallel --env autosomes -C' ' '
    qctool -g ${autosomes}/impute_{}_interval.bgen -s ${autosomes}/interval.samples \
           -incl-samples work/INTERVAL.samples -bgen-bits 8 -og work/INTERVAL-{}.bgen -os work/INTERVAL-{}.samples
    bgenix -g work/INTERVAL-{}.bgen -index -clobber
  '
# HLA region
  plink2 --bfile ${merged_imputation} --chr 6 --from-bp 25392021 --to-bp 33392022 --make-bed --out work/INTERVAL-HLA
# Chromosome X
  bcftools query -l ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz | head
  awk '{print $1 "_" $1}' work/INTERVAL.samples | \
  bcftools view -S - ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -O v --force-samples > work/INTERVAL-X.vcf
  export idno=$(awk '/POS/ && /QUAL/ {print NR} ' work/INTERVAL-X.vcf)
  awk -v idno=${idno} 'NR==idno{print} ' work/INTERVAL-X.vcf > work/INTERVAL.idline
  (
    cat work/INTERVAL.samples | \
    parallel --dry-run -C' ' "
      export s={}_{};
      export t={};
      sed -i 's/'\"\${s}\"'/'\"\${t}\"'/g' work/INTERVAL.idline
    "
  ) | bash
  echo 110001440667 > work/INTERVAL-X.excl-samples
  (
    echo X:2699676_G_A
    echo X:2699683_G_A
  ) > work/INTERVAL-X.excl-rsids
  grep -v 110001440667 work/INTERVAL.samples > work/INTERVAL-X.samples
  (
    awk -v idno=${idno} 'NR<idno' work/INTERVAL-X.vcf
    cat work/INTERVAL.idline
    awk -v idno=${idno} 'NR>idno' work/INTERVAL-X.vcf
  ) | \
  qctool -filetype vcf -g - -bgen-bits 8 -s work/INTERVAL-22.samples -excl-samples work/INTERVAL-X.excl-samples -og work/INTERVAL-X.bgen
  bgenix -g work/INTERVAL-X.bgen -index -clobber
}

# Analysis 5-susceptibility (phenotype name: ANA5):

module load gcc/6

# Single-variant association tests
## Autosomal

step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL \
   --phenoFile=work/INTERVAL-covid.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=output/INTERVAL \
   --nThreads=4

echo $(seq 22) | \
tr ' ' '\n' | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=work/INTERVAL-{}.bgen \
   --bgenFileIndex=work/INTERVAL-{}.bgen.bgi \
   --minMAF=0.0001 \
   --minMAC=1 \
   --sampleFile=work/INTERVAL.samples \
   --GMMATmodelFile=output/INTERVAL.rda \
   --varianceRatioFile=output/INTERVAL.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.txt \
   --IsOutputNinCaseCtrl=TRUE \
   --IsOutputHetHomCountsinCaseCtrl=TRUE \
   --IsOutputAFinCaseCtrl=TRUE
'
## X

grep -v -f work/INTERVAL-X.excl-samples work/INTERVAL-covid.txt > work/INTERVAL-covid-X.txt
step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL \
   --phenoFile=work/INTERVAL-covid-X.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=output/INTERVAL-X \
   --nThreads=4

echo X | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=work/INTERVAL-{}.bgen \
   --bgenFileIndex=work/INTERVAL-{}.bgen.bgi \
   --idstoExcludeFile=work/INTERVAL-X.excl-rsids \
   --minMAF=0.0001 \
   --minMAC=1 \
   --sampleFile=work/INTERVAL-X.samples \
   --GMMATmodelFile=output/INTERVAL-X.rda \
   --varianceRatioFile=output/INTERVAL-X.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.txt \
   --IsOutputNinCaseCtrl=TRUE \
   --IsOutputHetHomCountsinCaseCtrl=TRUE \
   --IsOutputAFinCaseCtrl=TRUE
'

# Gene-based association tests

createSparseGRM.R \
   --plinkFile=work/INTERVAL \
   --minMAF=0.0001 \
   --nThreads=8 \
   --outputPrefix=output/INTERVAL.sparseGRM \
   --numRandomMarkerforSparseKin=2000 \
   --relatednessCutoff=0.125

echo $(seq 22) | \
tr ' ' '\n' | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=work/INTERVAL-{}.bgen \
   --bgenFileIndex=work/INTERVAL-{}.bgen.bgi \
   --chrom={} \
   --minMAF=0 \
   --minMAC=0.5 \
   --maxMAFforGroupTest=0.01 \
   --sampleFile=work/INTERVAL.samples \
   --GMMATmodelFile=output/INTERVAL.rda \
   --varianceRatioFile=output/INTERVAL.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.SAIGE.gene.txt \
   --numLinesOutput=1 \
   --sparseSigmaFile=output/INTERVAL.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   --idstoIncludeFile=work/INTERVAL-{}.incl \
   --groupFile=work/INTERVAL-{}.gene \
   --IsOutputAFinCaseCtrl=TRUE \
   --IsSingleVarinGroupTest=TRUE \
   --IsOutputPvalueNAinGroupTestforBinary=TRUE \
   --IsAccountforCasecontrolImbalanceinGroupTest=TRUE
'

## X

paste -d' ' work/INTERVAL-X.samples work/INTERVAL-X.samples > work/INTERVAL-X.samples2
plink2 --bfile work/INTERVAL --make-bed -keep work/INTERVAL-X.samples2 --out work/INTERVAL-X

createSparseGRM.R \
   --plinkFile=work/INTERVAL-X \
   --minMAF=0.0001 \
   --nThreads=8 \
   --outputPrefix=output/INTERVAL-X.sparseGRM \
   --numRandomMarkerforSparseKin=2000 \
   --relatednessCutoff=0.125

step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL-X \
   --phenoFile=work/INTERVAL-covid-X.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --invNormalize=TRUE \
   --outputPrefix=output/INTERVAL-X \
   --outputPrefix_varRatio=output/INTERVAL-X_varRatio \
   --sparseGRMFile=output/INTERVAL-X.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   --sparseGRMSampleIDFile=output/INTERVAL-X.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
   --nThreads=8 \
   --LOCO=FALSE \
   --skipModelFitting=FALSE \
   --IsSparseKin=TRUE \
   --isCateVarianceRatio=TRUE \
   --IsOverwriteVarianceRatioFile=TRUE

echo X | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=work/INTERVAL-{}.bgen \
   --bgenFileIndex=work/INTERVAL-{}.bgen.bgi \
   --minMAF=0 \
   --minMAC=0.5 \
   --maxMAFforGroupTest=0.01 \
   --sampleFile=work/INTERVAL-X.samples \
   --GMMATmodelFile=output/INTERVAL-X.rda \
   --varianceRatioFile=output/INTERVAL.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.SAIGE.gene.txt \
   --numLinesOutput=1 \
   --sparseSigmaFile=output/INTERVAL-X.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   --idstoIncludeFile=work/INTERVAL-X.incl \
   --groupFile=work/INTERVAL-X.gene \
   --IsOutputAFinCaseCtrl=TRUE \
   --IsSingleVarinGroupTest=TRUE \
   --IsOutputPvalueNAinGroupTestforBinary=TRUE \
   --IsAccountforCasecontrolImbalanceinGroupTest=TRUE
'

# Aggregation of results

(
  cut -d' ' -f1-3,5-9,11-14,17,21-26 output/INTERVAL-*.txt | head -1
  echo $(seq 22) X | \
  tr ' ' '\n' | \
  parallel -C' ' '
    sed '1d' output/INTERVAL-{}.txt | \
    cut -d" " -f1-3,5-9,11-14,17,21-26
  '
) | gzip -f > output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.20200515.txt.gz

(
  echo second part

) | gzip -f > gene-based.results

# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[gwas software].[YYYYMMDD].txt.gz
# sex=M/MALE, F/FEMALE, ALL

function gsutil_install()
# https://cloud.google.com/storage/docs/gsutil_install#linux
{
  wget https://storage.googleapis.com/pub/gsutil.tar.gz
  tar xvfz gsutil.tar.gz -C ..
  cd ../gsutil
  pip install pyasn1==0.4.8  --user
  python setup.py install --prefix=$HPC_WORK
}

# Request an account
# https://docs.google.com/forms/d/1eAaf-4XNYkplBo5Appbf8LHl2KHJyks9R4t0E3h0jII/viewform?edit_requested=true

gsutil cp UKBB.Doe.ANA5.1.ALL.EUR.154.1341.SAIGE.20200515.txt.gz gs://covid19-hg-upload-INTERVAL/

# Fill the form here,
# https://airtable.com/shrdJDwSHRZKXv45H
# Download the data from here,
# https://tinyurl.com/y97u49jz

# SBATCH -A CARDIO-SL0-CPU -p cardio_intr --qos=cardio_intr
# createSparseGRM.R --help
# step1_fitNULLGLMM.R --help
# step2_SPAtests.R --help

cd 06-05-2020/INTERVAL
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVALdata_06MAY2020.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVALdata_P3_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
cd -
