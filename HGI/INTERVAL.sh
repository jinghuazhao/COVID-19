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
# GRM
  module load plink/2.00-alpha
  plink2 --bfile ${merged_imputation} --indep-pairwise 1000kb 1 0.1 --out work/INTERVAL
  plink2 --bfile ${merged_imputation} --make-bed --extract work/INTERVAL.prune.in --out work/INTERVAL
# Autosomes
  sed '1d' work/INTERVAL-covid.txt | \
  cut -d' ' -f1 > work/INTERVAL.samples
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
  (
    awk -v idno=${idno} 'NR<idno' work/INTERVAL-X.vcf
    cat work/INTERVAL.idline
    awk -v idno=${idno} 'NR>idno' work/INTERVAL-X.vcf
  ) | \
  qctool -filetype vcf -g - -bgen-bits 8 -og work/INTERVAL-X.bgen
  bgenix -g work/INTERVAL-X.bgen -index -clobber
}

# Analysis 5-susceptibility (phenotype name: ANA5):

module load gcc/6

# Single-variant association tests

step1_fitNULLGLMM.R \
   --plinkFile=INTERVAL \
   --phenoFile=INTERVAL-covid.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=INTERVAL \
   --nThreads=4

echo $(seq 22) X | \
tr ' ' '\n' | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=INTERVAL-{}.bgen \
   --bgenFileIndex=INTERVAL-{}.bgen.bgi \
   --minMAF=0.0001 \
   --minMAC=1 \
   --sampleFile=INTERVAL.samples \
   --GMMATmodelFile=INTERVAL.rda \
   --varianceRatioFile=INTERVAL.varianceRatio.txt \
   --SAIGEOutputFile=INTERVAL-{}.txt \
   --IsOutputNinCaseCtrl=TRUE \
   --IsOutputHetHomCountsinCaseCtrl=TRUE \
   --IsOutputAFinCaseCtrl=TRUE
'

# Gene-based association tests

createSparseGRM.R \
   --plinkFile=INTERVAL \
   --minMAF=0.0001 \
   --nThreads=4 \
   --outputPrefix=INTERVAL.sparseGRM \
   --numRandomMarkerforSparseKin=2000 \
   --relatednessCutoff=0.125

echo $(seq 22) X | \
tr ' ' '\n' | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=INTERVAL-{}.bgen \
   --bgenFileIndex=INTERVAL-{}.bgen.bgi \
   --chrom={} \
   --minMAF=0 \
   --minMAC=0.5 \
   --maxMAFforGroupTest=0.01 \
   --GMMATmodelFile=INTERVAL.rda \
   --varianceRatioFile=INTERVAL.varianceRatio.txt \
   --SAIGEOutputFile=INTERVAL-{}.SAIGE.gene.txt \
   --numLinesOutput=1 \
   --groupFile=INTERVAL-{}.eneBasedtest.txt \
   --sparseSigmaFile=INTERVAL-{}.sparseSigma.mtx \
   --IsOutputAFinCaseCtrl=TRUE \
   --IsSingleVarinGroupTest=TRUE \
   --IsOutputPvalueNAinGroupTestforBinary=TRUE \
   --IsAccountforCasecontrolImbalanceinGroupTest=TRUE
'
##

# SBATCH -A CARDIO-SL0-CPU -p cardio_intr --qos=cardio_intr
# createSparseGRM.R --help
# step1_fitNULLGLMM.R --help
# step2_SPAtests.R --help

(
  cut -f1,7,8,15,18,19 $ref/impute_*_interval.snpstats | \
  head -1
  parallel --env ref -C' ' 'sed '1d' $ref/impute_${chr}_interval.snpstats | cut -f1,7,8,15,18,19 $ref/impute_1_interval.snpstats'
) | gzip -f INTERVAL.snpstats.gz

cd 06-05-2020/INTERVAL
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVALdata_06MAY2020.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVALdata_P3_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
cd -
