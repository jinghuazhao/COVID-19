# 15-5-2020 JHZ

export dir=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
export merged_imputation=/home/jhz22/rds/post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/merged_imputation
export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export ev=$ref/annot_INT_50PCs_pcs.txt

module load gcc/6

function phenofile()
{
  module load ceuadmin/stata
  stata -b do INTERVAL.do
  rm INTERVAL.log
}

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
           -incl-samples work/INTERVAL.samples -bgen-bits 8 -og work/INTERVAL-{}-rsid.bgen -os work/INTERVAL-{}-rsid.samples
    bgenix -g work/INTERVAL-{}-rsid.bgen -index -clobber
  '
  seq 22 | \
  parallel --env INF -C' ' '
  ( 
    cat ${INF}/csd3/ld.hdr
    bgenix -g work/INTERVAL-{}-rsid.bgen -list 2>&1 | \
    sed "1,9d" | \
    awk "
    {
      if (\$3<10) CHR= \"0\" \$3; else CHR=\$3
      SNP=\$2
      POS=\$4
      a1=toupper(\$6)
      a2=toupper(\$7)
      if (a1>a2) {t=a1; a1=a2; a2=t}
      snpid=\$3 \":\" POS \"_\" a1 \"/\" a2
      if (NF==7) print SNP, SNP, CHR, POS, \$6, \$7, SNP, snpid, \$3, POS, a1, a2
    }"
  ) > work/INTERVAL-{}.map
  qctool -g work/INTERVAL-{}-rsid.bgen -map-id-data work/INTERVAL-{}.map -og work/INTERVAL-{}.bgen
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
  grep -v 110001440667 work/INTERVAL.samples > work/INTERVAL-X.samples
  (
    awk -v idno=${idno} 'NR<idno' work/INTERVAL-X.vcf
    cat work/INTERVAL.idline
    awk -v idno=${idno} 'NR>idno' work/INTERVAL-X.vcf
  ) | \
  qctool -filetype vcf -g - -bgen-bits 8 -s work/INTERVAL-22.samples -excl-samples work/INTERVAL-X.excl-samples -og work/INTERVAL-X.bgen
  bgenix -g work/INTERVAL-X.bgen -index -clobber
}

function aggregate()
{
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
    cut -d' ' -f1,2,11-14 output/INTERVAL-*.gene.txt | head -1
    echo $(seq 22) X | \
    tr ' ' '\n' | \
    parallel -C' ' '
      sed '1d' output/INTERVAL-{}.gene.txt | \
      cut -d" " -f1,2,11-14
    '
  ) | gzip -f > output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.gene.20200516.txt.gz

  ) | gzip -f > gene-based.results
}

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

# gsutil cp UKBB.Doe.ANA5.1.ALL.EUR.154.1341.SAIGE.20200515.txt.gz gs://covid19-hg-upload-INTERVAL/

# Fill the form here,
# https://airtable.com/shrdJDwSHRZKXv45H
# Download the data from here,
# https://tinyurl.com/y97u49jz

# SBATCH -A CARDIO-SL0-CPU -p cardio_intr --qos=cardio_intr
# createSparseGRM.R --help
# step1_fitNULLGLMM.R --help
# step2_SPAtests.R --help

function init()
{
  cd 06-05-2020/INTERVAL
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVALdata_06MAY2020.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_20200506.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVALdata_P3_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
  cd -
}

(
  echo X:2699676_G_A
  echo X:2699683_G_A
  echo X:2699739_G_C
  echo X:2699898_C_CT
  echo X:2699906_T_G
  echo X:2699934_C_A
  echo X:2699968_A_G
  echo X:2700027_T_C
  echo X:2700040_T_C
  echo X:2700089_A_G
  echo X:2700147_C_T
  echo X:2700157_G_A
  echo X:2700185_C_T
  echo X:2700189_A_G
  echo X:2700194_C_G
  echo X:2700202_A_G
  echo X:74273085_AATCTT_A
) > work/INTERVAL-X.excl-rsids