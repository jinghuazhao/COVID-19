# 18-5-2020 JHZ

export dir=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
export merged_imputation=/home/jhz22/rds/post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/merged_imputation
export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export ev=$ref/annot_INT_50PCs_pcs.txt
export snvResults=output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.20200515.txt.gz
export geneResults=output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.gene.20200515.txt.gz

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
  awk '{print $1, $1}' work/INTERVAL.samples > work/INTERAL.samples2
  awk '!/110001440667/' work/INTERVAL-covid.txt > work/INTERVAL-covid-X.txt
  awk '!/110001440667/' work/INTERVAL.samples > work/INTERVAL-X.samples
  awk '{print $1 "_" $1}' work/INTERVAL-X.samples > work/INTERVAL-X.samples2
# GRM
  module load plink/2.00-alpha
  plink2 --bfile ${merged_imputation} --indep-pairwise 1000kb 1 0.1 --out work/INTERVAL
  plink2 --bfile ${merged_imputation} --make-bed --extract work/INTERVAL.prune.in -keep work/INTERVAL.samples2 --out work/INTERVAL
  plink2 --bfile work/INTERVAL --make-bed -keep work/INTERVAL-X.samples2 --out work/INTERVAL-X
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
      CHR=\$3+0
      POS=\$4
      a1=toupper(\$6)
      a2=toupper(\$7)
      snpid=CHR \":\" POS \"_\" a1 \"/\" a2
      if (NF==7) print \$1, \$2, \$3, POS, \$6, \$7, \$1, snpid, CHR, POS, a1, a2
    }"
  ) | \
  awk "a[\$2]++==0" > work/INTERVAL-{}.map
  cut -d" " -f2 work/INTERVAL-{}.map > work/INTERVAL-{}.nodup
  qctool -g work/INTERVAL-{}-rsid.bgen -incl-rsids work/INTERVAL-{}.nodup -map-id-data work/INTERVAL-{}.map -bgen-bits 8 -og work/INTERVAL-{}.bgen
  bgenix -g work/INTERVAL-{}.bgen -index -clobber
  '
# Chromosome X -- NOTE that work/INTERVAL-X.samples has a single column
  awk '{print $1 "_" $1}' work/INTERVAL-X.samples | \
  bcftools view -S - ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -O v > work/INTERVAL-X.vcf
  export idno=$(awk '/POS/ && /QUAL/ {print NR} ' work/INTERVAL-X.vcf)
  awk -v idno=${idno} 'NR==idno{print} ' work/INTERVAL-X.vcf > work/INTERVAL.idline
# HLA region for future imputation
  plink2 --bfile ${merged_imputation} --chr 6 --from-bp 25392021 --to-bp 33392022 --make-bed --out work/INTERVAL-HLA
}

function X_vcf()
# Ploidy: GRCh37
# The pseudoautosomal regions 60,001-2,699,520 and 154,931,044-155,270,560 with the ploidy 2
# PAR1	X	60,001 	2,699,520	{ from=>1, to=>60_000, M=>1 }
# PAR2	X 	154,931,044 	155,260,560	{ from=>2_699_521, to=>154_931_043, M=>1 }
# https://en.wikipedia.org/wiki/Pseudoautosomal_region
{
  export idno=$(awk '/POS/ && /QUAL/ {print NR} ' work/INTERVAL-X.vcf)
  (
    cat work/INTERVAL-X.samples | \
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
  vcf-fix-ploidy -s work/INTERVAL-X.FM | \
  bgzip -f > work/INTERVAL-X-ploidy.vcf.gz
  bcftools index -tf work/INTERVAL-X-ploidy.vcf.gz
}

function aggregate()
{
  (
    cut -d' ' -f1-3,5-9,11-14,17,21-26 output/INTERVAL-*.txt | head -1
    echo $(seq 22) X | \
    tr ' ' '\n' | \
    parallel -j1 -C' ' '
      sed '1d' output/INTERVAL-{}.txt | \
      cut -d" " -f1-3,5-9,11-14,17,21-26
    '
  ) | gzip -f > ${snvResults}
  (
    cut -d' ' -f1,2,11-14 output/INTERVAL-*.gene.txt | head -1
    echo $(seq 22) X | \
    tr ' ' '\n' | \
    parallel -j1 -C' ' '
      sed '1d' output/INTERVAL-{}.gene.txt | \
      cut -d" " -f1,2,11-14
    '
  ) | gzip -f > ${geneResults}
}

# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[gwas software].[YYYYMMDD].txt.gz
# sex=M/MALE, F/FEMALE, ALL

function annotate()
{
  (
    awk 'BEGIN{print "##fileformat=VCFv4.0"}'
    awk 'BEGIN{print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"}'
    gunzip -c ${snvResults} | \
    awk -v OFS="\t" 'NR>1{print $1,$2,$1 ":" $2 "_" $4 "/" $5, $4, $5, ".", ".", $8}'
  ) > work/INTERVAL.vepinput
  export ANNOVAR=${HPC_WORK}/annovar
  export LEFTEE=${HPC_WORK}/loftee
  export POLYPHEN=$HPC_WORK/polyphen-2.2.2
  export VEP=${HPC_WORK}/ensembl-vep
  export TMPDIR=/rds/user/jhz22/hpc-work/work
  export dbNSFP_1=clinvar_id,clinvar_clnsig,clinvar_review,clinvar_trait,1000Gp3_EUR_AF,CADD_phred,Eigen-PC-phred_coding,ExAC_NFE_AF,LRT_pred
  export dbNSFP_2=FATHMM_pred,GERP++_RS,GTEx_V7_tissue,MutPred_protID,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,SIFT4G_pred,fathmm-MKL_coding_pred
  export dbNSFP_3=rs_dbSNP151,fathmm-MKL_coding_pred,gnomAD_exomes_NFE_AF,gnomAD_genomes_NFE_AF
  export dbNSFP_fields=${dbNSFP_1},${dbNSFP_2},${dbNSFP_3}
  vep -i work/INTERVAL.vepinput -o work/INTERVAL.dbNSFP --cache --distance 500000 --force --offline --pick --tab \
      --plugin LoF,loftee_path:${LOFTEE},human_ancestor_fa:human_ancestor.fa.gz \
      --plugin dbNSFP,${VEP}/dbNSFP4.0a/dbNSFP4.0a.gz,${dbNSFP_fields}
}

function init()
{
  cd 06-05-2020/INTERVAL
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVALdata_06MAY2020.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_20200506.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVALdata_P3_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
  cd -
}

function X_bgen_defunct()
{
  (
    cat work/INTERVAL-X.samples | \
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
  awk '!/110001440667/' work/INTERVAL-22.samples > work/INTERVAL-X.samples2
  qctool -filetype vcf -g - -bgen-bits 8 -s work/INTERVAL-X.samples2 -og work/INTERVAL-X.bgen
  bgenix -g work/INTERVAL-X.bgen -index -clobber
# requires the following SAIGE flags:
# --idstoExcludeFile=work/INTERVAL-X.excl-rsids
# Examples
  (
    echo X:2699676_G_A;  echo X:2699683_G_A;  echo X:2699739_G_C
    echo X:2699898_C_CT; echo X:2699906_T_G;  echo X:2699934_C_A
    echo X:2699968_A_G;  echo X:2700027_T_C;  echo X:2700040_T_C
    echo X:2700089_A_G;  echo X:2700147_C_T;  echo X:2700157_G_A
    echo X:2700185_C_T;  echo X:2700189_A_G;  echo X:2700194_C_G
    echo X:2700202_A_G;  echo X:74273085_AATCTT_A
  ) > work/INTERVAL-X.excl-rsids
}

# createSparseGRM.R --help
# step1_fitNULLGLMM.R --help
# step2_SPAtests.R --help
# SBATCH -A CARDIO-SL0-CPU -p cardio_intr --qos=cardio_intr

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

# gsutil cp ${svResults} ${geneResults} gs://covid19-hg-upload-INTERVAL/

# Fill the form here,
# https://airtable.com/shrdJDwSHRZKXv45H
# Download the data from here,
# https://tinyurl.com/y97u49jz
