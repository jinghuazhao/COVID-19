#!/usr/bin/bash

module load gcc/6
module load ceuadmin/stata

export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data

function phenotype()
{
  export dir=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
  export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
  export ev=$ref/annot_INT_50PCs_pcs.txt
  export merged_imputation=/home/jhz22/rds/post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/merged_imputation
  sed '1,2d' $autosomes/interval.samples | \
  cut -d' ' -f1 > work/INTERVAL.samples
  stata -b do INTERVAL.do
}

function genotype()
{
# Chromosome X
  vcf-fix-ploidy -s work/INTERVAL-X.FM | \
  bgzip -cf > work/INTERVAL-X-ploidy.vcf.gz
  bcftools index -tf work/INTERVAL-X-ploidy.vcf.gz
# HLA region for future imputation
  plink2 --bfile ${merged_imputation} --chr 6 --from-bp 25392021 --to-bp 33392022 --make-bed --out work/INTERVAL-HLA
}

function GRM()
{
  sed '1d' work/INTERVAL-covid.txt | \
  cut -d' ' -f1 > work/INTERVAL-covid.samples
  awk '{print $1, $1}' work/INTERVAL-covid.samples > work/INTERVAL-covid.samples2
  awk '{print $1 "_" $1}' work/INTERVAL-covid.samples | awk '{print $1,$1}' > work/INTERVAL-covid.id2
  awk '!/110001440667/' work/INTERVAL-covid.txt > work/INTERVAL-covid-X.txt
  awk '!/110001440667/' work/INTERVAL-covid.samples > work/INTERVAL-covid-X.samples
  awk '!/110001440667/' work/INTERVAL-covid.id2 > work/INTERVAL-covid-X.id2
  awk '{print $1 "_" $1}' work/INTERVAL-covid.id2 | awk '!/110001440667/' > work/INTERVAL-covid-X.samples2
# GRM
  module load plink/2.00-alpha
  plink2 --bfile ${merged_imputation} --indep-pairwise 1000kb 1 0.1 --out work/INTERVAL-covid
  plink2 --bfile ${merged_imputation} --make-bed --extract work/INTERVAL-covid.prune.in -keep work/INTERVAL-covid.samples2 --out work/INTERVAL-covid
  awk '{$1=$1 "_"_ $1; $2=$2 "_" $2;};1' work/INTERVAL-covid.fam > work/INTERVAL-covid.fam2
  plink --bed work/INTERVAL-covid.bed --bim work/INTERVAL-covid.bim --fam work/INTERVAL-covid.fam2 \
        --make-bed -keep work/INTERVAL-covid-X.id2 --out work/INTERVAL-covid-X
}

function aggregate()
{
# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[gwas software].[YYYYMMDD].txt.gz
# sex=M/MALE, F/FEMALE, ALL
  export snvResults=output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.20200526.txt.gz
  export geneResults=output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.gene.20200526.txt.gz
  (
    cut -d' ' -f1-3,5-9,11-14,17,21-26 output/INTERVAL-*.txt | head -1
    seq 22 | \
    tr ' ' '\n' | \
    parallel -j1 -C' ' '
      sed "1d" output/INTERVAL-{}.txt | \
      cut -d" " -f1-3,5-9,11-14,17,21-26
    '
    echo X | \
    parallel -j1 -C' ' '
      sed "1d" output/INTERVAL-{}.txt | \
      cut -d" " -f1-3,4-8,10-13,16,20-25 | \
      sed "s/X/23/"
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
# qqman
  export results=${snvResults}
  R --no-save -q <<\ \ END
    analysis <- "ANA5"
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

phenotype
