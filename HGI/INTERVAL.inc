#!/usr/bin/bash

export dir=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
export merged_imputation=/home/jhz22/rds/post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/merged_imputation
export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export ev=$ref/annot_INT_50PCs_pcs.txt
export snvResults=output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.20200526.txt.gz
export geneResults=output/INTERVAL.Zhao.ANA5.1.ALL.EUR.97.237.SAIGE.gene.20200526.txt.gz

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
# Chromosome X -- NOTE that work/INTERVAL-X.samples has a single column
  awk '{print $1 "_" $1}' work/INTERVAL-X.samples | \
  bcftools view -S - ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -O v > work/INTERVAL-X.vcf
  export idno=$(awk '/POS/ && /QUAL/ {print NR} ' work/INTERVAL-X.vcf)
  awk -v idno=${idno} 'NR==idno{print} ' work/INTERVAL-X.vcf > work/INTERVAL.idline
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
  bgzip -cf > work/INTERVAL-X-ploidy.vcf.gz
  bcftools index -tf work/INTERVAL-X-ploidy.vcf.gz
# sbatch INTERVAL.sb
# Autosomes
  seq 22 | \
  parallel --env autosomes -C' ' '
    (
      echo alternate_ids rsid chromosome position allele1 allele2 rsid SNPID chromosome position allele1 allele2
      bgenix -g ${autosomes}/imputed/impute_{}_interval.bgen -list 2>&1 | \
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
    qctool -g ${autosomes}/imputed/impute_{}_interval.bgen -s ${autosomes}/imputed/interval.samples \
           -incl-samples work/INTERVAL.samples \
           -incl-rsids work/INTERVAL-{}.nodup -map-id-data work/INTERVAL-{}.map \
           -bgen-bits 8 \
           -og work/INTERVAL-{}.bgen -os work/INTERVAL-{}.samples
  '
    seq 22 | \
    parallel -j2 -C 'bgenix -g work/INTERVAL-{}.bgen -index -clobber'
  '
# HLA region for future imputation
  plink2 --bfile ${merged_imputation} --chr 6 --from-bp 25392021 --to-bp 33392022 --make-bed --out work/INTERVAL-HLA
}

function aggregate()
{
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
}

# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# [dataset].[last name].[analysis_name].[freeze_number].[sex].[ancestry].[gwas software].[YYYYMMDD].txt.gz
# sex=M/MALE, F/FEMALE, ALL

function qqman()
{
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

function annotate()
{
  export TMPDIR=/rds/user/jhz22/hpc-work/work
  echo $(seq 22 -1 1) X | \
  tr ' ' '\n' | \
  parallel -j2 --env autosomes -C' ' '
    cd work
    (
      awk "BEGIN{print \"##fileformat=VCFv4.0\"}"
      awk -vOFS="\t" "BEGIN{print \"#CHROM\",\"POS\",\"ID\",\"REF\",\"ALT\",\"QUAL\",\"FILTER\",\"INFO\"}"
      if [ "{}" == "X" ]; then
         export vcfgz=${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz; else export vcfgz=${autosomes}/{}.pbwt_reference_impute.vcf.gz
      fi
      bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/INFO\n" ${vcfgz}  | \
      awk -v OFS="\t" "NR>1{print \$1,\$2,\$1 \":\" \$2 \"_\" \$3 \"/\" \$4, \$3, \$4, \$5, \$6, \$7}"
    ) | \
    gzip -f > INTERVAL-{}.vepinput.gz
  # Split large chromosomes into two chunks (at most comparable to chromosome 7)
    if [ {} -le 6 ] && [ "{}" != "X" ]; then
      gunzip -c INTERVAL-{}.vepinput.gz | \
      split -l 5000000 --numeric-suffixes=1 --additional-suffix=.vepinput - INTERVAL-{}.
      gzip -f INTERVAL-{}.01.vepinput
      (
        gunzip -c INTERVAL-{}.vepinput.gz | \
        awk "NR<3{print}"
        cat INTERVAL-{}.02.vepinput
        rm INTERVAL-{}.02.vepinput
      ) | \
      gzip -f > INTERVAL-{}.02.vepinput.gz
    fi
  '
  export ANNOVAR=${HPC_WORK}/annovar
  export LOFTEE=${HPC_WORK}/loftee
  export POLYPHEN=$HPC_WORK/polyphen-2.2.2
  export VEP=${HPC_WORK}/ensembl-vep
  export APV=VEP
  echo $(seq 22 -1 7) $(seq 6 -1 1 | parallel 'echo {}.01 {}.02') X | \
  tr ' ' '\n' | \
  parallel -j2 --env autosomes --env ANNOVAR --env POLYPHEN --env LOFTEE --env VEP --env APV -C' ' '
    cd work
    if [ ${APV} == "VEP" ]; then
   # VEP
      vep --input_file INTERVAL-{}.vepinput.gz --output_file INTERVAL-{}.tsv --cache --dir_cache ${VEP}/.vep --dir_plugins ${LOFTEE} --offline \
          --pick --force_overwrite --species homo_sapiens --assembly GRCh37 --tab
      #   --plugin LoF,loftee_path:${LOFTEE},human_ancestor_fa:human_ancestor.fa.gz,conservation_file:phylocsf_gerp.sql.gz
    elif [ ${APV} == "ANNOVAR" ]; then
   # ANNOVAR
      gunzip -c INTERVAL-${s}.vepinput.gz | \
      awk -v OFS="\t" "NR>2{print \$1,\$2,\$2,\$4,\$5}" > INTERVAL-${s}.avinput
      annotate_variation.pl -buildver hg19 INTERVAL-${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile INTERVAL-${s}
      table_annovar.pl INTERVAL-${s}.avinput $ANNOVAR/test -buildver hg19 -out ${s} \
           -protocol ensGene,refGene,ccdsGene,wgEncodeGencodeBasicV19,cytoBand,exac03,avsnp147,dbnsfp30a,gwasCatalog \
           -operation g,g,g,g,r,f,f,f,r \
           -remove -nastring . -csvout -polish -xref $ANNOVA/example/gene_xref.txt
    else 
   # Polyphen-2
      gunzip -c INTERVAL-${s}.vepintput | \
      awk "NR>2{gsub(/_/,\" \",\$3);print \"chr\" \$3}" | \
      sort -k1,1 | \
      uniq > INTERVAL-${s}.pph.list
      mapsnps.pl -g hg19 -m -U -y INTERVAL-${s}.pph.input INTERVAL-${s}.pph.list 1>INTERVAL-${s}.pph.features 2>INTERVAL-${s}.log
      run_pph.pl INTERVAL-${s}.pph.input 1>$INTERVAL-{s}.pph.output 2>INTERVAL-${s}.pph.log
      run_weka.pl INTERVAL-${s}.pph.output >INTERVAL-${s}.pph.humdiv.output
      run_weka.pl -l $POLYPHEN/models/HumVar.UniRef100.NBd.f11.model INTERVAL-${s}.pph.output >INTERVAL-${s}.pph.humvar.output
    fi
    cd -
  '
}

function init()
# ~/COVID-19/HGI/rds-asb38-ceu-restricted/datasets/interval/covid19/version_history/20200520
{
  cd 06-05-2020/INTERVAL
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVALdata_06MAY2020.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_20200506.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
  sed '1d' INTERVALdata_P3_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
  cd -
}

# srun -A CARDIO-SL0-CPU -p cardio_intr --qos=cardio_intr -N1 -n1 -c8 -t 12:0:0 --pty bash -i
# sintr -A CARDIO-SL0-CPU -p cardio_intr --qos=cardio_intr
# sbatch -A CARDIO-SL0-CPU -p cardio_long --qos=cardio_long
# sbatch -A CARDIO-SL0-CPU -p cardio_short --qos=cardio_short

# Request an account
# https://docs.google.com/forms/d/1eAaf-4XNYkplBo5Appbf8LHl2KHJyks9R4t0E3h0jII/viewform?edit_requested=true

# gsutil cp ${snvResults} gs://covid19-hg-upload-UK-Blood-Donors-Cohort/
# gsutil cp ${geneResults} gs://covid19-hg-upload-UK-Blood-Donors-Cohort/

# Fill the form here,
# https://airtable.com/shrdJDwSHRZKXv45H
# Download the data from here,
# https://tinyurl.com/y97u49jz
