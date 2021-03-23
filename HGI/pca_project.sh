#!/usr/bash/bin

export TMPDIR=${HPC_WORK}/work
# genotype
export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/
# location of PCs
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
# HGI working directory
export prefix=~/rds/rds-asb38-ceu-restricted/projects/covid/HGI
export dir=20201201-male-ANA_C2_V2

export PCA_projection=pca_projection
export PCA_loadings=hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.tsv
export PCA_af=hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh38.plink.afreq
export sscore=hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz

function snpid_rsid()
{
  cut -f1 ${PCA_projection}/${PCA_loadings} | tail -n +2 > variants.extract
  (
    cat variants.extract
    awk '{split($1,a,":");print "chr"a[1]":"a[2]"_"a[4]"_"a[3]}' variants.extract
  ) > variants.extract2
  awk '{split($1,a,":");if(a[3]<a[4]) print "chr"a[1]":"a[2]"_"a[3]"_"a[4];else print "chr"a[1]":"a[2]"_"a[4]"_"a[3]}' variants.extract > snpid.extract
  join <(sort -k1,1 snpid.extract) <(gunzip -c ${HPC_WORK}/SUMSTATS/snp150.snpid_rsid.gz) > snpid-rsid.extract
  cut -d' ' -f2 snpid-rsid.extract > rsid.extract
}

function bgen_only_as_documented()
{
  seq 22 | \
  parallel -C' ' --env prefix --env dir ' bgenix -g ${prefix}/${dir}/output/INTERVAL-{}.bgen -incl-rsids rsid.extract > ${prefix}/work/INTERVAL-{}.bgen'
  cat-bgen -g $(echo ${prefix}/work/INTERVAL-{1..22}.bgen) -og ${prefix}/work/INTERVAL.bgen -clobber
  bgenix -g ${prefix}/work/INTERVAL.bgen -index
  qctool -g ${prefix}/work/INTERVAL.bgen -og ${prefix}/work/INTERVAL.vcf.gz
}

function bgen_sample()
{
  seq 22 | \
  parallel -C' ' --env prefix --env dir '
  qctool -g ${prefix}/${dir}/output/INTERVAL-{}.bgen -incl-rsids rsid.extract \
         -s ${prefix}/${dir}/output/INTERVAL-{}.samples \
         -og ${prefix}/work/INTERVAL-{}.bgen -os ${prefix}/work/INTERVAL.samples
  bgenix -g ${prefix}/work/INTERVAL-{}.bgen -index -clobber
  '
}

function map()
{
  seq 22 | \
  parallel --env prefix -C' ' '
    (
      echo alternate_ids rsid chromosome position allele1 allele2 rsid SNPID chromosome position allele1 allele2
      bgenix -g ${prefix}/work/INTERVAL-{}.bgen -list 2>&1 | \
      sed "1,9d" | \
      awk "
      {
        CHR=\$3+0
        POS=\$4
        a1=toupper(\$6)
        a2=toupper(\$7)
        snpid=CHR \":\" POS \":\" a1 \":\" a2
        if (NF==7) print \$1, \$2, \$3, POS, \$6, \$7, \$1, snpid, CHR, POS, a1, a2
      }"
    ) | \
    awk "a[\$2]++==0" > ${prefix}/work/INTERVAL-{}.map
    cut -d" " -f2 ${prefix}/work/INTERVAL-{}.map > ${prefix}/work/INTERVAL-{}.nodup
  '
}

function snpid()
{
  seq 22 | \
  parallel --env prefix -C' ' '
    qctool -g ${prefix}/work/INTERVAL-{}.bgen -s ${prefix}/work/INTERVAL.samples \
           -incl-rsids ${prefix}/work/INTERVAL-{}.nodup -map-id-data ${prefix}/work/INTERVAL-{}.map \
           -bgen-bits 8 \
           -og ${prefix}/work/snpid-{}.bgen -os ${prefix}/work/snpid-{}.samples
    bgenix -g ${prefix}/work/snpid-{}.bgen -index -clobber
  '
  cat-bgen -g $(echo ${prefix}/work/snpid-{1..22}.bgen) -og ${prefix}/work/snpid.bgen -clobber
  bgenix -g ${prefix}/work/snpid.bgen -index
}

function install_R_packages()
{
  Rscript -e "install.packages(c("data.table", "hexbin", "optparse", "patchwork", "R.utils", "tidyverse"))"
}

snpid_rsid;bgen_sample;map;snpid

module load plink/2.00-alpha
plink2 --bgen ${prefix}/work/snpid.bgen ref-first --make-pfile --out ${prefix}/work/snpid

# The "chr" prefix in the provided file is unnecessary
sed 's/chr//g' ${PCA_projection}/${PCA_af} > ${prefix}/work/PCA_af.tsv
plink2 --pfile ${prefix}/work/snpid \
       --score ${PCA_projection}/${PCA_loadings} center cols=-scoreavgs,+scoresums list-variants header-read \
       --score-col-nums 3-22 --read-freq ${prefix}/work/PCA_af.tsv --out ${prefix}/work/init
# The score file does not have FID (=0)
awk '
{
  if (NR==1) $1="#FID IID"
  else $1=0" "$1
  print
}' ${prefix}/work/init.sscore > ${prefix}/work/snpid.sscore
ln -sf ${prefix}/work/init.sscore.vars ${prefix}/work/snpid.sscore.vars

# Our PCs are PC_# rather than PC#
# The phenotype and covariate file missed FID (=0) column
awk '
{
  if (NR==1)
  {
    $1="FID I"$1
    gsub(/PC_/,"PC",$0)
  }
  else $1=0" " $1
};1' ${dir}/work/INTERVAL-covid.txt | \
tr ' ' '\t' > ${prefix}/work/snpid.txt
cut -f1-3 ${prefix}/work/snpid.txt > ${prefix}/work/snpid.pheno
cut -f3 --complement ${prefix}/work/snpid.txt > ${prefix}/work/snpid.covars
Rscript ${PCA_projection}/plot_projected_pc.R \
  --sscore ${prefix}/work/snpid.sscore \
  --phenotype-file ${prefix}/work/snpid.pheno \
  --phenotype-col SARS_CoV \
  --covariate-file ${prefix}/work/snpid.covars \
  --pc-prefix PC \
  --pc-num 20 \
  --ancestry EUR \
  --study snpid \
  --out ${prefix}/work/snpid
