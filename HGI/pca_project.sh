#!/usr/bash/bin

export TMPDIR=${HPC_WORK}/work
# genotype
export autosomes=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/
# location of PCs
export ref=~/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
# HGI working directory
export prefix=~/rds/rds-asb38-ceu-restricted/projects/covid/HGI
export dir=20201201-ANA_C2_V2

export PCA_projection=pca_projection
export PCA_loadings=hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.tsv
export PCA_af=hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.afreq
export sscore=hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz

cut -f1 ${PCA_projection}/${PCA_loadings} | tail -n +2 > variants.extract
(
  cat variants.extract
  awk '{split($1,a,":");print "chr"a[1]":"a[2]"_"a[4]"_"a[3]}' variants.extract
) > variants.extract2

function extract_data()
{
  cp ${prefix}/${dir}/output/INTERVAL-*.bgen.bgi ${TMPDIR}
  seq 22 | \
  parallel -C' ' --env prefix --env dir '
    ln -sf ${prefix}/${dir}/output/INTERVAL-{}.bgen ${TMPDIR}/INTERVAL-{}.bgen
    python update_bgi.py --bgi ${TMPDIR}/INTERVAL-{}.bgen.bgi
    bgenix -g ${TMPDIR}/INTERVAL-{}.bgen -incl-rsids variants.extract2 > ${prefix}/work/INTERVAL-{}.bgen
  '
  cat-bgen -g $(echo ${prefix}/work/INTERVAL-{1..22}.bgen) -og ${prefix}/work/INTERVAL.bgen -clobber
  bgenix -g ${prefix}/work/INTERVAL.bgen -index -clobber
}

module load plink/2.00-alpha
plink2 --bgen ${prefix}/work/INTERVAL.bgen ref-first --make-pfile --out INTERVAL

export csvfile=INTERVAL.csv
python update_bgi.py --bgi INTERVAL.bgen.bgi
(
  head -1 INTERVAL.pvar
  paste <(sed '1d' INTERVAL.pvar | cut -f1,2) \
        <(sed '1d' INTERVAL.csv | cut -d, -f3) | \
  paste - <(sed '1d' INTERVAL.pvar | cut -f4,5)
) > ${prefix}/work/INTERVAL.pvar
cp INTEVAL.p?? ${prefix}/work

plink2 --pfile ${prefix}/work/INTERVAL \
       --score ${PCA_projection}/${PCA_loadings} center cols=-scoreavgs,+scoresums list-variants header-read \
       --score-col-nums 3-22 --read-freq ${PCA_projection}/${PCA_af} --out ${prefix}/work/init
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
  --study INTERVAL \
  --out ${prefix}/work/INTERVAL
