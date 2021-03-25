#!/usr/bin/bash

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

module load plink/2.00-alpha ceuadmin/stata

function extract_data()
{
  cut -f1 ${PCA_projection}/${PCA_loadings} | tail -n +2 > variants.extract
  (
    cat variants.extract
    awk '{split($1,a,":");print "chr"a[1]":"a[2]"_"a[4]"_"a[3]}' variants.extract
  ) > variants.extract2

  cp ${prefix}/${dir}/output/INTERVAL-*.bgen.bgi ${TMPDIR}
  seq 22 | \
  parallel -C' ' --env prefix --env dir '
    ln -sf ${prefix}/${dir}/output/INTERVAL-{}.bgen ${TMPDIR}/INTERVAL-{}.bgen
    python update_bgi.py --bgi ${TMPDIR}/INTERVAL-{}.bgen.bgi
    bgenix -g ${TMPDIR}/INTERVAL-{}.bgen -incl-rsids variants.extract2 > ${prefix}/work/INTERVAL-{}.bgen
  '
}

function twist()
{
  cat-bgen -g $(echo ${prefix}/work/INTERVAL-{1..22}.bgen) -og ${prefix}/work/INTERVAL.bgen -clobber
  bgenix -g ${prefix}/work/INTERVAL.bgen -index -clobber

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
}

function project_pc()
{
#!/bin/bash

  set -eu
################################################################################
# Please fill in the below variables
################################################################################
# Metadata
  STUDY_NAME="INTERVAL"
  ANALYST_LAST_NAME="ZHAO"
  DATE="$(date +'%Y%m%d')"
  OUTNAME="${prefix}/work/${STUDY_NAME}.${ANALYST_LAST_NAME}.${DATE}"
################################################################################
# Location of downloaded input files
  PCA_LOADINGS="${PCA_projection}/${PCA_loadings}"
  PCA_AF="${PCA_projection}/${PCA_af}"
################################################################################
# Location of imputed genotype files
# [Recommended]
# PLINK 2 binary format: a prefix (with directories) of .pgen/.pvar/.psam files
  PFILE="${prefix}/work/INTERVAL"
# [Acceptable]
# PLINK 1 binary format: a prefix of .bed/.bim/.fam files
  BFILE=""
################################################################################

  function error_exit() {
    echo "${1:-"Unknown Error"}" 1>&2
    exit 1
  }

# Input checks
  if [[ -z "${STUDY_NAME}" ]]; then
    error_exit "Please specify \$STUDY_NAME."
  fi

  if [[ -z "${ANALYST_LAST_NAME}" ]]; then
    error_exit "Please specify \$ANALYST_LAST_NAME."
  fi

  if [[ -z "${PCA_LOADINGS}" ]]; then
    error_exit "Please specify \$PCA_LOADINGS."
  fi

  if [[ -z "${PCA_AF}" ]]; then
    error_exit "Please specify \$PCA_AF."
  fi

  if [[ -n "${PFILE}" ]]; then
    input_command="--pfile ${PFILE}"
  elif [[ -n "${BFILE}" ]]; then
    input_command="--bfile ${BFILE}"
  else
    error_exit "Either \$PFILE or \$BFILE should be specified"
  fi

# Run plink2 --score
  plink2 \
    ${input_command} \
    --score ${PCA_LOADINGS} \
            variance-standardize \
            cols=-scoreavgs,+scoresums \
            list-variants \
            header-read \
    --score-col-nums 3-12 \
    --read-freq ${PCA_AF} \
    --out ${OUTNAME}

# The score file does not have FID (=0)
  awk '
  {
    if (NR==1) $1="#FID IID"
    else $1=0" "$1
    print
  }' ${OUTNAME}.sscore > ${prefix}/work/snpid.sscore
  ln -sf ${OUTNAME}.sscore.vars ${prefix}/work/snpid.sscore.vars

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

  stata -b do ethnic.do

  Rscript ${PCA_projection}/plot_projected_pc.R \
        --sscore ${prefix}/work/snpid.sscore \
        --phenotype-file ${prefix}/work/snpid.pheno \
        --phenotype-col SARS_CoV \
        --covariate-file ${prefix}/work/snpid.covars \
        --pc-prefix PC \
        --pc-num 20 \
        --ancestry-file ethnic.txt \
        --ancestry-col ethnic \
        --study ${STUDY_NAME} \
        --out ${OUTNAME}
}

twist;project_pc
