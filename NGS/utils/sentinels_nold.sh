#!/usr/bin/bash

export tag=_nold
export pval=1e-6
export src=plink2

if [ ! -d ${pval} ]; then mkdir ${pval}; fi

function pgz()
# 1. extract all significant SNPs
{
  ls ${src}/*.gz | \
  sed 's|'"$src"'/||g;s/.gz//g' | \
  parallel -j10 -C' ' '
  (
    zcat ${src}/{}.gz | head -1
    zcat ${src}/{}.gz | awk -v p=${pval} "NR>1 && \$12 <= p" | sort -k1,1n -k2,2n
  ) | gzip -f > ${pval}/{}.p.gz'
}

function _HLA()
# 2. handling HLA
{
  for p in $(ls ${src}/*.gz | sed 's|'"$src"'/||g;s/.gz//g')
  do
    (
      zcat ${src}/${p}.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      (
        zcat ${pval}/${p}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk 'NR>1 && !($1 == "6" && $3 >= 25392021 && $3 < 33392022)'
        zcat ${pval}/${p}.p.gz | \
        awk -vOFS="\t" '{start=$2-1;$2=start "\t" $2};1' | \
        awk '$1 == "6" && $3 >= 25392021 && $3 < 33392022' | \
        sort -k13,13g | \
        awk 'NR==1'
      ) | \
      sort -k1,1n -k2,2n -k3,3n | \
      awk -v OFS="\t" '{$1="chr" $1};1'
    ) > ${pval}/${p}${tag}.p
    export lines=$(wc -l ${pval}/${p}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}${tag} with $lines lines
      rm ${pval}/${p}${tag}.p
    fi
  done
}

for cmd in pgz _HLA; do $cmd; done
