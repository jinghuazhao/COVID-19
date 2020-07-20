#!/usr/bin/bash

export results=${HOME}/rds/results/public/gwas/covid19/hgi/covid19-hg-analysis/20200619/results
export harmonized=/home/jhz22/rds/results/public/gwas/covid19/hgi/covid19-hg-analysis/20200619/harmonized
export variants=L10RB_IFNAR2_variants.txt

export chr=($(sed '1d' ${variants} | cut -f1 | tr '\n' ' '))
export pos=($(sed '1d' ${variants} | cut -f2 | tr '\n' ' '))
export a1=($(sed '1d' ${variants} | cut -f3 | tr '\n' ' '))
export a2=($(sed '1d' ${variants} | cut -f4 | tr '\n' ' '))

if [ ! -d MR ]; then mkdir MR; fi

ls $harmonized/*gz | \
parallel -j10 -C' ' '
  export f=$(basename -s .gz {})
  (
    zcat {} | \
    head -1
    for i in $(seq 0 5)
    do
      awk -vOFS="\t" -vchr=${chr[$i]} -vpos=${pos[$i]} -va1=${a1[$i]} -va2=${a2[$i]} "
          NR==1 || (\$1==chr && \$2==pos && (\$3==a1 && \$4==a2 || \$3==a2 && \$4==a1))
      "
    done
  ) > MR/${f}.txt
'
