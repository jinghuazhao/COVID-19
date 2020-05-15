#!/usr/bin/bash

export glist=${INF}/csd3/glist-hg19
(
  awk '$1!="X" && $1!="Y" && $1!="XY"' ${glist} | sort -k1,1n -k2,2n
  awk '$1=="X"' ${glist} | sort -k2,2n
  awk '$1=="XY"' ${glist} | sort -k2,2n
  awk '$1=="Y"' ${glist} | sort -k2,2n
) | \
tr ' ' '\t' > work/glist-hg19.bed

seq 22 | \
parallel -C' ' '
  qctool -g work/INTERVAL-{}.bgen -annotate-bed4 work/glist-hg19.bed -osnp work/INTERVAL-{}.annotate
'

# SNP information
(
  cut -f1,7,8,15,18,19 $ref/impute_*_interval.snpstats | \
  head -1
  seq 22 | \
  parallel -j1 --env ref -C' ' '
    sed "1d" $ref/impute_{}_interval.snpstats | \
    cut -f1,3,4,7,8,15,18,19
  '
) | gzip -f > work/INTERVAL.snpstats.gz

# SNP information
(
  cut -f1,7,8,15,18,19 $ref/impute_*_interval.snpstats | \
  head -1
  seq 22 | \
  parallel -j1 --env ref -C' ' '
    sed "1d" $ref/impute_{}_interval.snpstats | \
    cut -f1,3,4,7,8,15,18,19
  '
) | gzip -f > work/INTERVAL.snpstats.gz

