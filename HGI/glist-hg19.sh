#!/usr/bin/bash

export glist=${INF}/csd3/glist-hg19
(
  awk '$1!="X" && $1!="Y" && $1!="XY"' ${glist} | sort -k1,1n -k2,2n | awk '{if ($1<10) $1="0"$1};1'
  awk '$1=="X"' ${glist} | sort -k2,2n
  awk '$1=="XY"' ${glist} | sort -k2,2n
  awk '$1=="Y"' ${glist} | sort -k2,2n
) | \
tr ' ' '\t' > work/glist-hg19.bed

echo $(seq 22) X | \
tr ' ' '\n' | \
parallel -C' ' '
  qctool -g work/INTERVAL-{}.bgen -annotate-bed4 work/glist-hg19.bed -osnp work/INTERVAL-{}.annotate
'

echo $(seq 22) X | \
tr ' ' '\n' | \
parallel -C' ' '
   awk "NR>9 && \$8!=\"NA\" && \$1!=\".\" && \$1!=\"#\"{print \$1}" work/INTERVAL-{}.annotate > work/INTERVAL-{}.incl
   export list=($(awk "NR>8 && \$8!=\"NA\"" work/INTERVAL-{}.annotate | cut -f8 | sort | uniq))
   (
     for g in ${list[@]}
     do
        awk -v g=${g} "\$8==g" work/INTERVAL-{}.annotate | \
        awk -vOFS="\t" "\$1!=\".\" {printf OFS \$1}" | \
        awk -v g=${g} -v OFS="\t" "{print g,\$0}"
     done
   ) > work/INTERVAL-{}.gene
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
