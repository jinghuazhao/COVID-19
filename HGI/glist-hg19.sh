#!/usr/bin/bash

export glist=${INF}/csd3/glist-hg19
(
  echo \#chrom Start End Gene
  awk '$1!="X" && $1!="Y" && $1!="XY"' ${glist} | sort -k1,1n -k2,2n | awk '{$1="chr" $1; print}'
  awk '$1=="X"' ${glist} | sort -k2,2n | awk '{$1="chr" $1; print}'
  awk '$1=="XY"' ${glist} | sort -k2,2n | awk '{$1="chr" $1; print}'
  awk '$1=="Y"' ${glist} | sort -k2,2n | awk '{$1="chr" $1; print}'
) > work/glist-hg19.bed
