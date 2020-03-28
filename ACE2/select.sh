#!/usr/bin/bash

export sumstats=2805-6_2-1.tbl.gz
function cmd()
{
  gunzip -c ${sumstats}
}

(
  cmd | head -1
  cmd | awk 'NR>1 && $12 < -3 && $1!="M" && $1!="X"' | sort -k1,1n -k2,2n
  cmd | awk 'NR>1 && $12 < -3 && $1=="X"'            | sort -k2,2n
) > 2805-6_2.001.txt
