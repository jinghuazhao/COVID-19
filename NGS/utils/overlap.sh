#!/usr/bin/bash

# TNF IL6 CXCL8
for uniprot in P01375 P05231 P10145
do
  echo ${uniprot}
  for pval in 1e-5 1e-6 1e-7
  do
    echo ${pval}
    grep ${uniprot} ${pval}/NGS.sentinels | wc -l
    grep -v -e P01375 -e P05231 -e P10145 ${pval}/NGS.sentinels > NGS.${pval}
    R --no-save <pQTLtools.R > ${pval}.out
  done
done

