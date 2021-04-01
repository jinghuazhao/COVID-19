# Mapping data to reference population

## GIANT reference (N=437, 79,399 variants, dated 8/8/2014)

* refdata.*
* EPIC-Omics.sh
* HGI.sh
* MDSplot.R

## 1000Genomes (N=2,504, 117,220 variants)

Data were downloaded from files distributed with VEGAS2 and extracted from via the following command

```bash
#!/usr/bin/bash

for pop in AFR AMR EAS EUR SAS
do
  echo ${pop}
  plink --bfile g1000p3_${pop} --extract keep.dat --make-bed --out ${pop}
  awk -v pop=${pop} '{$1=pop};1' ${pop}.fam > ${pop}-id.fam
done

echo AFR AMR EAS EUR SAS | tr ' ' '\n' | awk '{print $1".bed", $1".bim", $1"-id.fam"}' > keep.list
plink --merge-list keep.list --out keep
```

where `keep.dat` contains variants derived from HGI variants.

A call is made with `keep.sh`.
