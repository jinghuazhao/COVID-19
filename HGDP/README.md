# Mapping data to reference population

## GIANT reference (dated 8/8/2014)

* refdata.* (N=437, # of variants=79,399)
* MDSplot.R

and examples are in `EPIC-Omics.sh` and `HGI.sh` produces `HGI.C1-C2.png`.

## 1000Genomes phase 3

* keep.dat contains a list of RSid's derived from HGI variants (# of variants=117,220).
* keep.bed/bim/fam contains genotypes for N=2,504 samples (too big to upload here).
* keep.sh is the pipeline to produce `INTERVAL.C1-C2.png`.

The reference data were downloaded from files distributed with VEGAS2 

[https://vegas2.qimrberghofer.edu.au/](https://vegas2.qimrberghofer.edu.au/)

and extracted from via the following script,

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
