#!/usr/bin/bash

export harmonized=/home/jhz22/rds/results/public/gwas/covid19/hgi/covid19-hg-analysis/20200619/harmonized_n
export results=${HOME}/rds/results/public/gwas/covid19/hgi/covid19-hg-analysis/20200619/results
export variants=L10RB_IFNAR2_variants

function hg19()
{
  export chr=($(sed '1d' ${variants}.txt | cut -f2 | tr '\n' ' '))
  export pos=($(sed '1d' ${variants}.txt | cut -f3 | tr '\n' ' '))
  export a1=($(sed '1d' ${variants}.txt | cut -f4 | tr '\n' ' '))
  export a2=($(sed '1d' ${variants}.txt | cut -f5 | tr '\n' ' '))
}

function hg19Tohg38()
{
  awk 'NR>1{print "chr" $2,$3-1,$3,$2":"$3"_"toupper($4)"/"toupper($5)"-"$1}' OFS='\t' ${variants}.txt > ${variants}.bed
  liftOver ${variants}.bed ${HPC_WORK}/bin/hg19ToHg38.over.chain.gz ${variants}.hg38 ${variants}.unMapped
  export chr=($(awk '{sub(/chr/,"") $1;print $1}' $variants.hg38))
  export pos=($(awk '{print $3}' $variants.hg38))
  export a1=($(awk '{gsub("_|/|-"," ",$4);split($4,a);print a[2]}' $variants.hg38))
  export a2=($(awk '{gsub("_|/|-"," ",$4);split($4,a);print a[3]}' $variants.hg38))
}

hg19Tohg38

function viaParallel()
{
  if [ ! -d MR ]; then mkdir MR; fi
  ls $harmonized/*gz | \
  env_parallel --env chr --env pos --env a1 --env a2 -j10 -C' ' '
  export f=$(basename -s .gz {})
  echo ${f}
  (
    zcat {} | \
    head -1
    for i in $(seq 0 5)
    do
      zcat {} | \
      awk -vchr=${chr[$i]} -vpos=${pos[$i]} -va1=${a1[$i]} -va2=${a2[$i]} \
          "\$1==chr && \$2==pos && (\$3==a1 && \$4==a2 || \$3==a2 && \$4==a1)"
    done
  ) > MR/${f}.txt
  '
}

function viaLoop()
{
  export outdir=selected_n
  if [ ! -d ${outdir} ]; then mkdir ${outdir}; fi
  for f in $(ls $harmonized/*gz)
  do
  export r=$(basename -s .gz ${f})
  (
    zcat ${f} | \
    head -1
    for i in $(seq 0 5)
    do
      zcat ${f} | \
      awk -vchr=${chr[$i]} -vpos=${pos[$i]} -va1=${a1[$i]} -va2=${a2[$i]} '$1==chr && $2==pos && ($3==a1 && $4==a2 || $3==a2 && $4==a1)'
    done
  ) > ${outdir}/${r}.dat
  done
}

# viaLoop
# export loc=hgi/covid19-hg-analysis/20200619/L10RB_IFNAR2_variants
# pandoc ${loc}/README.md -o ${loc}/README.pdf
