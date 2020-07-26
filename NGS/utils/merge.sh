#!/usr/bin/bash

export TMPDIR=$HPC_WORK/work
export prefix=1e-6
export tag=_nold

if [ ! -d ${prefix}/work ]; then mkdir ${prefix}/work; fi

for p in $(ls ${prefix}/*${tag}.p | sed 's|'"$prefix"'/||g;s|'"$tag"'.p||g'); do

echo $p
export p=${p}
(
  mergeBed -i ${prefix}/${p}_nold.p -d 1000000 -c 13 -o min | \
  awk -v OFS="\t" -v prot=${p} '
  {
    if(NR==1) print "Chrom", "Start", "End", "P", "prot"
    print $0, prot
  }'
) > ${prefix}/work/${p}.merged
(
  cut -f1-4,13 ${prefix}/${p}_nold.p | \
  bedtools intersect -a ${prefix}/work/${p}.merged -b - -wa -wb | \
  awk '$4==$10' | \
  cut -f1-6,8-10 | \
  awk -v OFS="\t" '
  {
    if(NR==1) print "Chrom", "Start", "End", "P", "prot", "MarkerName", "CHR", "POS", "SNP", "P_check"
    $5=$5 OFS $6 ":" $7
    gsub(/chr/,"",$6)
    print
  }'
) | uniq > ${prefix}/work/${p}.sentinels

done

(
  cat ${prefix}/work/*sentinels | head -1
  for p in $(ls ${prefix}/*${tag}.p | sed 's|'"$prefix"'/||g;s|'"$tag"'.p||g'); do awk 'NR>1' ${prefix}/work/${p}.sentinels; done
) > ${prefix}/NGS.tmp
R --no-save -q <<END
  prefix <- Sys.getenv("prefix")
  f <- paste(prefix,"NGS.tmp",sep="/")
  ngs <- read.table(f,header=TRUE,as.is=TRUE)
  dim(ngs)
  head(ngs)
  library(dplyr)
  t <- ngs %>% group_by(prot,Chrom,Start,End) %>% slice(which.min(P))
  t
  p <- table(ngs$P)[table(ngs$P)>1]
  write.table(t,file=paste(prefix,"NGS.merge",sep="/"),quote=FALSE,row.names=FALSE,sep='\t')
  print(p)
END
cut -f5 ${prefix}/NGS.merge | sed '1d' | sort | uniq > ${prefix}/NGS.merge.prot

R --no-save -q <<END
  prefix <- Sys.getenv("prefix")
  merge <- read.delim(paste(prefix,"NGS.merge",sep="/"),as.is=TRUE)
  m <- subset(merge,MarkerName!=".")
  write.table(m[,1:6],file=paste(prefix,"NGS.sentinels",sep="/"),row.names=FALSE,quote=FALSE)
END

cut -d' ' -f5 ${prefix}/NGS.sentinels | sed '1d' | sort | uniq > ${prefix}/NGS.sentinels.prot
gunzip -c hgTables.gz | awk 'length($1)<=5' | grep -f ${prefix}/NGS.sentinels.prot -
