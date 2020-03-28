#!/usr/bin/bash

function cmd()
{
  gunzip -c ${sumstats}
}

function cvd2()
{
  export sumstats=2019-04-25_meta_ACE2.txt.gz
  (
    cmd | head -1
    cmd | awk 'NR>1 && $10 < 0.001 && $18!="M" && $18!="X"' | sort -k18,18n -k19,19n
    cmd | awk 'NR>1 && $10 < 0.001 && $18=="M"'             | sort -k19,19n
    cmd | awk 'NR>1 && $10 < 0.001 && $18=="X"'             | sort -k19,19n
  ) > 2019-04-25_meta_ACE2_0.001.txt
  R --no-save <<\ \ END
    ACE2 <- read.delim("2019-04-25_meta_ACE2_0.001.txt",as.is=TRUE)
    ord <- with(ACE2,order(PVAL))
    ACE2[ord,]
    write.table(ACE2[ord,],file="2019-04-25_meta_ACE2.001",quote=FALSE,row.names=FALSE,sep="\t")
  END
}

function ACE2()
{
  export sumstats=ACE2-1.tbl.gz
  (
    cmd | head -1
    cmd | awk 'NR>1 && 10^$12 < 0.001 && $1!="M" && $1!="X"' | sort -k1,1n -k2,2n
    cmd | awk 'NR>1 && 10^$12 < 0.001 && $1=="M"'            | sort -k2,2n
    cmd | awk 'NR>1 && 10^$12 < 0.001 && $1=="X"'            | sort -k2,2n
  ) > ACE2_0.001.txt
  R --no-save <<\ \ END
    ACE2 <- read.delim("ACE2_0.001.txt",as.is=TRUE)
    ord <- with(ACE2,order(log.P.))
    write.table(ACE2[ord,],file="ACE2.001",quote=FALSE,row.names=FALSE,sep="\t")
  END
}

$1
