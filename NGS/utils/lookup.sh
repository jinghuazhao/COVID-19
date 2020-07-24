#!/usr/bin/bash

export prefix=1e-6

function check()
{
  export prot=$(grep -w $1 ${prefix}/NGS.merge | cut -f5 | sed 's/_invn//g')
  if [ "${prot}" == "" ]; then
    echo Empty
  else
    grep -H -w ${prot} $INF/doc/hgTables.tsv
    grep -H -w utils/INTERVAL-box.tsv
  fi
}

function Sun()
{
  awk 'NR>1{sub(/_/," ");print $5,$6,$7}' ${prefix}/NGS.merge | \
  parallel -C' ' '
    echo {1} {2} {3}
    grep -w {2} ${INF}/pQTL.Sun-B_pQTL_EUR_2017 | grep {3}
    grep -H -w {2} ${INF}/INTERVAL_box.tsv
  '
}

Sun

function Olink()
{
  export OLINK=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/jp549/olink-merged-output
  ls $OLINK | \
  grep -v check > olink.list
  join -j2 \
       <(awk 'NR>1{sub(/_/," ");print $5,$6,$7}' ${prefix}/NGS.merge | sort -k2,2) \
       <(ls $OLINK/*gz  | sed 's/___/ /g;s/_chr_merged.gz//g' | sort -k2,2) | \
  awk '{
     gsub(/chr/,"",$3);
     split($3,a,":");
     chr=a[1];
     pos=a[2];
     print $1,pos,$4,chr
  }' | \
  parallel -C' ' '
    echo {1} {2}
    zgrep -H -w {2} {3}___{1}_chr_merged.gz | \
    awk -vchr={4} "(\$3==chr)"
  '
}

Olink
