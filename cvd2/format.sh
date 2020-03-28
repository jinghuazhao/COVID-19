# 27-3-2020 JHZ

export SCALLOP=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop
export INTERVAL=${SCALLOP}/jp549/olink-merged-output/

gunzip -c ${INTERVAL}/INTERVAL_cvd2_ACE2___Q9BYF1_chr_merged.gz | \
awk -f ${INF}/tryggve/INTERVAL.awk | \
awk -f ${INF}/tryggve/order.awk | \
gzip -f > sumstats/INTERVAL.ACE2.gz

ls tryggve | \
sed 's/.ACE2//g;s/.txt.gz\*//g;s/\./ /g' | \
parallel -C' ' '
  gunzip -c tryggve/{1}.ACE2.{2}.{3}.txt.gz | \
  awk -v OFS="\t" "{
      if(NR>1)
      {
        CHR=\$3;POS=\$4;A1=toupper(\$7);A2=toupper(\$8);
        if(A1<A2) A1A2=\"_\" A1 \"_\" A2;else A1A2=\"_\" A2 \"_\" A1
        SNPID=\"chr\" CHR \":\" POS A1A2
        \$1=SNPID
      }
  };1" | gzip -f > sumstats/{1}.ACE2.gz
'
