# 6-4-2020 JHZ

export TMPDIR=$HPC_WORK/work
export tag=_nold

export p=ACE2
(
  mergeBed -i ${p}_nold.p -d 1000000 -c 13 -o min | \
  awk -v OFS="\t" -v prot=${p} '
  {
    if(NR==1) print "Chrom", "Start", "End", "log10P", "prot"
    print $0, prot
  }'
) > ${p}.merged
(
  cut -f1-4,13 ${p}_nold.p | \
  bedtools intersect -a ${p}.merged -b - -wa -wb | \
  awk '$4==$10' | \
  cut -f1-6,8-10 | \
  awk -v OFS="\t" '
  {
    if(NR==1) print "Chrom", "Start", "End", "log10P", "prot", "MarkerName", "CHR", "POS", "SNP", "log10P_check"
    $5=$5 OFS $6 ":" $7
    gsub(/chr/,"",$6)
    print
  }'
) | uniq > ${p}.tmp

(
  cat *tmp | head -1
  for p in $(ls *${tag}.p | sed 's|'"$tag"'.p||g'); do awk 'NR>1' ${p}.tmp; done
) > ACE2.merge

R --no-save -q <<END
  merge <- read.delim("ACE2.merge",as.is=TRUE)
  m <- subset(merge,SNP!=".")
  write.table(m[,1:6],file="ACE2.sentinels",row.names=FALSE,quote=FALSE)
END
