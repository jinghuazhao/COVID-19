# 8-4-2020 JHZ

export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export TMPDIR=/rds/user/jhz22/hpc-work/work

# extract refGene via ANNOVAR annotation
R --no-save <<END
  csv <- read.csv("ACE2.proxy.hg19_multianno.csv",as.is=TRUE)
  csv <- within(csv,{
      Chr <- paste0("chr",Chr); 
      Start <- Start-1; 
      chrpos <- paste0(Chr,":",End); 
      snpid <- paste0(chrpos,"_",Ref,"_",Alt)
  })
  vars <- c("Chr","Start","End","chrpos","snpid","avsnp147","Gene.refGene","Gene.ensGene")
  write.table(csv[vars],file="ACE2.proxy.hg19_multianno.gene",quote=FALSE,row.names=FALSE,sep="\t")
END
# double check
cut -f7 ACE2.proxy.hg19_multianno.gene | tr ';' '\n' | grep -w -f - $INF/csd3/glist-hg19 | grep -v AS1
(
  head -1 $INF/doc/hgTables.tsv | cut -f1-4,7
  cut -f7 ACE2.proxy.hg19_multianno.gene | \
  tr ';' '\n' | \
  grep -w -f - $INF/doc/hgTables.tsv | \
  cut -f1-4,7 | \
  awk -vOFS="\t" -vd=1000000 '{$2=$2-d;$3=$3+d};1'
) > st.tmp

(
  awk -vOFS="\t" 'BEGIN{print "Chrom","Start","End","Name","Gene","SNPID"}'
  cut -f1,8,9 ACE2.merge | \
  awk -vOFS="\t" '{if(NR>1){$2=$2=$2-1 OFS $2;print}}' | \
  bedtools intersect -a st.tmp -b - -wa -wb | \
  cut -f1-5,9
) > st.bed

cut -f5 st.bed | sed '1d' | \
parallel -j4 -C' ' '
  (
     echo -e "MarkerName\tP-value\tWeight"
     grep -w {} st.bed > st.{}
     read chrom start end gene prot < st.{}
     rm st.{}
     gunzip -c ACE2-1.tbl.gz | \
     awk -vOFS="\t" -vchr=$chrom -vstart=$start -vend=$end \
         "(\$1 == chr && \$2 >= start && \$2 <= end){split(\$3,a,\"_\");print a[1],10^\$12,\$18}" | \
     sort -k1,1 | \
     join -12 -21 ${INF}/work/snp_pos - | \
     awk -vOFS="\t" "{print \$2, \$3, \$4}"
  )  > {}.lz'

for g in $(cut -f5 st.bed | sed '1d')
do
  parallel -j1 -C' ' '
     grep -w {} st.bed > st.tmp
     read chrom start end gene prot < st.tmp
     rm -f ld_cache.db
     locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal {}.lz \
               --plotonly --chr $chrom --start $start --end $end --no-date --rundir .
     mv chr${chrom}_${start}-${end}.pdf {}.lz.pdf
     pdftopng -r 300 {}.lz.pdf {}
     mv {}-000001.png {}.lz-1.png
     mv {}-000002.png {}.lz-2.png
     cd -
  '
done
