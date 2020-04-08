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

R --no-save -q <<END
  protein <- "ACE2";
  print(protein);
  gz <- gzfile(paste0(protein,"-1.tbl.gz"));
  require(qqman);
  tbl <- read.delim(gz,as.is=TRUE);
  tbl <- within(tbl,{
     SNP <- MarkerName
     CHR <- as.numeric(Chromosome)
     BP <- Position
  })
  tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
  qq <- paste0(protein,"_qq.png");
  png(qq,width=12,height=10,units="in",pointsize=4,res=300)
  qq(with(tbl,P))
  dev.off()
  manhattan <- paste0(protein,"_manhattan.png");
  png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
  manhattan(tbl,main=protein,genomewideline=-log10(8.210181e-12),suggestiveline=FALSE,ylim=c(0,25));
  dev.off();

# Q9BYF1 -- ACE2
END

cut -f5 st.bed | sed '1d' | \
parallel -j4 -C' ' '
  (
     echo -e "MarkerName\tP-value\tWeight"
     grep -w {} st.bed | sed 's/chr//g' > st.{}
     read chrom start end prot gene snpid < st.{}
     rm st.{}
     gunzip -c ACE2-1.tbl.gz | \
     awk -vOFS="\t" -vchr=$chrom -vstart=$start -vend=$end \
         "(\$1 == chr && \$2 >= start && \$2 <= end){split(\$3,a,\"_\");print a[1],10^\$12,\$18}" | \
     sort -k1,1 | \
     join -12 -21 ${INF}/work/snp_pos - | \
     awk -vOFS="\t" "{print \$2, \$3, \$4}"
  )  > {}.lz'

cut -f5 st.bed | sed '1d' | \
  parallel -j1 -C' ' '
     grep -w {} st.bed | sed 's/chr//g' > st.tmp
     read chrom start end prot gene snpid < st.tmp
     rm -f ld_cache.db
     locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal {}.lz \
               --plotonly --chr $chrom --start $start --end $end --no-date --rundir .
     mv chr${chrom}_${start}-${end}.pdf {}.lz.pdf
     pdftopng -r 300 {}.lz.pdf {}
     mv {}-000001.png {}.lz-1.png
     mv {}-000002.png {}.lz-2.png
     cd -
  '
