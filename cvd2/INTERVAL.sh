# 8-4-2020 JHZ

# On Tryggve call R.3.5.3

R --no-save -q <<END
  protein <- "ACE2";
  print(protein);
  gz <- gzfile(paste0("sumstats/INTERVAL.",protein,".gz"));
  require(qqman);
  tbl <- read.delim(gz,as.is=TRUE);
  tbl <- within(tbl,{
     SNP <- SNPID
     BP <- PVAL
  })
  tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
  qq <- "INTERVAL_qq.png";
  png(qq,width=12,height=10,units="in",pointsize=4,res=300)
  qq(with(tbl,P))
  dev.off()
  manhattan <- "INTERVAL_manhattan.png";
  png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
  manhattan(tbl,main=protein,genomewideline=-log10(8.210181e-12),suggestiveline=FALSE,ylim=c(0,25));
  dev.off();
END

cut -f5 st.bed | sed '1d' | \
parallel -j4 -C' ' '
  (
     echo -e "MarkerName\tP-value\tWeight"
     grep -w {} st.bed | sed 's/chr//g' > st.{}
     read chrom start end prot gene snpid < st.{}
     rm st.{}
     gunzip -c sumstats/INTERVAL.ACE2.gz | \
     awk -vOFS="\t" -vchr=$chrom -vstart=$start -vend=$end \
         "(\$2 == chr && \$3 >= start && \$3 <= end){split(\$1,a,\"_\");print a[1],10^\$11,\$5}" | \
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
