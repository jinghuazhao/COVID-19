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

function fp()
{
  (
    cat ${INF}/work/METAL.hdr
  # 1. PLINK --clumping. The default.
  # 3. GCTA --cojo.  Change Chr, Pos to CHR, BP as in 1.
  # 2. R/gap/sentinels output. use  SNPID or $4 below
  # awk 'NR>1 {print $1,$4}' work/INF1_nold.sentinels | \
  # parallel -j4 -C' ' 'zgrep -w -H {2} METAL/{1}-1.tbl.gz'
    awk 'NR>1 {print $5,$6}' ACE2.merge | \
    parallel -j4 -C' ' 'zgrep -w -H {2} {1}-1.tbl.gz'
  ) | \
  sed 's/-1.tbl.gz//g' > ACE2.tbl
  cut -f3 ACE2.tbl | \
  awk 'NR>1' | \
  sort -k1,1 | \
  uniq | \
  join -a2 -e "NA" ${INF}/work/INTERVAL.rsid - -o2.1,1.2> ACE2.rsid
  (
    cat ${INF}/work/sumstats.hdr
    awk 'NR>1' ACE2.tbl | \
    cut -f1,3,13 | \
    awk '{split($1,a,":");print a[1],$2,$3}' | \
    parallel -j4 -C' ' '
      export direction=$(zgrep -w {2} {1}-1.tbl.gz | cut -f13)
      let j=1
      for i in $(grep "Input File" {1}-1.tbl.info | cut -d" " -f7)
      do
         export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
         if [ "$n" != "?" ]; then zgrep -H -w {2} $i; fi
         let j=$j+1
      done
  '
  ) | \
  sed 's|sumstats||g;s/.gz//g' > ACE2.all
  if [ -f ACE2.fp.log ]; then rm ACE2.fp.log; fi
  (
  R -q --no-save <<\ \ END
    require(gap)
    t <- read.delim("ACE2.tbl",as.is=TRUE)
    tbl <- within(t, {
      prot <- sapply(strsplit(Chromosome,":"),"[",1)
      Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
    })
    a <- read.table("ACE2.all",as.is=TRUE, header=TRUE)
    all <- within(a, {
      dir.study.prot <- sapply(strsplit(SNPID,":"),"[",1)
      p1 <- sapply(strsplit(SNPID,":"),"[",2)
      p2 <- sapply(strsplit(SNPID,":"),"[",3)
      MarkerName <- paste(p1,p2,sep=":")
      study <- sapply(strsplit(dir.study.prot,"/"),"[",2)
      study.prot <- sapply(strsplit(dir.study.prot,"/"),"[",3)
      substudy <- sapply(strsplit(study.prot,"[.]"),"[",1)
      pos <- unlist(lapply(gregexpr("[.]",study.prot),"[",1))
      prot <- substring(study.prot,pos+1)
    })
    droplist <- c("SNPID","dir.study.prot","p1","p2","pos","study.prot","substudy")
    all <- all[setdiff(names(all),droplist)]
    rsid <- read.table("ACE2.rsid",as.is=TRUE,col.names=c("MarkerName","rsid"))
    save(tbl,all,rsid,file="ACE2.rda",version=2)
    METAL_forestplot(tbl,all,rsid,"ACE2.fp.pdf",width=8.75,height=5)
  END
  ) > 2 > &1 | tee > ACE2.fp.log
  (
    echo  prot MarkerName Q df p I2 lower.I2 upper.I2
    grep I2 ACE2.fp.log | \
    awk '{gsub(/prot|=|MarkerName|Q|df|p|lower.I2|upper.I2|I2/,"");print}' | \
    sed '1d'
  ) | sed 's/  //1;s/   / /g' > ACE2.Q
}
