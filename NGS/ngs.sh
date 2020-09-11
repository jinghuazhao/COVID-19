#!/usr/bin/bash

export opanels=(cvd2 cvd3 inf neu)
export qc_opanels=(qc_cvd2 qc_cvd3 qc_inf neu_qc)
export panels=(CARDIOMETABOLIC CARDIOMETABOLIC INFLAMMATION NEUROLOGY)

module load ceuadmin/stata

function init()
{
  export src=tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv
  sed '1d' ${src} | \
  cut -d';' -f4,5,7 | \
  sort | \
  uniq | \
  tr ';' ' ' > work/ids.dat
  R --no-save -q <<\ \ END
    ids <- read.table("work/ids.dat",col.names=c("UniProt","Prot","Panel"),as.is=TRUE)
    save(ids,file="work/ids.rda")
  END
  grep -v -e CSA ${src} | \
  sed 's/NaN/NA/g' > work/NPX.csv
  ln -sf $HOME/rds/post_qc_data/interval/phenotype/olink_proteomics
  stata <<\ \ END
    insheet using INTERVALdata_19JUN2020.csv, case clear
    save work/data, replace
    insheet using INTERVAL_OmicsMap_20200619.csv, case clear
    format Affymetrix_gwasQC_bl %15.0g
    format Olink_* %15.0g
    save work/omics, replace
    insheet using Olink_NGS.csv, case clear
    rename OLINK_NGS SampleID
    merge 1:1 identifier using work/data, gen(data)
    merge 1:1 identifier using work/omics,gen(data_omics)
    keep if SampleID!=""
    save work/dataomics, replace
  END
  for i in $(seq 0 3)
  do
    export opanel=${opanels[$i]}
    export qc_opanel=${qc_opanels[$i]}
    head -1 olink_proteomics/qc/olink_${qc_opanel}.csv | \
    tr ',' '\n' | \
    sed '1,2d' | \
    awk '{print toupper($1)}' | \
    grep -f - work/NPX.csv | \
    cut -d';' -f4 | \
    sort | \
    uniq > work/NPX-${opanel}.id
    stata <<\ \ \ \ END
       local opanel: env opanel
       insheet UniProt using work/NPX-`opanel'.id, case
       sort UniProt
       save work/NPX-`opanel', replace
    END
  done
}

init

function ngs()
{
for i in $(seq 0 3)
do
  export opanel=${opanels[$i]}
  export qc_opanel=${qc_opanels[$i]}
  export panel=${panels[$i]}
  for opt in raw LOD QC col1
  do
      export opt=${opt}
      echo ${panel} - ${opanel} - ${opt}
      if [ ${opt} == "raw" ]; then
         cat work/NPX.csv > work/NPX-${opt}.csv
      elif [ ${opt} == "LOD" ]; then
         awk -vFS=';' -vOFS=';' '$11 > $12 {$12="NA"};1' work/NPX.csv > work/NPX-${opt}.csv
      elif [ ${opt} == "QC" ]; then
         awk -vFS=';' -vOFS=';' '$10 =="WARN" {$12="NA"};1' work/NPX.csv > work/NPX-${opt}.csv
      else
         awk -vFS=';' -vOFS=';' '/01;/{$12="NA"};1' work/NPX.csv > work/NPX-${opt}.csv
      fi
      stata -b do ngs.do
      R --no-save -q < ngs.R
  done
done
}

ngs

function correlogram()
{
R --no-save -q <<END
# IL4, IL12B (INF), CD28 (ONC)
  correlogram <- function(opt)
  {
    for (i in 1:4)
    {
      opanel <- opanels[i]
      panel <- panels[i]
      rt <- paste0(panel,"-",opanel,"-", opt)
      f <- paste0("work/",rt,".dat")
      print(f)
      t <- subset(read.table(f,as.is=TRUE,header=TRUE),!is.na(r))
      if(nrow(t)<=1) break
      print(t)
      N <- nrow(t)
      d <- density(with(t,r))
      if (i%%2==1) {
        plot(d,ann=FALSE,axes=FALSE,main="",ylim=c(0,2))
        axis(2)
      } else plot(d,ann=FALSE,axes=FALSE,ylim=c(0,2))
      if(i>2) axis(1)
      title(rt)
    }
  }
  opanels <- c("cvd2","cvd3","inf","neu")
  qc_opanels <- c("qc_cvd2","qc_cvd3","qc_inf","neu_qc")
  panels <- c("CARDIOMETABOLIC","CARDIOMETABOLIC","INFLAMMATION","NEUROLOGY")
  pdf("work/correlogram.pdf")
  par(mfrow=c(2,2))
  for (opt in c("raw","LOD","QC","col1")) correlogram(opt)
  dev.off()
END
}

# NPX.csv
# OID20426;P05112;IL4;INFLAMMATION
# OID20666;P29460;IL12B;INFLAMMATION
# OID21265;P10747;CD28;ONCOLOCY

dos2unix Olink_NGS.csv
dos2unix work/NPX-QC.csv
dos2unix INTERVAL_OmicsMap_20200619.csv

function check()
{
  export f=olink_proteomics/qc/olink_qc_inf.csv
  export col=$(awk 'BEGIN {RS = ","}/^p29460$/{print $0"\t"NR}' ${f} | cut -f2)

  join  -12 -21 \
       <(awk -v FS=',' 'NR>1{print $1 " " $2}' Olink_NGS.csv | sort -k2,2) \
       <(awk -v FS=';' '/OID20666/{gsub(/;/," ");print}' work/NPX-raw.csv | sort -k1,1) | \
  sort -k2,2 | \
  join -22 <(cut -d, -f1,4,9 INTERVAL_OmicsMap_20200619.csv | sed '1d;s/,/\t/g' | sort -k1,1) - | \
  sort -k3,3 | \
  join -23 <(cut -d, -f1,${col} ${f} | sed '1d;s/,/ /g;s/\r//g' | sort -k1,1) -
}

function IL12B()
{
(
  echo Olink_inf_QC_24m qc identifier Affymetrix_gwasQC_bl SampleID Index OlinkID UniProt Assay MissingFreq Panel Panel_Version PlateID QC_Warning LOD NPX
  check
) > work/IL12B.txt

R --no-save -q <<END
  check <- read.table("work/IL12B.txt",as.is=TRUE,header=TRUE)
  method <- "spearman"
# raw
  raw <- subset(check,!is.na(NPX))
  with(raw,cor(qc,NPX,method=method,use="everything"))
# LOD
  lod <- subset(raw,NPX>=LOD)
  with(lod,cor(qc,NPX,method=method,use="everything"))
# QC
  qc <- subset(raw,QC_Warning=="PASS")
  with(qc,cor(qc,NPX,method=method,use="everything"))
# col1
  col1 <- subset(raw,substr(SampleID,11,12)!="01")
  with(col1,cor(qc,NPX,method=method,use="everything"))
END
}

function bgen()
{
  export dat=olink_ngs_proteomics/gwasqc/olink_ngs_gwasqc.txt
  export ids=INTERVAL_OmicsMap_20200619.csv
  export sam=olink_ngs_proteomics/gwasqc/sample_info.txt
  awk 'NF<1473{print $1,NF}' $dat | \
  parallel --env dat -C' ' '
    awk -vFS="\t" -vid={1} "{if(\$1==id)for(i=1;i<=NF;++i) if(\$i==\"\") print \$1,NR,i}" ${dat}' | \
    parallel -C' ' 'head -1 $dat | awk -vid={1} -vline={2} -vcol={3} "{print id,line,col,\$col}"'
  # -f7-10 ${ids} ==> 24m, line, column, ID:
  # 110003567640 177 782 NEUROLOGY_O95994
  # 110004126595 186 595 INFLAMMATION_Q13291
  # 110014074648 357 441 INFLAMMATION_P05112
  (
    head -1 ${dat} | \
    awk -v OFS='\t' '{$1="FID\tIID";print}'
    join -12 -21 <(cut -d, -f4,7 ${ids} | tr ',' '\t'| sort -k2,2) \
                 <(sed '1d' ${dat} | \
                   awk -vFS='\t' -vOFS='\t' '{if(NR==176) $782=-9; else if(NR==185) $595=-9; else if(NR==356) $441=-9};1') | \
    sort -k1,1 | \
    awk -v OFS='\t' '{$1=$2};1'
  ) > work/ngs.pheno
  cut -f2 work/ngs.pheno | \
  sed '1d' > work/affymetrix.id

  export interval=/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
  export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/
  export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
  export TMPDIR=/rds/user/jhz22/hpc-work/work

  gunzip -c ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz | \
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/INFO\n' | \
  awk -v OFS="\t" 'NR>1{print $1,$2,$1 ":" $2 "_" $3 "/" $4, $3, $4, $5, $6, $7}' | \
  awk '$8 >= 0.8 {print $1":"$2}' > work/chrX.incl-positions
  paste work/affymetrix.id -d_ work/affymetrix.id > work/chrX.incl-samples
  qctool -g ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz -og work/ngs-X.bgen -os work/ngs-X.samples \
         -incl-positions work/chrX.incl-positions -incl-samples work/chrX.incl-samples
  # MAF cutoff 0.05
  module load plink/2.00-alpha
  plink2 --bgen work/ngs-X.bgen --sample work/ngs-X.samples --maf 0.05 --make-bed --out work/ngs-X.05
  cut -f2 work/ngs-X.05.bim > work/ngs-X.05.snpids
  awk '{if(NR<=2) print; else {split($1,a,"_"); print a[1],a[2],$3}}' work/ngs-X.samples > work/chrX.samples
  qctool -g work/ngs-X.bgen -og work/chrX.bgen -bgen-bits 8 -incl-snpids work/ngs-X.05.snpids
  bgenix -g work/chrX.bgen -index -clobber
  seq 22 | \
  parallel -j5 --env interval --env ref -C' ' '
    sed "1d" ${ref}/impute_{}_interval.snpstats | \
    awk -v OFS="\t" "\$15>=0.05{if(\$1==\".\") \$1=\$3+0 \":\" \$4 \"_\" \$5 \"/\" \$6; print \$3,\$4,\$1,\$5,\$6,\".\",\".\",\$19}" | \
    awk "\$8 >= 0.8 {print \$1\":\"\$2}" > work/chr{}.incl-positions
    # NGS samples
    qctool -g ${interval}/impute_{}_interval.bgen -s ${interval}/interval.samples \
           -incl-samples work/affymetrix.id -incl-positions work/chr{}.incl-positions \
           -og work/chr{}.bgen -os work/chr{}.samples
  '
}

function qPCR()
{
  export OLINK=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/jp549/olink-merged-output
  ls $OLINK/*gz | \
  xargs -l basename -s _chr_merged.gz | \
  awk '{gsub(/INTERVAL_/,"");sub(/_/," ");gsub(/___/," ")};1'> qPCR.list
  (
    gunzip -c ${OLINK}/INTERVAL_cvd3_SELP___P16109_chr_merged.gz | \
    awk 'NR==1{print "UniProt","prot","chr","pos",$2,$5,$6,$24,$25,$22,"NGS"}'
    join -12 -23 <(awk 'NR>1&&!/NEUROLOGY/&&!/ONCOLOGY/{split($5,a,"_");print a[1],a[2],$6}' ${pval}/NGS.sentinels | sort -k2,2) \
                 <(sort -k3,3 qPCR.list) | \
    awk '{
       gsub(/chr/,"",$3);
       split($3,a,":");
       chr=a[1];
       pos=a[2];
       print $0,chr,pos
    }' | \
    parallel --env OLINK -C' ' '
      zgrep -H -w {7} ${OLINK}/INTERVAL_{4}_{5}___{1}_chr_merged.gz | \
      awk -vchr={6} "chr==\$3+0" | \
      awk -vuniprot={1} -v prot={5} -v chr={6} -v pos={7} -vpanel={2}-{4} "
         {print uniprot, prot, chr, pos, \$2,\$5,\$6,\$24,\$25,\$22,panel}"
    '
  ) > qPCR-${pval}.lookup
}

export pval=1e-5
qPCR
export pval=1e-6
qPCR
export pval=1e-7
qPCR

function cvt()
{
R --no-save <<END
  require(gap)
  require(pQTLtools)
  SCALLOP <- Sys.getenv("SCALLOP")
  pval <- Sys.getenv("pval")
  f <- file.path(SCALLOP,"NGS",paste0("qPCR-",pval,".lookup"))
  clumped <- read.table(f,as.is=TRUE,header=TRUE)[c("chr","pos","rsid","prot","UniProt")]
  names(clumped) <- c("Chr","bp","SNP","prot","UniProt")
  hg19 <- within(hg19Tables[c("X.chrom", "chromStart", "chromEnd","hgncSym","uniprotName","acc")],
          {chr <- sub("chr","",X.chrom);prot <- sub("_HUMAN","",uniprotName);uniprot <- acc})
  names(hg19) <- c("chrom","start","end","gene","uniprotName","UniProt","uniprot","prot","chr")
  hits <- merge(clumped,hg19,by="UniProt")
  cistrans <- cis.vs.trans.classification(hits[c("Chr","bp","SNP","UniProt")],hg19[c("chr","start","end","gene","prot","UniProt")],"UniProt")
  cis.vs.trans <- with(cistrans,data)
  write.table(cis.vs.trans,file=file.path(SCALLOP,"NGS",paste0("qPCR-",pval,".cis.vs.trans")),
              quote=FALSE,row.names=FALSE)
  sink(file.path(SCALLOP,"NGS",paste0("qPCR-",pval,".table")))
  with(cistrans,table)
  sink()
  names(clumped) <- c("CHR","BP","SNP","prot","UniProt")
  write.table(clumped,file=file.path(SCALLOP,"NGS",paste0("qPCR-",pval,".hits")),row.names=FALSE,quote=FALSE)
  pdf(file.path(SCALLOP,"NGS",paste0("qPCR",pval,".circos.pdf")))
  circos.cis.vs.trans.plot(hits=file.path(SCALLOP,"NGS",paste0("qPCR-",pval,".hits")),subset(hg19,prot %in% with(clumped,prot)),"uniprot")
  dev.off()
END
}
# UniProt prot chr pos rsid alleleA alleleB frequentist_add_beta_1 frequentist_add_se_1 frequentist_add_pvalue NGS

export pval=1e-5
cvt
export pval=1e-6
cvt
export pval=1e-7
cvt

function qPCR_ngs()
{
  cd work; ln -sf INFLAMMATION-inf-QC.dat INFLAMMATION-inf1-QC.dat; cd -
  export lookup=qPCR-${pval}.lookup
  (
    head -1 ${lookup} | \
    awk '{print $0,"P","r"}'
    cat ${lookup} | \
    parallel -C' ' '
      export MarkerName="chr{3}:{4}"
      export all="{1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}"
      export panel=$(echo {11} | cut -d"-" -f1)
      grep -e ${panel}_{1} ${pval}/NGS.sentinels | \
      grep -e ${MarkerName} | \
      awk -vall="${all}" "!/NEUROLOGY/&&!/ONCOLOGY/{sub(/_/,\" \",\$5);print all,\$4,\$5}"
    ' | \
    parallel -C' ' '
      export all="{1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}"
      if [ -f work/{11}-QC.dat ]; then
         grep -e {13} work/{11}-QC.dat | \
         grep -e {14} | \
         awk -vall="${all}" "{print all,\$4}"
      else
         echo ${all} NA
      fi
    '
  ) > qPCR-${pval}.Pr
# .lookup
# UniProt prot chr pos rsid alleleA alleleB frequentist_add_beta_1 frequentist_add_se_1 frequentist_add_pvalue NGS
# .sentinels
# Chrom Start End P prot MarkerName
# .dat
# NGS UniProt Prot r
# benchmarks
# P27930 - INFLAMMATION, cvd3
# P13500 - NEUROLOGY, cvd3/inf1
}

export pval=1e-5
qPCR_ngs
export pval=1e-6
qPCR_ngs
export pval=1e-7
qPCR_ngs
