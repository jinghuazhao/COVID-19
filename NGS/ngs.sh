#!/usr/bin/bash

export opanels=(cvd2 cvd3 inf neu)
export qc_opanels=(qc_cvd2 qc_cvd3 qc_inf neu_qc)
export panels=(CARDIOMETABOLIC CARDIOMETABOLIC INFLAMMATION NEUROLOGY)

function init()
{
  export src=tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv
  sed '1d' ${src} | \
  cut -d';' -f4,5 | \
  sort | \
  uniq | \
  tr ';' ' ' > work/ids.dat
  R --no-save -q <<\ \ END
    ids <- read.table("work/ids.dat",col.names=c("UniProt","Prot"),as.is=TRUE)
    save(ids,file="work/ids.rda")
  END
  grep -v -e CSA ${src} | \
  sed 's/NaN/NA/g' > work/NPX.csv
  ln -sf $HOME/rds/post_qc_data/interval/phenotype/olink_proteomics
  module load ceuadmin/stata
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

function check()
{
for i in $(seq 0 3)
do
  export opanel=${opanels[$i]}
  export qc_opanel=${qc_opanels[$i]}
  export panel=${panels[$i]}
  echo ${qc_opanel}
  wc -l work/NPX-${opanel}.id
  (
    head -1 work/NPX.csv
    grep -f work/NPX-${opanel}.id work/NPX.csv
  ) > work/NPX-${opanel}.csv
  for opt in raw LOD QC col1
  do
      export opt=${opt}
      echo ${panel} - ${opanel} - ${opt}
      if [ ${opt} == "raw" ]; then
         cat work/NPX-${opanel}.csv > work/NPX-${opt}.csv
      elif [ ${opt} == "LOD" ]; then
         awk -vFS=';' -vOFS=';' '$11 > $12 {$12="NA"};1' work/NPX-${opanel}.csv > work/NPX-${opt}.csv
      elif [ ${opt} == "QC" ]; then
         awk -vFS=';' -vOFS=';' '$10 =="WARN" {$12="NA"};1' work/NPX-${opanel}.csv > work/NPX-${opt}.csv
      else
         awk -vFS=';' -vOFS=';' '/01;/{$12="NA"};1' work/NPX-${opanel}.csv > work/NPX-${opt}.csv
      fi
      stata -b do ngs.do
      R --no-save -q < ngs.R
  done
done
}

ngs

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
      plot(d,axes=FALSE,main="",ylim=c(0,2))
      axis(2)
    } else plot(d,ann=FALSE,axes=FALSE,ylim=c(0,2))
    title(rt)
  }
}
opanels <- c("cvd2","cvd3","inf","neu")
qc_opanels <- c("qc_cvd2","qc_cvd3","qc_inf","neu_qc")
panels <- c("CARDIOMETABOLIC","CARDIOMETABOLIC","INFLAMMATION","NEUROLOGY")
pdf("correlogram.pdf")
par(mfrow=c(2,2))
for (opt in c("LOD","QC","col1")) correlogram(opt)
dev.off()
END
