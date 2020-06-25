#!/usr/bin/bash

function ngs()
{
  grep -v -e CSA tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv > work/NPX.csv
  module load ceuadmin/stata
  stata <<\ \ END
  local panel : env panel
  insheet using work/NPX.csv, case clear delim(;)
  save work/ngs, replace
  insheet using INTERVALdata_19JUN2020.csv, case clear
  save work/data, replace
  insheet using INTERVAL_OmicsMap_20200619.csv, case clear
  format Olink_* %15.0g
  save work/omics, replace
  insheet using Olink_NGS.csv, case clear
  rename OLINK_NGS SampleID
  merge 1:1 identifier using work/data, gen(data)
  merge 1:1 identifier using work/omics,gen(data_omics)
  keep if SampleID!=""
  save work/dataomics, replace
  use work/ngs if Panel=="`panel'"
  merge m:1 SampleID using work/dataomics, gen(olink_ngs)
  save work/`panel', replace
  d
  keep SampleID UniProt LOD NPX
  sort SampleID UniProt
  by SampleID: gen j=_n
  reshape wide UniProt LOD NPX, i(SampleID) j(j)
  outsheet using work/`panel', noquote replace
  END
}

# ln -sf $HOME/rds/post_qc_data/interval/phenotype/olink_proteomics

function corr()
{
  R --no-save -q <<\ \ END
# old Olink panels
  panels <- c("NEUROLOGY", "CARDIOMETABOLIC", "CARDIOMETABOLIC", "INFLAMMATION")
  qc <- c("neu_qc", "qc_cvd2", "qc_cvd3", "qc_inf")
  panel <- Sys.getenv("panel")
  sel <- which(panels==panel)
  od <- read.csv(paste0("olink_proteomics/qc/olink_",qc[sel],".csv"),as.is=TRUE)
  names(od) <- toupper(names(od))
# new data from NGS
  d <-read.delim(paste0("work/",panel,".out"),as.is=TRUE)
  ncols <- length(d[1,-1])/3
  enum <- (1:ncols)*3-1
  uniprot <- d[1,enum]
  names(d)[enum+1] <- paste0(uniprot,"_lod")
  names(d)[enum+2] <- uniprot
  d <- d[,-enum]
  write.table(d,file=paste0("work/",panel,".csv"),quote=FALSE,row.names=FALSE,sep=",")
# correlation
  for(i in intersect(names(d),names(od)))
  {
    x <- d[i]
    y <- od[i]
    cor(x,y,use="everything")
  }
  END
}

for panel in CARDIOMETABOLIC INFLAMMATION NEUROLOGY ONCOLOGY
do
  export panel=${panel}
# ngs
  corr
done
