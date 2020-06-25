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
  format Affymetrix_gwasQC_bl %15.0g
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
  keep Affymetrix_gwasQC_bl UniProt LOD NPX
  sort Affymetrix_gwasQC_bl UniProt
  by Affymetrix_gwasQC_bl: gen j=_n
  reshape wide UniProt LOD NPX, i(Affymetrix_gwasQC_bl) j(j)
  outsheet using work/`panel', noquote replace
  END
}

# ln -sf $HOME/rds/post_qc_data/interval/phenotype/olink_proteomics

function corr()
{
  R --no-save -q <<\ \ END
  corr <- function()
  {
    od <- read.csv(paste0("olink_proteomics/qc/olink_",qc[sel],".csv"),as.is=TRUE)
    names(od) <- toupper(names(od))
    Olink_id <- paste0("Olink_",ids[sel],"_QC_24m")
    omics <- read.csv("INTERVAL_OmicsMap_20200619.csv",as.is=TRUE)
    omics_id <- subset(omics,!is.na(omics[Olink_id]))[c(Olink_id,"Affymetrix_gwasQC_bl")]
    od <- merge(omics_id,od,by.x=Olink_id,by.y="ALIQUOT_ID")
    odd <- merge(d,od,by="Affymetrix_gwasQC_bl")
  # correlation
    for(i in setdiff(intersect(names(d),names(od)),"Affymetrix_gwasQC_bl"))
    {
      with(odd, {
      x <- odd[[paste0(i,".x")]]
      y <- odd[[paste0(i,".y")]]
      cat(panel, "-", i, ":", cor(x,y,use="everything"), "\n")
      })
    }
  }
# new data from NGS
  panels <- c("NEUROLOGY", "CARDIOMETABOLIC", "INFLAMMATION")
  qc <- c("neu_qc", "qc_cvd2", "qc_cvd3", "qc_inf")
  ids <- c("neu","cvd2","cvd3","inf")
  panel <- Sys.getenv("panel")
  sel <- which(panel==panels)
  d <-read.delim(paste0("work/",panel,".out"),as.is=TRUE)
  ncols <- length(d[1,-1])/3
  enum <- (1:ncols)*3-1
  uniprot <- d[1,enum]
  names(d)[enum+1] <- paste0(uniprot,"_lod")
  names(d)[enum+2] <- uniprot
  d <- d[,-enum]
  write.table(d,file=paste0("work/",panel,".csv"),quote=FALSE,row.names=FALSE,sep=",")
# old Olink panels
  cat("Panel =",panel,sel,"\n")
  corr()
  if (sel==2) {sel=3; corr()}
  END
}

for panel in CARDIOMETABOLIC INFLAMMATION NEUROLOGY ONCOLOGY
do
  export panel=${panel}
# ngs
  if [ "${panel}" != "ONCOLOGY" ]; then corr; fi
done
