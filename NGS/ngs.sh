#!/usr/bin/bash

function init()
{
  grep -v -e CSA tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv > work/NPX.csv
  ln -sf $HOME/rds/post_qc_data/interval/phenotype/olink_proteomics
  module load ceuadmin/stata
  staya <<\ \ END
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
}

# init

for panel in CARDIOMETABOLIC INFLAMMATION NEUROLOGY ONCOLOGY
do
  export panel=${panel}
  for opt in LOD QC column1
      export opt=${opt}
      if [ ${opt} == "LOD" ]; then
         awk -vFS=';' '$11 > $12 {$12="NA")};1' work/NPX.csv > work/NPX-${opt}.csv
      elif [ ${opt} == "QC" ]; then
         grep -v -e WARN work/NPX.csv > work/NPX-${opt}.csv
      else
         grep -v -e '01;' work/NPX.csv > work/NPX-${opt}.csv
      fi
      stata -b do ngs.do
      if [ "${panel}" != "ONCOLOGY" ]; then R --no-save -q <ngs.R; fi
  done
done
