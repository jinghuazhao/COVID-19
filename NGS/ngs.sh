#!/usr/bin/bash

function init()
{
  export src=tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv
  sed '1d' ${src} | \
  cut -d';' -f4,5 | \
  sort | \
  uniq | \
  tr ';' ' ' > work/ids.dat
  grep -v -e CSA ${src} | \
  sed 's/NaN/./g' > work/NPX.csv
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
}

init

export opanels=(cvd2 cvd3 inf neu)
export qc_opanels=(qc_cvd2 qc_cvd3 qc_inf neu_qc)
export panels=(CARDIOMETABOLIC CARDIOMETABOLIC INFLAMMATION NEUROLOGY)

for i in $(seq 0 3)
do
  export opanel=${opanels[$i]}
  export qc_opanel=${qc_opanels[$i]}
  export panel=${panels[$i]}
  for opt in LOD QC col1
  do
      export opt=${opt}
      echo ${panel} - ${opanel} - ${opt}
      if [ ${opt} == "LOD" ]; then
         awk -vFS=';' '$11 > $12 {$12="."};1' work/NPX.csv > work/NPX-${opt}.csv
      elif [ ${opt} == "QC" ]; then
         grep -v -e WARN work/NPX.csv > work/NPX-${opt}.csv
      else
         grep -v -e '01;' work/NPX.csv > work/NPX-${opt}.csv
      fi
      stata -b do ngs.do
      R --no-save -q < ngs.R
  done
done
