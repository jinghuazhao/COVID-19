#!/usr/bin/bash

function init()
{
  grep -v -e CSA tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv > work/NPX.csv
  ln -sf $HOME/rds/post_qc_data/interval/phenotype/olink_proteomics
  module load ceuadmin/stata
}

# init

for panel in CARDIOMETABOLIC INFLAMMATION NEUROLOGY ONCOLOGY
do
  export panel=${panel}
  for opt in LOD QC column1
      export opt=${opt}
#     stata -b do ngs.do
      if [ "${panel}" != "ONCOLOGY" ]; then R --no-save -q <ngs.R; fi
  done
done
