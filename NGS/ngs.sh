#!/usr/bin/bash

# cols 1-12
grep -v CSA tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv > NPX.csv

module load ceuadmin/stata
stata <<END
  insheet using tmp_raw_data/20200292_Danesh_NPX_2020-06-03.csv, case clear delim(;)
  save work/ngs, replace
  insheet using INTERVALdata_19JUN2020.csv, case clear
  save work/data, replace
  insheet using INTERVAL_OmicsMap_20200619.csv, case clear
  save work/omics, replace
  insheet using Olink_NGS.csv, case clear
  rename OLINK_NGS SampleID
  merge 1:1 identifier using work/data, gen(data)
  merge 1:1 identifier using work/omics,gen(data_omics)
  keep if SampleID!=""
  save work/dataomics, replace
  use work/ngs if Panel=="INFLAMMATION"
  merge m:1 SampleID using work/dataomics, gen(olink_ngs)
  d
END

// CARDIOMETABOLIC, INFLAMMATION, NEUROLOGY, ONCOLOGY
