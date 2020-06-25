
local panel : env panel
local opt : env opt
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
