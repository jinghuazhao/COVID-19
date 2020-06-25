
local panel : env panel
local opt : env opt
insheet using work/NPX-`opt'.csv, case clear delim(;)
save work/ngs-`opt', replace
use work/ngs-`opt' if Panel=="`panel'"
merge m:1 SampleID using work/dataomics, gen(olink_ngs)
save work/`panel', replace
d
keep Affymetrix_gwasQC_bl UniProt LOD NPX
sort Affymetrix_gwasQC_bl UniProt
by Affymetrix_gwasQC_bl: gen j=_n
reshape wide UniProt LOD NPX, i(Affymetrix_gwasQC_bl) j(j)
outsheet using work/`panel'-`opt', noquote replace
