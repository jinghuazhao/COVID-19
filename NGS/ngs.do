set maxvar 50000
local opanel: env opanel
local panel : env panel
local opt : env opt

insheet using work/NPX-`opt'.csv, case clear delim(;)
sort UniProt
merge m:1 UniProt using work/NPX-`opanel', gen(npx)
sort SampleID
keep if npx==3
// keep if Panel=="`panel'"
merge m:1 SampleID using work/dataomics, gen(olink_ngs)
keep Olink_`opanel'_QC_24m UniProt LOD NPX
sort Olink_`opanel'_QC_24m UniProt
by Olink_`opanel'_QC_24m: gen j=_n
reshape wide UniProt LOD NPX, i(Olink_`opanel'_QC_24m) j(j)
outsheet using work/`panel'-`opanel'-`opt', noquote replace
