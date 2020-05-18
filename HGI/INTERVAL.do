// 18-5-2020 JHZ

local dir : env dir
local ev : env ev
insheet using "`ev'", case clear delim(" ")
keep ID Batch PC_1-PC_20
sort ID
format ID %15.0g
save work/INTERVAL-pca, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_Covid_06MAY2020.csv", case clear
sort identifier
save work/covid, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVALdata_06MAY2020.csv", case clear
sort identifier
save work/INTERVAL, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_OmicsMap_20200506.csv", case clear
sort identifier
merge 1:1 identifier using work/INTERVAL, gen(dataid)
keep if data==3
rename Affymetrix_gwasQC_bl ID
format ID %15.0g
duplicates tag ID, gen(dup)
keep if dup==0
drop dup
merge 1:1 ID using work/INTERVAL-pca, gen(dataidpca)
keep if dataidpca==3
merge 1:1 identifier using work/covid, gen(intervalcovid)
keep if intervalcovid==3
rename agePulse age
rename sexPulse sex
egen SARS_CoV=rowtotal(SARS_CoV2_1 SARS_CoV2_2 SARS_CoV2_3 SARS_CoV2_4 SARS_CoV2_5 SARS_CoV2_6)
replace SARS_CoV=1 if SARS_CoV>0
outsheet SARS_CoV age sex PC_1-PC_20 ID using work/INTERVAL-covid.txt, delim(" ") noquote replace
gzsave work/INTERVAL-covid, replace
label define sexFM 1 "M" 2 "F"
label values sex sexFM
outsheet ID sex using work/INTERVAL.FM, noname noquote replace
