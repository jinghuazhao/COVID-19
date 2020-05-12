// 12-5-2020 JHZ

local dir : env dir
local ev : env ev
insheet using "`ev'", case clear delim(" ")
keep ID Batch PC_1-PC_20
sort ID
rename ID Affymetrix_gwasQC_bl
save INTERVAL-pca, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_Covid_06MAY2020.csv", case clear
sort identifier
save covid, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVALdata_06MAY2020.csv", case clear
sort identifier
save INTERVAL, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_OmicsMap_20200506.csv", case clear
sort identifier
merge 1:1 identifier using INTERVAL, gen(dataid)
keep if data==3
duplicates tag Affymetrix_gwasQC_bl, gen(dup)
keep if dup==0
drop dup
merge 1:1 Affymetrix_gwasQC_bl using INTERVAL-pca, gen(dataidpca)
keep if dataidpca==3
merge 1:1 identifier using covid, gen(intervalcovid)
keep if intervalcovid==3
outsheet Affymetrix_gwasQC_bl agePulse sexPulse PC_1-PC_20 using INTERVAL-covid.txt, noquote replace
gzsave INTERVAL-covid, replace
