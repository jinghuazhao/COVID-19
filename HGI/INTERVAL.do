// 26-5-2020 JHZ

// 1. PCA
local dir : env dir
local ev : env ev
insheet using "`ev'", case clear delim(" ")
keep ID Batch PC_1-PC_20
sort ID
format ID %15.0g
save work/INTERVAL-pca, replace

// 2. Data
insheet using "`dir'/06-05-2020/INTERVAL/INTERVALdata_06MAY2020.csv", case clear
sort identifier
save work/INTERVAL-data, replace

// 3. Omics
insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_OmicsMap_20200506.csv", case clear
sort identifier
merge 1:1 identifier using work/INTERVAL, gen(dataid)
rename Affymetrix_gwasQC_bl ID
format ID %15.0g
duplicates tag ID, gen(dup)
keep if dup==0
drop dup
merge 1:1 ID using work/INTERVAL-pca, gen(dataidpca)
gzsave work/INTERVAL, replace

// 4. COVID-19
insheet using "06-05-2020/INTERVAL/INTERVAL_Covid_06MAY2020.csv", case clear
sort identifier
save work/covid, replace

// 5. INTERVAL-COVID
gzuse work/INTERVAL
drop if identifier==.
merge 1:1 identifier using work/covid
keep if _merge==3
drop _merge
save work/INTERVAL-covid, replace

insheet ID using work/INTERVAL.samples, case clear
gen idn=_n
save work/INTERVAL-omics, replace

gzuse work/INTERVAL
merge 1:1 ID using work/INTERVAL-covid
drop _merge
merge 1:1 ID using work/INTERVAL-omics
rename agePulse age
rename sexPulse sex
outsheet ID SARS_CoV age sex PC_1-PC_20 using work/INTERVAL-covid.txt, delim(" ") noquote replace
tostring ID,gen(IDS) format(%15.0g)
gen str31 ID2=IDS + "_" + IDS
label define sexFM 1 "M" 2 "F"
label values sex sexFM
drop if idn==.
outsheet ID2 sex using work/INTERVAL-X.FM if sex!=., noname noquote replace
egen SARS_CoV=rowtotal(SARS_CoV2_1 SARS_CoV2_2 SARS_CoV2_3 SARS_CoV2_4 SARS_CoV2_5 SARS_CoV2_6)
replace SARS_CoV=1 if SARS_CoV>0
tab sex if SARS_Cov!=.
tabstat age if SARS_Cov!=., stat(mean sd) by(sex)
gzsave work/INTERVAL-covid, replace
