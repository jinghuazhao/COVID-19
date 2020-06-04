// Phenotype

// 1. PCA
local dir : env dir
local ev : env ev
insheet using "`ev'", case clear delim(" ")
keep ID Batch PC_1-PC_20
sort ID
format ID %15.0g
save work/INTERVAL-pca, replace

// 2. Data
// insheet using "`dir'/06-05-2020/INTERVAL/INTERVALdata_06MAY2020.csv", case clear
insheet using "20200520/INTERVALdata_20MAY2020.csv", case clear
sort identifier
save work/INTERVAL-data, replace

// 3. Omics
// insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_OmicsMap_20200506.csv", case clear
 insheet using "20200520/INTERVAL_OmicsMap_20200520.csv", case clear
sort identifier
merge 1:1 identifier using work/INTERVAL-data, gen(dataid)
rename Affymetrix_gwasQC_bl ID
format ID %15.0g
duplicates tag ID, gen(dup)
keep if dup==0
drop dup
merge 1:1 ID using work/INTERVAL-pca
drop _merge
rename agePulse age
gen age2=age*age
rename sexPulse sex
gen sex1=sex-1
gen sexage=sex*age
keep ID identifier sex1 sex age age2 sexage PC_1-PC_20
gzsave work/INTERVAL, replace

// 4. COVID-19
// insheet using "06-05-2020/INTERVAL/INTERVAL_Covid_06MAY2020.csv", case clear
insheet using "20200520/INTERVAL_Covid_20MAY2020.csv", case clear
sort identifier
egen SARS_CoV=rowtotal(SARS_CoV2_1 SARS_CoV2_2 SARS_CoV2_3 SARS_CoV2_4 SARS_CoV2_5 SARS_CoV2_6)
replace SARS_CoV=1 if SARS_CoV>0
drop SARS_CoV2_1-SARS_CoV2_6
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
insheet ID2 using work/INTERVAL-X.samples, case clear
gen idx=_n
save work/INTERVAL-omics-X, replace

gzuse work/INTERVAL
merge 1:1 ID using work/INTERVAL-covid
drop _merge
merge 1:1 ID using work/INTERVAL-omics
keep if _merge==3
drop _merge
tab SARS_CoV
tab sex if SARS_CoV!=.
tabstat age if SARS_CoV!=., stat(mean sd) by(sex)
outsheet ID if sex==. | age==. using work/INTERVAL.excl-samples, noname replace

program single_imputation
// https://stats.idre.ucla.edu/stata/seminars/mi_in_stata_pt1_new/
mi set wide
mi register imputed age sex1
mi register regular PC_1-PC_20
mi impute chained (regress) age (logit) sex1 = PC_1-PC_20, add(1) rseed(123456)
replace SARS_CoV=0 if SARS_CoV!=1
replace age=_1_age
replace age2=age*age
replace sex=_1_sex1+1
replace sexage=sex*age
end

// drop if SARS_CoV==.
// single_imputation
drop if sex==. | age==.
outsheet ID SARS_CoV sex age age2 sexage PC_1-PC_20 using work/INTERVAL-covid.txt, delim(" ") noquote replace
tostring ID,gen(IDS) format(%15.0g)
gen str31 ID2=IDS + "_" + IDS
label define sexFM 1 "M" 2 "F"
label values sex sexFM
drop if idn==.
gzsave work/INTERVAL-covid, replace
merge 1:1 ID2 using work/INTERVAL-omics-X
keep if _merge==3
drop _merge
outsheet ID2 sex using work/INTERVAL-X.FM if idx!=., noname noquote replace
