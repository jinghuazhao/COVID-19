// Phenotype

program single_imputation
// https://stats.idre.ucla.edu/stata/seminars/mi_in_stata_pt1_new/
// less useful after conversion of the whole cohort into 8-bit format
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

local dir : env dir
local ev : env ev

// 1. PCA
insheet using "`ev'", case clear delim(" ")
keep ID Batch PC_1-PC_20
format ID %15.0g
sort ID
save work/INTERVAL-pca, replace

// 2. Data
// insheet using "`dir'/06-05-2020/INTERVAL/INTERVALdata_06MAY2020.csv", case clear
// insheet using "20200520/INTERVALdata_20MAY2020.csv", case clear
insheet using "20200603/INTERVALdata_03JUN2020.csv", case clear
// insheet using "20200617/INTERVALdata_17JUN2020.csv", case clear
sort identifier
keep identifier sexPulse agePulse
rename agePulse age
gen age2=age*age
rename sexPulse sex
gen sexage=sex*age
save work/INTERVAL-data, replace

// 3. Omics
// insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_OmicsMap_20200506.csv", case clear
// insheet using "20200520/INTERVAL_OmicsMap_20200520.csv", case clear
// insheet using "20200603/INTERVAL_OmicsMap_20200603.csv", case clear
insheet using "20200617/INTERVAL_OmicsMap_20200617.csv", case clear
keep identifier Affymetrix_QC_bl Affymetrix_gwasQC_bl
format Affymetrix_QC_bl %15.0g
format Affymetrix_gwasQC_bl %15.0g
rename Affymetrix_gwasQC_bl ID
merge 1:1 identifier using work/INTERVAL-data, gen(omics_data)
count if ID==.
drop if ID==.
sum sex age
merge 1:1 ID using work/INTERVAL-pca, gen(omics_data_pca)
gen sex1=sex-1
sort identifier
drop if sex==. | age ==. | PC_1==.
gzsave work/INTERVAL, replace

insheet ID using work/INTERVAL.samples, case clear
format ID %15.0g
gen idn=_n
save work/INTERVAL-omics, replace
insheet ID2 using work/INTERVAL-X.samples, case clear
gen idx=_n
save work/INTERVAL-omics-X, replace

// 4. COVID-19
// insheet using "06-05-2020/INTERVAL/INTERVAL_Covid_06MAY2020.csv", case clear
// insheet using "20200520/INTERVAL_Covid_20MAY2020.csv", case clear
// insheet using "20200603/INTERVAL_Covid_03JUN2020.csv", case clear
insheet using "20200617/INTERVAL_Covid_17JUN2020.csv", case clear
sort identifier
egen SARS_CoV=rowtotal(SARS_CoV2_1 SARS_CoV2_2 SARS_CoV2_3 SARS_CoV2_4 SARS_CoV2_5 SARS_CoV2_6 SARS_CoV2_7)
replace SARS_CoV=1 if SARS_CoV>0
keep identifier SARS_CoV
save work/covid, replace

// 5. INTERVAL-COVID
gzuse work/INTERVAL, clear
drop if identifier==.
merge 1:1 identifier using work/covid, gen(omics_data_pca_covid)
keep if omics_data_pca_covid==3
save work/INTERVAL-covid, replace

gzuse work/INTERVAL,clear
merge 1:1 ID using work/INTERVAL-covid
drop _merge
merge 1:1 ID using work/INTERVAL-omics
keep if _merge==3
drop _merge
tab SARS_CoV
tab sex if SARS_CoV!=.
tabstat age if SARS_CoV!=., stat(mean sd) by(sex)
outsheet ID if sex==. | age==. using work/INTERVAL.excl-samples, noname replace
drop if sex==. | age==.
// replace SARS_CoV=0 if SARS_CoV==.
tab SARS_CoV
tab sex if SARS_CoV!=.
tabstat age if SARS_CoV!=., stat(mean sd) by(sex)
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
