// Phenotype preparation
// with focus on check for test dates

local dir : env HGI
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
// insheet using "20200603/INTERVALdata_03JUN2020.csv", case clear
// insheet using "20200731/INTERVALdata_31JUL2020.csv", case clear
insheet using "Data_INTERVAL_latest/INTERVALdata_01DEC2020.csv", case clear
sort identifier
keep identifier sexPulse agePulse attendanceDate
rename sexPulse sex
gen age_at_test=agePulse+(date("01DEC2020","DMY")-date(attendanceDate,"DMY"))/365.25
sum sex age_at_test
save work/INTERVAL-data, replace

// 3. Omics
// insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_OmicsMap_20200506.csv", case clear
// insheet using "20200520/INTERVAL_OmicsMap_20200520.csv", case clear
// insheet using "20200603/INTERVAL_OmicsMap_20200603.csv", case clear
// insheet using "20200731/INTERVAL_OmicsMap_20200731.csv", case clear
insheet using "Data_INTERVAL_latest/INTERVAL_OmicsMap_20201201.csv", case clear
keep identifier Affymetrix_QC_bl Affymetrix_gwasQC_bl
format Affymetrix_QC_bl %15.0g
format Affymetrix_gwasQC_bl %15.0g
rename Affymetrix_gwasQC_bl ID
merge 1:1 identifier using work/INTERVAL-data, gen(omics_data)
count if ID==.
drop if ID==.
merge 1:1 ID using work/INTERVAL-pca, gen(omics_data_pca)
sort identifier
drop if sex==. | agePulse ==. | PC_1==.
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
// insheet using "20200731/INTERVAL_Covid_31JUL2020.csv", case clear
insheet using "Data_INTERVAL_latest/INTERVAL_Covid_01DEC2020.csv", case clear
sort identifier
#delimit ;
egen SARS_CoV2=rowtotal(
               SARS_CoV2_1 SARS_CoV2_2 SARS_CoV2_3 SARS_CoV2_4 SARS_CoV2_5 SARS_CoV2_6 SARS_CoV2_7 SARS_CoV2_8 SARS_CoV2_9 SARS_CoV2_10
               SARS_CoV2_11 SARS_CoV2_12 SARS_CoV2_13 SARS_CoV2_14 SARS_CoV2_15 SARS_CoV2_16 SARS_CoV2_17 SARS_CoV2_18 SARS_CoV2_19 SARS_CoV2_20
               SARS_CoV2_21 SARS_CoV2_22 SARS_CoV2_23 SARS_CoV2_24 SARS_CoV2_25 SARS_CoV2_26 SARS_CoV2_27 SARS_CoV2_28 SARS_CoV2_29 SARS_CoV2_30
               SARS_CoV2_31 SARS_CoV2_32 SARS_CoV2_33 SARS_CoV2_34 SARS_CoV2_35
             );
#delimit cr
gen SARS_CoV=SARS_CoV2
replace SARS_CoV=1 if SARS_CoV2>0
save work/covid, replace

// 5. INTERVAL-COVID
gzuse work/INTERVAL, clear
drop if identifier==.
merge 1:1 identifier using work/covid, gen(omics_data_pca_covid)
drop if ID==.
tab SARS_CoV2
forval x=1/35 {
  gen age_`x'=agePulse+(date(specimenDate_`x',"DMY")-date(attendanceDate,"DMY"))/365.25
}
#delimit ;
replace age_at_test=max(
                     age_1,age_2,age_3,age_4,age_5,age_6,age_7,age_8,age_9,age_10,
                     age_11,age_12,age_13,age_14,age_15,age_16,age_17,age_18,age_19,age_20,
                     age_21,age_22,age_23,age_24,age_25,age_26,age_27,age_28,age_29,age_30,
                     age_31,age_32,age_33,age_34,age_35
                   ) if SARS_CoV2==0;
#delimit cr
forval x=1/35 {
  replace age_at_test=age_`x' if SARS_CoV2!=0 & SARS_CoV2_`x'==1
}
format age_* %6.4g
save work/INTERVAL-covid, replace

gzuse work/INTERVAL,clear
merge 1:1 ID using work/INTERVAL-covid
rename age_at_test age
gen agegroup=.
replace agegroup=1 if age<=60
replace agegroup=2 if age>60
tab agegroup SARS_CoV
tab sex SARS_CoV
gen age2=age*age
gen sexage=sex*age
drop _merge
merge 1:1 ID using work/INTERVAL-omics
keep if _merge==3
drop _merge
replace SARS_CoV=0 if SARS_CoV==.
tab agegroup SARS_CoV
tab sex SARS_CoV
gzsave work/INTERVAL-omics-covid, replace

program CxV2
//1=1/2 for C1/C2; 2=working directory
  gzuse work/INTERVAL-omics-covid, clear
  tab SARS_CoV
  tab sex if SARS_CoV!=.
  tabstat age if SARS_CoV!=., stat(mean sd n) by(sex) format
  outsheet ID if sex==. | age==. using `2'/work/INTERVAL.excl-samples, noname replace
  drop if sex==. | age==.
  replace SARS_CoV=0 if SARS_CoV==. & `1'==2
  tab SARS_CoV
  tab sex if SARS_CoV!=.
  tab sex SARS_CoV
  tab sex SARS_CoV if age<=60
  tab sex SARS_CoV if age>60
  tabstat age if SARS_CoV!=., stat(mean sd n) by(sex) format
  tabstat age if SARS_CoV!=., stat(mean sd n) by(SARS_CoV) format
  drop if SARS_CoV==.
  outsheet ID SARS_CoV `3' PC_1-PC_20 using `2'/work/INTERVAL-covid.txt, delim(" ") noquote replace
  tostring ID,gen(IDS) format(%15.0g)
  gen str31 ID2=IDS + "_" + IDS
  label define sexFM 1 "M" 2 "F"
  label values sex sexFM
  drop if idn==.
  gzsave `2'/INTERVAL-covid, replace
  merge 1:1 ID2 using work/INTERVAL-omics-X
  keep if _merge==3
  drop _merge
  outsheet ID2 sex using `2'/work/INTERVAL-X.FM if idx!=., noname noquote replace
end

local covlist sex age age2 sexage
CxV2 2 20201201-ANA_C2_V2 "`covlist'"

program C2sex
  // 1=working directory; 2=sex
  gzuse work/INTERVAL-omics-covid,clear
  outsheet ID if sex!=`2' | age==. using `1'/work/INTERVAL.excl-samples, noname replace
  drop if sex==. | age==.
  replace SARS_CoV=0 if SARS_CoV==.
  drop if SARS_CoV==. | sex!=`2'
  outsheet ID SARS_CoV age age2 PC_1-PC_20 using `1'/work/INTERVAL-covid.txt, delim(" ") noquote replace
  tostring ID,gen(IDS) format(%15.0g)
  gen str31 ID2=IDS + "_" + IDS
  label define sexFM 1 "M" 2 "F"
  label values sex sexFM
  drop if idn==.
  gzsave `1'/work/INTERVAL-covid, replace
  merge 1:1 ID2 using work/INTERVAL-omics-X
  keep if _merge==3
  drop _merge
  outsheet ID2 age age2 using `1'/work/INTERVAL-X.FM if idx!=., noname noquote replace
end

C2sex 20201201-male-ANA_C2_V2 1
C2sex 20201201-female-ANA_C2_V2 2

program C2age
// 1=working directory; 2=age
  gzuse work/INTERVAL-omics-covid, clear
  outsheet ID if sex==. | agegroup!=`2' using `1'/work/INTERVAL.excl-samples, noname replace
  drop if sex==. | age==.
  replace SARS_CoV=0 if SARS_CoV==.
  drop if SARS_CoV==. | agegroup!=`2'
  outsheet ID SARS_CoV sex PC_1-PC_20 using `1'/work/INTERVAL-covid.txt, delim(" ") noquote replace
  tostring ID,gen(IDS) format(%15.0g)
  gen str31 ID2=IDS + "_" + IDS
  label define sexFM 1 "M" 2 "F"
  label values sex sexFM
  drop if idn==.
  gzsave `1'/work/INTERVAL-covid, replace
  merge 1:1 ID2 using work/INTERVAL-omics-X
  keep if _merge==3
  drop _merge
  outsheet ID2 sex using `1'/work/INTERVAL-X.FM if idx!=., noname noquote replace
end

C2age 20201201-le_60-ANA_C2_V2 1
C2age 20201201-gt_60-ANA_C2_V2 2

exit,clear
// 31JUL2020 version
program CxV2
//1=1/2 for C1/C2; 2=working directory
  gzuse work/INTERVAL-omics-covid, clear
  tab SARS_CoV
  tab sex if SARS_CoV!=.
  tabstat age if SARS_CoV!=., stat(mean sd n) by(sex) format
  outsheet ID if sex==. | age==. using `2'/work/INTERVAL.excl-samples, noname replace
  drop if sex==. | age==.
  replace SARS_CoV=0 if SARS_CoV==. & `1'==2
  tab SARS_CoV
  tab sex if SARS_CoV!=.
  tab sex SARS_CoV
  tab sex SARS_CoV if age<=60
  tab sex SARS_CoV if age>60
  tabstat age if SARS_CoV!=., stat(mean sd n) by(sex) format
  tabstat age if SARS_CoV!=., stat(mean sd n) by(SARS_CoV) format
  drop if SARS_CoV==.
  outsheet ID SARS_CoV `3' PC_1-PC_20 using `2'/work/INTERVAL-covid.txt, delim(" ") noquote replace
  tostring ID,gen(IDS) format(%15.0g)
  gen str31 ID2=IDS + "_" + IDS
  label define sexFM 1 "M" 2 "F"
  label values sex sexFM
  drop if idn==.
  gzsave `2'/INTERVAL-covid, replace
  merge 1:1 ID2 using work/INTERVAL-omics-X
  keep if _merge==3
  drop _merge
  outsheet ID2 sex using `2'/work/INTERVAL-X.FM if idx!=., noname noquote replace
end

local covlist sex age age2 sexage
CxV2 1 20200731-ANA_C1_V2 "`covlist'"
CxV2 2 20200731-ANA_C2_V2 "`covlist'"

program sexCxV2
  // 1=1/2 for C1/C2; 2=working directory; 3=sex 4=covlist
  gzuse work/INTERVAL-omics-covid,clear
  outsheet ID if sex==. | age==. using `2'/work/INTERVAL.excl-samples, noname replace
  drop if sex==. | age==.
  replace SARS_CoV=0 if SARS_CoV==. & `1'==2
  drop if SARS_CoV==. | sex!=`3'
  outsheet ID SARS_CoV `4' PC_1-PC_20 using `2'/work/INTERVAL-covid.txt, delim(" ") noquote replace
  tostring ID,gen(IDS) format(%15.0g)
  gen str31 ID2=IDS + "_" + IDS
  label define sexFM 1 "M" 2 "F"
  label values sex sexFM
  drop if idn==.
  gzsave `2'/work/INTERVAL-covid, replace
  merge 1:1 ID2 using work/INTERVAL-omics-X
  keep if _merge==3
  drop _merge
  outsheet ID2 sex using `2'/work/INTERVAL-X.FM if idx!=., noname noquote replace
end

local covlist age age2
sexCxV2 1 20200731-male-ANA_C1_V2 1 "`covlist'"
sexCxV2 2 20200731-male-ANA_C2_V2 1 "`covlist'"
sexCxV2 1 20200731-female-ANA_C1_V2 2 "`covlist'"
sexCxV2 2 20200731-female-ANA_C2_V2 2 "`covlist'"

program sex60CxV2
  // 1=1/2 for C1/C2; 2=working directory; 3=sex
  gzuse work/INTERVAL-omics-covid, clear
  outsheet ID if sex==. | age==. using `2'/work/INTERVAL.excl-samples, noname replace
  drop if sex==. | age==.
  replace SARS_CoV=0 if SARS_CoV==. & `1'==2
  drop if SARS_CoV==. | sex!=`3' | age>60
  outsheet ID SARS_CoV PC_1-PC_20 using `2'/work/INTERVAL-covid.txt, delim(" ") noquote replace
  tostring ID,gen(IDS) format(%15.0g)
  gen str31 ID2=IDS + "_" + IDS
  label define sexFM 1 "M" 2 "F"
  label values sex sexFM
  drop if idn==.
  gzsave `2'/work/INTERVAL-covid, replace
  merge 1:1 ID2 using work/INTERVAL-omics-X
  keep if _merge==3
  drop _merge
  outsheet ID2 sex using `2'/work/INTERVAL-X.FM if idx!=., noname noquote replace
end

sex60CxV2 1 20200731-male-60-ANA_C1_V2 1
sex60CxV2 2 20200731-male-60-ANA_C2_V2 1
sex60CxV2 1 20200731-female-60-ANA_C1_V2 2
sex60CxV2 2 20200731-female-60-ANA_C2_V2 2

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
