//usr/local/Cluster-Apps/ceuadmin/stata/15/stata

insheet using INTERVALdata_28FEB2020.csv, case clear comma
sort identifier
save work/INTERVALdata_28FEB2020, replace
insheet using omicsMap.csv, case clear comma
sort identifier
merge 1:1 identifier using work/INTERVALdata_28FEB2020
d,f
keep identifier Wes_gwasQC_bl Wgs_gwasQC_bl Wgs_QC_bl Wgs_QC_24m Wgs_gwasQC_24m Olink_* ethnicPulse agePulse sexPulse
count if ethnicPulse=="Arab"
drop if ethnicPulse=="Arab"
drop ethnicPulse
rename agePulse age
rename sexPulse sex
gen age2=age*age
format Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl Wgs_QC_24m %15s
format Wgs_gwasQC_24m Olink_*_QC_24m Olink_*gwasQC_24m %15.0g
drop if Wes_gwasQC_bl=="" & Wgs_QC_bl=="" & Wgs_gwasQC_bl=="" & Wgs_QC_24m=="" & Wgs_gwasQC_24m==.
drop if Olink_cvd2_QC_24m==. & Olink_cvd3_QC_24m==. & Olink_inf_QC_24m==. & Olink_neu_QC_24m==.
count if Wes_gwasQC_bl!=""
count if Wgs_gwasQC_bl!=""
count if Wgs_QC_bl!=""
count if Wgs_QC_24m!=""
count if Wgs_gwasQC_24m!=.
drop Wgs_gwasQC_24m
count if Wgs_gwasQC_bl!="" & Wgs_QC_bl=="" & Wgs_QC_24m==""
count if Olink_cvd2_QC_24!=.
count if Olink_cvd3_QC_24!=.
count if Olink_inf_QC_24!=.
count if Olink_neu_QC_24!=.
count if !(Wgs_gwasQC_bl==Wgs_QC_bl | Wgs_gwasQC_bl==Wgs_QC_24m)
drop Wgs_gwasQC_bl
gen Wgs=Wgs_QC_bl+Wgs_QC_24m
mvencode _all, mv(-999)
replace Wes_gwasQC_bl="-999" if Wes_gwasQC_bl==""
replace Wgs_QC_bl="-999" if Wgs_QC_bl==""
replace Wgs_QC_24m="-999" if Wgs_QC_24m==""
outsheet identifier Wes_gwasQC_bl Olink_*QC_24m if Wes_gwasQC_bl!="-999" & Wgs_QC_bl=="-999" /*
*/       & Wgs_QC_24m=="-999" using work/wes.txt, noquote replace
outsheet identifier Wgs Olink_*QC_24m if Wes_gwasQC_bl=="-999" | (Wes_gwasQC_bl!="-999" & Wes_gwasQC_bl==Wgs_QC_bl)/*
*/       using work/wgs.txt, noquote replace
outsheet identifier Wes_gwasQC_bl Wgs Wgs_QC_bl Wgs_QC_24m Olink_*QC_24m sex age age2 using work/weswgs.txt, noquote replace
exit, clear

/* The following is now furnised in R instead */
program pheno
  insheet using `1', case clear comma
  format aliquot_id %15.0g
  rename aliquot_id Olink_`2'_QC_24m
  sort Olink_`2'_QC_24m
  merge 1:m Olink_`2'_QC_24m using work/weswgs
  drop plate* Olink_*_QC_24m Olink_*_gwasQC_24m
  keep if _merge==3
  tab centre, gen(centre)
  drop identifier Affymetrix_QC_bl _merge centre centre25
  save work/`2', replace
  outsheet using work/`2'.txt, comma noquote replace
  outsheet Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl sex age centre* using work/`2'-covariates.txt,comma noquote replace
  drop Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl sex age centre*
  outsheet using work/`2'-protein.txt,comma noquote replace
end

pheno high_dimensional_data/Olink_proteomics_cvd2/qc/olink_qc_cvd2.csv cvd2
pheno high_dimensional_data/Olink_proteomics_cvd3/qc/olink_qc_cvd3.csv cvd3
pheno high_dimensional_data/Olink_proteomics_inf/qc/olink_qc_inf.csv inf
pheno high_dimensional_data/Olink_proteomics_neurology/qc/olink_neu_qc.csv neu
