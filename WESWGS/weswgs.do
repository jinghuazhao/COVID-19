//usr/local/Cluster-Apps/ceuadmin/stata/15/stata

insheet using INTERVALdata_28FEB2020.csv, case clear comma
sort identifier
save work/INTERVALdata_28FEB2020, replace
insheet using omicsMap.csv, case clear comma
sort identifier
merge 1:1 identifier using work/INTERVALdata_28FEB2020
d,f
keep identifier Affymetrix_QC_bl Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl Olink_* sexPulse agePulse ethnicPulse centre ht_bl wt_bl
mvencode _all, mv(-999)
format Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl %15s
format Affymetrix_QC_bl Olink_*_QC_24m Olink_*gwasQC_24m %15.0g
drop if Wes_gwasQC_bl=="" & Wgs_QC_bl=="" & Wgs_gwasQC_bl==""
replace Wes_gwasQC_bl="-999" if Wes_gwasQC_bl==""
replace Wgs_QC_bl="-999" if Wgs_QC_bl==""
replace Wgs_gwasQC_bl="-999" if Wgs_gwasQC_bl==""
outsheet Olink_*_QC_24m Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl sexPulse agePulse using work/weswgs.txt, noquote replace 
save work/weswgs, replace

program test
insheet using `1', case clear comma
format aliquot_id %15.0g
rename aliquot_id Olink_`2'_QC_24m
sort Olink_`2'_QC_24m
merge 1:m Olink_`2'_QC_24m using work/weswgs
drop plate* Olink_*_QC_24m Olink_*_gwasQC_24m
keep if _merge==3 & ethnicPulse!="Arab"
tab centre, gen(centre)
drop identifier Affymetrix_QC_bl _merge ethnicPulse centre
save work/`2', replace
outsheet using work/`2'.txt, comma noquote replace
outsheet Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl sexPulse agePulse ht_bl wt_bl centre* using work/`2'-covariates.txt,comma noquote replace
drop Wes_gwasQC_bl Wgs_QC_bl Wgs_gwasQC_bl sexPulse agePulse ht_bl wt_bl centre*
outsheet using work/`2'-protein.txt,comma noquote replace
end

test high_dimensional_data/Olink_proteomics_cvd2/qc/olink_qc_cvd2.csv cvd2
test high_dimensional_data/Olink_proteomics_cvd3/qc/olink_qc_cvd3.csv cvd3
test high_dimensional_data/Olink_proteomics_inf/qc/olink_qc_inf.csv inf
test high_dimensional_data/Olink_proteomics_neurology/qc/olink_neu_qc.csv neu
