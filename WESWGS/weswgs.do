//usr/local/Cluster-Apps/ceuadmin/stata/15/stata

insheet using INTERVALdata_28FEB2020.csv, case clear comma
sort identifier
save work/INTERVALdata_28FEB2020, replace
insheet using omicsMap.csv, case clear comma
sort identifier
merge 1:1 identifier using work/INTERVALdata_28FEB2020
d,f
keep identifier Affymetrix_QC_bl Wes_gwasQC_bl Wgs_gwasQC_bl Olink_* sexPulse agePulse ethnicPulse centre ht_bl wt_bl
mvencode _all, mv(-999)
format Wes_gwasQC_bl Wgs_gwasQC_bl %15s
format Affymetrix_QC_bl Olink_*_QC_24m Olink_*gwasQC_24m %15.0g
drop if Wes_gwasQC_bl=="" & Wgs_gwasQC_bl==""
replace Wes_gwasQC_bl="-999" if Wes_gwasQC_bl==""
replace Wgs_gwasQC_bl="-999" if Wgs_gwasQC_bl==""
outsheet Olink_*_QC_24m Wes_gwasQC_bl Wgs_gwasQC_bl sexPulse agePulse using work/weswgs.txt, noquote replace 
save work/weswgs, replace
