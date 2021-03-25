/*
insheet using work/snpid.txt
save work/snpid
*/
gzuse work/INTERVAL-omics-covid.dta.gz
sort identifier
gzsave INTERVAL-omics-covid, replace
insheet using 20201201/INTERVALdata_01DEC2020.csv, case clear
sort identifier
gzmerge identifier using INTERVAL-omics-covid.dta.gz
table ethnicPulse
table ethnic_bl
table ethnicPulse ethnic_bl
gen str3 ethnic="NA"
replace ethnic="EUR" if inlist(ethnicPulse,"Eng/W/Scot/NI/Brit","White Irish")==1
replace ethnic="EAS" if inlist(ethnicPulse,"Asian- Bangladeshi","Asian- Indian","Asian- Pakistan","Chinese")==1
replace ethnic="MID" if ethnicPulse=="Arab"
gen FID=0
rename ID IID
outsheet FID IID ethnic using ethnic.txt if IID!=., noquote replace
rm INTERVAL-omics-covid.dta.gz
