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
tabulate ethnicPulse SARS_CoV, all
gen str20 ethnic=ethnicPulse
replace ethnic="EUR" if inlist(ethnicPulse,"Eng/W/Scot/NI/Brit","White Irish")==1
replace ethnic="EAS" if inlist(ethnicPulse,"Asian- Bangladeshi","Asian- Indian","Asian- Pakistani","Chinese")==1
replace ethnic="MID" if ethnicPulse=="Arab"
gen ethnic_NA=ethnic
replace ethnic_NA="NA" if inlist(ethnic,"EUR","EAS","MID")==0
tab ethnic SARS_CoV, all
tab ethnic_NA SARS_CoV, all
gen FID=0
rename ID IID
outsheet FID IID ethnic using ethnic.txt if IID!=., noquote replace
outsheet FID IID ethnic_NA using ethnic_NA.txt if IID!=., noquote replace
rm INTERVAL-omics-covid.dta.gz
