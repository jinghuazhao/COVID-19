// 12-5-2020 JHZ

local dir "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI"
insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_Covid_06MAY2020.csv", clear
sort identifier
save covid, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVALdata_06MAY2020.csv", clear
sort identifier
save INTERVAL, replace
insheet using "`dir'/06-05-2020/INTERVAL/INTERVAL_OmicsMap_20200506.csv", clear
sort identifier
merge 1:1 identifier using INTERVAL, gen(dataid)
merge 1:1 identifier using covid, gen(intervalcovid)
gzsave INTERVAL-covid, replace
