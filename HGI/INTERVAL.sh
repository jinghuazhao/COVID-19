# 12-5-2020 JHZ

cd 06-05-2020/INTERVAL
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVALdata_06MAY2020.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVALdata_P3_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
cd -
export dir=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export ev=$ref/annot_INT_50PCs_pcs.txt
export snpstats=$ref/impute_${chr}_interval.snpstats

module load ceuadmin/stata
stata -b do INTERVAL.do
