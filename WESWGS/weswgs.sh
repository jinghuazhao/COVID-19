#!/usr/bin/bash

module load ceuadmin/stata

stata -b do weswgs.do

(
for panel in cvd2 cvd3 inf neurology
do
  echo ${panel} | \
  tr '\n' '\t'
  export f=$(ls high_dimensional_data/Olink_proteomics_${panel}/qc)
  sed '1d' high_dimensional_data/Olink_proteomics_${panel}/qc/${f} | \
  cut -d, -f1 | \
  grep -f - work/weswgs.txt | \
  wc -l
done
) | xsel -i

# cvd2	1297
# cvd3	1305
# inf	1290
# neurology	1264

