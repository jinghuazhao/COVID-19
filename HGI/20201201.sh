#!/usr/bin/bash

function HGI20201201()
{
  for d in \
  20201201-male-ANA_C2_V2 \
  20201201-female-ANA_C2_V2 \
  20201201-le_60-ANA_C2_V2 \
  20201201-gt_60-ANA_C2_V2
  do
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}
    ln -s ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/work
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/output
  done
}

function rename()
# [dataset].[last name].[analysis_name].[freeze_number].[age].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# INTERVAL.Zhao.ANA_A2_V2.9.ALL.F.EUR.460.20698.SAIGE.20201212.txt.gz
{
  export d=${1}
  export dataset=INTERVAL
  export lastname=Zhao
  export analysisname=ANA_${2}
  export freezenumber=8
  export age=${3}
  export sex=${4}
  export ancestry=EUR
  export ncases=${5}
  export ncontrols=${6}
  export gwassoftware=SAIGE
  export YYYYMMDD=20201212
  export filename=${dataset}.${lastname}.${analysisname}.${freezenumber}.${age}.${sex}.${ancestry}.${ncases}.${ncontrols}.${gwassoftware}.${YYYYMMDD}
  mv ${SCALLOP}/HGI/${d}.txt.gz ${SCALLOP}/HGI/${d}/output/${filename}.txt.gz
  echo ${filename} | awk -F. '{print NF}'
}

rename 20201201-female-ANA_C2_V2    C2_V2 ALL   F    460 20698
rename 20201201-male-ANA_C2_V2      C2_V2 ALL   M    378 20296
rename 20201201-le_60-ANA_C2_V2     C2_V2 LE_60 ALL  676 28144
rename 20201201-gt_60-ANA_C2_V2     C2_V2 GT_60 ALL  162 12850

. tab agegroup SARS_CoV

           |       SARS_CoV
  agegroup |         0          1 |     Total
-----------+----------------------+----------
         1 |    28,144        676 |    28,820
         2 |    12,850        162 |    13,012
-----------+----------------------+----------
     Total |    40,994        838 |    41,832


. tab sex SARS_CoV

           |       SARS_CoV
       sex |         0          1 |     Total
-----------+----------------------+----------
         1 |    20,296        378 |    20,674
         2 |    20,698        460 |    21,158
-----------+----------------------+----------
     Total |    40,994        838 |    41,832



