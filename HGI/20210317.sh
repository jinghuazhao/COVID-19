#!/usr/bin/bash

function HGI20210317()
{
  for d in \
  20210317-ANA_C2_V2 \
  20210317-male-ANA_C2_V2 \
  20210317-female-ANA_C2_V2 \
  20210317-le_60-ANA_C2_V2 \
  20210317-gt_60-ANA_C2_V2
  do
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/work
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/output
    ln -s ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}
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
  export freezenumber=9
  export age=${3}
  export sex=${4}
  export ancestry=EUR
  export ncases=${5}
  export ncontrols=${6}
  export gwassoftware=SAIGE
  export YYYYMMDD=20210419
  export filename=${dataset}.${lastname}.${analysisname}.${freezenumber}.${age}.${sex}.${ancestry}.${ncases}.${ncontrols}.${gwassoftware}.${YYYYMMDD}
  mv ${SCALLOP}/HGI/${d}.txt.gz ${SCALLOP}/HGI/${d}/output/${filename}.txt.gz
  echo ${filename} | awk -F. '{print NF}'
  echo ${filename}.txt.gz
}

rename 20210317-ANA_C2_V2           C2_V2 ALL   ALL  2098 39733
rename 20210317-female-ANA_C2_V2    C2_V2 ALL   F    1168 19989
rename 20210317-male-ANA_C2_V2      C2_V2 ALL   M     930 19744
rename 20210317-le_60-ANA_C2_V2     C2_V2 LE_60 ALL  1713 27106
rename 20210317-gt_60-ANA_C2_V2     C2_V2 GT_60 ALL   385 12627

#   SARS_CoV |      Freq.     Percent        Cum.
#------------+-----------------------------------
#          0 |     39,733       94.98       94.98
#          1 |      2,098        5.02      100.00
#------------+-----------------------------------
#      Total |     41,831      100.00
#
#           |       SARS_CoV
#       sex |         0          1 |     Total
#-----------+----------------------+----------
#         1 |    19,744        930 |    20,674
#         2 |    19,989      1,168 |    21,157
#-----------+----------------------+----------
#     Total |    39,733      2,098 |    41,831
#
#           |       SARS_CoV
#  agegroup |         0          1 |     Total
#-----------+----------------------+----------
#         1 |    27,106      1,713 |    28,819
#         2 |    12,627        385 |    13,012
#-----------+----------------------+----------
#     Total |    39,733      2,098 |    41,831
