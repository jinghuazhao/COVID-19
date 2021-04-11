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
  export freezenumber=8
  export age=${3}
  export sex=${4}
  export ancestry=EUR
  export ncases=${5}
  export ncontrols=${6}
  export gwassoftware=SAIGE
  export YYYYMMDD=20201214
  export filename=${dataset}.${lastname}.${analysisname}.${freezenumber}.${age}.${sex}.${ancestry}.${ncases}.${ncontrols}.${gwassoftware}.${YYYYMMDD}
  mv ${SCALLOP}/HGI/${d}.txt.gz ${SCALLOP}/HGI/${d}/output/${filename}.txt.gz
  echo ${filename} | awk -F. '{print NF}'
  echo ${filename}.txt.gz
}

rename 20210317-ANA_C2_V2           C2_V2 ALL   ALL  838 40994
rename 20210317-female-ANA_C2_V2    C2_V2 ALL   F    460 20698
rename 20210317-male-ANA_C2_V2      C2_V2 ALL   M    378 20296
rename 20210317-le_60-ANA_C2_V2     C2_V2 LE_60 ALL  676 28144
rename 20210317-gt_60-ANA_C2_V2     C2_V2 GT_60 ALL  162 12850
