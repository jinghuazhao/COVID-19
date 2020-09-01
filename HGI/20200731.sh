#!/usr/bin/bash

function mvupload()
{
  export d=${1}
  export dataset=INTERVAL
  export lastname=Zhao
  export analysisname=ANA_${2}
  export freezenumber=5
  export age=${3]
  export sex=${4}
  export ancestry=EUR
  export ncases=${5}
  export ncontrols=${6}
  export gwassoftware=SAIGE
  export YYYYMMDD=20200901
  export filename=${dataset}.${lastname}.${analysisname}.${freezenumber}.${age}.${sex}.${ancestry}.${ncases}.${ncontrols}.${gwassoftware}.${YYYYMMDD}
  mv ${SCALLOP}/HGI/${d}.txt.gz ${SCALLOP}/HGI/${d}/output/${filename}.txt.gz
}

function HGI20200731()
{
  for d in \
  20200731-female-60-ANA_C1_V2 \
  20200731-female-60-ANA_C2_V2 \
  20200731-female-ANA_C1_V2 \
  20200731-female-ANA_C2_V2 \
  20200731-male-60-ANA_C1_V2 \
  20200731-male-60-ANA_C2_V2 \
  20200731-male-ANA_C1_V2 \
  20200731-male-ANA_C2_V2 \
  20200731-ANA_C1_V2 \
  20200731-ANA_C2_V2
  do
#   ln -s /rds/user/jhz22/rds-asb38-ceu-restricted/projects/covid/HGI/${d}
#   mkdir /rds/user/jhz22/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/work
#   mkdir /rds/user/jhz22/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/output
  cd ${d}
# [dataset].[last name].[analysis_name].[freeze_number].[age].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
# export snvResults=output/INTERVAL.Zhao.ANA_C2_V2.5.ALL.EUR.144.612.SAIGE.20200617.txt.gz
  mvupload 20200731-female-60-ANA_C1_V2 C1_V2 LE_60 FEMALE 98  702
  mvupload 20200731-female-60-ANA_C2_V2 C2_V2 LE_60 FEMALE 98  18665
  mvupload 20200731-female-ANA_C1_V2    C1_V2 ALL   FEMALE 100 730
  mvupload 20200731-female-ANA_C2_V2    C2_V2 ALL   FEMALE 100 21059
  mvupload 20200731-male-60-ANA_C1_V2   C1_V2 LE_60 MALE   57  287
  mvupload 20200731-male-60-ANA_C2_V2   C2_V2 LE_60 MALE   57  16925
  mvupload 20200731-male-ANA_C1_V2      C1_V2 ALL   MALE   61  380
  mvupload 20200731-male-ANA_C2_V2      C2_V2 ALL   MALE   61  20615
  mvupload 20200731-ANA_C1_V2           C1_V2 ALL   ALL    161 1119
  mvupload 20200731-ANA_C2_V2           C2_V2 ALL   ALL    161 41674
  done
}
