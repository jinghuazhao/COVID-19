#!/usr/bin/bash

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
  done
}

function rename()
# [dataset].[last name].[analysis_name].[freeze_number].[age].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
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

rename 20200731-female-60-ANA_C1_V2 C1_V2 LE_60 F    98  702
rename 20200731-female-60-ANA_C2_V2 C2_V2 LE_60 F    98  18665
rename 20200731-female-ANA_C1_V2    C1_V2 ALL   F    100 730
rename 20200731-female-ANA_C2_V2    C2_V2 ALL   F    100 21059
rename 20200731-male-60-ANA_C1_V2   C1_V2 LE_60 M    57  287
rename 20200731-male-60-ANA_C2_V2   C2_V2 LE_60 M    57  16925
rename 20200731-male-ANA_C1_V2      C1_V2 ALL   M    61  380
rename 20200731-male-ANA_C2_V2      C2_V2 ALL   M    61  20615
rename 20200731-ANA_C1_V2           C1_V2 ALL   ALL  161 1119
rename 20200731-ANA_C2_V2           C2_V2 ALL   ALL  161 41674
