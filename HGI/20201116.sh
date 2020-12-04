#!/usr/bin/bash

function HGI20201116()
{
  for d in \
  20201116-male-ANA_C2_V2 \
  20201116-female-ANA_C2_V2 \
  20201116-le_60-ANA_C2_V2 \
  20201116-gt_60-ANA_C2_V2
  do
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}
    ln -s ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/work
    mkdir ~/rds/rds-asb38-ceu-restricted/projects/covid/HGI/${d}/output
  done
}

function rename()
# [dataset].[last name].[analysis_name].[freeze_number].[age].[sex].[ancestry].[n_cases].[n_controls].[gwas software].[YYYYMMDD].txt.gz
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
  export YYYYMMDD=20201123
  export filename=${dataset}.${lastname}.${analysisname}.${freezenumber}.${age}.${sex}.${ancestry}.${ncases}.${ncontrols}.${gwassoftware}.${YYYYMMDD}
  mv ${SCALLOP}/HGI/${d}.txt.gz ${SCALLOP}/HGI/${d}/output/${filename}.txt.gz
}

rename 20201116-female-ANA_C2_V2    C2_V2 ALL   F    161 41674
rename 20201116-male-ANA_C2_V2      C2_V2 ALL   M    161 41674
rename 20201116-le_60-ANA_C2_V2     C2_V2 LE_60 ALL  161 41674
rename 20201116-gt_60-ANA_C2_V2     C2_V2 GT_60 ALL  161 41674
