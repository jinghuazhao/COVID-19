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
  cd ${d}
  sbatch ${SCALLOP}/HGI/20200731-autosomes.sb
  cd -
  done
}
