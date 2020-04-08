#!/usr/bin/bash

## INTERVAL data reformat
format.sh

## Meta-analysis
metal ACE2.metal 2>&1 | \
tee ACE2-1.tbl.log

## P values for 0.001 threshold and sentinels
select.sh ACE2
sentinels_nold.sh
merge.sh

## annotation
annotation.sh
for catalogue in GWAS eQTL
do
   export catalogue=${catalogue}
   R --no-save -q < ps.R
done

# QQ/Manhattan/LocusZoom/Forest plots
R --no-save -q < qml.R

INTERVAL.sh
