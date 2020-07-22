## The working directory

/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/NGS 

1. Phenotype information

Results at work/ are according to the combination of,

* NGS panels -- CARDIOMETABOLIC, INFLAMMATION, NEUROLOGY.
* Olink panels -- cvd2, cvd3, inf, neu.
* QC types -- raw: raw measurements, QC: NGS QC=PASS, LOD: set to be NA when < LOD, and col1: set to NA when 01.

for a total of 16 NGS-Olink-QC combinations, with suffexes .dat for correlation and .pdf for scatter plots.

Group | Combinations
--------|-------------
1-4 | CARDIOMETABOLIC-cvd2-col1, CARDIOMETABOLIC-cvd2-LOD, CARDIOMETABOLIC-cvd2-QC, CARDIOMETABOLIC-cvd2-raw
5-8 | CARDIOMETABOLIC-cvd3-col1, CARDIOMETABOLIC-cvd3-LOD, CARDIOMETABOLIC-cvd3-QC, CARDIOMETABOLIC-cvd3-raw
9-12 | INFLAMMATION-inf-col1, INFLAMMATION-inf-LOD, INFLAMMATION-inf-QC, INFLAMMATION-inf-raw
13-16 | NEUROLOGY-neu-col1, NEUROLOGY-neu-LOD, NEUROLOGY-neu-QC, NEUROLOGY-neu-raw

The Olink NGS QC=PASS appears to have the highest correlation.

The corresponding density plots are shown in correlogram.pdf. 

2. Genotype-trait association

**plink2/** contains genotype-protein association results for MAF>=0.05, INFO>0.8.

NGS.merge is contained for each of the following p value thresholds,

Directory | Description | sentinels
----------|-------------|---------:
1e-5/ | sentinel identification at p=1e-5 | 23,034
1e-6/ | sentinel identification at p=1e-6 |  2,836
1e-7/ | sentinel identification at p 1e-7 |    598
