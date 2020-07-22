## The working directory

/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/NGS 

### Phenotype information

Results at work/ are according to the combination of,

* NGS panels -- CARDIOMETABOLIC, INFLAMMATION, NEUROLOGY.
* qPCR panels -- cvd2, cvd3, inf, neu.
* QC types -- raw: raw measurements, QC: NGS QC=PASS, LOD: set to be NA when < LOD, and col1: set to NA when 01.

for a total of 16 NGS-Olink-QC combinations, with suffexes .dat for correlation and .pdf for scatter plots.

qPCR panel | Combinations
--------|----------------
cvd2 | CARDIOMETABOLIC-cvd2-col1[-LOD|-QC|-raw]
cvd3 | CARDIOMETABOLIC-cvd3-col1[-LOD|-QC|-raw]
inf | INFLAMMATION-inf-col1[-LOD|-QC|-raw]
neu | NEUROLOGY-neu-col1[-LOD|-QC|-raw]

The Olink NGS QC=PASS appears to have the highest correlation.

The corresponding density plots are shown in correlogram.pdf. 

### Genotype-trait association

**plink2/** contains genotype-protein association results for MAF>=0.05, INFO>0.8.

NGS.merge is contained for each of the following p value thresholds,

Directory | Description | sentinels
----------|-------------|---------:
1e-5/ | sentinel identification at p=1e-5 | 23,034
1e-6/ | sentinel identification at p=1e-6 |  2,836
1e-7/ | sentinel identification at p 1e-7 |    598
