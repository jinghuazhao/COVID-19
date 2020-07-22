## The working directory

/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/NGS 

## Phenotype information

Results at work/ is according to

NGS (CARDIOMETABOLIC, INFLAMMATION, NEUROLOGY) -- Olink (cvd2, cvd3, inf, neu) panels -- QC (raw, QC, LOD, col1 for raw measurements, Olink NGS QC=PASS, LOD (set to be NA when < LOD), and column 1 (set to NA when 01;).

There are 16 NGS - Olink - QC combinations shown in correlogram.pdf for NGS -- Olink correlation (.dat) and scatter plots (.pdf),

* CARDIOMETABOLIC-cvd2-col1, CARDIOMETABOLIC-cvd2-LOD, CARDIOMETABOLIC-cvd2-QC, CARDIOMETABOLIC-cvd2-raw
* CARDIOMETABOLIC-cvd3-col1, CARDIOMETABOLIC-cvd3-LOD, CARDIOMETABOLIC-cvd3-QC, CARDIOMETABOLIC-cvd3-raw
* INFLAMMATION-inf-col1, INFLAMMATION-inf-LOD, INFLAMMATION-inf-QC, INFLAMMATION-inf-raw
* NEUROLOGY-neu-col1, NEUROLOGY-neu-LOD, NEUROLOGY-neu-QC, NEUROLOGY-neu-raw

The Olink NGS QC=PASS appears to have the highest correlation.

2. Genotype-trait association

**plink2/** contains genotype-protein association results for MAF>=0.05, INFO>0.8.

NGS.merge is contained for each of the following p value thresholds,

Directory | Description | sentinels
----------|-------------|---------:
1e-5/ | sentinel identification at p=1e-5 | 23,034
1e-6/ | sentinel identification at p=1e-6 |  2,836
1e-7/ | sentinel identification at p 1e-7 |    598
