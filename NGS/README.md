## The working directory

/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/NGS 

## Phenotype information

Restuls at work/ is according to

* NGS -- Olink (cvd2, cvd3, inf, neu) panels
* raw, QC, LOD, col1 for raw measurements, Olink NGS QC=PASS, LOD (set to be NA when < LOD), and column 1 (set to NA when 01;).

correlogram.pdf shows NGS -- Olink correlation

NGS-Olink-QC combination
------------------------
CARDIOMETABOLIC-cvd2-col1.pdf
CARDIOMETABOLIC-cvd2-LOD.pdf
CARDIOMETABOLIC-cvd2-QC.pdf
CARDIOMETABOLIC-cvd2-raw.pdf
CARDIOMETABOLIC-cvd3-col1.pdf
CARDIOMETABOLIC-cvd3-LOD.pdf
CARDIOMETABOLIC-cvd3-QC.pdf
CARDIOMETABOLIC-cvd3-raw.pdf
INFLAMMATION-inf-col1.pdf
INFLAMMATION-inf-LOD.pdf
INFLAMMATION-inf-QC.pdf
INFLAMMATION-inf-raw.pdf
NEUROLOGY-neu-col1.pdf
NEUROLOGY-neu-LOD.pdf
NEUROLOGY-neu-QC.pdf
NEUROLOGY-neu-raw.pdf
---------------------

The Olink NGS QC=PASS appears to have the highest correlation.

2. Genotype-trait association

plink2/ contains genotype-protein association results for MAF>=0.05, INFO>0.8.

NGS.merge is contained for each of the following p value thresholds,

Directory | Description
----------|------------
1e-5 | sentinel identification at p=1e-5
1e-6 | sentinel identification at p=1e-6
1e-7 | sentinel identification at p 1e-7
