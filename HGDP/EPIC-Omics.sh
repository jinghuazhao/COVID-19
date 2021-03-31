#15-4-2015 MRC-Epid JHZ

# STEP1: Combine your dataset with reference dataset

mydata=/genetics/data/omics/EPICNorfolk/Axiom_UKB_EPICN_release_04Dec2014

plink2 --bfile $mydata --bmerge refdata.bed refdata.bim refdata.fam --keep test.dat --make-bed --out mydata.refdata
# ** If there is a strand issue, please flip (--flip in Plink) your SNPs based on the file of mydata.refdata.missnp created in STEP2 and redo this step. **

# STEP2: Select variants with geno 0.05 and MAF > 0.05

plink2 --bfile mydata.refdata --maf 0.05 --geno 0.01 --make-bed --out mydata.refdata.QCed

# STEP3: Make genome file

plink2 --bfile mydata.refdata.QCed --Z-genome --out mydata.refdata.QCed.Z
# ** this step may take very LONG! **

# STEP4: MDS

plink2 --bfile mydata.refdata.QCed --read-genome mydata.refdata.QCed.Z.genome.gz \
--cluster --mds-plot 10 --out mydata.refdata.QCed.MDS

# STEP5: MDS-PLOT

Rscript --slave mydata.refdata.QCed.MDS.mds refdata.fam MyOutput < MDSplot.R
