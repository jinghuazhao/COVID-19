#!/usr/bin/bash

export HGI=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
export TMPDIR=${HPC_WORK}/work

module load plink/2.00-alpha

# STEP1: Combine your dataset with reference dataset

mydata=${HGI}/INTERVAL
refdata=${HPC_WORK}/VEGAS2/keep

plink2 --pfile $mydata --make-bed --out mydata
plink --bfile mydata --bmerge ${refdata}.bed ${refdata}.bim ${refdata}.fam --make-bed --out mydata.refdata
# ** If there is a strand issue, please flip (--flip in Plink) your SNPs based on the file of mydata.refdata.missnp created in STEP2 and redo this step. **

# STEP2: Select variants with geno 0.05 and MAF > 0.05

plink --bfile mydata.refdata --maf 0.05 --geno 0.01 --make-bed --out mydata.refdata.QCed

# STEP3: Make genome file

sbatch --job-name=_Z-genome --account CARDIO-SL0-CPU --partition cardio --qos=cardio --mem=40800 --time=5-00:00:00 --export ALL \
       --output=${TMPDIR}/_Z-genome_%A_%a.out --error=${TMPDIR}/_Z-genome_%A_%a.err --wait --wrap ". keep.sb"

# plink --bfile mydata.refdata.QCed --Z-genome --out mydata.refdata.QCed.Z
# ** this step may take very LONG! **

# STEP4: MDS

plink --bfile mydata.refdata.QCed --read-genome mydata.refdata.QCed.Z.genome.gz \
--cluster --mds-plot 10 --out mydata.refdata.QCed.MDS

# STEP5: MDS-PLOT

Rscript MDSplot.R mydata.refdata.QCed.MDS.mds ${refdata}.fam MyOutput
