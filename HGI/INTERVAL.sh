#!/usr/bin/bash

source INTERVAL.inc

# Analysis 5-susceptibility (phenotype name: ANA5):

# Single-variant association tests

function svat()
{
step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL-covid \
   --phenoFile=work/INTERVAL-covid.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=output/INTERVAL-covid \
   --nThreads=8 \
   --IsOverwriteVarianceRatioFile=TRUE

seq 22 | \
parallel -j1 --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=$autosomes/imputed/impute_{}_interval.bgen \
   --bgenFileIndex=$autosomes/imputed/impute_{}_interval.bgen.bgi \
   --chrom={} \
   --minMAF=0.0001 \
   --minMAC=1 \
   --sampleFile=work/INTERVAL.samples \
   --idstoIncludeFile=work/INTERVAL-covid.samples \
   --GMMATmodelFile=output/INTERVAL-covid.rda \
   --varianceRatioFile=output/INTERVAL-covid.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.txt \
   --IsOutputNinCaseCtrl=TRUE \
   --IsOutputHetHomCountsinCaseCtrl=TRUE \
   --IsOutputAFinCaseCtrl=TRUE
'
}

# Gene-based association tests

function gbat()
{
createSparseGRM.R \
   --plinkFile=work/INTERVAL \
   --minMAF=0.0001 \
   --nThreads=8 \
   --outputPrefix=output/INTERVAL.sparseGRM \
   --numRandomMarkerforSparseKin=2000 \
   --relatednessCutoff=0.125

step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL \
   --phenoFile=work/INTERVAL-covid.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=output/INTERVAL \
   --outputPrefix_varRatio=output/INTERVAL \
   --nThreads=8 \
   --LOCO=FALSE \
   --skipModelFitting=FALSE \
   --IsSparseKin=TRUE \
   --isCateVarianceRatio=TRUE \
   --IsOverwriteVarianceRatioFile=TRUE

seq 22 | \
parallel -j1 --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=work/INTERVAL-{}.bgen \
   --bgenFileIndex=work/INTERVAL-{}.bgen.bgi \
   --chrom={} \
   --minMAF=0 \
   --minMAC=1 \
   --maxMAFforGroupTest=0.001 \
   --sampleFile=work/INTERVAL.samples \
   --idstoIncludeFile=work/INTERVAL-autosomes.samples \
   --GMMATmodelFile=output/INTERVAL.rda \
   --varianceRatioFile=output/INTERVAL.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.gene.txt \
   --numLinesOutput=1 \
   --sparseSigmaFile=output/INTERVAL.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   --idstoIncludeFile=work/INTERVAL-{}.incl \
   --groupFile=work/INTERVAL-{}.gene \
   --IsOutputAFinCaseCtrl=TRUE \
   --IsSingleVarinGroupTest=TRUE \
   --IsOutputPvalueNAinGroupTestforBinary=TRUE \
   --IsAccountforCasecontrolImbalanceinGroupTest=TRUE
'
}

svat
