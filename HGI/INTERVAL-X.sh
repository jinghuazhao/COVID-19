#!/usr/bin/bash

# step 4. association analysis on chromosome X variants

source INTERVAL.inc

step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL-covid-X \
   --phenoFile=work/INTERVAL-covid-X.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,age2,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=output/INTERVAL-covid-X \
   --nThreads=8 \
   --IsOverwriteVarianceRatioFile=TRUE 

step2_SPAtests.R \
   --vcfFile=work/INTERVAL-X-ploidy.vcf.gz \
   --vcfFileIndex=work/INTERVAL-X-ploidy.vcf.gz.tbi \
   --chrom=X \
   --minMAF=0.0001 \
   --minMAC=1 \
   --sampleFile=work/INTERVAL-X.samples \
   --GMMATmodelFile=output/INTERVAL-covid-X.rda \
   --varianceRatioFile=output/INTERVAL-covid-X.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-X.txt \
   --IsOutputNinCaseCtrl=TRUE \
   --IsOutputHetHomCountsinCaseCtrl=TRUE \
   --IsOutputAFinCaseCtrl=TRUE

# Gene-based association tests

createSparseGRM.R \
   --plinkFile=work/INTERVAL-X \
   --minMAF=0.0001 \
   --nThreads=8 \
   --outputPrefix=output/INTERVAL-X.sparseGRM \
   --numRandomMarkerforSparseKin=2000 \
   --relatednessCutoff=0.125

step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL-X \
   --phenoFile=work/INTERVAL-covid-X.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,age2,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=output/INTERVAL-X \
   --outputPrefix_varRatio=output/INTERVAL-X \
   --nThreads=8 \
   --LOCO=FALSE \
   --skipModelFitting=FALSE \
   --IsSparseKin=TRUE \
   --isCateVarianceRatio=TRUE \
   --IsOverwriteVarianceRatioFile=TRUE

step2_SPAtests.R \
   --vcfFile=work/INTERVAL-X-ploidy.vcf.gz \
   --vcfFileIndex=work/INTERVAL-X-ploidy.vcf.gz.tbi \
   --chrom=X \
   --minMAF=0 \
   --minMAC=1 \
   --maxMAFforGroupTest=0.001 \
   --sampleFile=work/INTERVAL-X.samples \
   --GMMATmodelFile=output/INTERVAL-X.rda \
   --varianceRatioFile=output/INTERVAL-X.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-X.gene.txt \
   --numLinesOutput=1 \
   --sparseSigmaFile=output/INTERVAL-X.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   --idstoIncludeFile=work/INTERVAL-X.incl \
   --groupFile=work/INTERVAL-X.gene \
   --IsOutputAFinCaseCtrl=TRUE \
   --IsSingleVarinGroupTest=TRUE \
   --IsOutputPvalueNAinGroupTestforBinary=TRUE \
   --IsAccountforCasecontrolImbalanceinGroupTest=TRUE
