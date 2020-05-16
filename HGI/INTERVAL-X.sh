# 15-5-2020 JHZ

source INTERVAL.inc

grep -v -f work/INTERVAL-X.excl-samples work/INTERVAL-covid.txt > work/INTERVAL-covid-X.txt
step1_fitNULLGLMM.R \
   --plinkFile=work/INTERVAL \
   --phenoFile=work/INTERVAL-covid-X.txt \
   --phenoCol=SARS_CoV \
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --outputPrefix=output/INTERVAL-X \
   --nThreads=4

echo X | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=work/INTERVAL-{}.bgen \
   --bgenFileIndex=work/INTERVAL-{}.bgen.bgi \
   --idstoExcludeFile=work/INTERVAL-X.excl-rsids \
   --minMAF=0.0001 \
   --minMAC=1 \
   --sampleFile=work/INTERVAL-X.samples \
   --GMMATmodelFile=output/INTERVAL-X.rda \
   --varianceRatioFile=output/INTERVAL-X.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.txt \
   --IsOutputNinCaseCtrl=TRUE \
   --IsOutputHetHomCountsinCaseCtrl=TRUE \
   --IsOutputAFinCaseCtrl=TRUE
'

# Gene-based association tests

paste -d' ' work/INTERVAL-X.samples work/INTERVAL-X.samples > work/INTERVAL-X.samples2
plink2 --bfile work/INTERVAL --make-bed -keep work/INTERVAL-X.samples2 --out work/INTERVAL-X

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
   --covarColList=age,sex,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
   --sampleIDColinphenoFile=ID \
   --traitType=binary \
   --invNormalize=TRUE \
   --outputPrefix=output/INTERVAL-X \
   --outputPrefix_varRatio=output/INTERVAL-X_varRatio \
   --sparseGRMFile=output/INTERVAL-X.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   --sparseGRMSampleIDFile=output/INTERVAL-X.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
   --nThreads=8 \
   --LOCO=FALSE \
   --skipModelFitting=FALSE \
   --IsSparseKin=TRUE \
   --isCateVarianceRatio=TRUE \
   --IsOverwriteVarianceRatioFile=TRUE

echo X | \
parallel --env autosomes -C' ' '
step2_SPAtests.R \
   --bgenFile=work/INTERVAL-{}.bgen \
   --bgenFileIndex=work/INTERVAL-{}.bgen.bgi \
   --minMAF=0 \
   --minMAC=0.5 \
   --maxMAFforGroupTest=0.01 \
   --sampleFile=work/INTERVAL-X.samples \
   --GMMATmodelFile=output/INTERVAL-X.rda \
   --varianceRatioFile=output/INTERVAL-X.varianceRatio.txt \
   --SAIGEOutputFile=output/INTERVAL-{}.SAIGE.gene.txt \
   --numLinesOutput=1 \
   --sparseSigmaFile=output/INTERVAL-X.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
   --idstoIncludeFile=work/INTERVAL-X.incl \
   --groupFile=work/INTERVAL-X.gene \
   --IsOutputAFinCaseCtrl=TRUE \
   --IsSingleVarinGroupTest=TRUE \
   --IsOutputPvalueNAinGroupTestforBinary=TRUE \
   --IsAccountforCasecontrolImbalanceinGroupTest=TRUE
'
