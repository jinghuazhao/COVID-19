#!/usr/bin/bash

# step 4. susceptibility association

module load gcc/6

function Cx_V2_X()
{
  step1_fitNULLGLMM.R \
     --plinkFile=${dir}/work/INTERVAL-covid-X \
     --phenoFile=${dir}/work/INTERVAL-covid-X.txt \
     --phenoCol=SARS_CoV \
     --covarColList=${covlist}
     --sampleIDColinphenoFile=ID \
     --traitType=binary \
     --outputPrefix=${dir}/output/INTERVAL-covid-X \
     --nThreads=8 \
     --IsOverwriteVarianceRatioFile=TRUE 

  step2_SPAtests.R \
     --vcfFile=${dir}/work/INTERVAL-X-ploidy.vcf.gz \
     --vcfFileIndex=${dir}/work/INTERVAL-X-ploidy.vcf.gz.tbi \
     --chrom=X \
     --minMAF=0.0001 \
     --minMAC=1 \
     --sampleFile=${dir}/work/INTERVAL-X.samples \
     --GMMATmodelFile=${dir}/output/INTERVAL-covid-X.rda \
     --varianceRatioFile=${dir}/output/INTERVAL-covid-X.varianceRatio.txt \
     --SAIGEOutputFile=${dir}/output/INTERVAL-X.txt \
     --IsOutputNinCaseCtrl=TRUE \
     --IsOutputHetHomCountsinCaseCtrl=TRUE \
     --IsOutputAFinCaseCtrl=TRUE
  gzip -f ${dir}/output/INTERVAL-X.txt
}

export d=20200731
export covlist=sex,age,age2,sexage,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20
for dir in ${d}-ANA_C1_V2 ${d}-ANA_C2_V2 \
           ${d}-male-ANA_C1_V2 ${d}-male-ANA_C2_V2 ${d}-female-ANA_C1_V2 ${d}-female-ANA_C2_V2 \
           ${d}-male-60-ANA_C1_V2 ${d}-male-60-ANA_C2_V2 ${d}-female-60-ANA_C1_V2 ${d}-female-60-ANA_C2_V2
do
  export dir=${dir}
  sbatch --wait ${SCALLOP}/HGI/autosomes.sb
  Cx_V2_X
done
