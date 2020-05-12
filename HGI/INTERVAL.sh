# 12-5-2020 JHZ

cd 06-05-2020/INTERVAL
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVALdata_06MAY2020.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVAL_Covid_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
sed '1d' INTERVALdata_P3_06MAY2020.csv | cut -d',' -f1 | sort | join - <(sed '1d' INTERVAL_OmicsMap_p3_20200506.csv | cut -d',' -f1) | wc -l
cd -
export dir=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HGI
export autosomes=/home/jhz22/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/imputed
export X=/rds/project/jmmh2/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data
export ref=/home/jhz22/rds/post_qc_data/interval/reference_files/genetic/interval
export ev=$ref/annot_INT_50PCs_pcs.txt
export snpstats=$ref/impute_${chr}_interval.snpstats

module load ceuadmin/stata
stata -b do INTERVAL.do

# ${autosomes}/impute_${chr}_interval.bgen.bgi
# ${X}/INTERVAL_X_imp_ann_filt_v2.vcf.gz

# Rscript createSparseGRM.R --help
# Rscript step1_fitNULLGLMM.R --help
# Rscript step2_SPAtests.R --help

# single-variant association tests

Rscript step1_fitNULLGLMM.R \
   --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
   --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
   --phenoCol=y_binary \
   --covarColList=x1,x2 \
   --sampleIDColinphenoFile=IID \
   --traitType=binary \
   --outputPrefix=./output/example_binary \
   --nThreads=4

Rscript step2_SPAtests.R \
   --bgenFile=./input/genotype_100markers.bgen \
   --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
   --minMAF=0.0001 \
   --minMAC=1 \
   --sampleFile=./input/samplefileforbgen_10000samples.txt \
   --GMMATmodelFile=./output/example.rda \
   --varianceRatioFile=./output/example.varianceRatio.txt \
   --SAIGEOutputFile=./output/example.SAIGE.bgen.txt \
   --numLinesOutput=2 \
   --IsOutputNinCaseCtrl=TRUE \
   --IsOutputHetHomCountsinCaseCtrl=TRUE \
   --IsOutputAFinCaseCtrl=TRUE

# Gene-based association tests

createSparseGRM.R \
   --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
   --nThreads=4 \
   --outputPrefix=./output/sparseGRM \
   --numRandomMarkerforSparseKin=2000 \
   --relatednessCutoff=0.125

Rscript step2_SPAtests.R \
   --vcfFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav \
   --vcfFileIndex=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav.s1r \
   --vcfField=DS \
   --chrom=chr1 \
   --minMAF=0 \
   --minMAC=0.5 \
   --maxMAFforGroupTest=0.01       \
   --sampleFile=./input/samplelist.txt \
   --GMMATmodelFile=./output/example_binary.rda \
   --varianceRatioFile=./output/example_binary_cate_v2.varianceRatio.txt \
   --SAIGEOutputFile=./output/example_binary.SAIGE.gene.txt \
   --numLinesOutput=1 \
   --groupFile=./input/groupFile_geneBasedtest.txt \
   --sparseSigmaFile=./output/example_binary_cate_v2.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx \
   --IsOutputAFinCaseCtrl=TRUE \
   --IsSingleVarinGroupTest=TRUE \
   --IsOutputPvalueNAinGroupTestforBinary=TRUE \
   --IsAccountforCasecontrolImbalanceinGroupTest=TRUE
