#!/bin/R

rm(list=ls())
closeAllConnections()

library(openxlsx); library(data.table)
wgs.ref <- fread("~/rds/rds-jmmh2-projects/interval_flagship/reference_files/klaudia_sanger_farm/INTERVAL_WGS_Sample_QC_09-12-2019.txt", header=T, colClasses="character", showProgress=F)
omics.map <- fread("~/rds/rds-jmmh2-projects/interval_flagship/data_management_files/batch4_updated/omicsMap.csv", header=T, sep=",", colClasses="character", showProgress=F)
dat.ref <- fread("~/rds/rds-jmmh2-projects/interval_flagship/data_management_files/batch4_updated/INTERVALdata_28FEB2020.csv", header=T, sep=",", colClasses="character", showProgress=F)	
pheno.ref <- omics.map[,.(identifier)]

##Add WGS variables

pheno.ref[,wgs_bl:=omics.map$Wgs_RAW_bl]
pheno.ref[,wgs_24m:=omics.map$Wgs_RAW_24m]

pheno.ref[wgs_bl!="" & wgs_24m=="", wgs_id:=wgs_bl]
pheno.ref[wgs_bl=="" & wgs_24m!="", wgs_id:=wgs_24m]
	
pheno.ref[,wgs_bl_pass:=wgs.ref[match(pheno.ref$wgs_bl, ID), PASS]]
pheno.ref[,wgs_24m_pass:=wgs.ref[match(pheno.ref$wgs_24m, ID), PASS]]

pheno.ref[wgs_bl!=wgs_24m & wgs_bl!="" & wgs_24m!="" & wgs_bl_pass=="1", wgs_id:=wgs_bl]
pheno.ref[wgs_bl!=wgs_24m & wgs_bl!="" & wgs_24m!="" & wgs_24m_pass=="1", wgs_id:=wgs_24m]

pheno.ref <- pheno.ref[,.(identifier, wgs_id)]
pheno.ref[wgs_id=="", wgs_id:=NA]

##Run this bit to ensure that the right matching variable is selected for merging

datname.check <- NULL
id.check <- pheno.ref$wgs_id
for(col.name in colnames(wgs.ref)) datname.check <- rbind(datname.check, data.table(name.ref=col.name, n.match=length(intersect(pheno.ref$wgs_id, as.character(unlist(wgs.ref[,col.name,with=F]))))))
datname.check <- datname.check[order(-n.match),]

pheno.ref <- cbind(pheno.ref, wgs.ref[match(pheno.ref$wgs_id, ID), .(phase=Study, qc.stat=PASS)])

cat("\nNotes:")
cat("\n\nThere are 35 IDs in the omicsMap 'Wgs_RAW_bl' column not present in Sanger manifest - 'INTERVAL_WGS_Sample_QC_09-12-2019.txt'\n\n")
tmp.list <- NULL; tmp.list <- setdiff(pheno.ref$wgs_id, wgs.ref$ID); tmp.list <- tmp.list[is.na(tmp.list)==F]
print(tmp.list)

pheno.ref[qc.stat!="1", wgs_id:=NA]
pheno.ref[is.na(qc.stat), wgs_id:=NA]
pheno.ref[,qc.stat:=NULL]

cat("\nProcessing data...")
cat("\n\tOlink_CVD2")
##Add Olink variables

olink <- NULL; olink <- fread("~/rds/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_qcgwas_cvd2.csv", header=T, sep=",", colClasses="character", showProgress=F)
tmp.dat <- NULL; tmp.dat <- omics.map[,.(identifier, id.merge=Olink_cvd2_gwasQC_24m)]
tmp.dat[id.merge=="", id.merge:=NA]

if(tmp.dat[!is.na(id.merge), length(unique(id.merge))]!=length(intersect(olink$aliquot_id, tmp.dat$id.merge))) stop("Some problem with IDs!!!")

tmp.dat <- cbind(tmp.dat, olink[match(tmp.dat$id.merge, aliquot_id), c(3:ncol(olink)), with=F]); tmp.dat[ tmp.dat=="" ] <- NA 
setnames(tmp.dat, paste0("Olink_CVD2_", colnames(tmp.dat)))

pheno.ref <- cbind(pheno.ref, tmp.dat[,3:ncol(tmp.dat), with=F])

cat("\n\tOlink_CVD3")
##Add Olink variables

olink <- NULL; olink <- fread("~/rds/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_qcgwas_cvd3.csv", header=T, sep=",", colClasses="character", showProgress=F)
tmp.dat <- NULL; tmp.dat <- omics.map[,.(identifier, id.merge=Olink_cvd3_gwasQC_24m)]
tmp.dat[id.merge=="", id.merge:=NA]

if(tmp.dat[!is.na(id.merge), length(unique(id.merge))]!=length(intersect(olink$aliquot_id, tmp.dat$id.merge))) stop("Some problem with IDs!!!")

tmp.dat <- cbind(tmp.dat, olink[match(tmp.dat$id.merge, aliquot_id), c(3:ncol(olink)), with=F]); tmp.dat[ tmp.dat=="" ] <- NA
setnames(tmp.dat, paste0("Olink_CVD3_", colnames(tmp.dat)))

pheno.ref <- cbind(pheno.ref, tmp.dat[,3:ncol(tmp.dat), with=F])

cat("\n\tOlink_CVD-INF")
##Add Olink variables

olink <- NULL; olink <- fread("~/rds/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_qcgwas_inf.csv", header=T, sep=",", colClasses="character", showProgress=F)
tmp.dat <- NULL; tmp.dat <- omics.map[,.(identifier, id.merge=Olink_inf_gwasQC_24m)]
tmp.dat[id.merge=="", id.merge:=NA]

if(tmp.dat[!is.na(id.merge), length(unique(id.merge))]!=length(intersect(olink$aliquot_id, tmp.dat$id.merge))) stop("Some problem with IDs!!!")

tmp.dat <- cbind(tmp.dat, olink[match(tmp.dat$id.merge, aliquot_id), c(3:ncol(olink)), with=F]); tmp.dat[ tmp.dat=="" ] <- NA
setnames(tmp.dat, paste0("Olink_INF_", colnames(tmp.dat)))

pheno.ref <- cbind(pheno.ref, tmp.dat[,3:ncol(tmp.dat), with=F])

cat("\n\tOlink_NEU")
##Add Olink variables

olink <- NULL; olink <- fread("~/rds/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_neu_qcgwas.csv", header=T, sep=",", colClasses="character", showProgress=F)
tmp.dat <- NULL; tmp.dat <- omics.map[,.(identifier, id.merge=Olink_neu_gwasQC_24m)]
tmp.dat[id.merge=="", id.merge:=NA]
	
if(tmp.dat[!is.na(id.merge), length(unique(id.merge))]!=length(intersect(olink$aliquot_id, tmp.dat$id.merge))) stop("Some problem with IDs!!!")

tmp.dat <- cbind(tmp.dat, olink[match(tmp.dat$id.merge, aliquot_id), c(3:ncol(olink)), with=F]); tmp.dat[ tmp.dat=="" ] <- NA
setnames(tmp.dat, paste0("Olink_NEU_", colnames(tmp.dat)))

pheno.ref <- cbind(pheno.ref, tmp.dat[,3:ncol(tmp.dat), with=F])

###Generate summary of numbers including the number of individuals with missing phenotypes
cat("\n\nGenerating a summary of the number of samples available for analyses...\n\n")	

##Add WES QC stat IDs directly from omicsMap file
pheno.ref[,wes_id:=omics.map$Wes_gwasQC_bl]
pheno.ref[wes_id=="", wes_id:=NA]

cat("\nCalculating mean/sd for all, wgs and wes data\n\n")

pheno.ref <- pheno.ref[!is.na(wgs_id) | !is.na(wes_id),]
summary.numbers <- NULL

summary.numbers <- rbind(summary.numbers,
	data.table(col.variable='wes_wgs_phenotype', pheno.ref[!is.na(wes_id) & !is.na(wgs_id), lapply(.SD, function(x) length(which(is.na(x)==F))), .SDcols=4:ncol(pheno.ref)]),
	data.table(col.variable='wes_phenotype', pheno.ref[!is.na(wes_id), lapply(.SD, function(x) length(which(is.na(x)==F))), .SDcols=4:ncol(pheno.ref)]),
	data.table(col.variable='wgs_phenotype', pheno.ref[!is.na(wgs_id), lapply(.SD, function(x) length(which(is.na(x)==F))), .SDcols=4:ncol(pheno.ref)]),
	data.table(col.variable='wgs_phenotype_phase1', pheno.ref[!is.na(wgs_id) & phase=="Phase1", lapply(.SD, function(x) length(which(is.na(x)==F))), .SDcols=4:ncol(pheno.ref)]),
	data.table(col.variable='wgs_phenotype_phase2', pheno.ref[!is.na(wgs_id) & phase=="Phase2", lapply(.SD, function(x) length(which(is.na(x)==F))), .SDcols=4:ncol(pheno.ref)]),
	data.table(col.variable='wgs_phenotype_phase3', pheno.ref[!is.na(wgs_id) & phase=="Phase3", lapply(.SD, function(x) length(which(is.na(x)==F))), .SDcols=4:ncol(pheno.ref)]))

summary.out <- NULL
summary.out <- dcast(melt(summary.numbers, id.vars = "col.variable"), variable ~ col.variable)
summary.out[,c('assay.tmp','tmp.var'):=tstrsplit(variable, split="_", fixed=T)[1:2]]

summary.out[,Assay:=assay.tmp]
summary.out[!grepl('NMR|Somalogic', assay.tmp), Assay:=paste(assay.tmp, tmp.var, sep="_")]
summary.out <- summary.out[Assay!="wes_id",]

summary.out[, wgs_all:=max(wgs_phenotype), by=Assay]
summary.out[, wes_all:=max(wes_phenotype), by=Assay]
summary.out[, wes_wgs_all:=max(wes_wgs_phenotype), by=Assay]

summary.out[, wgs_phenotype_missing:=wgs_all-wgs_phenotype]
summary.out[, wes_phenotype_missing:=wes_all-wes_phenotype]
summary.out[,Phenotype:=gsub(paste0(Assay, "_") ,"", variable),by=1:nrow(summary.out)]
summary.out <- summary.out[,.(Assay, Phenotype, wes_wgs_all, wes_all, wgs_all, wes_wgs_phenotype, wes_phenotype, wes_phenotype_missing, wgs_phenotype, wgs_phenotype_missing, wgs_phenotype_phase1, wgs_phenotype_phase2, wgs_phenotype_phase3)]

cat("Total number of samples with WGS data:", nrow(pheno.ref[!is.na(wgs_id),]))
cat("\nTotal number of samples with WES data:", nrow(pheno.ref[!is.na(wes_id),]))
cat("\nTotal number of samples with WGS and WES data:", nrow(pheno.ref[!is.na(wgs_id) & !is.na(wes_id),]))

cat("\n\nPath to sample summaries on CSD3: work/sample_size_summary_nmr_olink_somalogic.txt")
cat("\nPath to phenotype reference file on CSD3: work/interval_flagship_phenotype_nmr_olink_somalogic.txt")

pheno.ref <- pheno.ref[,c(grep('wes|wgs', colnames(pheno.ref)), grep('Somalogic|NMR|Olink', colnames(pheno.ref))), with=F]
colnames(pheno.ref)[match(c('wes_id','wgs_id'), colnames(pheno.ref))] <- c('ID_WES','ID_WGS')

write.table(summary.out, "work/sample_size_summary_olink.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(pheno.ref, "work/interval_flagship_phenotype_olink.txt", col.names=T, row.names=F, quote=F, sep="\t")

cat("\n\n")
