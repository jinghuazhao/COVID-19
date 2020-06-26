# NGS
options(digits=5, scipen=20, width=500)
opanel <- Sys.getenv("opanel")
panel <- Sys.getenv("panel")
qc_opanel <- Sys.getenv("qc_opanel")
opt <- Sys.getenv("opt")
rt <- paste0("work/",panel,"-",opanel,"-",opt)
d <-read.delim(paste0(rt,".out"),as.is=TRUE)
ncols <- length(d[1,-1])/3
enum <- (1:ncols)*3-1
uniprot <- d[1,enum]
names(d)[enum+1] <- paste0(uniprot,"_lod")
names(d)[enum+2] <- uniprot
d <- d[,-enum]
write.table(d,file=paste0(rt,".csv"),quote=FALSE,row.names=FALSE,sep=",")
# old Olink panels
od <- read.csv(paste0("olink_proteomics/qc/olink_",qc_opanel,".csv"),as.is=TRUE)
names(od) <- toupper(names(od))
Olink_id <- paste0("Olink_",opanel,"_QC_24m")
omics <- read.csv("INTERVAL_OmicsMap_20200619.csv",as.is=TRUE)
omics_id <- subset(omics,!is.na(omics[Olink_id]))[c(Olink_id,"Affymetrix_gwasQC_bl")]
od <- merge(omics_id,od,by.x=Olink_id,by.y="ALIQUOT_ID")
odd <- merge(d,od,by="Affymetrix_gwasQC_bl")
# correlations and scatter plots
overlaps <- setdiff(intersect(names(d),names(od)),"Affymetrix_gwasQC_bl")
load("work/ids.rda")
ids <- subset(ids, UniProt%in%overlaps)
pdf(paste0(rt,".pdf"))
cat("NGS","UniProt","Prot","r","\n",file=paste0(rt,".dat"))
for(i in overlaps) with(odd, {
  Prot <- ids[ids$UniProt==i,"Prot"]
  x <- odd[[paste0(i,".x")]]
  y <- odd[[paste0(i,".y")]]
  r <- round(cor(x,y,use="pairwise.complete.obs"),2)
  cat(i,r,"\n")
  print(cbind(x,y))
  cat(panel, i, Prot, r, "\n",append=TRUE,file=paste0(rt,".dat"))
  plot(x,y,main=paste0(i,"-",Prot," (r=",r,")"),xlab="Old panel",ylab="NGS")
})
dev.off()
