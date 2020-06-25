corr <- function()
{
  od <- read.csv(paste0("olink_proteomics/qc/olink_",qc[sel],".csv"),as.is=TRUE)
  names(od) <- toupper(names(od))
  Olink_id <- paste0("Olink_",ids[sel],"_QC_24m")
  omics <- read.csv("INTERVAL_OmicsMap_20200619.csv",as.is=TRUE)
  omics_id <- subset(omics,!is.na(omics[Olink_id]))[c(Olink_id,"Affymetrix_gwasQC_bl")]
  od <- merge(omics_id,od,by.x=Olink_id,by.y="ALIQUOT_ID")
  odd <- merge(d,od,by="Affymetrix_gwasQC_bl")
# correlations and scatter plots
  overlaps <- setdiff(intersect(names(d),names(od)),"Affymetrix_gwasQC_bl")
  pdf(paste0("work/",panel,"-",opt,".pdf"))
  for(i in overlaps)
  {
    with(odd, {
    x <- odd[[paste0(i,".x")]]
    y <- odd[[paste0(i,".y")]]
    cat(panel, "-", i, ":", cor(x,y,use="everything"), "\n")
    plot(x,y,main=paste0(i))
    })
  }
  dev.off()
}
# new data from NGS
panels <- c("NEUROLOGY", "CARDIOMETABOLIC", "INFLAMMATION")
qc <- c("neu_qc", "qc_cvd2", "qc_cvd3", "qc_inf")
ids <- c("neu","cvd2","cvd3","inf")
panel <- Sys.getenv("panel")
sel <- which(panel==panels)
d <-read.delim(paste0("work/",panel,".out"),as.is=TRUE)
ncols <- length(d[1,-1])/3
enum <- (1:ncols)*3-1
uniprot <- d[1,enum]
names(d)[enum+1] <- paste0(uniprot,"_lod")
names(d)[enum+2] <- uniprot
d <- d[,-enum]
write.table(d,file=paste0("work/",panel,".csv"),quote=FALSE,row.names=FALSE,sep=",")
# old Olink panels
cat("Panel =",panel,ids[sel],"\n")
corr()
if (sel==2) {sel=3; cat("Panel =",panel,ids[sel],"\n"); corr()}
