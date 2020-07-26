library(gap)

options(width=200)
pval <- Sys.getenv("pval")
f <- paste0("NGS.",pval)
# Chrom Start End P prot MarkerName
sentinels <- read.table(f,as.is=TRUE,header=TRUE)
hits <- within(sentinels,
{
  panel_uniprot <- strsplit(prot,"_")
  UniProt <- unlist(lapply(panel_uniprot,"[[",2))
  SNP <- MarkerName
  chrpos <- strsplit(MarkerName,":")
  chr <- unlist(lapply(chrpos,"[[",1))
  Chr <- sub("chr","",chr)
  bp <- unlist(lapply(chrpos,"[[",2))
})[ c("prot","Chr","bp","SNP","UniProt")]

library(pQTLtools)
names(Olink_NGS)[2] <- "gene"
p1 <- subset(Olink_NGS,Panel=="NEUROLOGY")
p2 <- subset(Olink_NGS,Panel=="CARDIOMETABOLIC")
p3 <- subset(Olink_NGS,Panel=="ONCOLOGY")
p4 <- subset(Olink_NGS,Panel=="INFLAMMATION")
nrow(p1)+nrow(p2)+nrow(p3)+nrow(p4)
intersect(p1$UniProt,p2$UniProt)
subset(p1,UniProt%in%c("P01375", "P05231", "P10145"))
subset(p2,UniProt%in%c("P01375", "P05231", "P10145"))
subset(p3,UniProt%in%c("P01375", "P05231", "P10145"))
subset(p4,UniProt%in%c("P01375", "P05231", "P10145"))
p <- list(p1$UniProt,p2$UniProt,p3$UniProt,p4$UniProt)
cnames <- c("NEUROLOGY","CARDIOMETABOLIC","ONCOLOGY","INFLAMMATION")
VennDiagram::venn.diagram(x = p, category.names=cnames, filename='ngs.png', imagetype="png", output=TRUE)

hg <- within(hgTables,
{
  chr <- sub("chr","",X.chrom)
  start <- chromStart
  end  <- chromEnd
  gene <- geneName
})[c("chr","start","end","gene","UniProt")]
ngs <- within(merge(Olink_NGS,hg,by.x=c("UniProt","gene"),by.y=c("UniProt","gene")),{prot <- paste0(Panel,"_",UniProt)})

cvt <- cis.vs.trans.classification(hits=hits, panel=ngs, id="UniProt")
cvt
