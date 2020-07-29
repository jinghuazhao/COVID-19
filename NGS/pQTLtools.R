bioc <- function()
{
  library(GenomicRanges)
  gr <- with(subset(hgTables,!grepl("hap",X.chrom)&!grepl("Un",X.chrom)&!grepl("random",X.chrom)&!grepl(";",geneName)&geneName!=""),
             GRanges(seqnames=X.chrom,ranges=IRanges(chromStart,chromEnd),strand=strand,gene=geneName,UniProt))
  rgr <- reduce(gr)
  m <- mergeByOverlaps(rgr,gr,type="within")
  ngs <- within(merge(Olink_NGS,as.data.frame(m),by.x=c("UniProt","gene"),by.y=c("UniProt","gene")),
  {
      prot <- paste0(Panel,"_",UniProt)
      chr <- sub("chr","",rgr.seqnames)
      start <- rgr.start
      end <- rgr.end
  })[c("prot","chr","start","end","UniProt","gene")]
  cvt <- cis.vs.trans.classification(hits=hits, panel=ngs, id="UniProt")
  d <- with(cvt, data)
  print(d[is.na(d["p.gene"]),c("prot","UniProt","Chr","bp","SNP","p.gene")])
  hr <- with(hits, GRanges(seqnames=paste0("chr",Chr),ranges=IRanges(as.numeric(bp)-1,as.numeric(bp))))
  radius <- 1000000
  fm <- flank(gr,radius,both=TRUE)
  findOverlapPairs(hr,fm)
  cvt
}

library(gap)
library(pQTLtools)
library(iBMQ)
options(width=220)

hg <- within(subset(hgTables,!grepl("hap",X.chrom)&!grepl("Un",X.chrom)&!grepl("random",X.chrom)&!grepl(";",geneName)&geneName!=""),
{
  chr <- sub("chr","",X.chrom)
  start <- chromStart
  end  <- chromEnd
})[c("chr","start","end","geneName","UniProt")]
subset(hg,grepl("Q14213",UniProt))
subset(hg,grepl("Q29980",UniProt))
subset(hg,grepl("P21217",UniProt))
subset(hg,grepl("O43521",UniProt))
# does not exist
subset(hg,grepl("P16284",UniProt))
hg <- rbind(hg,data.frame(chr="17", start=62396776, end=62407083, geneName="PECAM1", UniProt="P16284"))
names(Olink_NGS)[2] <- "gene"
subset(Olink_NGS,grepl("Q14213",UniProt))
subset(Olink_NGS,grepl("Q29980",UniProt))
subset(Olink_NGS,grepl("P21217",UniProt))
subset(Olink_NGS,grepl("O43521",UniProt))
Olink_NGS[with(Olink_NGS,UniProt=="Q14213_Q8NEV9"),"UniProt"] <- "Q14213"
Olink_NGS[with(Olink_NGS,UniProt=="Q14213"),"gene"] <- "EBI3"
Olink_NGS[with(Olink_NGS,UniProt=="Q29980_Q29983"),"UniProt"] <- "Q29980"
Olink_NGS[with(Olink_NGS,UniProt=="Q29980"),"gene"] <- "MICB"
Olink_NGS[with(Olink_NGS,UniProt=="P21217_Q11128"),"UniProt"] <- "P21217"
Olink_NGS[with(Olink_NGS,UniProt=="P21217"),"gene"] <- "PUT3"
Olink_NGS[with(Olink_NGS,UniProt=="O43521-2"),"UniProt"] <- "O43521"
Olink_NGS[with(Olink_NGS,UniProt=="O43521"),"gene"] <- "BCL2L11"
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

pval <- Sys.getenv("pval")
f <- paste0("NGS.",pval)
sentinels <- within(read.table(f,as.is=TRUE,header=TRUE),
{
# Chrom Start End P prot MarkerName
  panel_uniprot <- strsplit(prot,"_")
  Panel <- unlist(lapply(panel_uniprot,"[[",1))
  UniProt <- unlist(lapply(panel_uniprot,"[[",2))
  UniProt <- unlist(lapply(strsplit(UniProt,"-"),"[[",1))
  prot <- paste0(Panel,"_",UniProt)
  SNP <- MarkerName
  chrpos <- strsplit(MarkerName,":")
  chr <- unlist(lapply(chrpos,"[[",1))
  Chr <- sub("chr","",chr)
  bp <- as.integer(unlist(lapply(chrpos,"[[",2)))
})[ c("prot","Chr","bp","SNP","UniProt")]
hits <- merge(sentinels,Olink_NGS,by="UniProt")

ngs <- within(merge(hg,Olink_NGS,by="UniProt"),{prot <- paste0(Panel,"_",UniProt)})
pqtl <- hits[c("gene","SNP","prot")]
snp <- hits[c("SNP","Chr","bp")]
gene <- ngs[c("gene","chr","start","end")]
gene.type <- eqtlClassifier(pqtl, snp, gene, 1000000)
gene.table <- with(gene.type,table(gene,Type))
sink(paste(pval,"cis.vs.trans.out",sep="/"))
row.names(gene.type) <- 1:nrow(gene.type)
gene.type
gene.table
cat("Total",sum(gene.table[,1]),sum(gene.table[,2]),"\n")
sum(gene.table)
sink()
cvt <- cis.vs.trans.classification(hits=hits[c("Chr","bp","gene","prot")], panel=ngs[c("chr","start","end","gene","prot")], id="prot")
cvt
