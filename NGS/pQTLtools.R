library(gap)

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

pval <- Sys.getenv("pval")
f <- paste0("NGS.",pval)
# Chrom Start End P prot MarkerName
sentinels <- read.table(f,as.is=TRUE,header=TRUE)
hits <- within(sentinels,
{
  panel_uniprot <- strsplit(prot,"_")
  Panel <- unlist(lapply(panel_uniprot,"[[",1))
  UniProt <- unlist(lapply(panel_uniprot,"[[",2))
  UniProt <- unlist(lapply(strsplit(UniProt,"-"),"[[",1))
  prot <- paste(Panel,UniProt,sep="_")
  SNP <- MarkerName
  chrpos <- strsplit(MarkerName,":")
  chr <- unlist(lapply(chrpos,"[[",1))
  Chr <- sub("chr","",chr)
  bp <- unlist(lapply(chrpos,"[[",2))
})[ c("prot","Chr","bp","SNP","UniProt")]

options(width=200)

hg <- within(subset(hgTables,!grepl("hap",X.chrom)&!grepl("Un",X.chrom)&!grepl("random",X.chrom)&!grepl(";",geneName)&geneName!=""),
{
  chr <- sub("chr","",X.chrom)
  start <- chromStart
  end  <- chromEnd
  gene <- geneName
})[c("chr","start","end","gene","UniProt")]
ngs <- within(merge(Olink_NGS,hg,by.x=c("UniProt","gene"),by.y=c("UniProt","gene")),{prot <- paste0(Panel,"_",UniProt)})
pcvt <- cis.vs.trans.classification(hits=hits, panel=ngs, id="UniProt")

pcvt$table

bioc()$table

library(iBMQ)
# UniProt
pqtl <- hits[c("UniProt","SNP", "prot")]
snp <- hits[c("SNP","Chr","bp")]
gene <- ngs[c("UniProt","chr","start","end")]
UniProt.type <- eqtlClassifier(pqtl, snp, gene, 1000000)
UniProt.table <- with(UniProt.type,table(UniProt,Type))

# Gene
all <- merge(ngs,hits,by="UniProt")
pqtl <- all[c("gene","SNP","Panel")]
gene <- ngs[c("gene","chr","start","end")]
gene.type <- eqtlClassifier(pqtl, snp, gene, 1000000)
gene.table <- with(gene.type,table(gene,Type))
sink(paste(pval,"cis.vs.trans.out",sep="/"))
UniProt.table
cat("Total",sum(UniProt.table[,1]),sum(UniProt.table[,2]),sum(UniProt.table[,3]),"\n")
sum(UniProt.table)
gene.table
cat("Total",sum(gene.table[,1]),sum(gene.table[,2]),"\n")
sum(gene.table)
sink()
