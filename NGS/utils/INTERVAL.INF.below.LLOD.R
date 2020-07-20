rm(list=ls())

require(Biobase)

x <- readRDS("~/post-doc/o-link/esets/round2/post-qc/eset.inf1.flag.out.outlier.in.rds")

get.prop.below.LLOD <- function(eset, flagged = 'OUT'){
  
  ## A function to calculate no. of proteins i.e NA per sample (missing or <LLD per sample)
  # arguments 'eset' and 'flagged'
  # flagged indicates whether Flagged samples should be excluded (if they have not been already)
  
  if(!inherits(eset, what= "ExpressionSet")){
    stop("'eset' argument must inherit class ExpressionSet")
  }  
  
  if (!flagged %in% c('IN','OUT')){
    stop("'flagged' argument must be 'IN' or 'OUT")
  }  
  
  require(stringr)
  
  # best to cut flagged samples first at eset stage:
  # risk of messing up if cutting from matrix, and then dont edit pData
  
  ind.fl <- which(eset$Flagged == 'Flagged')
  
  if (flagged == "IN"){
    
    if (length(ind.fl) > 0){
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (flagged retained)")
    } else{
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (no. flagged samples = 0)")
    }
    
  } else if (flagged == "OUT"){
    # cut flagged samples
    
    if (length(ind.fl) > 0){
      eset <- eset[, -ind.fl] # nb annoying ESet behaviour: cols are samples
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (flagged removed)")
    } else{
      panel <- unique(eset$panel)
      mytit <- paste(toupper(panel), "panel \n (no. flagged samples = 0)")
    }
    
  }
  
  E <- t(exprs(eset))
  
  p.annot <- fData(eset)
  
  p.annot$pc.belowLOD.new <- NA
  
  # % proteins in each sample
  ##miss.by.prot <- apply(E, 2, FUN=function(x) 100*sum(is.na(x))/nrow(E) )
  
  for (i in 1:ncol(E)){
    m <- sum( E[,i] <= p.annot$lod.max[i], na.rm=T ) # number of samples <= LLOD
    t <- length( na.omit(E[,i] )) # denominator NB use this rather than just nrow(E) to make code robust in event of missing values
    p.annot$pc.belowLOD.new[i] <- 100*m/t
  }
  
  fData(eset) <- p.annot 
    
  eset
  #eof  
}

x <- get.prop.below.LLOD(x)

annot <- fData(x)

annot$MissDataProp <- as.numeric( gsub("\\%$", "", annot$MissDataProp) )

pdf("~/post-doc/o-link/scallop/plots/INF.interval.LLOD.orig.vs.post.QC.pdf")
plot(annot$MissDataProp, annot$pc.belowLOD.new, xlab="% <LLOD in Original", ylab="% <LLOD in post QC dataset")
dev.off()

# get the proteins with no pQTL
np <- read.table("~/post-doc/o-link/scallop/data-inputs/proteins.with.no.pqtl.txt", head=F,
                 col.names = c("prot", "uniprot.id"))

annot$pQTL <- rep(NA, nrow(annot))

no.pQTL.ind <- which(annot$uniprot.id %in% np$uniprot)

annot$pQTL[no.pQTL.ind] <- "red"
annot$pQTL[-no.pQTL.ind] <- "blue"

annot <- annot[order(annot$pc.belowLOD.new, decreasing = T),]

annot <- annot[-grep("^BDNF$", annot$ID),] 

pdf("~/post-doc/o-link/scallop/plots/INF.interval.LLOD.post.QC.vs.pQTL.pdf", height=5, width=7)
par(mar=c(5,4,1,1))
plot(annot$pc.belowLOD.new, col=annot$pQTL, las=1, ylab = "% samples with very low abundance per protein", xaxt= 'n', xlab = "ordered proteins")
axis(side=1, at = seq(from=0, to=nrow(annot), by=20))
legend("topright", legend= c("no pQTL", "pQTL"), col= c("red","blue"), pch = 1)
dev.off()
