function NPX_Inflammation()
{
R --no-save <<END
  xlsx <- "high_dimensional_data/Olink_proteomics_inf/rawData/20161007_Danesh_NPX_Inflammation.xlsx"
  library(openxlsx)
  readNPX <- function(sheet)
  {
    prop <- function(col)
    {
      prot <- data[,col]
      lod <- lodmax[col]
      below <- (prot < as.double(lod))
      plod[col-1] <<- length(prot[below&!is.na(prot)])/length(prot[!is.na(prot)])
      pnan[col-1] <<- length(prot[is.na(prot)])/length(prot)
    }
    data <- read.xlsx(xlsx, sheet=sheet, colNames=TRUE, skipEmptyRows=FALSE, cols=1:94, rows=3:5121)
    uniprot <- read.xlsx(xlsx, sheet=sheet, colNames=FALSE, skipEmptyRows=FALSE, cols=1:94, rows=1)
    prot <- read.xlsx(xlsx, sheet=sheet, colNames=FALSE, skipEmptyRows=FALSE, cols=1:94, rows=2)
    miprop <- read.xlsx(xlsx, sheet=sheet, colNames=FALSE, skipEmptyRows=FALSE, cols=1:94, rows=5122)
    lodmax <- read.xlsx(xlsx, sheet=sheet, colNames=FALSE, skipEmptyRows=FALSE, cols=1:94, rows=5123)
    names(data) <- uniprot
    data <- subset(data,substr(data[,1],1,7)!="Control")
    plod <- pnan <- numeric()
    sapply(seq(2, ncol(lodmax)), prop)
    names(plod) <- names(pnan) <- uniprot[c(2:93)]
    print(as.character(uniprot))
    print(as.character(miprop))
    print(substr(as.character(lodmax),1,4))
    print(plod)
    d <- data.frame(miprop=as.numeric(miprop[-1]),plod=plod,cols="black")
    d[names(plod)%in%with(nosig,uniprot),"cols"] <- "red"
    plodcol <- merge(data.frame(uniprot=substr(rownames(d),1,6),d),inf1,by="uniprot")
    xtick <- 1:91
    with(subset(plodcol,uniprot!="P23560"),{
      print(wilcox.test(plod~as.factor(cols)))
      plot(miprop,plod,cex=0.8,col=cols,main=sheet,pch=19,xlab="MissingDataFreq(%)",ylab="pLOD")
      plot(xtick,plod,xaxt="n",cex=0.8,col=cols,main=sheet,pch=19,xlab="",ylab="pLOD")
      axis(side=1,at=xtick,labels=prot,las=2,xpd=TRUE,cex.axis=0.3)
    })
    list(data=data,uniprot=uniprot,prot=prot,miprop=miprop,lod=lodmax,plod=plod,pnan=pnan)
  }
  options(width=200)
  INF <- Sys.getenv("INF")
  inf1 <- read.table(paste(INF,"work","inf1.tmp",sep="/"),col.names=c("prot","uniprot"),as.is=TRUE)
  nosig <- read.table(paste(INF,"work","INF1.merge.nosig",sep="/"),col.names=c("prot","uniprot"),as.is=TRUE)
  pdf("scallop-inf-lod.pdf")
  par(mfrow=c(3,2))
  NPX_NaN <- readNPX("NPX_NaN")
  NPX_LOD <- readNPX("NPX_LOD")
  NPX_Keep <- readNPX("NPX_Keep")
  dev.off()
  pdf("scallop-inf-NPX_Keep.pdf")
  NPX_Keep <- readNPX("NPX_Keep")
  dev.off()
END
}

function INTERVAL.INF.below.LLOD()
# INTERVAL.INF.below.LLOD.R from Jimmy
{
R --no-save -q <<\ \ END
  rm(list=ls())

  require(Biobase)

  x <- readRDS("olink_proteomics/post-qc/eset.inf1.flag.out.outlier.in.rds")

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

  pdf("INF.interval.LLOD.orig.vs.post.QC.pdf")
  plot(annot$MissDataProp, annot$pc.belowLOD.new, xlab="% <LLOD in Original", ylab="% <LLOD in post QC dataset", pch=19)
  dev.off()

  # get the proteins with no pQTL
  INF <- Sys.getenv("INF")
  np <- read.table(paste(INF, "work", "INF1.merge.nosig", sep="/"), head=F,
                   col.names = c("prot", "uniprot.id"))

  annot$pQTL <- rep(NA, nrow(annot))

  no.pQTL.ind <- which(annot$uniprot.id %in% np$uniprot)

  annot$pQTL[no.pQTL.ind] <- "red"
  annot$pQTL[-no.pQTL.ind] <- "blue"

  annot <- annot[order(annot$pc.belowLOD.new, decreasing = T),]

  annot <- annot[-grep("^BDNF$", annot$ID),] 

  pdf("INF.interval.LLOD.post.QC.vs.pQTL.pdf", height=5, width=7)
  par(mar=c(5,4,1,1))
  plot(annot$pc.belowLOD.new, col=annot$pQTL, las=1, ylab = "% samples with very low abundance per protein", xaxt= 'n', xlab = "ordered proteins", cex=0.8, pch=19)
  axis(side=1, at = seq(from=0, to=nrow(annot), by=20))
  legend("topright", legend= c("no pQTL", "pQTL"), col= c("red","blue"), pch = 1)
  dev.off()
  END
}

NPX_Inflammation
INTERVAL.INF.below.LLOD
