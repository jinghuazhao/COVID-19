R --no-save -q <<END
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
