# 11-4-2020 JHZ

library(GEOquery)
library(Biobase)
library(limma)

# GDS
gds <- getGEO('GDS2771')
m <- Meta(gds)
Columns(gds)
t <- Table(gds)
eset <- GDS2eSet(gds,do.log2=TRUE)
eset$y <- 2-as.numeric(eset$disease.state)
pData(eset)
gs <- c("1861_at","207187_at","203808_at","206254_at")
subset <- eset[gs,]
es <- exprs(subset)
lm(y ~ es["1861_at",]+es["207187_at",]+es["203808_at",]+es["206254_at",], data=subset)
gpl <- getGEO(filename="GPL96.annot.gz")
MA <- GDS2MA(gds,GPL=gpl)
lmFit(MA$M)

# GSM -- lung tissue from wild type mouse
gsm <- getGEO('GSM4115',GSEMatrix=TRUE)
Meta(gsm)
Columns(gsm)
Table(gsm)

# GSE
gse <- getGEO('GSE4115',GSEMatrix=FALSE)
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL96'},GSMList(gse))
probesets <- Table(GPLList(gse)[[1]])$ID
d <- do.call('cbind',lapply(gsmlist,function(x) {
              t <- Table(x)
              m <- match(probesets,t$ID_REF)
              return(t$VALUE[m])
     }))
d <- apply(d,2,function(x) {as.numeric(as.character(x))})
d <- log2(d)
rownames(d) <- probesets
colnames(d) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset <- new('ExpressionSet',exprs=d,phenoData=pheno)
pData(phenoData(eset))

# GSE -- the heart data
gse <- getGEO('GSE106118',GSEMatrix=FALSE)
Meta(gse)

# scRNA
gse <- getGEO("GSM3489195", GSEMatrix=FALSE)
gse

# SRA

library(SRAdb)
library(DBI)
srafile = getSRAdbFile()
con = dbConnect(RSQLite::SQLite(), srafile)
listSRAfile('SRP026197', con)
getSRAfile('SRP026197',con,fileType='sra')
