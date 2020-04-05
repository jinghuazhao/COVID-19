# 5-4-2020 JHZ

library(GEOquery)

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

# GSM
gsm <- getGEO('GSM4115',GSEMatrix=TRUE)
Meta(gsm)
Columns(gsm)
Table(gsm)

# GSE
gse <- getGEO('GSE4115',GSEMatrix=FALSE)
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL96'},GSMList(gse))
probesets <- Table(GPLList(gse)[[1]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x)
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
require(Biobase)
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
