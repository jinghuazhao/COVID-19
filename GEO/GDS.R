# 5-4-2020 JHZ

library(GEOquery)
gds <- getGEO('GDS27771')
m <- Meta(gds)
Columns(gds)
t <- Table(gds)
eset <- within(GDS2eSet(gds,do.log2=TRUE), {y=2-as.numeric(disease.state)})
pData(eset)
gs <- c("1861_at","207187_at","203808_at","206254_at")
subset <- eset[gs,]
es <- exprs(subset)
lm(y ~ es["1861_at",]+es["207187_at",]+es["203808_at",]+es["206254_at",], data=subset)

