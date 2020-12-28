olink <- function(panel)
{
  if (panel=="neu") x <- within(read.table("olink_proteomics/post_QC/interval_olink_neuro_post_qc_data_20180309.txt",as.is=TRUE,header=TRUE),
                               {Aliquot_Id <- as.character(Aliquot_Id)})
  else
  {
    require(Biobase)
    prefix <- "olink_proteomics/post-qc/eset"
    suffix <- "flag.out.outlier.in.rds"
    f <- paste(prefix,panel,suffix,sep=".")
    x <- as(readRDS(f),"data.frame")
  }
  renameList <- !names(x)=="Aliquot_Id"
  names(x)[renameList] <- paste(panel,names(x)[renameList],sep="_")
  x
}
cvd2 <- olink("cvd2")
cvd3 <- olink("cvd3")
inf1 <- olink("inf1")
neu <- olink("neu")

options(width=200)
require(plyr)
library(dplyr)
y <- cvd2[c(1:92,104)] %>%
     full_join(cvd3[c(1:92,104)],by="Aliquot_Id") %>%
     full_join(inf1[c(1:92,104)],by="Aliquot_Id") %>%
     full_join(neu[c(7,53:144)],by="Aliquot_Id")
dim(y)

f <- "WGS-WES-Olink_ID_map_INTERVAL_release_28FEB2020.txt"
idmap <- within(read.delim(f),
         {
           Olink_CVD2_id.merge <- as.character(Olink_CVD2_id.merge)
           Olink_CVD3_id.merge <- as.character(Olink_CVD3_id.merge)
           Olink_INF_id.merge <- as.character(Olink_INF_id.merge)
           Olink_NEU_id.merge <- as.character(Olink_NEU_id.merge)
         })
# WES/WGS samples
wes <- scan("work/wes.samples",what="")
wgs <- scan("work/wgs.samples",what="")
id_wes <- subset(idmap,wes_id%in%wes)
id_wgs <- subset(idmap,wgs_id%in%wgs)
overlap <- subset(id_wgs,wes_id==wgs_id)
subset(id_wes,wes_id%in%with(overlap,wes_id))
weswgs <- read.delim("work/weswgs.txt")

y_wes <- id_wes %>% select(-c(phase,wgs_id)) %>%
         rename(Aliquot_Id=Olink_CVD2_id.merge) %>% left_join(cvd2[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         rename(Aliquot_Id=Olink_CVD3_id.merge) %>% left_join(cvd3[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         rename(Aliquot_Id=Olink_INF_id.merge) %>% left_join(inf1[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         rename(Aliquot_Id=Olink_NEU_id.merge) %>% left_join(neu[c(7,53:144)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         left_join(weswgs[c("identifier","sexPulse","agePulse")])
rownames(y_wes) <- y_wes[["wes_id"]]

y_wgs <- id_wgs %>% select(-c(phase,wes_id)) %>%
         rename(Aliquot_Id=Olink_CVD2_id.merge) %>% left_join(cvd2[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         rename(Aliquot_Id=Olink_CVD3_id.merge) %>% left_join(cvd3[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         rename(Aliquot_Id=Olink_INF_id.merge) %>% left_join(inf1[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         rename(Aliquot_Id=Olink_NEU_id.merge) %>% left_join(neu[c(7,53:144)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
         left_join(weswgs[c("identifier","sexPulse","agePulse")])
rownames(y_wgs) <- y_wgs[["wgs_id"]]

library(gap)
normalize_sapply <- function(d,id)
{
  normfun <- function(col,verbose=FALSE)
  {
    if (verbose) cat(names(d[col]),col,"\n")
    y <- invnormal(d[[col]])
    l <- lm(y~sexPulse+agePulse,data=sexage)
    r <- y-predict(l,na.action=na.pass)
    invnormal(r)
  }
  proteins <- grep("cvd2|cvd3|inf1|neu",names(d))
  sexage <- d[grep("sex|age",names(d))]
  z <- sapply(names(d[proteins]), normfun)
  colnames(z) <- names(d[proteins])
  rownames(z) <- d[[id]]
  data.frame(id=d[[id]],z)
}
y_wes_sapply <- normalize_sapply(y_wes,"wes_id")
y_wgs_sapply <- normalize_sapply(y_wgs,"wgs_id")

require(doMC)
doMC::registerDoMC(cores = 14)
normalize_adply <- function(d)
{
  proteins <- names(d)[grep("cvd2|cvd3|inf1|neu",names(d))]
  sexage <- d[grep("sex|age",names(d))]
  z <- adply(d[proteins], 2, function(x)
       {
          y <- invnormal(x)
          l <- lm(y~sexPulse+agePulse,data=sexage)
          r <- y-predict(l,na.action=na.pass)
          invnormal(r)
       }, .progress = "none", .parallel = TRUE)
  tnames <- z[,1]
  z <- t(z[,-1])
  colnames(z) <- tnames
  rownames(z) <- rownames(d)
  data.frame(id=rownames(d),z)
}
y_wes_adply <- normalize_adply(y_wes)
y_wgs_adply <- normalize_adply(y_wgs)

identical(y_wgs_sapply,y_wgs_adply)
identical(y_wgs_sapply,y_wgs_adply)

a <- y_wes_sapply["cvd2_BMP.6"]
l <- lm(invnormal(cvd2_BMP.6)~sexPulse+agePulse,data=y_wes)
r <- y_wes["cvd2_BMP.6"]-predict(l,na.action=na.pass)
b <- invnormal(r)
plot(cbind(a,b))
head(cbind(a,b))

## overlaps

s <- subset(idmap[,c("wes_id","wgs_id")],wes_id==wgs_id)
o <- with(s,order(wes_id))
share <- s[o,]

f <- "work/interval_flagship_phenotype_nmr_olink_somalogic.txt"
d <- read.delim(f)
s <- subset(d[,c("ID_WGS","ID_WES")],ID_WGS==ID_WES)
o <- with(s,order(ID_WES))
overlap <- s[o,]
share==overlap

test <- function(d,id)
{
  normfun <- function(col,verbose=TRUE)
  {
    if (verbose) cat(col,names(y_wes[proteins])[col],"\n")
    y <- d[,col]
    l <- lm(y~sexPulse+agePulse,data=d[c("sexPulse","agePulse")])
    r <- y-predict(l,na.action=na.pass)
   invnormal(r)
  }
  z <- sapply(proteins, normfun)
  colnames(z) <- names(d)[proteins]
  rownames(z) <- d[[id]]
  z
}
y_wes_test <- test(y_wes[c(proteins,sexage)],"wes_id")
y_wgs_test <- test(y_wgs[c(proteins,sexage)],"wgs_id")
