olink <- function(panel)
{
  if (panel=="neu") {
     x <- within(read.table("olink_proteomics/post_QC/interval_olink_neuro_post_qc_data_20180309.txt",as.is=TRUE,header=TRUE),
                {Aliquot_Id <- as.character(Aliquot_Id)})
     renameList <- !names(x)=="Aliquot_Id"
     names(x)[renameList] <- paste(panel,names(x)[renameList],sep="_")
  } else {
     require(Biobase)
     prefix <- "olink_proteomics/post-qc/eset"
     suffix <- "flag.out.outlier.in.rds"
     f <- paste(prefix,panel,suffix,sep=".")
     eset <- readRDS(f)
     fd <- within(as(featureData(eset),"data.frame")[c("ID","uniprot.id")],{
                     ID <- as.character(ID)
                     fdi <- paste0(paste(panel,ID,sep="_"),"__",uniprot.id)
           }) %>% select(-uniprot.id)
    cb <-  cbind(rownames(eset),fd)
    rownames(eset) <- cb[,3]
    x <- as(eset,"data.frame")
  }
  x
}

options(width=200)
require(plyr)
library(dplyr)
cvd2 <- olink("cvd2")
cvd3 <- olink("cvd3")
inf1 <- olink("inf1")
neu <- olink("neu")

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
pca_wes <- read.table("work/wes.eigenvec",col.names=c("z","id",paste0("PC",1:20))) %>% select(-z)
pca_wgs <- read.table("work/wgs.eigenvec",col.names=c("z","id",paste0("PC",1:20))) %>% select(-z)
id_wes <- subset(idmap,wes_id%in%wes)
id_wgs <- subset(idmap,wgs_id%in%wgs)
overlap <- subset(id_wgs,wes_id==wgs_id)
subset(id_wes,wes_id%in%with(overlap,wes_id))
id_wes <- id_wes %>% select(-c(phase,wgs_id))
id_wgs <- id_wgs %>% select(-c(phase,wes_id))
weswgs <- read.delim("work/weswgs.txt")

panels <- function(d,weswgs_id,pca)
{
  d <- d %>%
  rename(Aliquot_Id=Olink_CVD2_id.merge) %>% left_join(cvd2[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(Aliquot_Id=Olink_CVD3_id.merge) %>% left_join(cvd3[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(Aliquot_Id=Olink_INF_id.merge) %>% left_join(inf1[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(Aliquot_Id=Olink_NEU_id.merge) %>% left_join(neu[c(7,53:144)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(id=weswgs_id) %>% mutate(average = rowMeans(select(., starts_with("cvd2_BMP.6__P22004")), na.rm = TRUE)) %>%
  left_join(weswgs[c("identifier","sexPulse","agePulse")]) %>% mutate(age2=agePulse*agePulse) %>%
  left_join(pca)
  rownames(d) <- d[["id"]]
  d
}

y_wes <- panels(id_wes,"wes_id",pca_wes)
y_wgs <- panels(id_wgs,"wgs_id",pca_wgs)

library(gap)
normalize_sapply <- function(d)
{
  normfun <- function(col,verbose=FALSE)
  {
    if (verbose) cat(names(d[col]),col,"\n")
    y <- invnormal(d[[col]])
    l <- lm(y~average+sexPulse+agePulse+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=d[covars])
    r <- y-predict(l,na.action=na.pass)
    invnormal(r)
  }
  proteins <- grep("cvd2|cvd3|inf1|neu",names(d))
  covars <- c(names(d)[grep("average|sex|age",names(d))],paste0("PC",1:20))
  z <- sapply(names(d[proteins]), normfun)
  colnames(z) <- names(d[proteins])
  rownames(z) <- d[["id"]]
  data.frame(id=d[["id"]],z)
}
y_wes_sapply <- normalize_sapply(y_wes)
y_wgs_sapply <- normalize_sapply(y_wgs)
write.table(y_wes_sapply,file="work/wes.pheno",quote=FALSE,row.names=FALSE,sep="\t")
write.table(y_wgs_sapply,file="work/wgs.pheno",quote=FALSE,row.names=FALSE,sep="\t")

require(doMC)
doMC::registerDoMC(cores = 14)
normalize_adply <- function(d)
{
  proteins <- names(d)[grep("cvd2|cvd3|inf1|neu",names(d))]
  covars <- c(names(d)[grep("average|sex|age",names(d))],paste0("PC",1:20))
  z <- adply(d[proteins], 2, function(x)
       {
          y <- invnormal(x)
          l <- lm(y~average+sexPulse+agePulse+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=d[covars])
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

examine <- function()
{
  identical(y_wgs_sapply,y_wgs_adply)
  identical(y_wgs_sapply,y_wgs_adply)

  check <- y_wes[grep("cvd2_BMP.6__P22004|average|sex|age|PC",names(y_wes))]
  check <- within(check,
           {
             y <- invnormal(cvd2_BMP.6__P22004)
             r <- y-predict(lm(y~average+sexPulse+agePulse+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20),
                            na.action=na.pass)
             b <- invnormal(r)
           }) %>% select(-names(check)[grep("sex|age|PC",names(check))])
  a <- y_wes_sapply["cvd2_BMP.6__P22004"]
  plot(cbind(a,check$b),pch=19)
  head(cbind(a,check$b))
# overlaps
  s <- subset(idmap[,c("wes_id","wgs_id")],wes_id==wgs_id)
  o <- with(s,order(wes_id))
  share <- s[o,]
# earlier version
  f <- "work/interval_flagship_phenotype_nmr_olink_somalogic.txt"
  d <- read.delim(f)
  s <- subset(d[,c("ID_WGS","ID_WES")],ID_WGS==ID_WES) %>% rename(id_wes=ID_WES,id_wgs=ID_WGS)
  o <- with(s,order(id_wes))
  overlap <- s[o,]
  cbind(share,overlap)
}

examine()

output <- function(weswgs,d)
{
  if(!dir.exists(weswgs)) dir.create(weswgs)
  n_vars <- ncol(d)
  invisible(sapply(2:n_vars,function(x)
            {
              print(names(d[x]))
              write.table(d[c(1,x)],file=file.path(weswgs,paste0(names(d[x]),".pheno")),quote=FALSE,row.names=FALSE,sep="\t")
            })
  )
}

output("work/wes",y_wes_sapply)
output("work/wgs",y_wgs_sapply)
