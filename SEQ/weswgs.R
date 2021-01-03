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

# WES/WGS samples
weswgs <- read.delim("work/weswgs.txt")
wes <- scan("work/wes.idmap",what="")
wgs <- scan("work/wgs.idmap",what="")
pca_wes <- read.table("work/wes.eigenvec",col.names=c("FID","id",paste0("PC",1:20))) %>% select(-FID)
pca_wgs <- read.table("work/wgs.eigenvec",col.names=c("FID","id",paste0("PC",1:20))) %>% select(-FID)

## Stata version
idmap <- function(f,weswgs_id)
# the upper bound is already set for all panels
{
  d <- within(read.delim(f), {
         Olink_cvd2_QC_24m <- as.character(Olink_cvd2_QC_24m)
         Olink_cvd3_QC_24m <- as.character(Olink_cvd3_QC_24m)
         Olink_inf_QC_24m <- as.character(Olink_inf_QC_24m)
         Olink_neu_QC_24m <- as.character(Olink_neu_QC_24m)
       })
  names(d)[2] <- "id"
  subset(d,id%in%weswgs_id)
}

id_wes <- idmap("work/wes.txt",wes) %>% select(-c(Olink_cvd2_gwasQC_24m, Olink_cvd3_gwasQC_24m, Olink_inf_gwasQC_24m, Olink_neu_gwasQC_24m))
id_wgs <- idmap("work/wgs.txt",wgs) %>% select(-c(Olink_cvd2_gwasQC_24m, Olink_cvd3_gwasQC_24m, Olink_inf_gwasQC_24m, Olink_neu_gwasQC_24m))

panels <- function(d,weswgs_id,pca)
{
  d <- d %>%
  rename(Aliquot_Id=Olink_cvd2_QC_24m) %>% left_join(cvd2[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(Aliquot_Id=Olink_cvd3_QC_24m) %>% left_join(cvd3[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(Aliquot_Id=Olink_inf_QC_24m) %>% left_join(inf1[c(1:92,104)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(Aliquot_Id=Olink_neu_QC_24m) %>% left_join(neu[c(7,53:144)],by="Aliquot_Id") %>% select(-Aliquot_Id) %>%
  rename(id=weswgs_id) %>% mutate(average = rowMeans(select(., starts_with("cvd2_BMP.6__P22004")), na.rm = TRUE)) %>%
  left_join(weswgs[c("identifier","sex","age","age2")]) %>% 
  left_join(pca)
  rownames(d) <- d[["id"]]
  d
}

y_wes <- panels(id_wes,"id",pca_wes)
y_wgs <- panels(id_wgs,"id",pca_wgs)

library(gap)
normalize_sapply <- function(d)
{
  normfun <- function(col,verbose=FALSE)
  {
    if (verbose) cat(names(d[col]),col,"\n")
    y <- invnormal(d[[col]])
    l <- lm(y~average+sex+age+age2, data=d[covars])
    r <- y-predict(l,na.action=na.pass)
    invnormal(r)
  }
  proteins <- grep("cvd2|cvd3|inf1|neu",names(d))
  covars <- names(d)[grep("average|sex|age",names(d))]
  z <- sapply(names(d[proteins]), normfun)
  colnames(z) <- names(d[proteins])
  rownames(z) <- d[["id"]]
  data.frame(id=d[["id"]],z)
}
y_wes_sapply <- normalize_sapply(y_wes)
y_wgs_sapply <- normalize_sapply(y_wgs)
wes.id <- y_wes_sapply %>% select(id) %>% rename(FID=id) %>% mutate(IID=FID)
wgs.id <- y_wgs_sapply %>% select(id) %>% rename(FID=id) %>% mutate(IID=FID)
wes.pheno <- y_wes_sapply %>% select(-id)
wgs.pheno <- y_wgs_sapply %>% select(-id)
write.table(data.frame(wes.id,wes.pheno),file="work/wes.pheno",quote=FALSE,row.names=FALSE,sep="\t")
write.table(data.frame(wgs.id,wgs.pheno),file="work/wgs.pheno",quote=FALSE,row.names=FALSE,sep="\t")

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

require(doMC)
doMC::registerDoMC(cores = 14)
normalize_adply <- function(d)
{
  proteins <- names(d)[grep("cvd2|cvd3|inf1|neu",names(d))]
  covars <- names(d)[grep("average|sex|age",names(d))]
  z <- adply(d[proteins], 2, function(x)
       {
          y <- invnormal(x)
          l <- lm(y~average+sex+age+age2, data=d[covars])
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

  check <- y_wes[grep("cvd2_BMP.6__P22004|average|sex|age",names(y_wes))]
  check <- within(check,
           {
             y <- invnormal(cvd2_BMP.6__P22004)
             r <- y-predict(lm(y~average+sex+age+age2), na.action=na.pass)
             b <- invnormal(r)
           }) %>% select(-names(check)[grep("sex|age",names(check))])
  a <- y_wes_sapply["cvd2_BMP.6__P22004"]
  identical(a[["cvd2_BMP.6__P22004"]],check$b)
}

examine()

y <- cvd2[c(1:92,104)] %>%
     full_join(cvd3[c(1:92,104)],by="Aliquot_Id") %>%
     full_join(inf1[c(1:92,104)],by="Aliquot_Id") %>%
     full_join(neu[c(7,53:144)],by="Aliquot_Id")
dim(y)

praveen <- function()
{
## Praveen's version, not fully incompatible with PCs
  f <- "WGS-WES-Olink_ID_map_INTERVAL_release_28FEB2020.txt"
  praveen <- within(read.delim(f),
             {
                Olink_CVD2_id.merge <- as.character(Olink_CVD2_id.merge)
                Olink_CVD3_id.merge <- as.character(Olink_CVD3_id.merge)
                Olink_INF_id.merge <- as.character(Olink_INF_id.merge)
                Olink_NEU_id.merge <- as.character(Olink_NEU_id.merge)
             })
  id_wes <- subset(praveen,wes_id%in%wes)
  id_wgs <- subset(praveen,wgs_id%in%wgs)
  overlap <- subset(id_wgs,wes_id==wgs_id)
  subset(id_wes,wes_id%in%with(overlap,wes_id))
  id_wes <- id_wes %>% select(-c(phase,wgs_id))
  id_wgs <- id_wgs %>% select(-c(phase,wes_id))
}

test <- function(d)
{
  proteins <- grep("cvd2|cvd3|inf1|neu",names(d))
  covars <- d[grep("average|sex|age",names(d))]
  normfun <- function(col,verbose=FALSE)
  {
    if (verbose) cat(col,names(d[col]),"\n")
    y <- d[,col]
    l <- lm(y~average+sex+age+age2, data=covars)
    r <- y-predict(l,na.action=na.pass)
    invnormal(r)
  }
  z <- sapply(proteins, normfun)
  colnames(z) <- names(d[proteins])
  rownames(z) <- rownames(d)
  data.frame(id=rownames(d),z)
}
y_wes_test <- test(y_wes)
y_wgs_test <- test(y_wgs)
