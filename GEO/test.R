# 29-3-2020 JHZ

source("test.ini")

library(dplyr)
library(Seurat)

labels <- c("Cryobiopsy","Donor","IPF","Myositis")
for (src in labels)
{
  lst <- scan(paste0(src,".list"),"")
  h5mm(lst)
  for(i in 1:length(lst)) investigate(src,lst[i])
}
