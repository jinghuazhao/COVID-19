# 29-3-2020 JHZ

source("test.ini")

library(dplyr)
library(Seurat)

for (src in c("Cryobiopsy","Donor","IPF","Myositis"))
{
  lst <- scan(paste0(src,".list"),"")
  h5mm(lst)
  for(i in 1:length(lst) investigate(lst[i])
}
