# 29-3-2020 JHZ

library(rhdf5)
library(data.table)
library(Matrix)
m <- function(rt,type="_filtered",suffix="_gene_bc_matrices_h5.h5")
{
  f <- paste0(rt,type,suffix)
  h5ls(f)
  h5readAttributes(f,"GRCh38")
  d <- h5read(f,name="GRCh38", read.attributes = TRUE)
  start <- h5read(f, "/GRCh38/indptr")
  dt <- data.table(
     row = h5read(f, "/GRCh38/indices", start = 1, count = tail(start, 1)) + 1,
     column = rep(seq_len(length(start) - 1), diff(start)),
     count = h5read(f, "/GRCh38/data", start = 1, count = tail(start, 1))
  )
  h5closeAll()
  list(d=d,dt=dt)
}


rt <-c(
"GSM3489182_Donor_01",
"GSM3489182_Donor_01",
"GSM3489183_IPF_01",
"GSM3489183_IPF_01",
"GSM3489184_IPF_02",
"GSM3489184_IPF_02",
"GSM3489185_Donor_02",
"GSM3489185_Donor_02",
"GSM3489186_Cryobiopsy_01",
"GSM3489186_Cryobiopsy_01",
"GSM3489187_Donor_03",
"GSM3489187_Donor_03",
"GSM3489188_IPF_03",
"GSM3489188_IPF_03",
"GSM3489189_Donor_04",
"GSM3489189_Donor_04",
"GSM3489190_IPF_04",
"GSM3489190_IPF_04",
"GSM3489191_Donor_05",
"GSM3489191_Donor_05",
"GSM3489192_HP_01",
"GSM3489192_HP_01",
"GSM3489193_Donor_06",
"GSM3489193_Donor_06",
"GSM3489194_SSc-ILD_01",
"GSM3489194_SSc-ILD_01",
"GSM3489195_Donor_07",
"GSM3489195_Donor_07",
"GSM3489196_Myositis-ILD_01",
"GSM3489196_Myositis-ILD_01",
"GSM3489197_Donor_08",
"GSM3489197_Donor_08",
"GSM3489198_SSc-ILD_02",
"GSM3489198_SSc-ILD_02")

donorID <- grep ("Donor",rt)
filtered <- rt[donorID][2*(1:8)]
library(DropletUtils)
for(i in 1:8)
{
  h5m <- m(filtered[i])
  d <- h5m$d
  dt <- as.matrix(h5m$dt)
  dat <- sparseMatrix(i=dt[,1],j=dt[,2],x=dt[,3],dims=d$shape)
  with(d, write10xCounts(filtered[i], dat, barcodes = barcodes, gene.id = genes,
                         gene.symbol = gene_names, gene.type = "Gene Expression", overwrite = TRUE)
  )
}
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
