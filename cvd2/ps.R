# 8-4-2020 JHZ

options(width=500)
require(phenoscanner)
catalogue <- Sys.getenv("catalogue")
rsid <- with(read.delim("ACE2.merge",as.is=TRUE),MarkerName)
batches <- split(rsid, ceiling(seq_along(rsid)/100))
s <- t <- list()
for(i in 1:length(batches))
{
  q <- phenoscanner(snpquery=batches[[i]], catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  s[[i]] <- with(q,snps)
  t[[i]] <- with(q,results)
}
snps <- do.call(rbind,s)
results <- do.call(rbind,t)
r <- list(snps=snps,results=results)
save(r,file=paste0("ACE2.",catalogue,".rda"))
results <- within(results,{
   a1 <- ref_a1
   a2 <- ref_a2
   swap <- ref_a1 > ref_a2
   a1[swap] <- ref_a2[swap]
   a2[swap] <- ref_a1[swap]
   ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
})
sink(paste0("ACE2.",catalogue))
s <- results[c("ref_rsid","ref_snpid","rsid","r2","p","trait","dataset","pmid")]
print(s)
sink()

print_by_study <- function()
{
  for(d in unique(with(results,dataset)))
  {
    cat(d,"\n")
    sink(paste(catalogue,d,sep="."))
    s <- subset(results[c("ref_rsid","ref_snpid","rsid","r2","p","trait","dataset","pmid")],dataset==d)
    print(s)
    sink()
  }
}
