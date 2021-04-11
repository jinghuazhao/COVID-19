# DEBUG
mdsfile <- "mydata.refdata.QCed.MDS.mds"
reffile <- "/rds/user/jhz22/hpc-work/VEGAS2/keep.fam"
cohort <- "INTERVAL"

cat(mdsfile,reffile,cohort,"\n")

mds <- read.table(mdsfile, header=TRUE)
rgl::plot3d(with(mds,cbind(C1,C2,C3)))

# read mds file
dat <- read.table(mdsfile, header=TRUE, sep="", as.is=T)
dat$name <- "mydata"
# read reference *.fam file for names of reference population
ref <- read.table(reffile, header=F, sep="", as.is=T)

refpos <- which(dat$FID %in% ref$V1)
dat$name[refpos] <- dat$FID[refpos]

# set 7 colors
col <- c(rgb(27,158,119,max=255),rgb(217,95,2,max=255),rgb(117,112,179,max=255),
    rgb(102,166,30,max=255),rgb(231,41,138,max=255),rgb(230,171,2,max=255),rgb(166,118,29,max=255))
dat$col <- col[as.numeric(as.factor(dat$name))]

# plot first 4 components
# outname <- paste(cohort,".C1-C4.png",sep="")
# png(outname, width=12, height=12, units="in", res=600)
# pairs(~C1+C2+C3+C4, data=dat, col=dat$col, pch=20)
# dev.off()

# plot first 2 components with legend
outname <- paste(cohort,".C1-C2.png",sep="")
png(outname, width=8, height=8, units="in", res=600)
plot(dat$C1,dat$C2,col=dat$col, pch=20)
legend("topleft", legend=unique(dat$name),col=unique(dat$col), pch=20); box()
dev.off()
