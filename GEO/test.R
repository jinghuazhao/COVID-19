# 29-3-2020 JHZ

filtered <- scan("filtered.list","")

source("test.ini")
# h5mm()

library(Seurat)
tx <- Read10X(filtered[1])
tx.seurat <- new("seurat", raw.data = tx)

tx.tx <- Setup(tx.seurat, min.cells = 3, min.genes = 200, project = "donor", do.scale = F, do.center = F, names.field = 2, names.delim = "\\-")
mito.genes <- grep("^MT-", rownames(tx@data), value = T)
percent.mito <- colSums(expm1(tx@data[mito.genes, ])) / colSums(expm1(tx@data))

# adds columns to object@data.info to stash QC stats
tx <- AddMetaData(tx, percent.mito, "percent.mito")
VlnPlot(tx, c("nGene", "nUMI", "percent.mito"), nCol = 3)

# - cells that have unique gene counts over 2,500 and under 500, and > 5% mitochondrial percentage
tx <- SubsetData(tx, subset.name = "nGene", accept.high = 2500)
tx <- SubsetData(tx, subset.name = "percent.mito", accept.high = 0.05)
tx <- SubsetData(tx, subset.name = "nGene", accept.low = 500)

# gene outliers on mean-variability plot
tx <- MeanVarPlot(tx, x.low.cutoff = 0, y.cutoff = 0.8)
length(tx@var.genes)

# nb regression on the variable genes, this sets their value in tx@scale.data, which is used for PCA/clustering
# with mitochondrial percentage, batch, and nUMI as confounding variables

tx <- RegressOut(tx,latent.vars = c("percent.mito", "orig.ident", "nUMI"), genes.regress = tx@var.genes, model.use = "negbinom")

# PCA with the IRLBA package (iteratively computes the top dimensions, dramatic increase in speed since we are throwing away most PCs anyway)
tx <- PCAFast(tx, pc.genes = tx@var.genes, pcs.compute = 40, pcs.print = 30)

PCElbowPlot(tx, num.pc = 40)
PrintPCA(tx,pcs.print = 1:36)
PCHeatmap(tx,pc.use = 1:12,100)
PCHeatmap(tx,pc.use = 13:24,100)
PCHeatmap(tx,pc.use = 25:36,100)

# select 25 PCs for downstream analysis
tx <- RunTSNE(tx, dims.use = 1:25, do.fast = T)

# save.SNN == re-run with different resolution values. Here we run with a few different res v
tx <- FindClusters(tx ,pc.use = 1:25, resolution = seq(2,4,0.5), save.SNN = T, do.sparse = T)

#Explore results with a resolution value of 2
tx <- SetAllIdent(tx, id = "res.2")
TSNEPlot(tx,do.label = T)

# slightly over-clusters the data
tx <- SetAllIdent(tx, id = "res.4")
TSNEPlot(tx,do.label = T)

# the 'multi-resolution' problem in graph-based clustering resolved by slightly over-clustering the data with a post-hoc merging step
#  Out-of-bag error (OOBE) from a random forest classifier as a test for merging

tx <- SetAllIdent(tx, id = "res.4")

#Build a classification hierarchy that places transcriptionally similar clusters adjacent on a tree
tx <- BuildClusterTree(tx, do.reorder = T, reorder.numeric = T)

#calculate the classification error for left/right cells on each branch of the tree
# For nodes with high OOBE, we cannot accurately tell the left/right children apart based on random forests, so the clusters may need to be merged
node.scores <- AssessNodes(tx)
node.scores[order(node.scores$oobe,decreasing = T),] -> node.scores

# choose the first eight splits to merge )
nodes.merge=node.scores[1:8,]
nodes.to.merge=sort(nodes.merge$node)
tx.merged <- tx
for (n in nodes.to.merge){
  tx.merged <- MergeNode(tx.merged, n)
}

# rebuild the classification hierarchy using all genes for interpretation
tx.merged <- BuildClusterTree(tx.merged, do.reorder = T, reorder.numeric = T,genes.use = rownames(tx.merged@data))
TSNEPlot(tx.merged, do.label = T)
PlotClusterTree(tx.merged)

# color TSNE based on a hierarchical split
ColorTSNESplit(tx.merged, node = 33)
ColorTSNESplit(tx.merged, node = 31,color1 = "red",color3="blue")
