# installation of DoubletFinder:
#devtools::install_github('chris-mcginnis-ucsf/doubletFinder')

#include required libraries:
library("doubletFinder")
library("Seurat")

#load mdsc.combined.cca:
load("~/Documents/kessenbrock_rotation/Doublet_Pipelines_Data/Combined_CCA/mdsc.combined.cca.Seurat.object.Rda")

#subset dataset into wt and pymt:
mdsc.wt <- SubsetData(object = mdsc.combined.cca, subset.name = "orig.ident", accept.value="wt")
mdsc.pymt <- SubsetData(object = mdsc.combined.cca, subset.name = "orig.ident", accept.value = "pymt")

#compute doublets with doubletFinder:
#doubletFinder(mdsc.wt, expected.doublets = 0.06)
doubletFinder(mdsc.combined.cca, expected.doublets=0.06) #this gives an error so let's process it from the beginning

########################################################################################################

#Separately label PYMT and WT
mdsc.combined.data <- Read10X("~/Documents/kessenbrock_rotation/MDSCs/filtered_gene_bc_matrices_mex/mm10/")
wt.ind<-grep("-1",colnames(mdsc.combined.data))
pymt.ind<-grep("-2",colnames(mdsc.combined.data))
colnames(mdsc.combined.data)[wt.ind]<-paste("wt",colnames(mdsc.combined.data)[wt.ind],sep="_")
colnames(mdsc.combined.data)[pymt.ind]<-paste("pymt",colnames(mdsc.combined.data)[pymt.ind],sep="_")

# Set up WT object
WT <- CreateSeuratObject(raw.data = mdsc.combined.data[,wt.ind], min.cells = 3,min.genes = 500, project = "WT")
WT <- NormalizeData(WT)
mito.genes <- grep("^mt-", rownames(WT@data), value = T)
percent.mito <- colSums(expm1(as.matrix(WT@data[mito.genes, ])))/colSums(expm1(as.matrix(WT@data)))
WT <- AddMetaData(WT, percent.mito, "percent.mito")
VlnPlot(WT, c("nGene", "nUMI", "percent.mito"), nCol = 3)
WT <- SubsetData(WT, subset.name = "percent.mito", accept.high = 0.08)
WT <- SubsetData(WT, subset.name = "nGene", accept.high = 5000)
WT <- ScaleData(WT, display.progress = F)

# Set up PYMTulated object
PYMT <- CreateSeuratObject(raw.data = mdsc.combined.data[,pymt.ind], min.cells = 3,min.genes = 500, project = "PYMT")
PYMT <- NormalizeData(PYMT)
mito.genes <- grep("^mt-", rownames(PYMT@data), value = T)
percent.mito <- colSums(expm1(as.matrix(PYMT@data[mito.genes, ])))/colSums(expm1(as.matrix(PYMT@data)))
PYMT <- AddMetaData(PYMT, percent.mito, "percent.mito")
VlnPlot(PYMT, c("nGene", "nUMI", "percent.mito"), nCol = 3)
PYMT <- SubsetData(PYMT, subset.name = "percent.mito", accept.high = 0.08)
PYMT <- SubsetData(PYMT, subset.name = "nGene", accept.high = 5000)
PYMT <- ScaleData(PYMT, display.progress = F)

# find variable genes for WT and PYMT
WT <- FindVariableGenes(WT, do.plot = T)
PYMT <- FindVariableGenes(PYMT, do.plot = T)

#run pca:
WT <- RunPCA(WT, pc.genes = WT@var.genes)
PYMT <- RunPCA(PYMT, pc.genes = PYMT@var.genes)

#figure out how many PCs to use for clustering:
PCElbowPlot(object = WT)

#find clusters:
WT <- FindClusters(object = WT, reduction.type = "pca", dims.use = 1:20, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
PYMT <- FindClusters(object = PYMT, reduction.type = "pca", dims.use = 1:20, 
                   resolution = 0.6, print.output = 0, save.SNN = TRUE)

# run TSNEs
WT <- RunTSNE(WT, dims.use = 1:20, do.fast = T)
PYMT <- RunTSNE(PYMT, dims.use = 1:20, do.fast = T)

#try doubletFinder:
WT <- doubletFinder(WT, expected.doublets = 250)
PYMT <- doubletFinder(PYMT, expected.doublets=500)

# merge back together:
merged_mdsc <- MergeSeurat(object1 = WT, object2 = PYMT) 
mdsc.combined.cca@meta.data$doubletFinder_prediction <- merged_mdsc@meta.data$pANNPredictions

# tsne of predicted doublets:
p1 <- TSNEPlot(mdsc.combined.cca, group.by="doublet_prediction")
p2 <- TSNEPlot(mdsc.combined.cca, group.by="doubletFinder_prediction", colors.use = ("red", "blue"))
plot_grid(p1,p2)
? TSNEPlot

#example just for myself of how loops and conditional statements work in R
singlet_counter = 0
doublet_counter = 0
for (i in 1:length(PYMT@meta.data$pANNPredictions)){
  if (PYMT@meta.data$pANNPredictions[i] == "Singlet"){
    singlet_counter = singlet_counter + 1
  } else if (PYMT@meta.data$pANNPredictions[i] == "Doublet"){
    doublet_counter = doublet_counter + 1
  }
}
print(doublet_counter)
print(singlet_counter)