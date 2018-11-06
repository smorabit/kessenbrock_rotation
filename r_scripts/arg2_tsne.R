library("Seurat")

#load mdsc.combined.cca and doublet data:
load("~/Documents/kessenbrock_rotation/Doublet_Pipelines_Data/Combined_CCA/mdsc.combined.cca.Seurat.object.Rda")
doublet_data <- read.table(file = "~/Documents/kessenbrock_rotation/detected_doublets.tsv", sep="\t", header=TRUE)

#set index to barcodes, remove barcodes column
rownames(doublet_data) <- doublet_data$barcode 
doublet_data$barcode <- NULL

#put doublet detection columns in the metadata df
mdsc.combined.cca@meta.data$doublet_score <- doublet_data$doublet_score
mdsc.combined.cca@meta.data$doublet_prediction <- doublet_data$doublet_prediction

#plot some stuff
FeaturePlot(mdsc.combined.cca, c("doublet_score"), cols.use = c("lightgrey","darkblue"), do.return=TRUE)
p2<-TSNEPlot(mdsc.combined.cca, group.by="doublet_prediction")

#save .rda
save(mdsc.combined.cca, file="~/Documents/kessenbrock_rotation/Doublet_Pipelines_Data/Combined_CCA/mdsc.combined.cca.Seurat.object.Rda")
