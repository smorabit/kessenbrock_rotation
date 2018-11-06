#load required packages:
library("Seurat")
library("DoubletDecon")

#load existing Seurat object:
load("/Volumes/shared/Sam/Doublet_Pipelines_Data/Combined_CCA/mdsc.combined.cca.Seurat.object.Rda")

#get marker genes
mdsc.combined.cca.markers <- FindAllMarkers(object = mdsc.combined.cca, only.pos = TRUE, min.pct = 0.25, hresh.use = 0.25)
saveRDS(mdsc.combined.cca.markers, file = "/Volumes/shared/Sam/doublet_detection/mdsc/mdsc.combined.cca.markers.rds")

#create required files for DoubletDecon from Seurat object
write.table(mdsc.combined.cca@scale.data, "/Volumes/shared/Sam/doublet_detection/mdsc/mdsc.combined.cca_scale.data.txt", sep="\t")
write.table(mdsc.combined.cca.markers,
            "/Volumes/shared/Sam/doublet_detection/mdsc/mdsc.combined.cca_allmarkers.txt", sep="\t")
#write.table(mdsc.combined.cca.markers %>% group_by(cluster) %>% top_n(50, avg_logFC),
            #"/Volumes/shared/Sam/doublet_detection/mdsc/mdsc.combined.cca_markers.txt", sep="\t")
write.table(mdsc.combined.cca@ident, 
            "/Volumes/shared/Sam/doublet_detection/mdsc/mdsc.combined.cca_clusters.txt", sep="\t")

#testing files:
#location="/Volumes/shared/Sam/doublet_detection/"
#expressionFile=paste0(location, "expressionMatrix.txt")
#genesFile=paste0(location, "markers.txt")
#clustersFile=paste0(location, "clusters.txt")

#actual data:
location="/Volumes/shared/Sam/doublet_detection/mdsc/"
expressionFile=paste0(location, "mdsc.combined.cca_scale.data.txt")
genesFile=paste0(location, "mdsc.combined.cca_allmarkers.txt")
clustersFile=paste0(location, "mdsc.combined.cca_clusters.txt")

newFiles=Seurat_Pre_Process(expressionFile, genesFile, clustersFile)

results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                           groupsFile=newFiles$newGroupsFile, 
                           filename="mdsc", 
                           location=location,
                           fullDataFile=NULL, 
                           removeCC=FALSE, 
                           species="mmu", 
                           rhop=1.1, 
                           write=TRUE, 
                           recluster="doublets_decon", 
                           PMF=TRUE, 
                           useFull=FALSE, 
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100)



