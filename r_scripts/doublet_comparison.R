# compare doublet scores
#
# Need to have already run doublet prediction from scrublet, and doubletFinder, 
# and incorporated into mdsc.combined.cca

library(Seurat)
library(nopaco)
library(DescTools)

load("~/Documents/kessenbrock_rotation/Doublet_Pipelines_Data/Combined_CCA/mdsc.combined.cca.Seurat.object.Rda")

mat = matrix(mdsc.combined.cca@meta.data$doublet_prediction)
mat = cbind(mdsc.combined.cca@meta.data$doubletFinder_prediction)
doublet_concordance = concordance.test(mat)

fisher = matrix(0, nrow=2, ncol=2)

for (i in 1:length(mdsc.combined.cca@meta.data$doublet_prediction)){
  scr = mdsc.combined.cca@meta.data$doublet_prediction[i]
  dubfind =  mdsc.combined.cca@meta.data$doubletFinder_prediction[i]
  # if they are both singlets:
  if (scr == 0 && dubfind == 0){
    fisher[1,1] = fisher[1,1] + 1
  # if scrublet is singlet and doubletFinder is doublet:
  } else if (scr == 0 && dubfind == 1) {
    fisher[1,2] = fisher[1,2] + 1
  }
  # if scrublet is doublet and doubletFinder is singlet:
  else if (scr == 1 && dubfind == 0){
    fisher[2,1] = fisher[2,1] + 1
  }
  # if both are predicted doublets:
  else if (scr == 1 && dubfind == 1){
    fisher[2,2] = fisher[2,2] + 1
  }
}

