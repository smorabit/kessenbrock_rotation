#install loomR
#install.packages("devtools")
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

library(loomR)
library(Seurat)

#convert WT, PYMT and mdsc.combined.cca to .loom files (those objects were made in the doublet_Finder.R script):
mdsc.combined.cca_loom <- Convert(from = mdsc.combined.cca, to = "loom", filename="~/Documents/kessenbrock_rotation/RNA_velocity/loom_files/mdsc.combined.cca.loom")

mdsc.WT_loom <- Convert(from = WT, to = "loom", filename="~/Documents/kessenbrock_rotation/RNA_velocity/loom_files/mdsc.WT.loom")

lfile = Convert(from = pbmc_small, to = 'loom', filename = '~/Documents/kessenbrock_rotation/RNA_velocity/loom_files/small.loom')