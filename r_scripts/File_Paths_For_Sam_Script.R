#### READ IN STROMAL DATASETS ####

#Read in ind2 BRCA stromal data set 

stromal.ind2 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind2_m1161181a_brca\\stromal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind2.stromal to retain the ind2 label throughout the analysis 

colnames(stromal.ind2)<-paste("ind2.stromal.BRCA",colnames(stromal.ind2),sep="_")

#Read in ind3 BRCA stromal data set

stromal.ind3 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind3_m1161055a_brca\\stromal\\filtered_gene_bc_matrices\\GRCh38")

colnames(stromal.ind3)<-paste("ind3.stromal.BRCA",colnames(stromal.ind3),sep="_")

#Read in the ind4 brca data set 

stromal.ind4 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind4_m1161343a_brca\\stromal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind4.stromal to retain the ind4.stromal label throughout the analysis 

colnames(stromal.ind4)<-paste("ind4.stromal.BRCA",colnames(stromal.ind4),sep="_")

#Read in the ind1 normal data set 

stromal.ind1 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind1_wd37026_normal\\stromal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind1.stromal to retain the ind1.stromal label throughout the analysis 

colnames(stromal.ind1)<-paste("ind1.stromal.NORMAL",colnames(stromal.ind1),sep="_")

#Read in the ind9 normal data set 

stromal.ind9 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind9_wd57922_normal\\stromal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind9.stromal to retain the ind9.stromal label throughout the analysis 

colnames(stromal.ind9)<-paste("ind9.stromal.NORMAL",colnames(stromal.ind9),sep="_")

#Read in the ind10 normal data set 

stromal.ind10 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind10_wd58673_normal\\stromal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind10.stromal to retain the ind10.stromal label throughout the analysis 

colnames(stromal.ind10)<-paste("ind10.stromal.NORMAL",colnames(stromal.ind10),sep="_")




#### READ IN EPITHELIAL DATASETS ####


#Read in the ind2 brca epithelial dataset 

epithelial.ind2 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind2_m1161181a_brca\\basal_luminal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind2.epithelial to retain the ind2 label throughout the analysis 

colnames(epithelial.ind2)<-paste("ind2.epithelial.BRCA",colnames(epithelial.ind2),sep="_")

#Read in ind3 brca epithelial

epithelial.ind3 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind3_m1161055a_brca\\basal_luminal\\filtered_gene_bc_matrices\\GRCh38")

colnames(epithelial.ind3)<-paste("ind3.epithelial.BRCA",colnames(epithelial.ind3),sep="_")

#Read in the ind4 brca data set 

epithelial.ind4 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind4_m1161343a_brca\\basal_luminal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind4.epithelial to retain the ind4.epithelial label throughout the analysis 

colnames(epithelial.ind4)<-paste("ind4.epithelial.BRCA",colnames(epithelial.ind4),sep="_")

#Read in the ind1 normal data set 

epithelial.ind1 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind1_wd37026_normal\\basal_luminal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind1.epithelial to retain the ind1.epithelial label throughout the analysis 

colnames(epithelial.ind1)<-paste("ind1.epithelial.NORMAL",colnames(epithelial.ind1),sep="_")

#Read in the ind9 normal data set 

epithelial.ind9 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind9_wd57922_normal\\basal_luminal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind9.epithelial to retain the ind9.epithelial label throughout the analysis 

colnames(epithelial.ind9)<-paste("ind9.epithelial.NORMAL",colnames(epithelial.ind9),sep="_")

#Read in the ind10 normal data set 

epithelial.ind10 <- Read10X(data.dir = "Z:\\Quy Alig\\10x_grch38_cr_210\\ind10_wd58673_normal\\basal_luminal\\filtered_gene_bc_matrices\\GRCh38")

#paste in ind10.epithelial to retain the ind10.epithelial label throughout the analysis 

colnames(epithelial.ind10)<-paste("ind10.epithelial.NORMAL",colnames(epithelial.ind10),sep="_")