library(BiocGenerics)
library(GenomeInfoDb)
library(Seurat)
#library(Signac)

library(BiocGenerics)
gbm  <- Load10X_Spatial(data.dir = "data/glioblastoma/",
                        filename = "CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
)


