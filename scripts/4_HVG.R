library(Seurat)

gbm <- FindVariableFeatures(gbm,
                            assay = "RNA",
                            selection.method = "vst",
                            nfeatures = 2000)
#ranked_genes <- VariableFeatures(gbm)

#length(VariableFeatures(gbm[["RNA"]]))
#ranked_genes

ranked_genes <- VariableFeatures(gbm[["RNA"]])
top_genes <- ranked_genes[1:15]
rm(ranked_genes)

gbm <- FindVariableFeatures(gbm,
                            assay = "Protein",
                            selection.method = "vst",
                            nfeatures = 2000)
ranked_proteins <- VariableFeatures(gbm[["Protein"]])
top_proteins <-ranked_proteins[1:15]
rm(ranked_proteins)

