library(Seurat)
#expression of the RNA higher - > lognormalize
gbm[["RNA"]] <- NormalizeData(
  gbm[["RNA"]],
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
gbm[["Protein"]] <- NormalizeData(
  gbm[["Protein"]],
  normalization.method = "CLR",
  margin = 2
)


rna_raw <- LayerData(gbm, assay = "RNA", layer = "counts")
rna_norm <- LayerData(gbm, assay = "RNA", layer = "data")


protein_raw <- LayerData(gbm, assay = "Protein", layer = "counts")
protein_norm <- LayerData(gbm, assay = "Protein", layer = "data")