library(Seurat)

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