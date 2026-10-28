
library(Seurat)

gbm[["RNA"]]<- ScaleData(gbm[["RNA"]])
gbm[["Protein"]] <- ScaleData(gbm[["Protein"]])


gbm <- RunPCA(
  gbm,
  assay = "RNA",
  reduction.name = "pca.RNA",
  npcs = 10
)
gbm <- RunPCA(
  gbm,
  assay = "Protein",
  reduction.name = "pca.Protein",
  npcs = 10
)

p1 <- DimPlot(gbm,
        group.by = "orig.ident",
        reduction = "pca.RNA")


p2 <- DimPlot(gbm,
              group.by = "orig.ident",
              reduction = "pca.Protein")

p3 <- VizDimLoadings(gbm, reduction = "pca.RNA", dims = 1:2)
p4<-VizDimLoadings(
  gbm,
  reduction = "pca.Protein",
  dims = 1:2,
  nfeatures = 10,
  balanced = TRUE
)

gbm <- RunUMAP(gbm, 
                reduction = "pca.RNA", 
                reduction.name = "umap.RNA", 
                return.model = TRUE, 
                dims = 1:10,
                seed.use = 12345)

p5 <- DimPlot(gbm,
  reduction = "umap.RNA",
  group.by = "orig.ident",
  label = TRUE)


gbm <- RunUMAP(gbm, 
               reduction = "pca.Protein", 
               reduction.name = "umap.Protein", 
               return.model = TRUE, 
               dims = 1:10,
               seed.use = 12345)

p6 <- DimPlot(gbm,
              reduction = "umap.Protein",
              group.by = "orig.ident",
              label = TRUE)



library(ggplot2)


output_dir <- "figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(output_dir, "RNA_PCA.png"),
  plot = p1,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(output_dir, "Protein_PCA.png"),
  plot = p2,
  width = 10, height = 7, dpi = 300
)

ggsave(
  file.path(output_dir, "RNA_PCA_loadings.png"),
  plot = p3,
  width = 14, height = 9, dpi = 300
)

ggsave(
  file.path(output_dir, "Protein_PCA_loadings.png"),
  plot = p4,
  width = 14, height = 9, dpi = 300
)

ggsave(
  file.path(output_dir, "RNA_UMAP.png"),
  plot = p5,
  width = 10, height = 7, dpi = 300
)
ggsave(
  file.path(output_dir, "Protein_UMAP.png"),
  plot = p6,
  width = 10, height = 7, dpi = 300
)