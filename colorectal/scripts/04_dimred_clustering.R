library(Seurat)
library(tidyverse)
set.seed(12345)

seurat_sketch <- readRDS("data/seurat_sketch.RDS")





seurat_processed <- FindVariableFeatures(seurat_sketch, 
                                         selection.method = "vst",
                                         nfeatures = 2000,
                                         assay = "sketch")


ranked_variable_genes <- VariableFeatures(seurat_processed,
                                          assay = "sketch")
top_genes <- ranked_variable_genes[1:15]
p <- VariableFeaturePlot(seurat_processed)
p <- LabelPoints(plot = p, points = top_genes, repel = TRUE)

ggsave(
  filename = paste0("figures/hvg", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

seurat_processed <- ScaleData(seurat_processed)

#########################################################################333


seurat_processed <- RunPCA(seurat_processed,
                           assay = "sketch",
                           reduction.name = "pca.sketch")


p <- DimPlot(seurat_processed,
        group.by = "orig.ident",
        reduction = "pca.sketch")

ggsave(
  filename = paste0("figures/binsPCA", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

p<-ElbowPlot(object = seurat_processed,
          reduction = "pca.sketch",
          ndims = 50)

ggsave(
  filename = paste0("figures/elbow_PCA", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

seurat_processed <- RunUMAP(seurat_processed, 
                            reduction = "pca.sketch", 
                            reduction.name = "umap.sketch", 
                            return.model = TRUE, 
                            dims = 1:50,
                            seed.use = 12345)


p<-DimPlot(seurat_processed, 
        group.by = "orig.ident",
        reduction = "umap.sketch")
ggsave(
  filename = "figures/umap.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


seurat_processed <- FindNeighbors(seurat_processed, 
                                  assay = "sketch", 
                                  reduction = "pca.sketch",
                                  dims = 1:50)

seurat_processed <- FindClusters(
  seurat_processed,
  graph.name = "sketch_snn",
  cluster.name = "seurat_cluster.sketched",
  algorithm = 4,
  resolution = 0.65,
  random.seed = 12345
)

# print(seurat_processed@graphs)
# print(seurat_processed@meta.data %>% 
#         head() %>% 
#         relocate("seurat_cluster.sketched"))
# print(table(
#   seurat_processed$seurat_cluster.sketched,
#   useNA = "ifany"
# ))

print(Idents(seurat_processed) %>% head())

p <- DimPlot(seurat_processed, 
             reduction = "umap.sketch", label = T) + 
  ggtitle("Sketched clustering")
p<-LabelClusters(p, id = "ident",  fontface = "bold", size = 5, 
              bg.colour = "white", bg.r = .2, force = 0)

ggsave(
  filename = "figures/umap_clusters.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

p<-SpatialDimPlot(seurat_processed,
               pt.size.factor = 16,
               image.alpha = 0)

ggsave(
  filename = "figures/clusters_on_tissue.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

options(future.globals.maxSize = 8000 * 1024^2)

# Project from sketch assay onto entire dataset
#preserving the clusters numbers
#umap.sketch model is used to move sketch pca.sketch -> full.pca_sketch with seurat_cluster_projetcted
seurat_processed <- ProjectData(
  object = seurat_processed,
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched"))
DefaultAssay(seurat_processed) <- "Spatial.008um"
Idents(seurat_processed) <- "seurat_cluster.projected"


p <- DimPlot(
  seurat_processed,
  reduction = "full.umap.sketch",
  group.by = "seurat_cluster.projected",
  label = TRUE,
  repel = TRUE
)

ggsave(
  filename = "figures/umap_clusters_whole_dataset.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

p <- SpatialDimPlot(
  seurat_processed,
  group.by = "seurat_cluster.projected",
  pt.size.factor = 16,
  image.alpha = 0
)

ggsave(
  filename = "figures/clusters_on_tissue_whole_dataset.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

saveRDS(seurat_processed, "data/seurat_processed.RDS")