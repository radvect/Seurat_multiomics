####TASK PORT DOWNSAMPLING HERE:

library(Seurat)
library(leidenbase)

library(ggplot2)
gbm <- FindNeighbors(gbm, 
                                  assay = "RNA", 
                                  reduction = "pca.RNA",
                                  dims = 1:10,
                                  graph.name = c("RNA_nn", "RNA_snn"))
gbm <- FindNeighbors(gbm, 
                     assay = "Protein", 
                     reduction = "pca.Protein",
                     dims = 1:10, 
                     graph.name = c("Protein_nn", "Protein_snn"))



gbm <- FindClusters(
  gbm,
  graph.name = "RNA_snn",
  cluster.name = "RNA_clusters",
  algorithm = 4,
  resolution = 0.65,
  random.seed = 12345
)

# Protein clustering
gbm <- FindClusters(
  gbm,
  graph.name = "Protein_snn",
  cluster.name = "Protein_clusters",
  algorithm = 4,
  resolution = 0.65,
  random.seed = 12345
)


p <- DimPlot(
  gbm,
  reduction = "umap.RNA",
  group.by = "RNA_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("RNA clustering")

# Visualize Protein clusters
p1 <- DimPlot(
  gbm,
  reduction = "umap.Protein",
  group.by = "Protein_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Protein clustering")



output_dir <- "figures/"
ggsave(
  file.path(output_dir, "RNA_UMAP_clusters.png"),
  plot = p,
  width = 10, height = 7, dpi = 300
)

output_dir <- "figures/"
ggsave(
  file.path(output_dir, "Proteins_UMAP_clusters.png"),
  plot = p1,
  width = 10, height = 7, dpi = 300
)

#Should I do the clustering together? 



options(future.globals.maxSize = 8000 * 1024^2)



p_rna <- SpatialDimPlot(
  gbm,
  images = "slice1",
  group.by = "RNA_clusters"
)

p_protein <- SpatialDimPlot(
  gbm,
  images = "slice1",
  group.by = "Protein_clusters"
# )

combined <- p_rna + p_protein



ggsave(
  file.path(output_dir, "clusters_on_gbm.png"),
  plot = combined,
  width = 10, height = 7, dpi = 300
)