library(Seurat)
library(tidyverse)
library(scales)


set.seed(12345)


seurat_filtered <- readRDS("data/seurat_filtered.RDS")
seurat_sketch <- NormalizeData(seurat_filtered)
seurat_sketch <- FindVariableFeatures(seurat_sketch, 
                                      selection.method = "vst",
                                      nfeatures = 2000)

ranked_variable_genes <- VariableFeatures(seurat_sketch)
top_genes <- ranked_variable_genes[1:15]

seurat_sketch <- SketchData(object = seurat_sketch,
                            assay = 'Spatial.008um',
                            ncells = 10000,
                            method = "LeverageScore",
                            sketched.assay = "sketch")


p <- ggplot(seurat_sketch@meta.data) +
  geom_histogram(aes(x=leverage.score, 
                     fill=orig.ident), 
                 alpha=0.5, bins=100) +
  theme_classic() +
  scale_x_log10(labels = scales::label_number())


ggsave(
  filename = paste0("figures/histogram_leverage", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


p <- SpatialFeaturePlot(seurat_sketch,
                   "leverage.score",
                   pt.size.factor = 15,
                   image.alpha = 0,
                   max.cutoff = 2)



ggsave(
  filename = paste0("figures/spatial_leverage", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

saveRDS(seurat_sketch, "data/seurat_sketch.RDS")

rm(seurat_filtered)
rm(seurat_sketch)
rm(p)

