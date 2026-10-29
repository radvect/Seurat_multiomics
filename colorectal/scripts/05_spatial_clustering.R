library(Seurat)
library(Banksy)
library(SeuratWrappers)
library(tidyverse)
library(future)
library(harmony)
library(patchwork)


#parallel computing and setting stuff
set.seed(12345)
options(future.globals.maxSize = 8000 * 1024^2)
plan(multisession, workers = parallel::detectCores() - 1)

seurat_processed <- readRDS("data/seurat_processed.RDS")


#getting coords from images
seurat_processed$cell <- rownames(seurat_processed@meta.data)
coords_crc <- GetTissueCoordinates(seurat_processed, 
                                   image = "P5CRC.008um")
coords_nat <- GetTissueCoordinates(seurat_processed, 
                                   image = "P5NAT.008um")
coords <- rbind(coords_crc, coords_nat)
print(coords %>% head())

seurat_processed@meta.data <- left_join(x = seurat_processed@meta.data, 
                                        y = coords, 
                                        by = "cell")
print(seurat_processed@meta.data %>% head())

print(rownames(seurat_processed@meta.data) <- seurat_processed$cell)


seurat_banksy <- readRDS("data/seurat_banksy.rds")
# seurat_banksy <- RunBanksy(seurat_processed, 
#                            lambda = 0.8, 
#                            k_geom = 15, 
#                            use_agf = TRUE,
#                            dimx = 'x',
#                            dimy = 'y', 
#                            group = "orig.ident", 
#                            split.scale = TRUE,
#                            assay = 'Spatial.008um', 
#                            slot = 'data',
#                            verbose = TRUE)

print(Features(seurat_banksy)[1:6])
print(Features(seurat_banksy)[2000:2006])
print(Features(seurat_banksy)[4000:4006])
#############################################Yu-Chen Cutoff

# seurat_banksy <- RunPCA(seurat_banksy, 
#                         assay = "BANKSY", 
#                         reduction.name = "pca.banksy", 
#                         features = Features(seurat_banksy), 
#                         npcs = 30)

p_orig <- DimPlot(seurat_banksy,
                  group.by = "orig.ident",
                  reduction = "pca.banksy") 


p_igkc <- FeaturePlot(seurat_banksy,
                      feature = "IGKC",
                      reduction = "pca.banksy")

ggsave(
  filename = paste0("figures/banksy_pca", ".png"),
  plot = p_orig+p_igkc,
  width = 6,
  height = 4,
  dpi = 300
)

# seurat_banksy <- RunHarmony(seurat_banksy,
#                             group.by.vars = "orig.ident",
#                             reduction = "pca.banksy",
#                             reduction.save = "harmony.banksy")


p_orig <- DimPlot(seurat_banksy,
                  group.by = "orig.ident",
                  reduction = "harmony.banksy") 

# Harmony BANKSY PCA colored by IGKC expression
p_igkc <- FeaturePlot(seurat_banksy,
                      feature = "IGKC",
                      reduction = "harmony.banksy")

ggsave(
  filename = paste0("figures/after_harmony_clusters", ".png"),
  plot = p_orig+p_igkc,
  width = 6,
  height = 4,
  dpi = 300
)

######################################################
#Clustering under harmony and without

# seurat_banksy <- FindNeighbors(seurat_banksy, 
#                                reduction = "harmony.banksy", 
#                                dims = 1:30)

# Leiden clustering
# seurat_banksy <- FindClusters(seurat_banksy, 
#                               algorithm = 4,
#                               cluster.name = "banksy_cluster",
#                               resolution = 0.65,
#                               random.seed = 12345)


DefaultAssay(seurat_banksy) <- "Spatial.008um"

# Delete BANKSY assay by setting it to NULL
#seurat_banksy[["BANKSY"]] <- NULL


##Spatial clutering change
p <-SpatialDimPlot(seurat_banksy, 
               group.by = c("banksy_cluster",
                            "seurat_cluster.projected"), 
               pt.size.factor = 15, 
               image.alpha = 0)

ggsave(
  filename = paste0("figures/banksy_space", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


p<-DimPlot(seurat_banksy, 
        group.by = c("banksy_cluster",
                     "seurat_cluster.projected"),
        reduction = "full.umap.sketch",
        label = TRUE)


ggsave(
  filename = paste0("figures/banksy_pca_clustering_comparison", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


saveRDS(seurat_banksy, "data/seurat_banksy.RDS")