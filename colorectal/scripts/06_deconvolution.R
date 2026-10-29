library(Seurat)
library(tidyverse)
library(spacexr)

plan("multicore", workers = parallel::detectCores() - 1)
options(future.globals.maxSize = 8000 * 1024^2)

seurat_banksy <- readRDS("data/seurat_banksy.rds")

crc <- subset(seurat_banksy, 
              subset = (orig.ident == "P5CRC"))
DefaultAssay(crc) <- 'Spatial.008um'
crc <- FindVariableFeatures(object = crc,
                            assay = "Spatial.008um")
crc <- SketchData(
  object = crc,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch_crc",
  var.name = "leverage.score_crc")

crc <- FindVariableFeatures(crc)
crc <- ScaleData(crc)
crc <- RunPCA(crc,
              reduction.name = "pca.sketch_crc")
DefaultAssay(crc) <- "sketch_crc"
counts_sketch <- crc[["sketch_crc"]]$counts
cells_sketch <- colnames(crc[["sketch_crc"]])
coords <- GetTissueCoordinates(crc)[cells_sketch, c("x", "y")]
nUMI_sketch <- crc@meta.data[cells_sketch, "nCounts_Spatial.008um"]
query <- SpatialRNA(coords = coords, 
                    counts = counts_sketch, 
                    nUMI = nUMI_sketch)



seurat_ref <- readRDS("data/crc_flex_ref_downsample.RDS")

p <- DimPlot(seurat_ref) +
  ggtitle("RCTD Reference Dataset")
# Add cluster labels
p <- LabelClusters(p,
              id = "ident",
              fontface = "bold",
              size = 3, 
              bg.colour = "white",
              bg.r = .2,
              force = 0)
ggsave(
  filename = paste0("figures/ref_dataset", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

counts <- seurat_ref[["RNA"]]$counts

# cell type annotation of reference dataset (must be a factor)
cluster <- seurat_ref$Level1
cluster <- as.factor(cluster)

# Total counts of each cell
nUMI <- seurat_ref$nCount_RNA

# Create the RCTD reference object
reference <- Reference(counts = counts, 
                       cell_types = cluster, 
                       nUMI = nUMI)
RCTD <- create.RCTD(spatialRNA = query, 
                    reference = reference, 
                    max_cores = parallel::detectCores() - 1)
# Run deconvolution
# This step may take 5-10 minutes to run
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

df_rctd_results <- RCTD@results$results_df
df_rctd_results %>%
  head()


ggplot(df_rctd_results,
       aes(x = first_type, 
           fill = first_type)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  ggtitle("P5CRC: RCTD First Type Results") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))


###########before projecting

# Add "sketch" to column names
colnames(df_rctd_results) <- paste0(colnames(df_rctd_results), 
                                    "_sketch")

# Create column "cell" to merge on
# Removing factors for easier merging later
df_rctd_results <- df_rctd_results %>%
  mutate(cell = rownames(df_rctd_results)) %>%
  mutate(first_type_sketch = as.character(first_type_sketch),
         second_type_sketch = as.character(second_type_sketch),
         spot_class_sketch = as.character(spot_class_sketch))

# Make sure columns were generated
colnames(df_rctd_results)


crc$cell <- rownames(crc@meta.data)
meta <- crc@meta.data

# Merge together RCTD results and metadata
meta <- left_join(x = meta, 
                  y = df_rctd_results,
                  by = "cell")
rownames(meta) <- meta$cell

# Put meta into crc Seurat object
crc@meta.data <- meta

cols <- c("spot_class_sketch", 
          "first_type_sketch", 
          "second_type_sketch")

# Iterate through each column to set "unassigned" to NA values
for (col in cols) {
  # Find rows that are in sketch_cells and are NA
  na_rows <- cells_sketch[is.na(crc@meta.data[cells_sketch, col])]
  # Assign the value to unassigned
  crc@meta.data[na_rows, col] <- "unassigned"
}

# Project back first_type to new column called first_type
crc <- ProjectData(
  object = crc,
  assay = "Spatial.008um",
  sketched.assay = "sketch_crc",
  full.reduction = "pca_crc",
  sketched.reduction = "pca.sketch_crc",
  dims = 1:30,
  refdata = list(spot_class = "spot_class_sketch",
                 first_type = "first_type_sketch",
                 second_type = "second_type_sketch"))



# Barplot of first_type labels after ProjectData()
p <- ggplot(crc@meta.data,
       aes(x = first_type, fill = first_type)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  ggtitle("P5CRC: Projected First Type RCTD") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Spatial plot of first_type labels after ProjectData()
p <- SpatialDimPlot(crc,
               group.by = "first_type",
               pt.size.factor = 15,
               image.alpha = 0.5)
ggsave(
  filename = paste0("figures/spatial_histogram", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


# UMAP of reference dataset
Idents(crc) <- "first_type"
p <- DimPlot(crc,
             reduction = "full.umap.sketch") +
  ggtitle("P5CRC RCTD First Type")
p<-LabelClusters(p, id = "ident",  fontface = "bold", size = 3, 
              bg.colour = "white", bg.r = .2, force = 0)
ggsave(
  filename = paste0("figures/spatial_cell_type", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)
