library(Seurat)
library(arrow)
library(dplyr)


#crc <- Load10X_Spatial(data.dir = "data/P5CRC_cropped/",
#                       bin.size = c(8, 16),
#                       slice = "P5CRC")



# how many bins belong to each bin size  (smaller bin -> more bins)
#ncol(crc[["Spatial.008um"]])
#ncol(crc[["Spatial.016um"]])

#ncount - number of UMIs(number of amplificated molecules of RNA), #nfeatures - number of genes 
#crc@meta.data %>% View()

#rm(crc)



samples <- c("P5CRC", "P5NAT")
list_seurat <- list()


for (sample in samples) {
  data_dir <- paste0("data/", sample, "_cropped")
  print(data_dir)
  
  # Create seurat object and set orig ident to be sample
  seurat <- Load10X_Spatial(data.dir = data_dir,
                            bin.size = 8,
                            slice = sample)
  
  seurat$orig.ident <- sample
  list_seurat[[sample]] <- seurat
  rm(seurat)
}

#merging two tissues

seurat_merged <- merge(x = list_seurat[["P5CRC"]],
                       y = list_seurat[["P5NAT"]],
                       add.cell.id = c("P5CRC", "P5NAT"))


# combining layers in a single one
seurat_merged <- JoinLayers(seurat_merged)
print(seurat_merged)

print(ncol(list_seurat[["P5CRC"]]) + ncol(list_seurat[["P5NAT"]]))
print(ncol(seurat_merged))
print(seurat_merged@meta.data %>% head())
print(seurat_merged@meta.data %>% tail())


# Set Idents to sample IDs - setting groups as default tissues. If we need to compare groups they will be compared by the tissue.
Idents(seurat_merged) <- "orig.ident"


saveRDS(seurat_merged, "data/01_seurat_merged.RDS")
rm(list_seurat)
rm(seurat_merged)
