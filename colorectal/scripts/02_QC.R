
library(ggplot2)


seurat_merged <- readRDS("data/01_seurat_merged.rds")

#speaking about identities
# print(ggplot(seurat_merged@meta.data) +
#   geom_bar(aes(x = orig.ident, fill = orig.ident),
#            color = "black") +
#   geom_text(aes(x = orig.ident, label=after_stat(count)), 
#             stat='count', vjust=-1) +
#   theme_classic()
# )

meta <- seurat_merged@meta.data


seurat_merged$log10GenesPerUMI <- log10(seurat_merged$nFeature_Spatial.008um) / log10(seurat_merged$nCount_Spatial.008um)
seurat_merged$log10GenesPerUMI[is.na(seurat_merged$log10GenesPerUMI)] <- 0


seurat_merged$mitoRatio <- PercentageFeatureSet(object = seurat_merged, 
                                                pattern = "^MT-")
seurat_merged$mitoRatio <- seurat_merged@meta.data$mitoRatio / 100
seurat_merged$mitoRatio[is.na(seurat_merged$mitoRatio)] <- 1.00


features_list <- c("nCount_Spatial.008um", "nFeature_Spatial.008um", "mitoRatio", "log10GenesPerUMI")

for (i in features_list){
 p <- SpatialFeaturePlot(seurat_merged,
                      i,
                      pt.size.factor = 15,
                      image.alpha = 0,
                      max.cutoff = "q90")
  
  ggsave(
    filename = paste0("figures/spatial_",i, ".png"),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  meta <- seurat_merged@meta.data
  p <- ggplot(meta) +
    geom_density(aes(x = .data[[i]],
                     fill = orig.ident),
                 alpha = 0.4,
                 color = "black")+
    scale_x_log10() +
    theme_classic()
  
  if (i %in% c("nCount_Spatial.008um",
               "nFeature_Spatial.008um")) {
    
    p <- p +
      scale_x_log10() +
      geom_vline(xintercept = 10, color = "lightblue") +
      geom_vline(xintercept = 30, color = "pink")
    
  } else if (i == "mitoRatio") {
    
    p <- p +
      geom_vline(xintercept = 0.25, color = "red")
    
  } else if (i == "log10GenesPerUMI") {
    
    p <- p +
      geom_vline(xintercept = 0.80, color = "red")
  }
  
  ggsave(
    filename = paste0("figures/before_filtration_",i, ".png"),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  
}


for (i in features_list){
  if (i %in% c("nCount_Spatial.008um",
               "nFeature_Spatial.008um")) {
  meta_filt <- subset(meta,
                      ((orig.ident == "P5CRC") & (meta[[i]] > 30)) |
                        ((orig.ident == "P5NAT") & (meta[[i]] > 10)))
  }
  else if (i == "mitoRatio") {
    
    meta_filt <- subset(meta, mitoRatio < 0.25)
    
  } else if (i == "log10GenesPerUMI") {
    
    meta_filt <- subset(meta, log10GenesPerUMI > 0.80)
  }
  
  p<-ggplot(meta_filt) +
    geom_density(aes(x = .data[[i]],
                     fill = orig.ident),
                 alpha = 0.4,
                 color = "black") +
    theme_classic()
  if (i %in% c("nCount_Spatial.008um",
               "nFeature_Spatial.008um")) {
    
    p <- p +
      scale_x_log10() +
      geom_vline(xintercept = 10, color = "lightblue") +
      geom_vline(xintercept = 30, color = "pink")
    
  } else if (i == "mitoRatio") {
    
    p <- p +
      geom_vline(xintercept = 0.25, color = "red")
    
  } else if (i == "log10GenesPerUMI") {
    
    p <- p +
      geom_vline(xintercept = 0.80, color = "red")
  }
  
  ggsave(
    filename = paste0("figures/after_filtration_",i, ".png"),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  
}



seurat_filtered <- subset(seurat_merged,
                          ((orig.ident == "P5CRC") & (nCount_Spatial.008um > 30)) |
                            ((orig.ident == "P5NAT") & (nCount_Spatial.008um > 10)))

seurat_filtered <- subset(seurat_filtered,
                          ((orig.ident == "P5CRC") & (nFeature_Spatial.008um > 30)) |
                            ((orig.ident == "P5NAT") & (nFeature_Spatial.008um > 10)))

seurat_filtered <- subset(seurat_filtered, mitoRatio < 0.25)
seurat_filtered <- subset(seurat_filtered, log10GenesPerUMI > 0.80)




p_ncount <- VlnPlot(seurat_filtered, 
                    features = "nCount_Spatial.008um", 
                    pt.size = 0, group.by = 'orig.ident') +
  NoLegend()

# Violin plot of number of genes
p_nfeats <- VlnPlot(seurat_filtered, 
                    features = "nFeature_Spatial.008um", 
                    pt.size = 0, group.by = 'orig.ident') + 
  NoLegend()

# Plot UMIs and gene count violin plots side-by-side
p <- p_ncount | p_nfeats

ggsave(
  filename = "figures/violins.png",
  plot = p,
  width = 12,
  height = 5,
  dpi = 300
)


p <- SpatialFeaturePlot(seurat_filtered, 
                   c("nFeature_Spatial.008um", 
                     "nCount_Spatial.008um"),
                   pt.size.factor = 16,
                   image.alpha = 0)

ggsave(
  filename = paste0("figures/filtered_spatial", ".png"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)

saveRDS(seurat_filtered, "data/seurat_filtered.RDS")

   
rm(seurat_merged)
rm(seurat_filtered)
rm(meta)
rm(meta_filt)
