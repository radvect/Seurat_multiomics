library(ggplot2)
library(Seurat)




spatial_bins <-SpatialFeaturePlot(gbm, 
                                  "nCount_Protein",
                                  pt.size.factor = 2,
                                  image.alpha = 0,
                                  max.cutoff = "q90")


ggsave("figures/proteins_spatial.png", plot = spatial_bins,
       width = 10, height = 7, dpi = 300)




number_of_genes <-SpatialFeaturePlot(gbm, 
                                     "nFeature_Protein",
                                     pt.size.factor = 2,
                                     image.alpha = 0,
                                     max.cutoff = "q90")


density_nprotein <-ggplot(gbm@meta.data) +
  geom_density(aes(x = nCount_Spatial, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  scale_x_log10() +
  theme_classic()
ggsave("figures/density_nprotein.png", plot = density_nprotein,
       width = 10, height = 7, dpi = 300)
