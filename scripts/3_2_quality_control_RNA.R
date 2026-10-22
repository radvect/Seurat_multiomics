library(ggplot2)
library(Seurat)



number_of_genes <-SpatialFeaturePlot(gbm, 
                                     "nFeature_Spatial",
                                     pt.size.factor = 2,
                                     image.alpha = 0,
                                     max.cutoff = "q90")


density_nfeatures <-ggplot(gbm@meta.data) +
  geom_density(aes(x = nCount_Spatial, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  scale_x_log10() +
  theme_classic()
ggsave("figures/density_ngenes.png", plot = density_nfeatures,
       width = 10, height = 7, dpi = 300)


#how to improve the pipelinex