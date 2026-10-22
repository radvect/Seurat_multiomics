#TODO
#CHECK THE TILE SIZE
library(ggplot2)
library(Seurat)

plot_bins <- ggplot(gbm@meta.data) +
  geom_bar(
    aes(x = orig.ident, fill = orig.ident),
    color = "black"
  ) +
  geom_text(
    aes(x = orig.ident, label = after_stat(count)),
    stat = "count",
    vjust = -1
  ) +
  theme_classic() +
  labs(
    title = "Number of bins per sample",
    x = "Sample",
    y = "Number of bins"
  )


ggsave("figures/bins.png", plot = plot_bins,
       width = 10, height = 7, dpi = 300)


spatial_bins <-SpatialFeaturePlot(gbm, 
                   "nCount_Spatial",
                   pt.size.factor = 2,
                   image.alpha = 0,
                   max.cutoff = "q90")


ggsave("figures/umi_number.png", plot = spatial_bins,
       width = 10, height = 7, dpi = 300)


density_ncount<-ggplot(gbm@meta.data) +
  geom_density(aes(x = nCount_Spatial, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  scale_x_log10() +
  theme_classic()

ggsave("figures/density_ncount.png", plot = density_ncount,
       width = 10, height = 7, dpi = 300)

