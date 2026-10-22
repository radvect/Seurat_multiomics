
#PCA the reduced number of cells based on the Leverage score. 
#SVG 
library(Seurat)


gbm <- SketchData(object = gbm,
                              assay = "RNA",
                              ncells = 2000,
                              method = "LeverageScore",
                              sketched.assay = "sketch"
                            )
#question: how can we accommodate proteins and RNA together.
#how can we include the spatial information as well? 

library(ggplot2)

ggplot(gbm@meta.data) +
  geom_histogram(aes(x=leverage.score, 
                     fill=orig.ident), 
                 alpha=0.5, bins=100) +
  theme_classic() +
  scale_x_log10(labels = scales::label_number())


