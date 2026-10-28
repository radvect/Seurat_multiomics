library(Seurat)

gbm <- FindVariableFeatures(gbm,
                            assay = "RNA",
                            selection.method = "vst",
                            nfeatures = 2000)
#ranked_genes <- VariableFeatures(gbm)

#length(VariableFeatures(gbm[["RNA"]]))
#ranked_genes

ranked_genes <- VariableFeatures(gbm[["RNA"]])
top_genes <- ranked_genes[1:15]
rm(ranked_genes)

gbm <- FindVariableFeatures(gbm,
                            assay = "Protein",
                            selection.method = "vst",
                            nfeatures = 2000)
ranked_proteins <- VariableFeatures(gbm[["Protein"]])
top_proteins <-ranked_proteins[1:15]
rm(ranked_proteins)


DefaultAssay(gbm) <- "RNA"

plot_RNA <- VariableFeaturePlot(gbm)

LabelPoints(
  plot = plot_RNA,
  points = top_genes,
  repel = TRUE
)

DefaultAssay(gbm) <- "Protein"

plot_proteins <- VariableFeaturePlot(gbm)

LabelPoints(
  plot = plot_proteins,
  points = top_proteins,
  repel = TRUE
)


plot_RNA <- LabelPoints(
  plot = plot_RNA,
  points = top_genes,
  repel = TRUE,
  size = 4
) + theme_classic()

plot_proteins <- LabelPoints(
  plot = plot_proteins,
  points = top_proteins,
  repel = TRUE,
  size = 4
)+ theme_classic()




ggsave("figures/RNA_HVG.png", plot = plot_RNA,
       width = 10, height = 7, dpi = 300)

ggsave("figures/Protein_HVG.png", plot = plot_proteins,
       width = 10, height = 7, dpi = 300)