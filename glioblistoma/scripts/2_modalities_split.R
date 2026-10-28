rna_counts <- LayerData(gbm, assay="Spatial", layer="counts.Gene Expression")
protein_counts  <- LayerData(gbm, assay="Spatial", layer="counts.Antibody Capture")

gbm[["RNA"]] <- CreateAssayObject(counts = rna_counts)
gbm[["Protein"]] <- CreateAssayObject(counts = protein_counts)

rm(rna_counts)
rm(protein_counts)
