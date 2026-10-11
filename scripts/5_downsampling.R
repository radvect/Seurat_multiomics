
#PCA the reduced number of cells based on the Leverage score. 
#SVG 
gbm <- SketchData(object = gbm,
                              assay = "RNA",
                              ncells = 2000,
                              method = "LeverageScore",
                              sketched.assay = "sketch"
                            )
#question: how can we accommodate proteins and RNA together.