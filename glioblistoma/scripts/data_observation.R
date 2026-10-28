gbm
head(rownames(gbm))
tail(rownames(gbm))
Layers(gbm[["Spatial"]])


#metadata load
#
gbm@meta.data %>% View()