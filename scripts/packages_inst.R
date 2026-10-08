install.packages("tidyverse")
library(tidyverse)
install.packages("Seurat")
library(Seurat)
install.packages("curl")
library("curl")
install.packages("Signac")
install.packages("hdf5r")
BiocManager::install(
  c(
    "GenomeInfoDb",
    "GenomicRanges",
    "IRanges",
    "S4Vectors"
  ),
  update = FALSE,
  ask = FALSE
)
install.packages("BiocManager")
