library(curl)

url_barcodes = "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_filtered_feature_bc_matrix.h5"#"https://cf.10xgenomics.com/samples/cell-arc/2.2.0/10k_Human_Brain_MO_gemx/10k_Human_Brain_MO_gemx_filtered_feature_bc_matrix.h5"
#url_atac_beds = "https://cf.10xgenomics.com/samples/cell-arc/2.2.0/10k_Human_Brain_MO_gemx/10k_Human_Brain_MO_gemx_atac_peaks.bed"
#url_atac_fragments_index = "https://cf.10xgenomics.com/samples/cell-arc/2.2.0/10k_Human_Brain_MO_gemx/10k_Human_Brain_MO_gemx_atac_fragments.tsv.gz.tbi"
#url_atac_fragments = "https://cf.10xgenomics.com/samples/cell-arc/2.2.0/10k_Human_Brain_MO_gemx/10k_Human_Brain_MO_gemx_atac_fragments.tsv.gz" 
url_spatial = "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma/CytAssist_FFPE_Protein_Expression_Human_Glioblastoma_spatial.tar.gz"


urls <- list(url_barcodes,url_spatial)#list(url_barcodes, url_atac_beds, url_atac_fragments_index, url_atac_fragments)

dir.create("data/glioblastoma/",recursive = TRUE)



for (URL in urls){
  dest <- file.path("data/glioblastoma", basename(URL))
  curl_download(
    url = URL,
    destfile = dest,
    mode = "wb")
  }
untar(
  file.path("data/glioblastoma", basename(url_spatial)),
  exdir = "data/glioblastoma/"
)