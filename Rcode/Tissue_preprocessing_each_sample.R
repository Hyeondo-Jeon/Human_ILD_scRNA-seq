## ------------------------------------------------------------
## library load
## ------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(umap)
library(ComplexHeatmap)
library(DGEobj.utils)
library(rbioapi)
library(ggplot2)
library(dplyr)
library("readxl")
library(scales)
library(scDblFinder)

## ------------------------------------------------------------
## Tissue data preprocessing
## ------------------------------------------------------------

sample_ids <- c("UIP_1", "UIP_2", "UIP_3", "UIP_4", "DIP_1", "NSIP_1")

set.seed(1234)

for (sid in sample_ids) {
  
  message("Processing ", sid)
  
  data_dir <- file.path(
    "data",
    "single_cell_tissue",
    sid,
    "CellRanger",
    "filtered_feature_bc_matrix"
  )
  
  seu <- Read10X(data.dir = data_dir)
  
  seu <- CreateSeuratObject(
    counts = seu,
    project = sid,
    min.cells = 3,
    min.features = 0
  )
  
  ## mitochondrial percentage
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  ## QC filtering
  seu <- subset(
    seu,
    subset =
      nFeature_RNA > 1000 &
      nFeature_RNA < 9000 &
      percent.mt < 5
  )
  
  ## normalization & clustering
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
  seu <- ScaleData(seu, features = rownames(seu))
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:10)
  
  ## doublet detection
  set.seed(1234)
  scdbl <- scDblFinder(
    GetAssayData(seu, slot = "counts"),
    clusters = Idents(seu)
  )
  
  seu$scDblFinder <- scdbl$scDblFinder.class
  
  ## doublet summary & visualization (Before filtering)
  dbl_summary <- table(seu$scDblFinder)
  print(dbl_summary)
  
  scDB <- DimPlot(
    seu,
    group.by = "scDblFinder",
    cols = c("darkslategray3", "salmon2")
  ) +
    ggtitle(
      paste0(
        "scDblFinder ", sid,
        " | Singlet: ", scdbl_summary["singlet"],
        " | Doublet: ", scdbl_summary["doublet"]
      )
    )
  print(scDB)
  
  ## retain singlets only
  seu <- subset(seu, scDblFinder == "singlet")
  
  ## GAPDH / ACTB filtering
  hk <- FetchData(seu, vars = c("GAPDH", "ACTB"))
  valid_cells <- rownames(hk)[!(hk$GAPDH == 0 & hk$ACTB == 0)]
  seu <- subset(seu, cells = valid_cells)
  
  ## seurat object creation
  assign(sid, seu, envir = .GlobalEnv)
  message(
      "Completed: ", sid,
      " | Cells retained: ", ncol(seu)
  )
}
