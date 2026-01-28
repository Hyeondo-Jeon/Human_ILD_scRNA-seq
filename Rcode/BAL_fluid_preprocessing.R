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
library(DoubletFinder)
library(scDblFinder)

## ------------------------------------------------------------
## PBMC data preprocessing
## ------------------------------------------------------------

sample_ids <- c('BALF_UIP_1','BALF_DIP_1','BALF_HP_1')

set.seed(1234)

for (sid in sample_ids) {
  
  message("Processing ", sid)
  
  data_dir <- file.path(
    "data",
    "single_cell_balf",
    sid,
    "CellRanger",
    "filtered_feature_bc_matrix"
  )
  
  seu <- Read10X(data.dir = data_dir)
  seu <- CreateSeuratObject(seu, project = sid, min.cells = 3, min.features = 0)
  
  ## mitochondrial percentage
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  ## QC filtering
  seu <- subset(
    seu, 
    subset = 
      nFeature_RNA > 200 &
      nFeature_RNA < 9000 &
      percent.mt < 5
  )
  
  ## normalization & clustering
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu, features = rownames(seu))
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:10)
  
  ## DoubletFinder
  set.seed(1234)
  sweep <- paramSweep(seu, PCs = 1:10, sct = FALSE)
  sweep <- summarizeSweep(sweep, GT = FALSE)
  bcmvn <- find.pK(sweep)
  
  homotypic.prop <- modelHomotypic(seu$seurat_clusters)
  
  expected_doublet_rate <- 0.076 * (ncol(seu) / 10000)
  nEXP_poi <- round(expected_doublet_rate * ncol(seu))
  nEXP_poi.adj <- round(nEXP_poi * (1 - homotypic.prop))
  
  pK <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]
  ))
  
  seu <- doubletFinder(
    seu,
    PCs = 1:10,
    pN = 0.25,
    pK = pK,
    nExp = nEXP_poi.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )
  
  seu@meta.data <- seu@meta.data %>%
    rename_with(~"DoubletFinder", matches("^DF"))
  
  ## DoubletFinder summary & visualization (Before filtering)
  df_summary <- table(seu$DoubletFinder)
  print(df_summary)
  
  DF <- DimPlot(
    seu,
    group.by = "DoubletFinder",
    cols = c("salmon2", "darkslategray3")
  ) +
    ggtitle(
      paste0(
        "DoubletFinder ", sid,
        " | Singlet: ", df_summary["Singlet"],
        " | Doublet: ", df_summary["Doublet"]
      )
    )
  print(DF)
  
  ## scDblFinder
  set.seed(1234)
  scdbl <- scDblFinder(
    GetAssayData(seu, slot = "counts"),
    clusters = Idents(seu)
  )
  seu$scDblFinder <- scdbl$scDblFinder.class
  
  ## scDblFinder summary & visualization (Before filtering)
  scdbl_summary <- table(seu$scDblFinder)
  print(scdbl_summary)
  
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
  
  ## Union doublet definition
  seu$Doublet <- "Doublet"
  seu$Doublet[
    seu$DoubletFinder == "Singlet" &
      seu$scDblFinder == "singlet"
  ] <- "Singlet"
  
  dbl_summary <- table(seu$Doublet)
  print(dbl_summary)
  
  unionDB <- DimPlot(
    seu,
    group.by = "Doublet",
    cols = c("salmon2", "darkslategray3")
  ) +
    ggtitle(
      paste0(
        "Doublet ", sid,
        " | Singlet: ", dbl_summary["Singlet"],
        " | Doublet: ", dbl_summary["Doublet"]
      )
    )
  print(unionDB)
  
  ## retain singlets only
  seu <- subset(seu, Doublet == "Singlet")
  
  ## GAPDH / ACTB filtering
  hk <- FetchData(seu, vars = c("GAPDH", "ACTB"))
  keep <- rownames(hk)[!(hk$GAPDH == 0 & hk$ACTB == 0)]
  seu <- subset(seu, cells = keep)
  
  ## seurat object creation
  assign(sid, seu, envir = .GlobalEnv)
  message(
    "Completed: ", sid,
    " | Cells retained: ", ncol(seu)
  )
}
