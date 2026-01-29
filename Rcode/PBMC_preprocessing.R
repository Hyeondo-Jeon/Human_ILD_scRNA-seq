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

sample_ids <- c('PBMC_UIP_1','PBMC_UIP_2','PBMC_UIP_3','PBMC_UIP_4','PBMC_UIP_5','PBMC_DIP_1')

set.seed(1234)

for (sid in sample_ids) {
  
  message("Processing ", sid)
  
  data_dir <- file.path(
    "data",
    "single_cell_pbmc",
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
      nFeature_RNA > 500 &
      nFeature_RNA < 6000 &
      percent.mt < 5
    )
  
  ## normalization & clustering
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu, features = rownames(seu))
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:15)
  
  ## DoubletFinder
  set.seed(1234)
  sweep <- paramSweep(seu, PCs = 1:15, sct = FALSE)
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
    PCs = 1:15,
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

## ------------------------------------------------------------
## PBMC data merge
## ------------------------------------------------------------
Pmerge <- merge(pbmc_DIP_1, 
                y = c(pbmc_UIP_1, pbmc_UIP_2, pbmc_UIP_3, pbmc_UIP_4, pbmc_UIP_5), 
                add.cell.ids = c("D1", "U1", "U2", "U3", "U4","U5")
               )
## subtype info
Pmerge$diagnosis <- "a"
Pmerge$diagnosis[Pmerge$orig.ident %in% c("PBMC_UIP_1","PBMC_UIP_2","PBMC_UIP_3","PBMC_UIP_4","PBMC_UIP_5")] <- "UIP"
Pmerge$diagnosis[Pmerge$orig.ident %in% c("PBMC_DIP_1")] <- "DIP"

## scRNA-seq QC metrics by sample
VlnPlot(Pmerge, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "orig.ident", 
        pt.size = 0
       )

## processing
Pmerge <- JoinLayers(Pmerge)
Pmerge <- NormalizeData(Pmerge)
Pmerge <- FindVariableFeatures(Pmerge, selection.method = "vst", nfeatures = 2800)
Pmerge <- ScaleData(Pmerge, features = rownames(Pmerge))
Pmerge <- RunPCA(Pmerge, features = VariableFeatures(Pmerge))
ElbowPlot(Pmerge, ndims = 50)
Pmerge <- FindNeighbors(Pmerge, dims = 1:15)
Pmerge <- FindClusters(Pmerge, resolution = 1.0)
Pmerge <- RunUMAP(Pmerge, dims = 1:15)
DimPlot(Pmerge, label = T)

## DEGs
pbmc <- FindAllMarkers(Pmerge, test.use = "wilcox", min.pct = 0.5, logfc.threshold = 1, group.by = "seurat_clusters")
FeaturePlot(Pmerge,  c("PTPRC","CD79A", "MS4A1","BANK1","FCER2","TCL1A","JCHAIN","CD3E","MKI67", "TOP2A","STMN1","CD4","CCR7","SELL","IL7R","FOXP3","CTLA4",
                       "CD8A", "TCF7","NELL2","GZMK","CD28","GZMB", "PRF1","NKG7", "NCAM1", "TYROBP","FCGR3A","LILRA1", "CX3CR1","CD14", "S100A8", "S100A9","CLEC10A","CD1C",
                       "CLEC4C","IRF7"))

## cell type annotation
Pmerge <- FindSubCluster(Pmerge, cluster = "24", graph.name = "RNA_snn", subcluster.name = "subcluster",resolution = 0.1)
DimPlot(Pmerge, label = T, group.by = "subcluster")
Pmerge$cell <- "a"
Pmerge$cell[Pmerge$seurat_clusters %in% c("13","21","4")] <- "Naive B"
Pmerge$cell[Pmerge$seurat_clusters %in% c("15","7")] <- "Activated B"
Pmerge$cell[Pmerge$seurat_clusters %in% c("25")] <- "Plasma"
Pmerge$cell[Pmerge$seurat_clusters %in% c("26")] <- "Cycling T"
Pmerge$cell[Pmerge$seurat_clusters %in% c("16","12","8","9")] <- "CD4+ Naive"
Pmerge$cell[Pmerge$seurat_clusters %in% c("22")] <- "CD4+ Treg"
Pmerge$cell[Pmerge$seurat_clusters %in% c("11")] <- "CD8+ Naive"
Pmerge$cell[Pmerge$seurat_clusters %in% c("6")] <- "CD8+ Memory"
Pmerge$cell[Pmerge$seurat_clusters %in% c("0","14","17")] <- "CD8+ CTL"
Pmerge$cell[Pmerge$seurat_clusters %in% c("2","18")] <- "NK-like T"
Pmerge$cell[Pmerge$subcluster %in% c("20","24_2")] <- "CD56 bright NK"
Pmerge$cell[Pmerge$seurat_clusters %in% c("1","3","10")] <- "CD56 dim NK"
Pmerge$cell[Pmerge$seurat_clusters %in% c("5","23","27")] <- "cMono"
Pmerge$cell[Pmerge$seurat_clusters %in% c("19")] <- "Non-cMono"
Pmerge$cell[Pmerge$subcluster %in% c("24_0")] <- "cDC2"
Pmerge$cell[Pmerge$subcluster %in% c("24_1")] <- "pDC"
Idents(Pmerge) <- Pmerge$cell
DimPlot(Pmerge, label = T, label.size = 5, repel = T)

## dot plot of marker genes
Idents(Pmerge) <- factor(Pmerge$cell,levels = c("Activated B","Naive B","Plasma","Cycling T","CD4+ Naive","CD4+ Treg","CD8+ Naive","CD8+ Memory",
                                                "CD8+ CTL","NK-like T","CD56 dim NK","CD56 bright NK","Non-cMono","cMono","cDC2","pDC"))
marker_list_pbmc <- c("PTPRC","CD79A", "MS4A1","BANK1", # Activated B
                      "FCER2","TCL1A", # Naive B
                      "JCHAIN", # Plasma
                      "CD3E","MKI67", "TOP2A", "STMN1", # Cycling T
                      "CD4","CCR7","SELL","IL7R", # CD4+ Naive
                      "FOXP3","CTLA4", # CD4+ Treg
                      "CD8A", "TCF7","NELL2", # CD8+ Memory
                      "GZMK","CD28", # CD8+ CTL
                      "GZMB", "PRF1","TYROBP", # NK-like T
                      "NKG7", "NCAM1", # CD56 dim NK & CD56 bright NK
                      "FCGR3A","LILRA1", "CX3CR1",  # Non-cMono
                      "CD14", "S100A8", "S100A9", # cMono
                      "CLEC10A","CD1C", # cDC2
                      "CLEC4C","IRF7" # pDC
                     )
DotPlot(Pmerge, features = marker_list_pbmc, cols = c("lightgray", "slateblue4")) + 
  labs(x = " ", y = " ") +
  theme(
    axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 16, face = "italic"), 
    axis.text.y = ggtext::element_markdown(hjust = 1, size = 18), 
    axis.title = element_text(size = 26),                        
    legend.title = element_text(size = 16),                     
    legend.text = element_text(size = 18)
  )
