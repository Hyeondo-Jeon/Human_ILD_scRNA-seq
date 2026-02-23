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
## BAL fluid data preprocessing
## ------------------------------------------------------------

sample_ids <- c('BALF_UIP_1','BALF_DIP_1','BALF_UIP_2')

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

## ------------------------------------------------------------
## BAL fluid data merge
## ------------------------------------------------------------
Bmerge <- merge(BAL_HP_1, 
                y = c(BAL_DIP_1, BAL_UIP_1), 
                add.cell.ids = c("D1", "N1", "U1")
               )

## scRNA-seq QC metrics by sample
VlnPlot(Bmerge, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3,
        group.by = "orig.ident",
        pt.size = 0
       )

## subtype info
Bmerge$diagnosis <- "a"
Bmerge$diagnosis[Bmerge$orig.ident %in% c("BAL_DIP_1")] <- "DIP"
Bmerge$diagnosis[Bmerge$orig.ident %in% c("BAL_UIP2")] <- "UIP"
Bmerge$diagnosis[Bmerge$orig.ident %in% c("BAL_UIP_1")] <- "UIP"

## processing
Bmerge <- JoinLayers(Bmerge)
Bmerge <- NormalizeData(Bmerge)
Bmerge <- FindVariableFeatures(Bmerge, selection.method = "vst", nfeatures = 3200) 
Bmerge <- ScaleData(Bmerge, features = rownames(Bmerge))
Bmerge <- RunPCA(Bmerge, features = VariableFeatures(Bmerge))
ElbowPlot(Bmerge, ndims = 50)
Bmerge <- FindNeighbors(Bmerge, dims = 1:15)
Bmerge <- FindClusters(Bmerge, resolution = 1.0)
Bmerge <- RunUMAP(Bmerge, dims = 1:15)
DimPlot(Bmerge, label = T)

## DEGs
balf <- FindAllMarkers(Bmerge, test.use = "wilcox", min.pct = 0.5, logfc.threshold = 1, group.by = "seurat_clusters")
FeaturePlot(Bmerge, c("PTPRC","CD79A", "MS4A1","BANK1","BLK","MZB1","JCHAIN","CD3E","MKI67", "TOP2A","STMN1","CD4","FOXP3","IL2RA", "CTLA4","CCR7", "IL7R",
                      "SELL","CXCR6","CD8A", "GZMB", "PRF1","NKG7", "NCAM1", "FCGR3A","CD14", "S100A8", "S100A9","SPP1","CLEC9A","XCR1",
                      "CLEC10A","CD1C","CLEC4C","IRF7","LAMP3","FSCN1","C1QA","CDK1","PCLAF","TREM2","MRC1","FABP4","FFAR4","TPSAB1","KIT", "CPA3","CSF3R","FCGR3B","NAMPT"))

## cell type annotation
Bmerge$cell <- "a"
Bmerge$cell[Bmerge$seurat_clusters %in% c("13")] <- "Neutrophil"
Bmerge$cell[Bmerge$seurat_clusters %in% c("25")] <- "Mast cell"
Bmerge$cell[Bmerge$seurat_clusters %in% c("24")] <- "Cycling MΦ"
Bmerge$cell[Bmerge$seurat_clusters %in% c("18")] <- "Activated B"
Bmerge$cell[Bmerge$seurat_clusters %in% c("23")] <- "Plasma"
Bmerge$cell[Bmerge$seurat_clusters %in% c("21")] <- "pDC"
Bmerge$cell[Bmerge$seurat_clusters %in% c("19")] <- "cDC1"
Bmerge$cell[Bmerge$seurat_clusters %in% c("17")] <- "cDC2"
Bmerge$cell[Bmerge$seurat_clusters %in% c("22")] <- "mDC"
Bmerge$cell[Bmerge$seurat_clusters %in% c("6")] <- "CD4+ Treg"
Bmerge$cell[Bmerge$seurat_clusters %in% c("7")] <- "CD4+ Naive"
Bmerge$cell[Bmerge$seurat_clusters %in% c("2")] <- "CD4+ Effector memory"
Bmerge$cell[Bmerge$seurat_clusters %in% c("1","15")] <- "CD8+ Memory"
Bmerge$cell[Bmerge$seurat_clusters %in% c("0")] <- "CD8+ CTL"
Bmerge$cell[Bmerge$seurat_clusters %in% c("14")] <- "Cycling T"
Bmerge$cell[Bmerge$seurat_clusters %in% c("12")] <- "CD56 bright NK"
Bmerge$cell[Bmerge$seurat_clusters %in% c("11")] <- "LAM"
Bmerge$cell[Bmerge$seurat_clusters %in% c("10","5","16","9")] <- "Alveolar MΦ"
Bmerge$cell[Bmerge$seurat_clusters %in% c("20")] <- "SPP1+ MΦ"
Bmerge$cell[Bmerge$seurat_clusters %in% c("3","4")] <- "SPP1+ Mono"
Bmerge$cell[Bmerge$seurat_clusters %in% c("8")] <- "iMono"
Idents(Bmerge) <- Bmerge$cell
DimPlot(Bmerge, label = T, label.size = 5, repel = T) 

## dot plot for marker genes
Idents(Bmerge) <- factor(Bmerge$cell,levels = c("Activated B","Plasma","Cycling T","CD4+ Treg","CD4+ Naive","CD4+ Effector memory","CD8+ Memory","CD8+ CTL",
                                                 "CD56 bright NK","SPP1+ Mono","iMono","cDC1","cDC2","pDC","mDC","Cycling MΦ","LAM","SPP1+ MΦ","Alveolar MΦ","Mast cell","Neutrophil"))
marker_list_bal <- c("PTPRC","CD79A", "MS4A1","BANK1","BLK", # Activated B
                     "MZB1","JCHAIN", # Plasma
                     "CD3E","MKI67", "TOP2A","STMN1", # Cycling T
                     "CD4","FOXP3","IL2RA", "CTLA4", # CD4+ Treg
                     "CCR7","IL7R","SELL", # CD4+ Naive
                     "CXCR6", # CD4+ Effector memory
                     "CD8A","GZMK","CD28", # CD8+ Memory
                     "GZMB", "PRF1", # CD8+ CTL
                     "NKG7","NCAM1", "FCGR3A", # CD56 bright NK
                     "CD14", "S100A8", "S100A9","SPP1", # SPP1+ Mono & iMono
                     "CLEC9A","XCR1", # cDC1
                     "CLEC10A","CD1C", # cDC2
                     "CLEC4C","IRF7", # pDC
                     "LAMP3","FSCN1", # mDC
                     "C1QA","CDK1","PCLAF", # Cycling MΦ
                     "TREM2","MRC1", # LAM
                     "FABP4","FFAR4", # Alveolar MΦ
                     "TPSAB1","KIT", "CPA3", # Mast cell
                     "CSF3R","FCGR3B","NAMPT" # Neutrophil
                    )
DotPlot(Bmerge, features = marker_list_bal, cols = c("lightgray", "slateblue4")) + 
  labs(x = " ", y = " ") +
  theme(
    axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 16, face = "italic"), 
    axis.text.y = ggtext::element_markdown(hjust = 1, size = 18), 
    axis.title = element_text(size = 26),                        
    legend.title = element_text(size = 16),                     
    legend.text = element_text(size = 18)
  )
