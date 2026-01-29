## ------------------------------------------------------------
## Tissue data merge
## ------------------------------------------------------------
Tmerge <- merge(UIP_1, 
                y = c(UIP_2, UIP_3, UIP_4, DIP_1, NSIP_1), 
                add.cell.ids = c("U1", "U2", "U3", "U4", "D1","N1")
               )

## subtype info
Tmerge$diagnosis <- "a"
Tmerge$diagnosis[Tmerge$orig.ident %in% c("UIP_1","UIP_2","UIP_3","UIP_4")] <- "UIP"
Tmerge$diagnosis[Tmerge$orig.ident %in% c("NSIP_1")] <- "NSIP"
Tmerge$diagnosis[Tmerge$orig.ident %in% c("DIP_1")] <- "DIP"
Tmerge$diagnosis <- as.factor(Tmerge$diagnosis)

## scRNA-seq QC metrics by sample
VlnPlot(Tmerge,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3,
        group.by = "orig.ident",
        pt.size = 0
       )

## processing
Tmerge <- JoinLayers(Tmerge)
Tmerge <- NormalizeData(Tmerge)
Tmerge <- FindVariableFeatures(Tmerge, selection.method = "vst", nfeatures = 3800)
Tmerge <- ScaleData(Tmerge, features = rownames(Tmerge))
Tmerge <- RunPCA(Tmerge, features = VariableFeatures(Tmerge))
ElbowPlot(Tmerge, ndims = 50)
Tmerge <- FindNeighbors(Tmerge, dims = 1:15)
Tmerge <- FindClusters(Tmerge, resolution = 1.5)
Tmerge <- RunUMAP(Tmerge, dims = 1:15)
DimPlot(Tmerge, label = T) + ggtitle("")

## DEGs
tissue <- FindAllMarkers(Tmerge, test.use = "wilcox", min.pct = 0.5, logfc.threshold = 1, group.by = "seurat_clusters")
FeaturePlot(Tmerge, c("CLEC4C","PTPRC","CD3E","MS4A1","JCHAIN","EPCAM","COL1A1","CDH5","NRCAM"))

## lineage cell type
Tmerge$population <- "a"
Tmerge$population[Tmerge$seurat_clusters %in% c("1","0","40","38","39","36","19","35","33","13")] <- "Lymphoid"
Tmerge$population[Tmerge$seurat_clusters %in% c("22","17","28","30","25","27","6","23","11","3","15","8","37")] <- "Myeloid"
Tmerge$population[Tmerge$seurat_clusters %in% c("20","26","9","41","21","29")] <- "Epithelial"
Tmerge$population[Tmerge$seurat_clusters %in% c("14","18","34","12","10","16","34")] <- "Stromal"
Tmerge$population[Tmerge$seurat_clusters %in% c("32","4","2","5","7","24","42")] <- "Endothelial"
Tmerge$population[Tmerge$seurat_clusters %in% c("31")] <- "Unknown"
DimPlot(Tmerge, label = T, group.by = "population") + ggtitle("")
Tmerge$cell <- "a"
Tmerge$cell[Tmerge$seurat_clusters %in% c("31")] <- "Unknown"

## dot plot for marker genes
marker_list <- c("CDH5","PECAM1","EPCAM","MUC1","PTPRC","CD3E","CD79A","CD14","ITGAX","COL1A1","COL3A1","MAP1B","NRCAM")
DotPlot(Tmerge, features = marker_list, cols = c("lightgray", "slateblue4"), group.by = "population") + 
  labs(x = " ", y = " ") +
  theme(
    axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 14), 
    axis.text.y = ggtext::element_markdown(hjust = 1, size = 16), 
    axis.title = element_text(size = 26),                        
    legend.title = element_text(size = 14),                     
    legend.text = element_text(size = 14)
  )

## tissue lymphoid sub clustering
Lym <- subset(Tmerge, population %in% "Lymphoid")
Lym <- NormalizeData(Lym)
Lym <- FindVariableFeatures(Lym, selection.method = "vst", nfeatures = 2600)
Lym <- ScaleData(Lym, features = rownames(Lym))
Lym <- RunPCA(Lym, features = VariableFeatures(Lym))
Lym <- FindNeighbors(Lym, dims = 1:10)
Lym <- FindClusters(Lym, resolution = 1.5)
Lym <- RunUMAP(Lym, dims = 1:10)
DimPlot(Lym, label = T)
FeaturePlot(Lym,  c("PTPRC","MS4A1","CD19","BANK1","BLK","FCER2",
                    "TCL1A","MME","DTX1","MKI67","TOP2A","MZB1","JCHAIN",
                    "CD3E","CD4","FOXP3","CTLA4","IL7R","SELL","CCR7",
                    "CD8A","CCL5","GZMK","CD28","GZMB","GZMK","PRF1",
                    "NKG7","NCAM1","FCGR3A","FGFBP2","GZMH","KIT","PCDH9"))

## cell type annotation
Lym$cell[Lym$seurat_clusters %in% c("17")] <- "Naive B"
Lym$cell[Lym$seurat_clusters %in% c("1","23")] <- "Activated B"
Lym$cell[Lym$seurat_clusters %in% c("21")] <- "GC B"
Lym$cell[Lym$seurat_clusters %in% c("20")] <- "Cycling B"
Lym$cell[Lym$seurat_clusters %in% c("9","22","14","24")] <- "Plasma"
Lym$cell[Lym$seurat_clusters %in% c("11")] <- "CD4+ Treg"
Lym$cell[Lym$seurat_clusters %in% c("8","15","7","6","2","5")] <- "CD4+ Naive"
Lym$cell[Lym$seurat_clusters %in% c("5","3")] <- "CD8+IL7R+ Memory"
Lym$cell[Lym$seurat_clusters %in% c("4","0","12")] <- "CD8+GZMK+ Memory"
Lym$cell[Lym$seurat_clusters %in% c("19")] <- "ILC"
Lym$cell[Lym$seurat_clusters %in% c("18")] <- "CD56 bright NK"
Lym$cell[Lym$seurat_clusters %in% c("16","13")] <- "CD56 dim NK"
Lym$cell[Lym$seurat_clusters %in% c("10")] <- "CD8+ CTL"
Idents(Lym) <- Lym$cell
DimPlot(Lym, label = T, label.size = 5, repel = T) 

## cell type labeling to merged tissue data
ind_lym <- match(rownames(Lym@meta.data), rownames(Tmerge@meta.data))
Tmerge@meta.data$cell[ind_lym] <- Lym@meta.data$cell

## dot plot for marker genes
Lym$cell <- factor(Lym$cell,levels = c("Activated B","Naive B","GC B","Cycling B","Plasma","CD4+ Treg","CD4+ Naive",
                                         "CD8+IL7R+ Memory","CD8+GZMK+ Memory","CD8+ CTL","CD56 dim NK","CD56 bright NK","ILC"))
marker_list_lym <- c("PTPRC","MS4A1","CD19","BANK1", # Activated B
                     "BLK","FCER2","TCL1A", # Naive B
                     "MME","DTX1", # GC B
                     "MKI67","TOP2A", # Cycling B
                     "MZB1","JCHAIN", # Plasma
                     "CD3E","CD4","FOXP3","CTLA4", # CD4+ Treg
                     "IL7R","SELL","CCR7", # CD4+ Naive
                     "CD8A","CCL5", # CD8+IL7R+ Memory
                     "GZMK","CD28", # CD8+GZMK+ Memory
                     "GZMB","PRF1", # CD8+ CTL
                     "NKG7","NCAM1","FCGR3A","FGFBP2","GZMH", # CD56 dim NK & CD56 bright NK
                     "KIT","PCDH9" # ILC
                    )
DotPlot(Lym, features = marker_list_lym, cols = c("lightgray", "slateblue4"), group.by = "cell") + 
  labs(x = " ", y = " ") +
  theme(axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 14),
        axis.text.y = ggtext::element_markdown(hjust = 1, size = 14, face = "italic"),
        legend.text = element_text(size = 14)) + 
  coord_flip() 

## cell population
Lym_diag <- Lym@meta.data %>% 
  group_by(diagnosis, cell) %>% 
  summarise(n = n()) %>% 
  ungroup()

ggplot(Lym_diag, aes(x = diagnosis, y = n, fill = cell)) +
  geom_bar(stat = "identity", position = "fill") + 
  ylab("Proportion") + xlab(" ") +
  theme_classic(base_size = 18) +  
  labs(fill = " ") +  
  theme(
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black")
  ) +
  scale_fill_manual(values = c("CD4+ Naive" = "#F8766D",
                               "CD4+ Treg" = "#E18A00",
                               "CD8+GZMK+ Memory" = "#BE9C00",
                               "CD56 bright NK"="#8CAB00",
                               "CD56 dim NK" = "#24B700",
                               "CD8+ CTL" = "#00BE70",
                               "CD8+IL7R+ Memory" ="#00C1AB",
                               "ILC" = "#00BBDA",
                               "Plasma" = "#00ACFC",
                               "Activated B" = "#8B93FF",
                               "Cycling B" ="#D575FE",
                               "Naive B"="#F962DD",
                               "GC B"="#FF65AC")
                   )

## tissue myeloid sub clustering
Mye <- subset(Tmerge, population %in% "Myeloid")
Mye <- NormalizeData(Mye)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2300)
Mye <- ScaleData(Mye, features = rownames(Mye))
Mye <- RunPCA(Mye, features = VariableFeatures(Mye))
Mye <- FindNeighbors(Mye, dims = 1:15)
Mye <- FindClusters(Mye, resolution = 1.5)
Mye <- RunUMAP(Mye, dims = 1:15)
DimPlot(Mye, label = T)
FeaturePlot(Mye, c("C1QA","FABP4","MKI67","SPP1", "TREM2", "GPNMB", "LPL",
                   "FCN1","FCGR3A","LYZ","S100A9","CLEC9A","CD1C","CLEC4C",
                   "MS4A2","CSF3R"), label = T)

## cell type annotation
Mye$cell[Mye$seurat_clusters %in% c("22","21","14")] <- "Cycling MΦ"
Mye$cell[Mye$seurat_clusters %in% c("1","12","24","25","10")] <- "Alveolar MΦ"
Mye$cell[Mye$seurat_clusters %in% c("8","7","4","5","0","6","3","16","20","27")] <- "SPP1+ MΦ"
Mye$cell[Mye$seurat_clusters %in% c("2","17")] <- "SPP1+ Mono"
Mye$cell[Mye$seurat_clusters %in% c("13")] <- "LAM"
Mye$cell[Mye$seurat_clusters %in% c("9","26")] <- "cMono"
Mye$cell[Mye$seurat_clusters %in% c("18")] <- "iMono"
Mye$cell[Mye$seurat_clusters %in% c("15","19")] <- "cDC2"
Mye$cell[Mye$seurat_clusters %in% c("23")] <- "pDC"
Mye$cell[Mye$seurat_clusters %in% c("28")] <- "cDC1"
Mye$cell[Mye$seurat_clusters %in% c("11","29")] <- "Mast cell"
Idents(Mye) <- Mye$cell
DimPlot(Mye, label = T, label.size = 5, repel = T) 

## cell type labeling to merged tissue data
ind_mye <- match(rownames(Mye@meta.data), rownames(Tmerge@meta.data))
Tmerge@meta.data$cell[ind_mye] <- Mye@meta.data$cell

## dot plot for marker genes
Mye$cell <- factor(Mye$cell, levels = c("Alveolar MΦ","Cycling MΦ","SPP1+ MΦ","SPP1+ Mono","LAM","cMono", "iMono", "cDC1", "cDC2", "pDC", "Mast cell"))
marker_list_mye <- c("C1QA","C1QB","FABP4","MARCO","PPARG", # Alveolar MΦ
                     "MKI67","TOP2A","PCNA", # Cycling MΦ
                     "SPP1","MERTK","LGALS3", # SPP1+ MΦ & SPP1+ Mono
                     "TREM2", "GPNMB", "LPL", "CD9", # LAM
                     "FCN1","CD14", # cMono
                     "FCGR3A","LYZ","S100A9", # iMono
                     "CLEC4C","IL3RA","TCF4", # cDC1
                     "CD1C","FCER1A","ITGAX", # cDC2
                     "CLEC9A","XCR1","BATF3", # pDC
                     "MS4A2","TPSAB1","KIT" # Mast cell
                    )
DotPlot(Mye, features = marker_list_mye, cols = c("lightgray", "slateblue4"), group.by = "cell") + 
  labs(x = " ", y = " ") +
  theme(axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 14),
        axis.text.y = ggtext::element_markdown(hjust = 1, size = 14, face = "italic"),
        legend.text = element_text(size = 14)
       ) + 
  coord_flip()

## cell population
Mye_diag <- Mye@meta.data %>% 
  group_by(diagnosis, cell) %>% 
  summarise(n = n()) %>% 
  ungroup()
ggplot(Mye_diag, aes(x = diagnosis, y = n, fill = cell)) + 
  geom_bar(stat = "identity", position = "fill") + 
  ylab("Proportion") + xlab(" ") +
  theme_classic(base_size = 18) + 
  theme(axis.title.x = element_text(colour = "black"), 
        axis.title.y = element_text(colour = "black")) +
  labs(fill = " ") + 
scale_fill_manual(values = c("SPP1+ MΦ" = "#F8766D",
                             "SPP1+ Mono" = "#DB8E00",
                             "cMono" = "#AEA200",
                             "cDC2"="#64B200",
                             "Cycling MΦ" = "#00BD5C",
                             "Alveolar MΦ" = "#00C1A7",
                             "pDC"="#00BADE",
                             "LAM" = "#00A6FF",
                             "Mast cell"="#B385FF",
                             "iMono"="#EF67EB",
                             "cDC1" = "#FF63B6")
                 )

## tissue stromal sub clustering
Stromal <- subset(Tmerge, population %in% "Stromal")
Stromal <- NormalizeData(Stromal)
Stromal <- FindVariableFeatures(Stromal, selection.method = "vst", nfeatures = 2300)
Stromal <- ScaleData(Stromal, features = rownames(Stromal))
Stromal <- RunPCA(Stromal, features = VariableFeatures(Stromal))
Stromal <- FindNeighbors(Stromal, dims = 1:15)
Stromal <- FindClusters(Stromal, resolution = 1.0)
Stromal <- RunUMAP(Stromal, dims = 1:15)
DimPlot(Stromal, label = T)
FeaturePlot(Stromal, c("COL1A1","WT1","PLIN2","HAS2","C3","SFRP1","MGST1","LEPR","BMP5","COL8A1","WNT2","FGFR4","MECOM","SFRP4",
                       "COL3A1","POSTN","MOXD1","PDGFRB","RGS5","HIGD1B","ACTG2","ACTA2","CNN1","MYH11","UPK3B","MSLN"))

## cell type annotation
Stromal$cell[Stromal$seurat_clusters %in% c("14")] <- "PLIN2+ fibroblast"
Stromal$cell[Stromal$seurat_clusters %in% c("11","1","18")] <- "Fibroblast"
Stromal$cell[Stromal$seurat_clusters %in% c("0","15")] <- "FGFR4+ myofibroblast"
Stromal$cell[Stromal$seurat_clusters %in% c("16","17","6","3","10")] <- "POSTN+ myofibroblast"
Stromal$cell[Stromal$seurat_clusters %in% c("5")] <- "SFRP4+ myofibroblast"
Stromal$cell[Stromal$seurat_clusters %in% c("7","2","13")] <- "Pericyte"
Stromal$cell[Stromal$seurat_clusters %in% c("9","8","4","12")] <- "SMC"
Stromal$cell[Stromal$seurat_clusters %in% c("19")] <- "Mesothelial cell"
Idents(Stromal) <- Stromal$cell
DimPlot(Stromal, label = T, label.size = 5, repel = T)

## cell type labeling to merged tissue data
ind_str <- match(rownames(Stromal@meta.data), rownames(Tmerge@meta.data))
Tmerge@meta.data$cell[ind_str] <- Stromal@meta.data$cell

## dot plot for marker genes
Stromal$cell <- factor(Stromal$cell,levels = c("PLIN2+ fibroblast","Fibroblast","FGFR4+ myofibroblast","SFRP4+ myofibroblast",
                                               "POSTN+ myofibroblast","Pericyte","SMC","Mesothelial cell"))
marker_list_str <- c("COL1A1","WT1","PLIN2","HAS2", # PLIN2+ fibroblast
                     "C3","SFRP1","MGST1","LEPR", # Fibroblast
                     "BMP5","COL8A1","WNT2","FGFR4", # FGFR4+ myofibroblast
                     "MECOM","SFRP4","COL3A1", # SFRP4+ myofibroblast
                     "POSTN","MOXD1", # POSTN+ myofibroblast
                     "PDGFRB","RGS5","HIGD1B", # Pericyte
                     "ACTG2","ACTA2","CNN1","MYH11", # SMC
                     "UPK3B","MSLN" # Mesothelial cell
                    )
DotPlot(Stromal, features = marker_list_str, cols = c("lightgray", "slateblue4"), group.by = "cell") + 
  labs(x = " ", y = " ") +
  theme(axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 14),
        axis.text.y = ggtext::element_markdown(hjust = 1, size = 14, face = "italic"),
        legend.text = element_text(size = 14)) + 
  coord_flip()

## cell population
Str_diag <- Stromal@meta.data %>% 
  group_by(diagnosis, cell) %>% 
  summarise(n = n()) %>% ungroup()
ggplot(Str_diag, aes(x = diagnosis, y = n, fill = cell)) + 
  geom_bar(stat = "identity", position = "fill") + 
  ylab("Proportion") + xlab(" ") +
  theme_classic(base_size = 18) + 
  theme(axis.title.x = element_text(colour = "black"), 
        axis.title.y = element_text(colour = "black")
       ) +
  labs(fill = " ") + 
  scale_fill_manual(values = c("PLIN2+ fibroblast" = "#F8766D",
                               "Fibroblast" = "#CD9600",
                               "FGFR4+ myofibroblast" = "#7CAE00",
                               "SFRP4+ myofibroblast" = "#00BE67",
                               "POSTN+ myofibroblast" = "#00BFC4",
                               "Pericyte" = "#00A9FF",
                               "SMC" = "#C77CFF",
                               "Mesothelial cell" = "#FF61CC")
                   )

##  tissue epithelial sub clustering
Epithelial <- subset(Tmerge, population %in% "Epithelial")
Epithelial <- NormalizeData(Epithelial)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2300)
Epithelial <- ScaleData(Epithelial, features = rownames(Epithelial))
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(Epithelial))
Epithelial <- FindNeighbors(Epithelial, dims = 1:15)
Epithelial <- FindClusters(Epithelial, resolution = 1.5)
Epithelial <- RunUMAP(Epithelial, dims = 1:15)
DimPlot(Epithelial, label = T)
FeaturePlot(Epithelial, c("EPCAM","SFTPC","NAPSA","SFTPA1","AGER","CAV1","PDPN","SCGB3A2","SFTPB","SCGB3A1","SCGB1A1","TFF3","MUC5B",
                          "BPIFB1","MUC5AC","KRT5","KRT15","KRT17","CDK1","TOP2A","MKI67","FOXJ1","PIFO","TPPP3"))

## cell type annotation
Epithelial$cell[Epithelial$seurat_clusters %in% c("6")] <- "AT2"
Epithelial$cell[Epithelial$seurat_clusters %in% c("21","10")] <- "AT1"
Epithelial$cell[Epithelial$seurat_clusters %in% c("7")] <- "Transitional AT2"
Epithelial$cell[Epithelial$seurat_clusters %in% c("8","2","13")] <- "SCGB3A2+ Secretory"
Epithelial$cell[Epithelial$seurat_clusters %in% c("9","11","12","16")] <- "Club"
Epithelial$cell[Epithelial$seurat_clusters %in% c("1","0")] <- "Basal"
Epithelial$cell[Epithelial$seurat_clusters %in% c("18")] <- "KRT5- Basal"
Epithelial$cell[Epithelial$seurat_clusters %in% c("5","14")] <- "Goblet"
Epithelial$cell[Epithelial$seurat_clusters %in% c("20")] <- "MUC5AC+ Goblet"
Epithelial$cell[Epithelial$seurat_clusters %in% c("17","3","15","4")] <- "Ciliated"
Epithelial$cell[Epithelial$seurat_clusters %in% c("19")] <- "Cycling Epithelial"
Idents(Epithelial) <- Epithelial$cell
DimPlot(Epithelial, label = T, label.size = 5, repel = T)

## cell type labeling to merged tissue data
ind_epi <- match(rownames(Epithelial@meta.data), rownames(Tmerge@meta.data))
Tmerge@meta.data$cell[ind_epi] <- Epithelial@meta.data$cell

## dot plot for marker genes
Epithelial$cell <- factor(Epithelial$cell,levels = c("AT2","Transitional AT2","AT1","SCGB3A2+ Secretory","Club","Goblet",
                                                     "MUC5AC+ Goblet","Basal","KRT5- Basal","Cycling Epithelial","Ciliated"))
marker_gene_epi <- c("EPCAM","SFTPC","NAPSA","SFTPA1", # AT2
                     "AGER","CAV1","PDPN", # AT1
                     "SCGB3A2","SFTPB", # SCGB3A2+ Secretory
                     "SCGB3A1","SCGB1A1", # Club
                     "TFF3","MUC5B","BPIFB1", # Goblet
                     "MUC5AC", # MUC5AC+ Goblet
                     "KRT5","KRT15","KRT17", # Basal
                     "CDK1","TOP2A","MKI67", # Cycling Epithelial
                     "FOXJ1","PIFO","TPPP3" # Ciliated
                    )
DotPlot(Epithelial, features = marker_gene_epi, cols = c("lightgray", "slateblue4"), group.by = "cell") + 
  labs(x = " ", y = " ") +
  theme(axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 12),
        axis.text.y = ggtext::element_markdown(hjust = 1, size = 14, face = "italic"),
        legend.text = element_text(size = 14)) + 
  coord_flip()

## cell population
Epi_diag <- Epithelial@meta.data %>% 
  group_by(diagnosis, cell) %>% 
  summarise(n = n()) %>% 
  ungroup()
ggplot(Epi_diag, aes(x = diagnosis, y = n, fill = cell)) + 
  geom_bar(stat = "identity", position = "fill") + 
  ylab("Proportion") + 
  xlab(" ") +
  theme_classic(base_size = 18) +  
  theme(axis.title.x = element_text(colour = "black"), 
        axis.title.y = element_text(colour = "black")
       ) + 
  labs(fill = " ") +
  scale_fill_manual(values = c("SCGB3A2+ Secretory" = "#F8766D",
                               "Ciliated" = "#DB8E00",
                               "Transitional AT2" = "#AEA200",
                               "Goblet" = "#64B200",
                               "Club" = "#00BD5C",
                               "Basal" = "#00C1A7",
                               "MUC5AC+ Goblet" = "#00BADE",
                               "AT2" = "#00A6FF",
                               "AT1" = "#B385FF",
                               "Cycling Epithelial" = "#EF67EB",
                               "KRT5- Basal" = "#FF63B6")
                   )

##  tissue endothelial sub clustering
Endothelial <- subset(Tmerge, population %in% "Endothelial")
Endothelial <- NormalizeData(Endothelial)
Endothelial <- FindVariableFeatures(Endothelial, selection.method = "vst", nfeatures = 2300)
Endothelial <- ScaleData(Endothelial, features = rownames(Endothelial))
Endothelial <- RunPCA(Endothelial, features = VariableFeatures(Endothelial))
Endothelial <- FindNeighbors(Endothelial, dims = 1:15)
Endothelial <- FindClusters(Endothelial, resolution = 0.4)
Endothelial <- RunUMAP(Endothelial, dims = 1:15)
DimPlot(Endothelial, label = T)
FeaturePlot(Endothelial,  c("HEY1","GJA4","GJA5","LYVE1","PROX1","GYPC","ACKR1","CCL14","SELE"), label = T)

## cell type annotation
Endothelial$cell[Endothelial$seurat_clusters %in% c("7")] <- "Lymphatic"
Endothelial$cell[Endothelial$seurat_clusters %in% c("5","6")] <- "Artery"
Endothelial$cell[Endothelial$seurat_clusters %in% c("0","1","2","3","4","8")] <- "Vein"
Idents(Endothelial) <- Endothelial$cell
DimPlot(Endothelial, label = T, label.size = 5, repel = T)

## cell type labeling to merged tissue data
ind_end <- match(rownames(Endothelial@meta.data), rownames(Tmerge@meta.data))
Tmerge@meta.data$cell[ind_end] <- Endothelial@meta.data$cell

## dot plot for marker genes
Endothelial$cell <- factor(Endothelial$cell,levels = c("Artery","Lymphatic","Vein"))
marker_gene_end <- c("HEY1","GJA4","GJA5", # Artery
                     "LYVE1","PROX1","GYPC", # Lymphatic
                     "ACKR1","CCL14","SELE" # Vein
                    )
DotPlot(Endothelial, features = marker_gene_end, cols = c("lightgray", "slateblue4"), group.by = "cell") + 
  labs(x = " ", y = " ") +
  theme(axis.text.x = ggtext::element_markdown(angle= 45,hjust = 1, size = 14),
        axis.text.y = ggtext::element_markdown(hjust = 1, size = 14, face = "italic"),
        legend.text = element_text(size = 14)) + 
  coord_flip()

## cell population
End_diag <- Endothelial@meta.data %>%  
  group_by(diagnosis, cell) %>% 
  summarise(n = n()) %>%  
  ungroup()
ggplot(End_diag, aes(x = diagnosis, y = n, fill = cell)) + 
  geom_bar(stat = "identity", position = "fill") + 
  ylab("Proportion") + 
  xlab(" ") +
  theme_classic(base_size = 18) +  
  theme(axis.title.x = element_text(colour = "black"), 
        axis.title.y = element_text(colour = "black")) + 
  labs(fill = " ") +
  scale_fill_manual(values = c("Vein" = "#F8766D",
                               "Artery" = "#00BA38",
                               "Lymphatic" = "#619CFF")
                   )

## merged tissue cell type
DimPlot(Tmerge, label = T, group.by = "cell", label.size = 5, repel = T) + ggtitle("")

