## ------------------------------------------------------------
## library load
## ------------------------------------------------------------
library(pheatmap)
library(viridis)

## immune cell subset
Tmerge$cor <- Tmerge$cell
Tmerge$cor[Tmerge$cell %in% c("CD8+GZMK+ Memory","CD8+IL7R+ Memory")] <- "CD8+ Memory"
Timmune <- subset(Tmerge, population %in% c("Myeloid","Lymphoid"))

## aaverage expression (UIP)
UIP_T <- subset(Timmune, diagnosis %in% c("UIP"))
UIP_P <- subset(Pmerge, diagnosis %in% c("UIP"))
UIP_B <- subset(Bmerge, diagnosis %in% c("UIP"))
avg_UIP_T <- AggregateExpression(UIP_T, assays = "RNA", slot = "data", group.by = "cell")$RNA
avg_UIP_P   <- AggregateExpression(UIP_P, assays = "RNA", slot = "data", group.by = "cell")$RNA
avg_UIP_B   <- AggregateExpression(UIP_B, assays = "RNA", slot = "data", group.by = "cell")$RNA

## shared genes across compartments
common_genes_UIP <- Reduce(intersect, list(rownames(avg_UIP_T),rownames(avg_UIP_P),rownames(avg_UIP_B)))
avg_UIP_T <- avg_UIP_T[common_genes_UIP, , drop = FALSE]
avg_UIP_P   <- avg_UIP_P[common_genes_UIP, , drop = FALSE]
avg_UIP_B   <- avg_UIP_B[common_genes_UIP, , drop = FALSE]
common_avg_UIP <- cbind(Tissue = avg_UIP_T[,1], PBMC = avg_UIP_P[,1], BALF = avg_UIP_B[,1])

## compartment-wise correlation (UIP)
cor_mat_UIP <- cor(common_avg_UIP, method = "spearman")
pheatmap(cor_mat_UIP, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_format = "%.2f", 
         number_color = "black",  
         fontsize = 15,
         color  = colorRampPalette(c("blue3", "white", "tomato3"))(length(seq(0, 1, by = 0.01)) - 1),
         breaks = seq(0, 1, by = 0.01),
         main = " "
        )

## average expression (DIP)
DIP_T <- subset(Timmune, diagnosis %in% c("DIP"))
DIP_P <- subset(Pmerge, diagnosis %in% c("DIP"))
DIP_B <- subset(Bmerge, diagnosis %in% c("DIP"))
avg_DIP_T <- AggregateExpression(DIP_T, assays = "RNA", slot = "data", group.by = "cell")$RNA
avg_DIP_P   <- AggregateExpression(DIP_P, assays = "RNA", slot = "data", group.by = "cell")$RNA
avg_DIP_B   <- AggregateExpression(DIP_B, assays = "RNA", slot = "data", group.by = "cell")$RNA

## shared genes across compartments
common_genes_DIP <- Reduce(intersect, list(rownames(avg_DIP_T),rownames(avg_DIP_P),rownames(avg_DIP_B)))
avg_DIP_T <- avg_DIP_T[common_genes_DIP, , drop = FALSE]
avg_DIP_P   <- avg_DIP_P[common_genes_DIP, , drop = FALSE]
avg_DIP_B   <- avg_DIP_B[common_genes_DIP, , drop = FALSE]
common_avg_DIP <- cbind(Tissue = avg_DIP_T[,1], PBMC = avg_DIP_P[,1], BALF = avg_DIP_B[,1])

## compartment-wise correlation (DIP)
cor_mat_DIP <- cor(common_avg_DIP, method = "spearman")
pheatmap(cor_mat_DIP, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE, 
         number_format = "%.2f", 
         number_color = "black",  
         fontsize = 15,
         color  = colorRampPalette(c("blue3", "white", "tomato3"))(length(seq(0, 1, by = 0.01)) - 1),
         breaks = seq(0, 1, by = 0.01),
         main = " "
        )

## overlapping cell types correlation (UIP)
common_ct_UIP <- Reduce(intersect, list(UIP_T$cor, UIP_P$cell, UIP_B$cell))
avg_UIP_T_cell <- AggregateExpression(UIP_T, assays = "RNA", slot = "data", group.by = "cor")$RNA
avg_UIP_P_cell   <- AggregateExpression(UIP_P, assays = "RNA", slot = "data", group.by = "cell")$RNA
avg_UIP_B_cell   <- AggregateExpression(UIP_B, assays = "RNA", slot = "data", group.by = "cell")$RNA
common_genes_UIP <- Reduce(intersect, list(rownames(avg_UIP_T),rownames(avg_UIP_P),rownames(avg_UIP_B)))
avg_UIP_T_cell <- avg_UIP_T_cell[common_genes_UIP, common_ct_UIP]
avg_UIP_P_cell <- avg_UIP_P_cell[common_genes_UIP, common_ct_UIP]
avg_UIP_B_cell <- avg_UIP_B_cell[common_genes_UIP, common_ct_UIP]
cor_mat_UIP <- sapply(common_ct_UIP, function(ct) {
  c(Tissue_vs_PBMC = cor(avg_UIP_T_cell[, ct], avg_UIP_P_cell[, ct], method = "spearman"),
    Tissue_vs_BALF = cor(avg_UIP_T_cell[, ct], avg_UIP_B_cell[, ct], method = "spearman"),
    PBMC_vs_BALF   = cor(avg_UIP_P_cell[, ct], avg_UIP_B_cell[, ct], method = "spearman"))
})
cor_mat_UIP <- t(cor_mat_UIP)
cor_mat_UIP <- cor_mat_UIP[order(rownames(cor_mat_UIP)), , drop = FALSE]
pheatmap(cor_mat_UIP,
         cluster_rows = F,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 15,
         color = colorRampPalette(c("blue3", "white", "tomato3"))(100),breaks = seq(0, 1, length.out = 101),
         main = " "
        )

## overlapping cell types correlation (DIP)
common_ct_DIP <- Reduce(intersect, list(DIP_T$cor, DIP_P$cell, DIP_B$cell))
avg_DIP_T_cell <- AggregateExpression(DIP_T, assays = "RNA", slot = "data", group.by = "cor")$RNA
avg_DIP_P_cell   <- AggregateExpression(DIP_P, assays = "RNA", slot = "data", group.by = "cell")$RNA
avg_DIP_B_cell   <- AggregateExpression(DIP_B, assays = "RNA", slot = "data", group.by = "cell")$RNA
common_genes_DIP <- Reduce(intersect, list(rownames(avg_DIP_T),rownames(avg_DIP_P),rownames(avg_DIP_B)))
avg_DIP_T_cell <- avg_DIP_T_cell[common_genes_DIP, common_ct_DIP]
avg_DIP_P_cell <- avg_DIP_P_cell[common_genes_DIP, common_ct_DIP]
avg_DIP_B_cell <- avg_DIP_B_cell[common_genes_DIP, common_ct_DIP]
cor_mat_DIP <- sapply(common_ct_DIP, function(ct) {
  c(Tissue_vs_PBMC = cor(avg_DIP_T_cell[, ct], avg_DIP_P_cell[, ct], method = "spearman"),
    Tissue_vs_BALF = cor(avg_DIP_T_cell[, ct], avg_DIP_B_cell[, ct], method = "spearman"),
    PBMC_vs_BALF   = cor(avg_DIP_P_cell[, ct], avg_DIP_B_cell[, ct], method = "spearman"))
  })
cor_mat_DIP <- t(cor_mat_DIP)
cor_mat_DIP <- cor_mat_DIP[order(rownames(cor_mat_DIP)), , drop = FALSE]
pheatmap(cor_mat_DIP,
         cluster_rows = F,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize = 15,
         color = colorRampPalette(c("blue3", "white", "tomato3"))(100),breaks = seq(0, 1, length.out = 101),
         main = " "
        )
