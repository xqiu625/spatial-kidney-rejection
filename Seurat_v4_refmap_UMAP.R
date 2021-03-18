.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

### UMAP already ref map data
library(Seurat)
library(stringr)
integrated <- readRDS("SeuratV4_mapped_0729.rds")
DefaultAssay(object = integrated) <- "RNA"

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Figures")
dpi = 300
png(file = "test.png", width = dpi * 18,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
p1 = DimPlot(integrated, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(integrated, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p <- p1 + p2
p
dev.off()

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Figures")
dpi = 300
png(file = "test.png", width = dpi * 18,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
p2 = DimPlot(integrated, 
             reduction = "ref.umap", 
             group.by = "predicted.celltype.l2", 
             split.by = "status",
             label = TRUE, 
             label.size = 3 ,
             repel = TRUE) + NoLegend()
p2
dev.off()

df <- integrated@assays[["RNA"]]
df[1:20, 1:20]


### Check prediction score for the Seraut ref map
library(dplyr)
score <- integrated@meta.data
score1 <- score %>% filter(predicted.celltype.l2.score > 0.75)
dim(score1)
score2 <- score %>% filter(predicted.celltype.l2.score < 0.75)
table(score2$predicted.celltype.l2)


### Save subset data for further analysis
Platelet <- subset(integrated, predicted.celltype.l2 == "Platelet")
Platelet
# An object of class Seurat
# 55672 features across 5446 samples within 5 assays
# Active assay: RNA (33538 features, 0 variable features)
# 4 other assays present: SCT, prediction.score.celltype.l1, prediction.score.celltype.l2, predicted_ADT
# 2 dimensional reductions calculated: ref.spca, ref.umap
saveRDS(Platelet, "/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Mapped_Platelets.rds")

Eryth <- subset(integrated, predicted.celltype.l2 == "Eryth")
Eryth
# An object of class Seurat
# 55672 features across 2312 samples within 5 assays
# Active assay: RNA (33538 features, 0 variable features)
# 4 other assays present: SCT, prediction.score.celltype.l1, prediction.score.celltype.l2, predicted_ADT
# 2 dimensional reductions calculated: ref.spca, ref.umap
saveRDS(Eryth, "/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Mapped_Eryth.rds")
