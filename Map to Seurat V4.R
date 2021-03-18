### need to conda Seurat v4 envir
# conda activate r4.0.3

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(stringr)

reference <- LoadH5Seurat("/bigdata/godziklab/shared/Xinru/v2_013020_sepsis/seurat_v4/pbmc_multimodal.h5seurat")
MK <- readRDS("Megakaryocytes_final.rds")
DefaultAssay(object = MK) <- "RNA"
df <- MK[["RNA"]]@counts
df[1:20, 1:20]

### Use counts for Seurat object
MK2 <- CreateSeuratObject(counts = df, project = "MKcounts", min.cells = 3, min.features = 200)
MK2
df2 <- MK2@assays[["RNA"]]
df2[1:20, 1:20]

### Map to Seurat ref
df <- SCTransform(MK2, verbose = FALSE)
anchors <- FindTransferAnchors(
  reference = reference,
  query = df,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
df <- MapQuery(
  anchorset = anchors,
  query = df,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Figures")
dpi = 300
png(file = "test.png", width = dpi * 18,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
p1 = DimPlot(df, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(df, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p <- p1 + p2
p
dev.off()


# Cells with high-confidence annotations (for example, prediction scores > 0.75) reflect predictions that are supported by mulitple consistent anchors.

