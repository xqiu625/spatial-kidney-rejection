### https://satijalab.org/seurat/articles/integration_rpca.html

.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(glmGamPoi)

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/cellType_data")
platetlet_sepsis <-  readRDS("Platelet.rds")

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19")
mk_covid <-  readRDS("Megakaryocytes_final.rds")

head(platetlet_sepsis@meta.data)
head(mk_covid@meta.data)

unique(platetlet_sepsis@meta.data$status)
unique(mk_covid@meta.data$Pseudotime_name)

platetlet_sepsis_ns <- subset(platetlet_sepsis, status == "Sepsis_NS")
mk_covid_severe <- subset(mk_covid, Pseudotime_name %in% c("Critical", "Complicated"))

df <- platetlet_sepsis_ns[["RNA"]]@counts
platetlet_sepsis_ns <- CreateSeuratObject(counts = df, project = "platetlet_sepsis_counts", min.cells = 3, min.features = 200)
df <- mk_covid_severe[["RNA"]]@counts
mk_covid_severe <- CreateSeuratObject(counts = df, project = "mk_covid_counts", min.cells = 3, min.features = 200)

platetlet_sepsis_ns@meta.data$Disease<- c("Sepsis")
mk_covid_severe@meta.data$Disease <- c("COVID-19")

### Integrate sepsis and covid samples
merged_platelet_severe <- merge(x = platetlet_sepsis_ns, y = mk_covid_severe)
merged_platelet_severe.list <- SplitObject(merged_platelet_severe, split.by = "Disease")
merged_platelet_severe.list <- lapply(X = merged_platelet_severe.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = merged_platelet_severe.list, nfeatures = 3000)
merged_platelet_severe.list <- PrepSCTIntegration(object.list = merged_platelet_severe.list, anchor.features = features)
merged_platelet_severe.list <- lapply(X = merged_platelet_severe.list, FUN = RunPCA, features = features)
immune.anchors <- FindIntegrationAnchors(object.list = merged_platelet_severe.list, normalization.method = "SCT", 
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)

# Save combined
saveRDS(immune.combined.sct, "/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/Severe_platelet_sepsis_covid19.rds")

# Visualization
setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Figures")
dpi = 300
png(file = "test.png", width = dpi * 9,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Disease")
p1
dev.off()

#---------------------------------Non-severe platelets Sepsis VS. COVID-19------------------
# Performing integration on datasets normalized with SCTransform

platetlet_sepsis_s <- subset(platetlet_sepsis, status == "Sepsis_S")
mk_covid_non_severe <- subset(mk_covid, Pseudotime_name %in% c("Moderate/early convalescence", "Late convalescence"))

df <- platetlet_sepsis_s[["RNA"]]@counts
platetlet_sepsis_s <- CreateSeuratObject(counts = df, project = "platetlet_sepsis_counts", min.cells = 3, min.features = 200)
df <- mk_covid_non_severe[["RNA"]]@counts
mk_covid_non_severe <- CreateSeuratObject(counts = df, project = "mk_covid_counts", min.cells = 3, min.features = 200)

platetlet_sepsis_s@meta.data$Disease<- c("Sepsis")
mk_covid_non_severe@meta.data$Disease <- c("COVID-19")

### Integrate sepsis and covid samples
merged_platelet_nonsevere <- merge(x = platetlet_sepsis_s, y = mk_covid_non_severe)
merged_platelet_nonsevere.list <- SplitObject(merged_platelet_nonsevere, split.by = "Disease")
merged_platelet_nonsevere.list <- lapply(X = merged_platelet_nonsevere.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = merged_platelet_nonsevere.list, nfeatures = 3000)
merged_platelet_nonsevere.list <- PrepSCTIntegration(object.list = merged_platelet_nonsevere.list, anchor.features = features)
merged_platelet_nonsevere.list <- lapply(X = merged_platelet_nonsevere.list, FUN = RunPCA, features = features)
immune.anchors <- FindIntegrationAnchors(object.list = merged_platelet_nonsevere.list, normalization.method = "SCT", 
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)

# Save combined
saveRDS(immune.combined.sct, "/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/NonSevere_platelet_sepsis_covid19.rds")

# Visualization
setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Figures")
dpi = 300
png(file = "test.png", width = dpi * 9,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Disease")
p1
dev.off()
