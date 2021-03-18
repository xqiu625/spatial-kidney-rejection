library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

for (file in c("HC1_filtered_feature_bc_matrix", 
               "HC2_filtered_feature_bc_matrix", 
               "NS1_T0_filtered_feature_bc_matrix", 
               "NS1_T6_filtered_feature_bc_matrix",
               "NS2_T0_filtered_feature_bc_matrix",
               "NS2_T6_filtered_feature_bc_matrix",
               "S1_T0_filtered_feature_bc_matrix",
               "S1_T6_filtered_feature_bc_matrix",
               "S2_T0_filtered_feature_bc_matrix",
               "S2_T6_filtered_feature_bc_matrix",
               "S3_T0_filtered_feature_bc_matrix",
               "S3_T6_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("~/Documents/Projects/scRNA/sepsis_2nd/data/filtered matrix/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

merged_seurat <- merge(x = HC1_filtered_feature_bc_matrix, 
                       y = c(HC2_filtered_feature_bc_matrix, 
                             NS1_T0_filtered_feature_bc_matrix, 
                             NS1_T6_filtered_feature_bc_matrix, 
                             NS2_T0_filtered_feature_bc_matrix,
                             NS2_T6_filtered_feature_bc_matrix,
                             S1_T0_filtered_feature_bc_matrix,
                             S1_T6_filtered_feature_bc_matrix,
                             S2_T0_filtered_feature_bc_matrix,
                             S2_T6_filtered_feature_bc_matrix,
                             S3_T0_filtered_feature_bc_matrix,
                             S3_T6_filtered_feature_bc_matrix
                             ), 
                       add.cell.id = c("HC1", 
                                       "HC2", 
                                       "NS1_T0", 
                                       "NS1_T6",
                                       "NS2_T0",
                                       "NS2_T6",
                                       "S1_T0",
                                       "S1_T6",
                                       "S2_T0",
                                       "S2_T6",
                                       "S3_T0",
                                       "S3_T6"))
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
View(merged_seurat@meta.data)


# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^HC1"))] <- "HC1"
metadata$sample[which(str_detect(metadata$cells, "^HC2"))] <- "HC2"
metadata$sample[which(str_detect(metadata$cells, "^NS1_T0"))] <- "NS1_T0"
metadata$sample[which(str_detect(metadata$cells, "^NS1_T6"))] <- "NS1_T6"
metadata$sample[which(str_detect(metadata$cells, "^NS2_T0"))] <- "NS2_T0"
metadata$sample[which(str_detect(metadata$cells, "^NS2_T6"))] <- "NS2_T6"
metadata$sample[which(str_detect(metadata$cells, "^S1_T0"))] <- "S1_T0"
metadata$sample[which(str_detect(metadata$cells, "^S1_T6"))] <- "S1_T6"
metadata$sample[which(str_detect(metadata$cells, "^S2_T0"))] <- "S2_T0"
metadata$sample[which(str_detect(metadata$cells, "^S2_T6"))] <- "S2_T6"
metadata$sample[which(str_detect(metadata$cells, "^S3_T0"))] <- "S3_T0"
metadata$sample[which(str_detect(metadata$cells, "^S3_T6"))] <- "S3_T6"

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
# Create .RData object to load at any time
save(merged_seurat, file="~/Documents/Projects/scRNA/sepsis_2nd/data/raw_merged_sepsis.RData")
# Visualize the number of cell counts per sample
p1 <- metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave(p1,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/cell counts per sample.png", 
       width = 4, height = 4)

# Visualize the number UMIs/transcripts per cell
p2 <- metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave(p2,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/UMIs_transcripts per cell.png", 
       width = 4, height = 4)

# Genes detected per cell
p3 <- metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
ggsave(p3,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/Genes detected per cell.png", 
       width = 4, height = 4)


# Visualize the distribution of genes detected per cell via boxplot
p4 <- metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave(p4,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/Genes detected per cell boxplot.png", 
       width = 4, height = 4)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
p5 <- metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
ggsave(p5,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/UMIs vs genes detected.png", 
       width = 4, height = 4)

# Visualize the distribution of mitochondrial gene expression detected per cell
p6 <-metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave(p6,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/Mitochondrial counts ratio.png", 
       width = 4, height = 4)

# Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI
p7 <-metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave(p7,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/novelty.png", 
       width = 4, height = 4)


# Filtering

# nGene >= 200
# nGene <= 6000
# mitoRatio < 0.20
# nUMI > 1000

# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nGene >= 200) & 
                            (nGene <= 6000) &
                            (mitoRatio < 0.20) & 
                            (nUMI > 1000))
filtered_seurat
table(filtered_seurat$sample)

# Create .RData object to load at any time
save(filtered_seurat, file="~/Documents/Projects/scRNA/sepsis_2nd/data/filtered_merged_sepsis.RData")


##### Below run in HPCC


load("~/Documents/Projects/scRNA/sepsis_2nd/data/filtered_merged_sepsis.RData")
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
# Load cell cycle markers
load("~/Documents/Projects/scRNA/sepsis_2nd/data/cycle.rda")
# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data) 

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

# SCTransform
options(future.globals.maxSize = 4000 * 1024^2)
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
split_seurat <- split_seurat[c("HC1", 
                               "HC2", 
                               "NS1_T0", 
                               "NS1_T6",
                               "NS2_T0",
                               "NS2_T6",
                               "S1_T0",
                               "S1_T6",
                               "S2_T0",
                               "S2_T6",
                               "S3_T0",
                               "S3_T6")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}


# Check which assays are stored in objects
split_seurat$HC1@assays

# Integrate samples using shared highly variable genes
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
# Save integrated seurat object
saveRDS(seurat_integrated, "~/Documents/Projects/scRNA/sepsis_2nd/data/integrated_sepsis_0729.RData")


#UMAP visualization
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
# Plot PCA
p8 <- PCAPlot(seurat_integrated,
        split.by = "sample")  
ggsave(p8,file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/pca_sample.png", 
       width = 4, height = 4)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:50)
p10 <- DimPlot(seurat_integrated,
        split.by = "sample",
        reduction = "umap", 
        label.size = 4, 
        pt.size = 0.5,
        plot.title = "UMAP")

ggsave(p10, file = "~/Documents/Projects/scRNA/sepsis_2nd/out/figures/QC/umap_sample.png", 
       width = 4, height = 4)

