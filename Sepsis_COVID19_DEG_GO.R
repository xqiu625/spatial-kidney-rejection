.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )

library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(scater)
library(MAST)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(rrvgo)

#---------------------------------Severe platelets Sepsis VS. COVID-19------------------

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19")
severe <- readRDS("Severe_platelet_sepsis_covid19.rds")
DefaultAssay(severe) <- "SCT"

### Sepsis vs COVID-19
celltype = 'Platelet'
Idents(severe) = 'Disease'
severe@misc$markers <- FindAllMarkers(object = severe, assay = 'SCT',only.pos = TRUE, test.use = 'MAST')
setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/DEG")
write.table(severe@misc$markers,file=paste0(celltype,'_Sepsis_vs_COVID19_Severe.txt'),row.names = FALSE,quote = FALSE,sep = '\t')

degs <- severe@misc$markers
groups = unique(degs$cluster)
for(group_ in groups){
  degs_ = degs %>% filter(.,cluster==group_ & p_val_adj<0.05)
  gene_ = degs_$gene
  B_fun.gene.1 <- bitr(gene_, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
  B_fun.gene = B_fun.gene.1$ENTREZID
  B_fun_bp = enrichGO(gene = B_fun.gene,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable= TRUE)
  if(!is.null(B_fun_bp)){
    B_fun_bp <- setReadable(B_fun_bp, 'org.Hs.eg.db', 'ENTREZID')
    setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/GO_analysis")
    write.table(B_fun_bp,file = paste(celltype, "_",group_,'_go-bp_severe.txt',sep = ''),sep='\t',quote = FALSE, row.names = F)
    
    dpi = 300
    png(file = paste(celltype, "_",group_,'_go-bp_severe.png',sep = ''),width = dpi * 12,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
    print(barplot(B_fun_bp,showCategory=20,drop=T)+theme(axis.text.y = element_text(size = 16),legend.text = element_text(size = 16),legend.title = element_text(size = 16)))
    dev.off()
  }
}

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/GO_analysis")
go_analysis <- read.delim("Platelet_Sepsis_go-bp_severe.txt")
simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Figures")
dpi = 300
png(file = "test.png", width = dpi * 9,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
p<- treemapPlot(reducedTerms)
p
dev.off()


#---------------------------------Convalescence platelets Sepsis VS. COVID-19------------------
setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19")
Convalescence <- readRDS("Convalescence_platelet_sepsis_covid19.rds")
DefaultAssay(Convalescence) <- "SCT"

### Sepsis vs COVID-19
celltype = 'Platelet'
Idents(Convalescence) = 'Disease'
Convalescence@misc$markers <- FindAllMarkers(object = Convalescence, assay = 'SCT',only.pos = TRUE, test.use = 'MAST')
setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/DEG")
write.table(Convalescence@misc$markers,file=paste0(celltype,'_Sepsis_vs_COVID19_Convalescence.txt'),row.names = FALSE,quote = FALSE,sep = '\t')

degs <- Convalescence@misc$markers
groups = unique(degs$cluster)
for(group_ in groups){
  degs_ = degs %>% filter(.,cluster==group_ & p_val_adj<0.05)
  gene_ = degs_$gene
  B_fun.gene.1 <- bitr(gene_, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
  B_fun.gene = B_fun.gene.1$ENTREZID
  B_fun_bp = enrichGO(gene = B_fun.gene,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable= TRUE)
  if(!is.null(B_fun_bp)){
    B_fun_bp <- setReadable(B_fun_bp, 'org.Hs.eg.db', 'ENTREZID')
    setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/GO_analysis")
    write.table(B_fun_bp,file = paste(celltype, "_",group_,'_go-bp_Convalescence.txt',sep = ''),sep='\t',quote = FALSE, row.names = F)
    
    dpi = 300
    png(file = paste(celltype, "_",group_,'_go-bp_Convalescence.png',sep = ''),width = dpi * 12,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
    print(barplot(B_fun_bp,showCategory=20,drop=T)+theme(axis.text.y = element_text(size = 16),legend.text = element_text(size = 16),legend.title = element_text(size = 16)))
    dev.off()
  }
}

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/Bernardes_2020_COVID19/GO_analysis")
go_analysis <- read.delim("Platelet_Sepsis_go-bp_Convalescence.txt")
simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

setwd("/bigdata/godziklab/shared/Xinru/072920_sepsis/seurat/Figures")
dpi = 300
png(file = "test.png", width = dpi * 9,height = dpi * 9,units = "px",res = dpi,type = 'cairo')
p<- treemapPlot(reducedTerms)
p
dev.off()





  
