library(CelliD)
library(Seurat)
library(tidyverse)
library(ggpubr)

library(msigdbr)
library(scater)
library(SeuratDisk)

routine <- function(Seuratz_dir,Nature_dir){
  
  nam2 <- Nature_dir
  data_dir <- Seuratz_dir
  list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
  data <- Read10X(data.dir = data_dir)
  seurat_object = CreateSeuratObject(counts = data)
  counts <- seurat_object@assays$RNA@counts
  counts <- counts[rownames(counts) %in% HgProteinCodingGenes,]
  Seurat_temp <- CreateSeuratObject(counts = counts ,min.cells = 5)
  mett1 <- readRDS(nam2)
  Seurat_temp<-AddMetaData(Seurat_temp,mett1)
  Seurat_temp <- ScaleData(Seurat_temp, features = rownames(Seurat_temp))
  Seurat_temp <- RunMCA(Seurat_temp)
  Seurat_temp <- FindVariableFeatures(Seurat_temp)
  Seurat_temp <- RunPCA(Seurat_temp)
  Seurat_temp <- RunUMAP(Seurat_temp,dims = 1:30)
  
  Seurat_temp
}



data("HgProteinCodingGenes")

signnature <-read.csv("glioSign.csv")
####################################################  MARKER LIST ################################################
pathwaysH <- list()
pathwaysH$MES2 <- signnature$MES2
pathwaysH$MES1 <- signnature$MES1
pathwaysH$AC <- signnature$AC
pathwaysH$OPC <- signnature$OPC
pathwaysH$NPC1 <- signnature$NPC1
pathwaysH$NPC2 <- signnature$NPC2
pathwaysH$G1S <- signnature$G1S
pathwaysH$G2M <- signnature$G2M

pathwaysH <- pathwaysH[sapply(pathwaysH, length) >= 4]
list.of.sample.Seurat1  <- c("BT333","BT346","BT368","BT389","BT400","BT402","BT409")
#list.of.sample.Seurat1 <- c("BT324-GSC","BT326-GSC","BT333-GSC")#"BT363-GSC","BT368-GSC"

dimnz<- c(50)
lzp <- list()
lzp1 <- list()

for( prp in 1:length(dimnz)){
  
  for( j in 1:length(list.of.sample.Seurat1)){
    
    nam <- paste("data/",list.of.sample.Seurat1[j],".filtered_gene_matrices",sep="")
    nam2 <- paste("data/GlioMeta/Meta",list.of.sample.Seurat1[j],".rds",sep="")
    Seurat_temp1 <- routine(nam,nam2)
    HGT_cell_gs <- RunCellHGT(Seurat_temp1, pathways = pathwaysH , reduction = "mca", dims = 1:dimnz[prp],minSize = 4)[[1]]
    pancreas_gs_prediction <- rownames(HGT_cell_gs)[apply(HGT_cell_gs, 2, which.max)]
    Seurat_temp1$celltypetot <- pancreas_gs_prediction
    Pl = apply(HGT_cell_gs, 2,max)
    Seurat_temp1$pvalue <- Pl
    pancreas_gs_prediction_signif <- ifelse(apply(HGT_cell_gs, 2, max)>2, yes = pancreas_gs_prediction, "unasigned")
    Seurat_temp1$celltypesign <- pancreas_gs_prediction_signif
    
    nam <- paste("DataPerPy/SeuratObj/",list.of.sample.Seurat1[j],".h5ad",sep="")
    nam1 <- paste("DataPerPy/SeuratObj/",list.of.sample.Seurat1[j],".h5Seurat",sep="")
    SaveH5Seurat(Seurat_temp1, filename = nam1)
    Convert(nam1, dest = "h5ad")
        
  }
    
}





