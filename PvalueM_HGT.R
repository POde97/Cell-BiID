library(CelliD)
library(tidyverse)
library(ggpubr)
#source("myFunctions.R")



data("HgProteinCodingGenes")


dimnz <- c(50)
prp <- 1

samplez1 <- c("BT322","BT333","BT346","BT368","BT389","BT390","BT400","BT402","BT407","BT409")
samplez2 <- c("BT324-GSC","BT326-GSC","BT333-GSC","BT363-GSC","BT368-GSC") #-GSC
samplez3 <- c("BT338","BT363","BT364","BT397") #_1of2 2of2


#samplez1Vssamplez1
list.of.sample.Seurat1 <- samplez2
list.of.sample.Seurat2 <- list.of.sample.Seurat1

for( j in 1:length(samplez1)){
  
  nam1 <- paste("data/",list.of.sample.Seurat1[j],".filtered_gene_matrices",sep="")
  
  data_dir <- nam1
  list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
  data <- Read10X(data.dir = data_dir)
  seurat_object = CreateSeuratObject(counts = data)
  counts <- seurat_object@assays$RNA@counts
  counts <- counts[rownames(counts) %in% HgProteinCodingGenes,]
  Seurat_temp1 <- CreateSeuratObject(counts = counts ,min.cells = 5)
  Seurat_temp1 <- ScaleData(Seurat_temp1, features = rownames(Seurat_temp1))
  Seurat_temp1 <- RunMCA(Seurat_temp1,nmcs=dimnz[prp])
  
  
  
  list.of.sample.Seurat2 <- list.of.sample.Seurat2[list.of.sample.Seurat2 != list.of.sample.Seurat1[j]]
  for( i in 1:length(list.of.sample.Seurat2)){#
    if(list.of.sample.Seurat1[j]!=list.of.sample.Seurat2[i]){
      
      nam <-paste(list.of.sample.Seurat2[i],list.of.sample.Seurat1[j],"dimz",dimnz[prp],sep="-")
      nam1 <- paste("matrixPvalue/",nam,".csv",sep="")
      ceck <- file.exists(nam1)
      #ceck <- FALSE
      if(ceck==FALSE){
        
        nam1 <- paste("data/",list.of.sample.Seurat2[i],".filtered_gene_matrices",sep="")
        
        data_dir <- nam1
        list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
        data <- Read10X(data.dir = data_dir)
        seurat_object = CreateSeuratObject(counts = data)
        counts <- seurat_object@assays$RNA@counts
        counts <- counts[rownames(counts) %in% HgProteinCodingGenes,]
        Seurat_temp2 <- CreateSeuratObject(counts = counts ,min.cells = 5)
        
        Seurat_temp2 <- ScaleData(Seurat_temp2, features = rownames(Seurat_temp2))
        Seurat_temp2 <- RunMCA(Seurat_temp2,nmcs = dimnz[prp])
        Seurat_temp2_cell_gs <- GetCellGeneSet(Seurat_temp2, dims = 1:dimnz[prp], n.features = 200)
        
        
        rm(Seurat_temp2)
        HGT_cell_gs <- RunCellHGT(Seurat_temp1, pathways = Seurat_temp2_cell_gs, reduction = "mca", dims = 1:dimnz[prp],minSize = 10)
        #PvalueM <- HGT_cell_gs[[1]] 
        #rm(Seurat_temp1)
        
        nam <-paste(list.of.sample.Seurat2[i],list.of.sample.Seurat1[j],"dimz",dimnz[prp],sep="-")
        nam1 <- paste("matrixPvalue/",nam,".csv",sep="")
        write.csv(PvalueM,file=nam1) #Per ricostruzione di comunità !!
        
        
      }
    }
  }
  
  
  
}

