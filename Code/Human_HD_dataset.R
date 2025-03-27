library(Seurat)
library(scCustomize)
library(BPCells)
library(tidyverse)
library(data.table)
library(ggrepel)
library(extrafont)

font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()

options(Seurat.object.assay.version = "v3")

seurat_list <- list()

#From compbio.mit.edu/scBBB/preprocessed/
BBB_HD_counts <- readRDS("brain.HD.snRNAseq.counts.rds")
BBB_HD_meta <- read.table("brain.HD.snRNAseq.metadata.txt", header = T)

###Create Seurat Object
BBB_HD_seurat <- CreateSeuratObject(counts = BBB_HD_counts, meta.data = BBB_HD_meta)

BBB_HD_seurat[["percent.mt"]] <- PercentageFeatureSet(BBB_HD_seurat, pattern = "^MT-|^mt-|^Mt-")

#Normalize and scale data
BBB_HD_seurat <- NormalizeData(BBB_HD_seurat)
BBB_HD_seurat <- FindVariableFeatures(BBB_HD_seurat)
BBB_HD_seurat <- ScaleData(BBB_HD_seurat)
BBB_HD_seurat <- RunPCA(BBB_HD_seurat, dims = 30)

#Dimensionality reduction and clustering
BBB_HD_seurat <- FindNeighbors(BBB_HD_seurat, reduction = "pca", dims = 1:30)
BBB_HD_seurat <- FindClusters(BBB_HD_seurat, resolution = 0.5)

BBB_HD_seurat <- RunUMAP(BBB_HD_seurat, reduction = "pca", dims = 1:30, min.dist = 0.01)

Idents(BBB_HD_seurat) <- BBB_HD_seurat$Condition

DimPlot(
  BBB_HD_seurat,
  reduction = "umap",
  label = T,
  repel = T
)

VlnPlot(BBB_HD_seurat, "RRM2B") +
  NoLegend() +
  xlab("") +
  ggtitle("") +
  scale_fill_manual(values = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey")) +
  theme(axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 45, face = "plain"),
        text=element_text(family="Arial", size = 16))
