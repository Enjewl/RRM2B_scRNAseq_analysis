library(tidyverse)
library(scCustomize)
library(extrafont)
library(Seurat)
library(tidyverse)

font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()

options(Seurat.object.assay.version = "v5")

#Downloaded feature/cell matrix (filtered) from https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-whole-transcriptome-analysis-1-standard-4-0-0

seurat <- Read10X(data.dir = "matrices/")

seurat <- CreateSeuratObject(counts = seurat, project = "Glioblastoma", min.cells = 3, min.features = 200)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

#Filter out cells with <200 genes and >10% mitochondrial reads
seurat <- subset(seurat, subset = nFeature_RNA >200 & percent.mt < 10)

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

#Dimensionality reduction and clustering
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), dims = 30)

seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat, resolution = 1.5)

seurat <- RunUMAP(seurat, dims = 1:30, min.dist = 0.03)

DimPlot(seurat, reduction = "umap", label = T, pt.size = 1)

#Labelling with sctype
library(sctype)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Brain"
gs_list <- gene_sets_prepare(db_, tissue)
scRNAseqData_scaled <- as.matrix(seurat[["RNA"]]@scale.data)

es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls <- do.call("rbind", lapply(unique(seurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat@meta.data[seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

seurat@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  seurat@meta.data$sctype_classification[seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(seurat, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')

Idents(seurat) <- seurat$sctype_classification

#Correlating cell barcodes with isoforms
Identifiers <- FetchData(object = seurat, vars = c("seurat_clusters", "orig.ident", "ident", "RRM2B"))
Isoform1 <- read.table("Isoform1.txt")
Isoform1 <- mutate(Isoform1, Isoform = "1")
Isoform2 <- read.table("Isoform2.txt")
Isoform2 <- mutate(Isoform2, Isoform = "2")

Isoform <- full_join(Isoform1, Isoform2)
colnames(Isoform) <- c("Barcode", "Isoform")

snRNA <- tibble::rownames_to_column(Identifiers, "Barcode")

Merged <- left_join(snRNA, Isoform)

Graph <- Merged %>%
  group_by(ident, Isoform) %>%
  summarise(Value = sum(n())) %>%
  mutate(total = sum(Value)) %>%
  mutate(Percentage = round(Value/total, 3)*100) %>%
  mutate(Isoform_Cells = paste0(ident, "_Isoform", Isoform))

Graph$Isoform <- as.character(Graph$Isoform)

ggplot(Graph, aes(y=`Value`, x=`ident`, fill = Isoform)) +
  geom_bar(stat="identity", color="black") +
  ylab("Number of cells expressing RRM2B isoforms") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c("grey", "white")) +
  theme(axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 45, face = "plain"),
        text=element_text(family="Arial", size = 16))

