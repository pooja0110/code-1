library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(tidyverse)
library(SeuratDisk)
library(ggplot2)
remotes::install_github("mojaveazure/seurat-disk")

install.packages("remotes")
library(patchwork)
#reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
setwd("D:\\wnn")
BiocManager::install(c('tidyverse'))
#Joint RNA-seq and ATAC-seq Data Analysis Using Weighted Nearest Neighbour (WNN) Analysis

# 1. Creating a object containing both RNA-seq and ATAC-seq data

# the 10x hdf5 file contains both data types. 
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")

fragpath <- "atac_fragments.tsv.gz"

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79 )
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA data
RNA_dt <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
# create ATAC assay and add it to the object
RNA_dt[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
memory.size() 
memory.limit(size=56000)
RNA_dt
DefaultAssay(RNA_dt) <- "ATAC"

RNA_dt <- NucleosomeSignal(RNA_dt)
RNA_dt <- TSSEnrichment(RNA_dt)
RNA_dt[["percent.mt"]] <- PercentageFeatureSet(RNA_dt, pattern = "^MT-")
library(dplyr)
VlnPlot(
  object = RNA_dt,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 4,
  pt.size = 0
)

RNA_dt <- subset(
  x = RNA_dt,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)
RNA_dt

# RNA analysis
DefaultAssay(RNA_dt) <- "RNA"
RNA_dt <- SCTransform(RNA_dt, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
install.packages("dplyr")
# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(RNA_dt) <- "ATAC"
RNA_dt <- RunTFIDF(RNA_dt)
RNA_dt <- FindTopFeatures(RNA_dt, min.cutoff = 'q0')
RNA_dt <- RunSVD(RNA_dt)
RNA_dt <- RunUMAP(RNA_dt, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
RNA_dt <- FindMultiModalNeighbors(RNA_dt, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
RNA_dt <- RunUMAP(RNA_dt, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
RNA_dt <- FindClusters(RNA_dt, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
RNA_dt
# perform sub-clustering on cluster 6 to find additional structure
RNA_dt <- FindSubCluster(RNA_dt, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(RNA_dt) <- "sub.cluster"
RNA_dt
RNA_dt$celltype <- Idents(RNA_dt)
p1 <- DimPlot(RNA_dt, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p1
p2 <- DimPlot(RNA_dt, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(RNA_dt, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# compute gene activities
gene.activities <- GeneActivity(RNA_dt)

# add the gene activity matrix to the Seurat object as a new assay
RNA_dt[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = RNA_dt,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(RNA_dt$nCount_RNA)
)