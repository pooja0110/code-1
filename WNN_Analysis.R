getwd()
install.packages('limma')
setwd("D:/wnn")
list.files()
library(Signac)
install.packages("Signac")
install.packages("Seurat")
library(Seurat)
install.packages("biovizBase")
BiocManager::install("biovizBase")
library(EnsDb.Mmusculus.v79)
install.packages("EnsDb.Mmusculus.v79")
BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
installed.packages("BiocManager")
#reading RNA+ATAC DATA
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
install.packages("hdf5r")
fragpath <- "atac_fragments.tsv.gz"

## get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
dtcells <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
dtcells

# create ATAC assay and add it to the object
dtcells[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
dtcells

DefaultAssay(dtcells) <- "ATAC"

dtcells <- NucleosomeSignal(dtcells)
dtcells <- TSSEnrichment(dtcells)

VlnPlot(
  object = dtcells,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

