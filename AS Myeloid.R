# Import library----------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Import data ------------------------------------------       -----------
AS.list<-list()
AR1LRN <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AR1_2/outs/per_sample_outs/AR1LNR/count/sample_feature_bc_matrix")
AR2LRN <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AR1_2/outs/per_sample_outs/AR2LNR/count/sample_feature_bc_matrix")
AR3L <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AR3/outs/per_sample_outs/AR3L/count/sample_feature_bc_matrix")
AR4L <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AR4/outs/per_sample_outs/AR4L/count/sample_feature_bc_matrix")
AR4R <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AR4/outs/per_sample_outs/AR4R/count/sample_feature_bc_matrix")
AR4N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AR4/outs/per_sample_outs/AR4N/count/sample_feature_bc_matrix")


#AR3R AS7L exclude due to no cells
AR3N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AR3/outs/per_sample_outs/AR3N/count/sample_feature_bc_matrix")
AS1 <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS1/outs/filtered_feature_bc_matrix/")
AS2LR <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS2_3/outs/per_sample_outs/AS2LR/count/sample_feature_bc_matrix")
AS2N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS2_3/outs/per_sample_outs/AS2N/count/sample_feature_bc_matrix")
AS3L <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS2_3/outs/per_sample_outs/AS3L/count/sample_feature_bc_matrix")
AS3R <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS2_3/outs/per_sample_outs/AS3R/count/sample_feature_bc_matrix")
AS3N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS2_3/outs/per_sample_outs/AS3N/count/sample_feature_bc_matrix")
AS4L <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS4/outs/per_sample_outs/AS4L/count/sample_feature_bc_matrix")
AS4NR <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS4/outs/per_sample_outs/AS4NR/count/sample_feature_bc_matrix")
AS5L <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS5/outs/per_sample_outs/AS5L/count/sample_feature_bc_matrix")
AS5N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS5/outs/per_sample_outs/AS5N/count/sample_feature_bc_matrix")
AS5R <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS5/outs/per_sample_outs/AS5R/count/sample_feature_bc_matrix")
AS6L <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS6/outs/per_sample_outs/AS6L/count/sample_feature_bc_matrix")
AS6N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS6/outs/per_sample_outs/AS6N/count/sample_feature_bc_matrix")
AS6R <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS6/outs/per_sample_outs/AS6R/count/sample_feature_bc_matrix")
AS7L <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS7_8/outs/per_sample_outs/AS7L/count/sample_feature_bc_matrix")
#AS7N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS7_8/outs/per_sample_outs/AS7N/count/sample_feature_bc_matrix")
AS7R <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS7_8/outs/per_sample_outs/AS7R/count/sample_feature_bc_matrix")
AS8LR <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS7_8/outs/per_sample_outs/AS8LR/count/sample_feature_bc_matrix")
AS8N <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS7_8/outs/per_sample_outs/AS8N/count/sample_feature_bc_matrix")
AS9 <- Read10X(data.dir = "~/Desktop/Human_AS/CellRanger/AS9/outs/filtered_feature_bc_matrix/")


AS.list$AR1LRN <- CreateSeuratObject(counts = AR1LRN$`Gene Expression`, project = "A", min.cells = 3, min.features = 200)
AS.list$AR1LRN <- AddMetaData(AS.list$AR1LRN, "AR1LRN", col.name = "stim")

AS.list$AR2LRN <- CreateSeuratObject(counts = AR2LRN$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AR2LRN <- AddMetaData(AS.list$AR2LRN , "AR2LRN", col.name = "stim")

AS.list$AR3L <- CreateSeuratObject(counts = AR3L$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AR3L <- AddMetaData(AS.list$AR3L , "AR3L", col.name = "stim")

AS.list$AR3N <- CreateSeuratObject(counts = AR3N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AR3N <- AddMetaData(AS.list$AR3N , "AR3N", col.name = "stim")

AS.list$AR4L <- CreateSeuratObject(counts = AR4L$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AR4L <- AddMetaData(AS.list$AR4L , "AR4L", col.name = "stim")

AS.list$AR4N <- CreateSeuratObject(counts = AR4N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AR4N <- AddMetaData(AS.list$AR4N , "AR4N", col.name = "stim")

AS.list$AR4R <- CreateSeuratObject(counts = AR4R$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AR4R <- AddMetaData(AS.list$AR4R , "AR4R", col.name = "stim")


AS.list$AS1<- CreateSeuratObject(counts = AS1, project = "A", min.cells = 3, min.features = 200)
AS.list$AS1 <- AddMetaData(AS.list$AS1 , "AS1", col.name = "stim")

AS.list$AS2LR <- CreateSeuratObject(counts = AS2LR$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS2LR<- AddMetaData(AS.list$AS2LR , "AS2LR", col.name = "stim")

AS.list$AS2N <- CreateSeuratObject(counts = AS2N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS2N <- AddMetaData(AS.list$AS2N , "AS2N", col.name = "stim")

AS.list$AS3N <- CreateSeuratObject(counts = AS3N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS3N <- AddMetaData(AS.list$AS3N , "AS3N", col.name = "stim")

AS.list$AS3L <- CreateSeuratObject(counts = AS3L$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS3L <- AddMetaData(AS.list$AS3L , "AS3L", col.name = "stim")

AS.list$AS3R <- CreateSeuratObject(counts = AS3R$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS3R <- AddMetaData(AS.list$AS3R , "AS3R", col.name = "stim")

AS.list$AS4L <- CreateSeuratObject(counts = AS4L$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS4L <- AddMetaData(AS.list$AS4L , "AS4L", col.name = "stim")

AS.list$AS4NR <- CreateSeuratObject(counts = AS4NR$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS4NR <- AddMetaData(AS.list$AS4NR , "AS4NR", col.name = "stim")

AS.list$AS5L <- CreateSeuratObject(counts = AS5L$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS5L <- AddMetaData(AS.list$AS5L , "AS5L", col.name = "stim")

AS.list$AS5N <- CreateSeuratObject(counts = AS5N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS5N <- AddMetaData(AS.list$AS5N , "AS5N", col.name = "stim")

AS.list$AS5R <- CreateSeuratObject(counts = AS5R$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS5R <- AddMetaData(AS.list$AS5R , "AS5R", col.name = "stim")

AS.list$AS6L <- CreateSeuratObject(counts = AS6L$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS6L <- AddMetaData(AS.list$AS6L , "AS6L", col.name = "stim")

AS.list$AS6N <- CreateSeuratObject(counts = AS6N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS6N <- AddMetaData(AS.list$AS6N , "AS6N", col.name = "stim")

AS.list$AS6R <- CreateSeuratObject(counts = AS6R$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS6R <- AddMetaData(AS.list$AS6R , "AS6R", col.name = "stim")

AS.list$AS7L <- CreateSeuratObject(counts = AS7L$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS7L <- AddMetaData(AS.list$AS7L , "AS7L", col.name = "stim")

#AS.list$AS7N <- CreateSeuratObject(counts = AS7N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
#AS.list$AS7N <- AddMetaData(AS.list$AS7N , "AS7N", col.name = "stim")

AS.list$AS7R <- CreateSeuratObject(counts = AS7R$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS7R <- AddMetaData(AS.list$AS7R , "AS7R", col.name = "stim")

AS.list$AS8LR <- CreateSeuratObject(counts = AS8LR$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS8LR <- AddMetaData(AS.list$AS8LR , "AS8LR", col.name = "stim")

AS.list$AS8N <- CreateSeuratObject(counts = AS8N$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
AS.list$AS8N <- AddMetaData(AS.list$AS8N , "AS8N", col.name = "stim")

AS.list$AS9 <- CreateSeuratObject(counts = AS9, project = "A", min.cells = 3, min.features = 200)
AS.list$AS9 <- AddMetaData(AS.list$AS9 , "AS9", col.name = "stim")

# Add sample-specific labels to the barcodes
for (j in 1:length(AS.list)) {
  for (i in 1:length(AS.list[[j]]@assays[["RNA"]]@data@Dimnames[[2]])){
    AS.list[[j]]@assays[["RNA"]]@data@Dimnames[[2]][i] <- paste0(AS.list[[j]]@assays[["RNA"]]@data@Dimnames[[2]][i],"_",j)
    AS.list[[j]]@assays[["RNA"]]@counts@Dimnames[[2]][i] <- paste0(AS.list[[j]]@assays[["RNA"]]@counts@Dimnames[[2]][i],"_",j)
  }
}

AS.list$AR1LRN[["percent.mt"]] <- PercentageFeatureSet(AS.list$AR1LRN, pattern = "^MT-")
AS.list$AR2LRN[["percent.mt"]] <- PercentageFeatureSet(AS.list$AR2LRN, pattern = "^MT-")
AS.list$AR3L[["percent.mt"]] <- PercentageFeatureSet(AS.list$AR3L, pattern = "^MT-")
AS.list$AR3N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AR3N, pattern = "^MT-")
AS.list$AR4L[["percent.mt"]] <- PercentageFeatureSet(AS.list$AR4L, pattern = "^MT-")
AS.list$AR4N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AR4N, pattern = "^MT-")
AS.list$AR4R[["percent.mt"]] <- PercentageFeatureSet(AS.list$AR4R, pattern = "^MT-")

AS.list$AS1[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS1, pattern = "^MT-")
AS.list$AS2LR[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS2LR, pattern = "^MT-")
AS.list$AS2N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS2N, pattern = "^MT-")
AS.list$AS3L[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS3L, pattern = "^MT-")
AS.list$AS3R[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS3R, pattern = "^MT-")
AS.list$AS3N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS3N, pattern = "^MT-")
AS.list$AS4L[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS4L, pattern = "^MT-")
AS.list$AS4NR[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS4NR, pattern = "^MT-")
AS.list$AS5L[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS5L, pattern = "^MT-")
AS.list$AS5N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS5N, pattern = "^MT-")
AS.list$AS5R[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS5R, pattern = "^MT-")
AS.list$AS6L[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS6L, pattern = "^MT-")
AS.list$AS6N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS6N, pattern = "^MT-")
AS.list$AS6R[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS6R, pattern = "^MT-")
AS.list$AS7L[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS7L, pattern = "^MT-")
#AS.list$AS7N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS7N, pattern = "^MT-")
AS.list$AS7R[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS7R, pattern = "^MT-")
AS.list$AS8LR[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS8LR, pattern = "^MT-")
AS.list$AS8N[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS8N, pattern = "^MT-")
AS.list$AS9[["percent.mt"]] <- PercentageFeatureSet(AS.list$AS9, pattern = "^MT-")

VlnPlot(AS.list$AS2LR,, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Cells with extremely high percent.mt (>12%) or very low (≤300) or very high (>4500) nFeature_RNA are likely to be dead cells or doublets.
plot1 <- FeatureScatter(AS.list$AS2LR, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AS.list$AS2LR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter low-quality cells ----------------------------------------------------
# Retain only cells with mitochondrial RNA content ≤12% and between 300 and 4,500 detected genes.
#AR
AS.list$AR1LRN<- subset(AS.list$AR1LRN, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AR2LRN <- subset(AS.list$AR2LRN, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AR3L<- subset(AS.list$AR3L, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AR3N<- subset(AS.list$AR3N, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AR4L<- subset(AS.list$AR4L, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AR4R<- subset(AS.list$AR4R, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AR4N<- subset(AS.list$AR4N, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
#AS
AS.list$AS1<- subset(AS.list$AS1, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS2LR<- subset(AS.list$AS2LR, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS2N<- subset(AS.list$AS2N, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS3L<- subset(AS.list$AS3L, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS3R<- subset(AS.list$AS3R, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS3N<- subset(AS.list$AS3N, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS4L<- subset(AS.list$AS4L, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS4NR<- subset(AS.list$AS4NR, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS5L<- subset(AS.list$AS5L, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS5N<- subset(AS.list$AS5N, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS5R<- subset(AS.list$AS5R, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS6L<- subset(AS.list$AS6L, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS6N<- subset(AS.list$AS6N, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS6R<- subset(AS.list$AS6R, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS7L<- subset(AS.list$AS7L, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
#AS.list$AS7N<- subset(AS.list$AS7N, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
AS.list$AS7R<- subset(AS.list$AS7R, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS8LR<- subset(AS.list$AS8LR, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS8N<- subset(AS.list$AS8N, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)
AS.list$AS9<- subset(AS.list$AS9, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & percent.mt < 12)

# Confirm that low-quality cells have been filtered out.
plot1 <- FeatureScatter(AS.list$AS2LR, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AS.list$AS2LR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Data normalization-----------------------------------------------
all.genes <- append(rownames(AS.list$AR1LRN),rownames(AS.list$AR2LRN))
all.genes <- append(all.genes,rownames(AS.list$AR3L))
all.genes <- append(all.genes,rownames(AS.list$AR3N))
all.genes <- append(all.genes,rownames(AS.list$AR4L))
all.genes <- append(all.genes,rownames(AS.list$AR4R))
all.genes <- append(all.genes,rownames(AS.list$AR4N))

all.genes <- append(all.genes,rownames(AS.list$AS1))
all.genes <- append(all.genes,rownames(AS.list$AS2LR))
all.genes <- append(all.genes,rownames(AS.list$AS2N))
all.genes <- append(all.genes,rownames(AS.list$AS3N))
all.genes <- append(all.genes,rownames(AS.list$AS3R))
all.genes <- append(all.genes,rownames(AS.list$AS3L))
all.genes <- append(all.genes,rownames(AS.list$AS4L))
all.genes <- append(all.genes,rownames(AS.list$AS4NR))
all.genes <- append(all.genes,rownames(AS.list$AS5L))
all.genes <- append(all.genes,rownames(AS.list$AS5N))
all.genes <- append(all.genes,rownames(AS.list$AS5R))
all.genes <- append(all.genes,rownames(AS.list$AS6L))
all.genes <- append(all.genes,rownames(AS.list$AS6N))
all.genes <- append(all.genes,rownames(AS.list$AS6R))
all.genes <- append(all.genes,rownames(AS.list$AS7L))
#all.genes <- append(all.genes,rownames(AS.list$AS7N))
all.genes <- append(all.genes,rownames(AS.list$AS7R))
all.genes <- append(all.genes,rownames(AS.list$AS8LR))
all.genes <- append(all.genes,rownames(AS.list$AS8N))
all.genes <- append(all.genes,rownames(AS.list$AS9))

AS.list <- lapply(X = AS.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- ScaleData(x, features = all.genes)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- RunPCA(x, features = VariableFeatures(object = x), verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = AS.list)

# Integration
AS.anchors <- FindIntegrationAnchors(object.list = AS.list, anchor.features = features, reduction = "rpca",k.filter = 35)

# This command creates an 'integrated' data assay
AS.combined <- IntegrateData(anchorset =AS.anchors, k.weight = 35)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(AS.combined) <- "integrated"

AS.combined <- ScaleData(AS.combined, verbose = FALSE)
AS.combined <- RunPCA(AS.combined, npcs = 30, verbose = FALSE)
AS.combined <- RunUMAP(AS.combined, reduction = "pca", dims = 1:30)
AS.combined <- FindNeighbors(AS.combined, reduction = "pca", dims = 1:30)
AS.combined <- FindClusters(AS.combined, resolution = 0.10)

DimPlot(AS.combined, reduction = "umap", split.by = "stim")
DimPlot(AS.combined, reduction = "umap", split.by = "ASAR", label = TRUE) 
DimPlot(AS.combined, reduction = "umap", split.by = "Sample", label = TRUE) 
DimPlot(AS.combined, reduction = "umap", label = TRUE, repel = TRUE)

levels(AS.combined) # [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9""10""11"
levels(AS.combined) <- c("0", "1", "3", "7", "8","10","11","2","9", "5","6","4")

new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11")
names(new.cluster.ids) <- levels(AS.combined)
AS.combined<- RenameIdents(AS.combined, new.cluster.ids)

# Extract the stim column
cell.stim <- AS.combined@meta.data[["stim"]]
# Assign disease and non-disease groups based on the "stim" column
newgroup <- c("AR","AR","AR","AR","AR","AR","AR","AS","AS","AS","AS", "AS","AS","AS","AS","AS","AS","AS","AS","AS","AS","AS","AS","AS","AS","AS")
names(newgroup) <- c("AR1LRN","AR2LRN","AR3L","AR3N","AR4L","AR4R","AR4N","AS1","AS2LR","AS2N","AS3N","AS3R","AS3L","AS4L","AS4NR","AS5L","AS5N","AS5R","AS6L","AS6N","AS6R","AS7L","AS7R","AS8LR","AS8N","AS9")
# Convert the "stim" column to the new group labels
cell.newgroup <- newgroup[cell.stim]
# Add the new group as metadata
AS.combined <- AddMetaData(AS.combined, cell.newgroup, col.name = "ASAR")

#newgroup2
newgroup2 
newgroup2 <- c("AR1","AR2","AR3","AR3","AR4","AR4","AR4","AS1","AS2","AS2","AS3","AS3", "AS3","AS4","AS4","AS5","AS5","AS5","AS6","AS6","AS6","AS7","AS7","AS8","AS8","AS9")
names(newgroup2) <- c("AR1LRN","AR2LRN","AR3L","AR3N", "AR4L", "AR4R", "AR4N","AS1","AS2LR","AS2N","AS3N","AS3R","AS3L","AS4L","AS4NR","AS5L","AS5N","AS5R","AS6L","AS6N","AS6R","AS7L","AS7R","AS8LR","AS8N","AS9")
cell.newgroup2 <- newgroup2[cell.stim]
AS.combined <- AddMetaData(AS.combined, cell.newgroup2, col.name = "Sample")


# Heatmap------------------------------------------------------------------------------
AS.combined.markers <- FindAllMarkers(AS.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AS.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AS.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AS.combined, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

# Save as RDS file
saveRDS(AS.combined, "~/Desktop/Human_AS/RDSfiles/AS_combined_saved_240122.rds")
AS.combined<- readRDS ("~/Desktop/Human_AS/RDSfiles/AS_combined_saved_240122.rds")

# Cell number count
AR.combined1LRN <- subset(AS.combined, subset = stim == c("AR1LRN"))#564sample
AR.combined2LRN <- subset(AS.combined, subset = stim == c("AR2LRN"))#5301sample
AR.combined3L <- subset(AS.combined, subset = stim == c("AR3L"))#245sample
AR.combined3N <- subset(AS.combined, subset = stim == c("AR3N"))#782sample
AR.combined4L <- subset(AS.combined, subset = stim == c("AR4L"))#401sample
AR.combined4R <- subset(AS.combined, subset = stim == c("AR4R"))#112sample
AR.combined4N <- subset(AS.combined, subset = stim == c("AR4N"))#253sample

AS.combined1 <- subset(AS.combined, subset = stim == c("AS1"))#1265sample
AS.combined2LR <- subset(AS.combined, subset = stim == c("AS2LR"))#196sample
AS.combined2N <- subset(AS.combined, subset = stim == c("AS2N"))#994sample
AS.combined3L <- subset(AS.combined, subset = stim == c("AS3L"))#3018sample
AS.combined3R <- subset(AS.combined, subset = stim == c("AS3R"))#281sample
AS.combined3N <- subset(AS.combined, subset = stim == c("AS3N"))#2979sample
AS.combined4L <- subset(AS.combined, subset = stim == c("AS4L"))#948sample
AS.combined4NR <- subset(AS.combined, subset = stim == c("AS4NR"))#2449sample
AS.combined5L <- subset(AS.combined, subset = stim == c("AS5L"))#418sample
AS.combined5N <- subset(AS.combined, subset = stim == c("AS5N"))#129sample
AS.combined5R <- subset(AS.combined, subset = stim == c("AS5R"))#202sample
AS.combined6L <- subset(AS.combined, subset = stim == c("AS6L"))#350sample
AS.combined6N <- subset(AS.combined, subset = stim == c("AS6N"))#257sample
AS.combined6R <- subset(AS.combined, subset = stim == c("AS6R"))#77sample
AS.combined7L <- subset(AS.combined, subset = stim == c("AS7L"))#73sample
AS.combined7R <- subset(AS.combined, subset = stim == c("AS7R"))#201sample
AS.combined8LR <- subset(AS.combined, subset = stim == c("AS8LR"))#365sample
AS.combined8N <- subset(AS.combined, subset = stim == c("AS8N"))#498sample
AS.combined9 <- subset(AS.combined, subset = stim == c("AS9"))#201sample

"Lymphoid"
VlnPlot(AS.combined,pt.size = 0, features = c("CD8A","IL7R","CD40LG") )
VlnPlot(AS.combined,pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))
VlnPlot(AS.combined,pt.size = 0, features = c("MKI67","CD79A","CD79B","COL1A1","IBSP"))

"Myeloid"
VlnPlot(AS.combined,pt.size = 0, features = c("CD14","C1QA","CD68","MKI67","CD79A","CD3E","KIT") )
VlnPlot(AS.combined,pt.size = 0, features = c("CD68","CD14","CSF1R","","S100A8", "S100A9") )
VlnPlot(AS.combined,pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","TNF","CXCL8") )
VlnPlot(AS.combined,pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AS.combined,pt.size = 0, features = c("CD3E","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1B","CXCR2") )
VlnPlot(AS.combined,pt.size = 0, features = c("SIRPA","CD47","SIRPB1","SIRPB2","TNF","SLAMF7"), split.by = "ASAR" )
VlnPlot(AS.combined,pt.size = 0, features = c("CSF1","CCL2","TNF","TNFRS1A","TNFRS1B","NRF2"), split.by = "ASAR")

# Heatmap
AS.combined.markers <- FindAllMarkers(AS.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AS.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AS.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AS.combined, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)
DoHeatmap(AS.combined, features = top10$gene) 

AS.combined_trans <- SCTransform (AS.combined)

cd_genes <- c("CD3E" ,"CD8A","GZMA","CD14","CD68","CD79B","IGHG1","ACTA2","COL1A1")
DotPlot(AS.combined_trans,features = cd_genes)+RotatedAxis()+coord_flip()

#Myeloid cell extraction-------
AS.combined_1 <- subset(AS.combined, idents = c(2,9))
AS.combined_1 <- ScaleData(AS.combined_1, vedsbose = FALSE)
AS.combined_1<- RunPCA(AS.combined_1, npcs = 30, verbose = FALSE)
AS.combined_1 <- RunUMAP(AS.combined_1 , reduction = "pca", dims = 1:30)
AS.combined_1 <- FindNeighbors(AS.combined_1 , reduction = "pca", dims = 1:30)
AS.combined_1 <- FindClusters(AS.combined_1, resolution = 0.10)

DimPlot(AS.combined_1, reduction = "umap", split.by = "ASAR", label = TRUE, repel = TRUE)
DimPlot(AS.combined_1, reduction = "umap", split.by = "stim")

VlnPlot(AS.combined_1, pt.size = 0, features = c("CD14","CD3E","KIT","GZMA","FCBR3B","S100A8", "S100A9","FCGR3A","FCGR3B","HLA-DQA1","CLEC10A","FCER1A","MKI67","CLEC9A") )

# Remove T cell and Mast cell
# Remove CD3 contamination-------
AS.combined_1_1 <- subset(AS.combined_1, idents = c(0,1,2,3,5,6,7,8,9,10,11))
AS.combined_1_1 <- ScaleData(AS.combined_1_1, verbose = FALSE)
AS.combined_1_1<- RunPCA(AS.combined_1_1, npcs = 50, verbose = FALSE)
AS.combined_1_1 <- RunUMAP(AS.combined_1_1, reduction = "pca", dims = 1:8)
AS.combined_1_1 <- FindNeighbors(AS.combined_1_1, reduction = "pca", dims = 1:8)
AS.combined_1_1 <- FindClusters(AS.combined_1_1, resolution = 0.40)

DimPlot(AS.combined_1_1, reduction = "umap", pt.size = 1.0, split.by = "Sample", label = TRUE, repel = TRUE)
DimPlot(AS.combined_1_1, reduction = "umap",split.by = "ASAR")

# Save as RDS file
saveRDS(AS.combined_1_1, "~/Desktop/Human_AS/RDSfiles/AS_combined_myeloid_240122.rds")
AS.combined_1_1<- readRDS ("~/Desktop/Human_AS/RDSfiles/AS_combined_myeloid_2312240122.rds")

# myeloid-CD3 contamination remove again----
AS.combined_1_2 <- subset(AS.combined_1_1, idents = c(0,1,2,3,4,5,6))
AS.combined_1_2 <- ScaleData(AS.combined_1_2, verbose = FALSE)
AS.combined_1_2<- RunPCA(AS.combined_1_2, npcs = 50, verbose = FALSE)
AS.combined_1_2 <- RunUMAP(AS.combined_1_2, reduction = "pca", dims = 1:6)
AS.combined_1_2 <- FindNeighbors(AS.combined_1_2, reduction = "pca", dims = 1:6)
AS.combined_1_2 <- FindClusters(AS.combined_1_2, resolution = 0.30)
DimPlot(AS.combined_1_2, reduction = "umap",split.by = "ASAR")

DimPlot(AS.combined_1_2, reduction = "umap", pt.size = 1.0, split.by = "Sample", label = TRUE, repel = TRUE)
DimPlot(AS.combined_1_2, reduction = "umap",split.by = "ASAR")

### Save as RDS file (!!This is the final RDS file!!)###
saveRDS(AS.combined_1_2, "~/Desktop/Human_AS/RDSfiles/AS_combined_myeloid_240125.rds")
AS.combined_1_2<- readRDS ("~/Desktop/Human_AS/RDSfiles/AS_combined_myeloid_240125.rds")

###Devide into AS and AR
AS.combined_1_2_AS <- subset(AS.combined_1_2, subset = ASAR =="AS")
AS.combined_1_2_AR <- subset(AS.combined_1_2, subset = ASAR =="AR")

###Count cell number of the cluster
AS.combined_1_2_AS_0 <- subset(AS.combined_1_2_AS, idents = 0)
AS.combined_1_2_AS_1 <- subset(AS.combined_1_2_AS, idents = 1)
AS.combined_1_2_AS_2 <- subset(AS.combined_1_2_AS, idents = 2)
AS.combined_1_2_AS_3 <- subset(AS.combined_1_2_AS, idents = 3)
AS.combined_1_2_AS_4 <- subset(AS.combined_1_2_AS, idents = 4)
AS.combined_1_2_AS_5 <- subset(AS.combined_1_2_AS, idents = 5)
AS.combined_1_2_AS_6 <- subset(AS.combined_1_2_AS, idents = 6)

AS.combined_1_2_AR_0 <- subset(AS.combined_1_2_AR, idents = 0)
AS.combined_1_2_AR_1 <- subset(AS.combined_1_2_AR, idents = 1)
AS.combined_1_2_AR_2 <- subset(AS.combined_1_2_AR, idents = 2)
AS.combined_1_2_AR_3 <- subset(AS.combined_1_2_AR, idents = 3)
AS.combined_1_2_AR_4 <- subset(AS.combined_1_2_AR, idents = 4)
AS.combined_1_2_AR_5 <- subset(AS.combined_1_2_AR, idents = 5)
AS.combined_1_2_AR_6 <- subset(AS.combined_1_2_AR, idents = 6)


# Extract Macrophage
AS.combined_1_3 <- subset(AS.combined_1_2, idents = c(0,1,2,3))
AS.combined_1_3 <- ScaleData(AS.combined_1_3, verbose = FALSE)
AS.combined_1_3<- RunPCA(AS.combined_1_3, npcs = 50, verbose = FALSE)
AS.combined_1_3 <- RunUMAP(AS.combined_1_3, reduction = "pca", dims = 1:6)
AS.combined_1_3 <- FindNeighbors(AS.combined_1_3, reduction = "pca", dims = 1:6)
AS.combined_1_3 <- FindClusters(AS.combined_1_3, resolution = 0.30)
DimPlot(AS.combined_1_3, reduction = "umap",split.by = "ASAR")

FeaturePlot(AS.combined_1_2, features=c('IL1B'), split.by = "ASAR",   min.cutoff=0.3, max.cutoff='q90' )
FeaturePlot(AS.combined_1_2, features=c('IL1B'),   min.cutoff=0.3, max.cutoff='q90' )
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("IL1B"),split.by = "ASAR")

FeaturePlot(AS.combined_1_2, features=c('NFKB1'), split.by = "ASAR",   min.cutoff=0.1, max.cutoff='q80' )
FeaturePlot(AS.combined_1_2, features=c('NFKB1'),   min.cutoff=0.1, max.cutoff='q80' )
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("NFKB1"),split.by = "ASAR")

VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("TNF"),split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("IFNG"),split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("HCST"),split.by = "ASAR")

# Cell death related genes
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("ZBP1","CASP8","RIPK1","RIPK3","FADD","AIM2","PYCARD","NLRP3"), split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("CASP1","CASP4","CASP5","GSDMD","IL1B","IL18","NLRP3","AIM2","PYCARD"), split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("CASP3","CASP7","CASP8","BAX","BAK1","BCL2L11","FAS","FASLG","TNFRSF10B"), split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("RIPK1","RIPK3","MLKL","TNF","TNFRSF1A","FADD","CASP8"), split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("ZBP1","CASP8","RIPK3","MLKL","GSDMD","CASP1","IL1B","FADD","TNF"), split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0.1,assay ="RNA", features = c("CASP3","CASP8","CASP1","NLRP3","IL1B","IL18","BAX","RIPK1","RIPK3"), split.by = "ASAR",cols = c("AR" = "blue", "AS" = "red"))

FeaturePlot(AS.combined_1_2, features=c('TNFR2'), split.by = "ASAR",   min.cutoff=0.1, max.cutoff='q80' )
FeaturePlot(AS.combined_1_2, features=c('KLRK1'), split.by = "ASAR",   min.cutoff=0.1, max.cutoff='q90' )
FeaturePlot(AS.combined_1_2, features=c('HCST'), split.by = "ASAR",   min.cutoff=0.1, max.cutoff='q90' )
FeaturePlot(AS.combined_1_2, features=c('TNFRSA'), split.by = "ASAR",  min.cutoff=0.1, max.cutoff='q90',keep.scale = 'all')

p <- FeaturePlot(
  AS.combined_1_2,
  features = "HCST",
  split.by = "ASAR",
  min.cutoff = 1,
  max.cutoff = "q90",
  keep.scale = "all"
)

p + theme(
  legend.position = "bottom",
  legend.direction = "horizontal"
)

VlnPlot(AS.combined_1_2,pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","IBSP") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("CXCL8","LYVE1","FOLR2","CD163","SPP1","MMP9") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("NFE2L2","HMOX1","TXNRD1","SPP1","CCR2","MARCO") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("IFIT1","IFIT2","IFIT3","VCAN","FCGR3A","FCGR3B") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("CD38","Il4l1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("CD38","Il4l1","IDO1","CD40","CD274","CD1C") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("CD47","SIRPA","COLEC12","MRC1","STAB1", "MARCO"))

VlnPlot(AS.combined_1_2,pt.size = 0, features = c("SIRPA","CD47","SIRPB1","SIRPB2","TNF","SLAMF7"))
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("FSCN1","CCR7","PDL1","CD86","TIM3","IL4R"))
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("NFKB1","NFKB2","TNF","IL1B","IL1R2","IL6") )
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("IL10","IL12","FTL","IL4","CD3E","TIMP1"))
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("MSR1","SCARB1","COLEC12","SCARA5","CD163", "MARCO"))
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("CSF1","CSF1R","CCL2","SCARA5","CD163", "CCR2"),split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0, assay ="RNA", features = c("LYVE1","FOLR2","MRC1","DAB2","COLEC12", "STAB1"),split.by = "ASAR", cols = c("AR" = "blue", "AS" = "red"))
VlnPlot(AS.combined_1_2,pt.size = 0, assay ="RNA", features = c("C1QA","IL1B","TET2","AREG","EREG","HBEGF"),split.by = "ASAR", cols = c("AR" = "blue", "AS" = "red"))
VlnPlot(AS.combined_1_2,pt.size = 0, assay ="RNA", features = c("KLRK1","HCST","CD226","TIGIT","PDCD1","FASLG","TNFSF10", "IFNG","ITGAL","ITGB2"),split.by = "ASAR", cols = c("AR" = "blue", "AS" = "red"))
VlnPlot(AS.combined_1_2,pt.size = 0, assay ="RNA", features = c("KLRK1","HCST"),split.by = "ASAR", cols = c("AR" = "blue", "AS" = "red"))

# scavengers related genes
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("MSR1","SCARB1","COLEC12","SCARA5","CD163", "MARCO"),split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("COLEC12","DAB2","MRC1", "MSR1", "SCARA5","CD163"),split.by = "ASAR",assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("SIRPA","CD47","MRC1", "MSR1", "SCARA5","CD163"),split.by = "ASAR")
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("COLEC12","DAB2","MRC1", "MSR1", "SCARA5","CD163"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("SQSTM1","NBR1","OPTN", "CALCOCO2", "TAX1BP1","BNIP3"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("LMAN1","LMAN2","SURF4", "TMED2", "TMED9"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("SORT1","SORL1","LRP1", "LDLR", "M2PR"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("CCPG1","FAM134B","TOLLIP", "WDFY3","FUNDC1"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("DAB2","MAF","STAB1", "SLCO2B1","CD163","RNASE1"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("TTN","SELENOP","FOLR2", "MRC1","COLEC12","SLC40A1"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("LGMN","PDK4","JUN", "HIF1A","MERTK","TIMD4"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("LYVE1","FOLR2","CD163", "MRC1","SLC40A1","LRP1"),assay = 'RNA')
VlnPlot(AS.combined_1_2,pt.size = 0, features = c("IGF1","STAB1","COLEC12", "DAB2","MERTK","PDK4"),assay = 'RNA')


cDC1---
  VlnPlot(AS.combined_1_2, pt.size = 0, features = c("CPVL","WDFY4","CADM1","CLEC9A","TNF") )
cDC2---
  VlnPlot(AS.combined_1_2, pt.size = 0, features = c("HLA-DQA1","CLEC10A","FCER1A","CD1C","PAK1") )
mREGDC---
  VlnPlot(AS.combined_1_2, pt.size = 0, features = c("FSCN1","CCR7","BIRC3","LAMP3","CCL22","EBI3") )
preDC-----
  VlnPlot(AS.combined_1_2, pt.size = 0, features = c("TCF4","CXCR3","IRF4","IFNA","CD86") )
CD16Mo----
  VlnPlot(AS.combined_1_2, pt.size = 0, features = c("FCGR3A","NR4A1","S100A4","LST1","IFITM2") )
cMo------
  VlnPlot(AS.combined_1_2, pt.size = 0, features = c("IL1B","VEGFA","CCL20","IL1RN","CXCL3","EREG") )
VlnPlot(AS.combined_1_2, pt.size = 0, features = c("PLAUR","CCL3","CCL4","TNFAIP3","CD44") )

cMo------
  VlnPlot(AS.combined_1_2, pt.size = 0, features = c("HSPA1","HSPA6","FOS","LYZ","","") )
VlnPlot(AS.combined_1_2, pt.size = 0, features = c("HBB","NPC2","ATOX1","STAT1","KLF6") )
VlnPlot(AS.combined_1_2, pt.size = 0, features = c("HSPA1","HSPA6","FOS","N","CD3E","CD69") )

# Dot plot myeloid---
AS.combined_1_2_trans <- SCTransform (AS.combined_1_2)
cd_genes <- c("CD68","CD14","FOSB","IL1B","TNF", "LYVE1","FOLR2","C1QA","CD9","TREM2","VCAN","LYZ","FCER1A","HLA-DQA1","CLEC9A","CLEC10A")
DotPlot(AS.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

# Dot plot myeloid---
AS.combined_1_2_trans <- SCTransform (AS.combined_1_2)
cd_genes <- c("CD68","CD14","LYVE1","FOLR2","CD163","MRC1","IL1B","TNF","CXCL2","CCL3L1","C1QA","FTL","CD9","TREM2","S100A8","S100A9","S100A12","VCAN","FCER1A","HLA-DQA1","CLEC10A","CCR7","CCL22")
DotPlot(AS.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

# volcano  Myeloid------
AS.combined_1_2$celltype <- Idents(AS.combined_1_2)
AS.combined_1_2$celltype.stim <- paste(Idents(AS.combined_1_2), AS.combined_1_2$ASAR, sep="_")
Idents(AS.combined_1_2) <- "celltype.stim" #Update the identities using the Idents() function
levels(AS.combined_1_2) #Check that the identities have been updated

AS.combined_1_2.table <- FindMarkers(AS.combined_1_2, ident.1 = c("1_AS", "4_AS"), 
                                     ident.2 =c("1_AR", "4_AR"), verbose = FALSE, logfc.threshold = 0)
AS.combined_1_2.table$logp <- -log10(AS.combined_1_2.table$p_val)
AS.combined_1_2_filtered_left = subset(AS.combined_1_2.table, logp>=2.0 & avg_log2FC <= -0.5)
AS.combined_1_2_filtered_right = subset(AS.combined_1_2.table, logp>=2.0 & avg_log2FC >= 0.5)

genes.to.label.left <- rownames(AS.combined_1_2_filtered_left)
genes.to.label.right <- rownames(AS.combined_1_2_filtered_right)

p1 <- ggplot(AS.combined_1_2.table, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

# Heatmap
AS.combined_1_2.markers <- FindAllMarkers(AS.combined_1_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AS.combined_1_2.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

top15 <- AS.combined_1_2.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(AS.combined_1_2, features = top15$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)
DoHeatmap(AS.combined_1_2, features = top10$gene) 

# Proportion
AS.combined_1_2_prop. <- prop.table(table(Idents(AS.combined_1_2),AS.combined_1_2$ASAR))
write.table(AS.combined_1_2_prop., "~/Desktop/AS.combined_1_2_prop_240123.txt", quote=F, col.names=F, append=T)
write.csv(AS.combined_1_2_prop., "AS.combined_1_2.csv")

AS.combined_1_2_prop. <- prop.table(table(Idents(AS.combined_1_2),AS.combined_1_2$Sample))
write.table(AS.combined_1_2_prop., "~/Desktop/AS.combined_1_2_prop_240123.txt", quote=F, col.names=F, append=T)
write.csv(AS.combined_1_2_prop., "AS.combined_1_2.csv")


# T cells
AS.combined_2 <- subset(AS.combined, idents = c(0,1,2,3,4,5,6))
AS.combined_2 <- ScaleData(AS.combined_2, verbose = FALSE)
AS.combined_2<- RunPCA(AS.combined_2, npcs = 30, verbose = FALSE)
AS.combined_2 <- RunUMAP(AS.combined_2 , reduction = "pca", dims = 1:10)
AS.combined_2 <- FindNeighbors(AS.combined_2 , reduction = "pca", dims = 1:10)
AS.combined_2 <- FindClusters(AS.combined_2, resolution = 0.20)

DimPlot(AS.combined_2, reduction = "umap", split.by = "ASAR", repel = TRUE)
DimPlot(AS.combined_2, reduction = "umap", split.by = "stim")

VlnPlot(AS.combined_2,pt.size = 0, assay = "RNA", features = c("CD3E","CD4","CD8A","GZMB","NKG7","TYROBP"))
VlnPlot(AS.combined_2,pt.size = 0, features = c("CD3E","CD4","CD8A","GZMB","TYROBP","MKI67"))

levels(AS.combined_2) <- c("0","2","3","1","4","6","7","5")
new.cluster.ids <- c("0","1","2","3","4","5","6","7")
names(new.cluster.ids) <- levels(AS.combined_2)
AS.combined_2<- RenameIdents(AS.combined_2, new.cluster.ids)

# Save as RDS file (T cell)
saveRDS(AS.combined_2, "~/Desktop/Human_AS/RDSfiles/AS_combined_Tcell_240125.rds")
AS.combined_2<- readRDS ("~/Desktop/Human_AS/RDSfiles/AS_combined_Tcell_240125.rds")

