# Import library----------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

##Data definition
ASFLEX_combined <- readRDS("~/Desktop/Human_AS/RDSfiles/ASFLEX_rpca_all_240519.rds")

DimPlot(ASFLEX_combined, reduction = "umap", label = TRUE, repel = TRUE) 

# Define AS and AR groups
AS_ASFLEX_combined <- subset(ASFLEX_combined, subset = stim %in% c("AS2", "AS2_2", "AS5", "AS11"))
AR_ASFLEX_combined <- subset(ASFLEX_combined, subset = stim %in% c("AR1", "AR2"))

### Perform random sampling while preserving the cell proportions of each cluster
### Create one sampling group for each cluster in your dataset
cluster_0 <- subset(AR_ASFLEX_combined, idents = c("0"))
cluster_1 <- subset(AR_ASFLEX_combined, idents = c("1"))
cluster_2 <- subset(AR_ASFLEX_combined, idents = c("2"))
cluster_3 <- subset(AR_ASFLEX_combined, idents = c("3"))
cluster_4 <- subset(AR_ASFLEX_combined, idents = c("4"))
cluster_5 <- subset(AR_ASFLEX_combined, idents = c("5"))
cluster_6 <- subset(AR_ASFLEX_combined, idents = c("6"))


## Set the total number of cells to sample (10,000 in this example; the actual number may vary slightly).
all_counts <- ncol(AR_ASFLEX_combined)

set.seed(12)
object_cluster_0 <- subset(AR_ASFLEX_combined, cells = Cells(cluster_0), downsample = (2299*ncol(cluster_0)/all_counts))
object_cluster_1 <- subset(AR_ASFLEX_combined, cells = Cells(cluster_1), downsample = (2299*ncol(cluster_1)/all_counts))
object_cluster_2 <- subset(AR_ASFLEX_combined, cells = Cells(cluster_2), downsample = (2299*ncol(cluster_2)/all_counts))
object_cluster_3 <- subset(AR_ASFLEX_combined, cells = Cells(cluster_3), downsample = (2299*ncol(cluster_3)/all_counts))
object_cluster_4 <- subset(AR_ASFLEX_combined, cells = Cells(cluster_4), downsample = (2299*ncol(cluster_4)/all_counts))
object_cluster_5 <- subset(AR_ASFLEX_combined, cells = Cells(cluster_5), downsample = (2299*ncol(cluster_5)/all_counts))
object_cluster_6 <- subset(AR_ASFLEX_combined, cells = Cells(cluster_6), downsample = (2299*ncol(cluster_6)/all_counts))
set.seed(NULL)
object_subset_2 <- merge(object_cluster_0, y = c(object_cluster_1, object_cluster_2,object_cluster_3,object_cluster_4,object_cluster_5,object_cluster_6,AS_ASFLEX_combined))

##Perform downsampling
ASFLEX_combined_downsampled <- subset(ASFLEX_combined, cells = Cells(object_subset_2))

ASFLEX_combined_downsampled <- ScaleData(ASFLEX_combined_downsampled, verbose = FALSE)
ASFLEX_combined_downsampled <- FindVariableFeatures(ASFLEX_combined_downsampled, selection.method = "vst", nfeatures = 2000)


DimPlot(ASFLEX_combined_downsampled, reduction = "umap", label = TRUE, repel = TRUE) 
DimPlot(ASFLEX_combined_downsampled, reduction = "umap", label = TRUE, repel = TRUE, split.by = "ASAR") 
DimPlot(ASFLEX_combined_downsampled, reduction = "umap", label = FALSE, repel = FALSE, split.by = "ASAR") 


#Check the cell count (AR)
ncol(AR_ASFLEX_combined_downsampled)


#Check the cell count ASAR（n=4549)
print(ncol(ASFLEX_combined_downsampled))


# Save as RDS file
saveRDS(ASFLEX_combined_downsampled, "~/Desktop/Human_AS/RDSfiles/ASFLEX_all_downsample_250804.rds")
ASFLEX_combined_downsampled<- readRDS ("~/Desktop/Human_AS/RDSfiles/ASFLEX_all_downsample_250804.rds")


# Dotplot-----
ASFLEX_combined_downsampled_trans <- SCTransform (ASFLEX_combined_downsampled)
cd_genes <- c("COL1A1","COL6A1","SPARC","BGN","ACTA2","MYH11","CD68","CD163","FOLR2","CD3E","IL7R","CDH5","VWF")
DotPlot(ASFLEX_combined_downsampled_trans,features = cd_genes)+RotatedAxis()+coord_flip()


# Heatmap----
ASFLEX_combined_downsampled.markers <- FindAllMarkers(ASFLEX_combined_downsampled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- ASFLEX_combined_downsampled.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(ASFLEX_combined_downsampled, features = top10$gene) 

###Split into ASAR
ASFLEX_combined_downsampled_AS <- subset(ASFLEX_combined_downsampled, subset = ASAR =="AS")
ASFLEX_combined_downsampled_AR <- subset(ASFLEX_combined_downsampled, subset = ASAR =="AR")

###Count cell number of the subcluster
ASFLEX_combined_downsampled_AS_0 <- subset(ASFLEX_combined_downsampled_AS, idents = 0)
ASFLEX_combined_downsampled_AS_1 <- subset(ASFLEX_combined_downsampled_AS, idents = 1)
ASFLEX_combined_downsampled_AS_2 <- subset(ASFLEX_combined_downsampled_AS, idents = 2)
ASFLEX_combined_downsampled_AS_3 <- subset(ASFLEX_combined_downsampled_AS, idents = 3)
ASFLEX_combined_downsampled_AS_4 <- subset(ASFLEX_combined_downsampled_AS, idents = 4)
ASFLEX_combined_downsampled_AS_5 <- subset(ASFLEX_combined_downsampled_AS, idents = 5)
ASFLEX_combined_downsampled_AS_6 <- subset(ASFLEX_combined_downsampled_AS, idents = 6)

###Count cell number of the subcluster
ASFLEX_combined_downsampled_AR_0 <- subset(ASFLEX_combined_downsampled_AR, idents = 0)
ASFLEX_combined_downsampled_AR_1 <- subset(ASFLEX_combined_downsampled_AR, idents = 1)
ASFLEX_combined_downsampled_AR_2 <- subset(ASFLEX_combined_downsampled_AR, idents = 2)
ASFLEX_combined_downsampled_AR_3 <- subset(ASFLEX_combined_downsampled_AR, idents = 3)
ASFLEX_combined_downsampled_AR_4 <- subset(ASFLEX_combined_downsampled_AR, idents = 4)
ASFLEX_combined_downsampled_AR_5 <- subset(ASFLEX_combined_downsampled_AR, idents = 5)
ASFLEX_combined_downsampled_AR_6 <- subset(ASFLEX_combined_downsampled_AR, idents = 6)


#VIC-------
ASFLEX_combined_downsampled_VIC <- subset(ASFLEX_combined_downsampled, idents = c(0,1,2))
ASFLEX_combined_downsampled_VIC <- ScaleData(ASFLEX_combined_downsampled_VIC, verbose = FALSE)
ASFLEX_combined_downsampled_VIC <- FindVariableFeatures(ASFLEX_combined_downsampled_VIC, selection.method = "vst", nfeatures = 2000)
ASFLEX_combined_downsampled_VIC <- RunPCA(ASFLEX_combined_downsampled_VIC, npcs = 30, verbose = FALSE)
ASFLEX_combined_downsampled_VIC <- RunUMAP(ASFLEX_combined_downsampled_VIC , reduction = "pca", dims = 1:15)
ASFLEX_combined_downsampled_VIC <- FindNeighbors(ASFLEX_combined_downsampled_VIC , reduction = "pca", dims = 1:15)
ASFLEX_combined_downsampled_VIC <- FindClusters(ASFLEX_combined_downsampled_VIC, resolution = 0.20)
DimPlot(ASFLEX_combined_downsampled_VIC, reduction = "umap", split.by = "ASAR")

#Reorder----
levels(ASFLEX_combined_downsampled_VIC) 
levels(ASFLEX_combined_downsampled_VIC) <- c("0", "1", "2", "4", "3","5")

new.cluster.ids <- c("0","1","2","3","4","5")
names(new.cluster.ids) <- levels(ASFLEX_combined_downsampled_VIC)
ASFLEX_combined_downsampled_VIC<- RenameIdents(ASFLEX_combined_downsampled_VIC, new.cluster.ids)


# Heatmap----
ASFLEX_combined_downsampled_VIC.markers <- FindAllMarkers(ASFLEX_combined_downsampled_VIC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- ASFLEX_combined_downsampled_VIC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(ASFLEX_combined_downsampled_VIC, features = top10$gene) 

write.table(top10, "table.txt", quote=F, col.names=F, append=T)
DoHeatmap(ASFLEX_combined_downsampled_VIC, features = top10$gene) 


# Dotplot-----
ASFLEX_combined_downsampled_VIC_trans <- SCTransform (ASFLEX_combined_downsampled_VIC)
cd_genes <- c("EGR1","BGN","VIM","COL1A1","COL6A1","COMP","SPARC","ACTA2","MYH11","CSF1R","CD68")
DotPlot(ASFLEX_combined_downsampled_VIC_trans,features = cd_genes)+RotatedAxis()+coord_flip()


VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("EGR1","COL1A1","COL6A1","VIM","SPARC","COMP"))
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("ACTA2","MYH11","TAGLN","CSF1R"))
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("EGR1","COMP","VIM", "SPARC","ACTA2","BMP4"))
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("VIM","DCN","COL1A1","COL1A2","FBN1","ELN"))
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("POSTN","CCL19","HSPG2","FRY","LAMA5","TNC","SULF1","TBX2"))
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("SPARC","COMP","EGR1","ACTA2","MYH11","TAGLN"), split.by = "ASAR")
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("MICA", "MICB", "ULBP1"))
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("MICA", "MICB", "ULBP1"), split.by = "ASAR")
VlnPlot(ASFLEX_combined_downsampled_VIC, assay ="RNA", pt.size = 0, features = c("STAT1","IRF1","CXCL9","CXCL10","CXCL11","HLA-DRA","HLA-DRB1","B2M","TAP1"), split.by = "ASAR")


FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('COMP'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('SPARC'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('COL1A1'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('EGR1'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('BGN'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('COMP'))
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('SPARC'))
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('COL1A1'))
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('EGR1'))
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('BGN'))

FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('TNFRSF1A'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('TNFRSF1B'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('TNFRSF1A'))
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('TNFRSF1B'))


FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('MICA'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('MICB'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('ULBP2'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('ULBP1'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('RAET1E'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('RAET1G'), split.by = "ASAR", min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('ICAM1'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('ICAM2'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('IFNGR1'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('IFNGR2'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('NFKB1'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('TNFRSF1A'), split.by = "ASAR")
FeaturePlot(ASFLEX_combined_downsampled_VIC, features=c('TNFRSF1B'), split.by = "ASAR")

p <- FeaturePlot(
  ASFLEX_combined_downsampled_VIC,
  features = "TNFRSF1A",
  split.by = "ASAR",
  min.cutoff = 0.1,
  max.cutoff = "q90",
  keep.scale = "all"
)

p + theme(
  legend.position = "right",
  legend.direction = "vertical",
)

p <- FeaturePlot(
  ASFLEX_combined_downsampled_VIC,
  features = "TNFRSF1B",
  split.by = "ASAR",
  min.cutoff = 0.1,
  max.cutoff = "q90",
  keep.scale = "all"
)

p + theme(
  legend.position = "right",
  legend.direction = "vertical",
)

p <- FeaturePlot(
  ASFLEX_combined_downsampled_VIC,
  features = "IFNGR1",
  split.by = "ASAR",
  min.cutoff = 0.1,
  max.cutoff = "q90",
  keep.scale = "all"
)

p + theme(
  legend.position = "right",
  legend.direction = "vertical",
)

p <- FeaturePlot(
  ASFLEX_combined_downsampled_VIC,
  features = "IFNGR2",
  split.by = "ASAR",
  min.cutoff = 0.1,
  max.cutoff = "q90",
  keep.scale = "all"
)

p + theme(
  legend.position = "right",
  legend.direction = "vertical",
)


p <- FeaturePlot(
  ASFLEX_combined_downsampled_VIC,
  features = "MICA",
  split.by = "ASAR",
  min.cutoff = 0.1,
  max.cutoff = "q90",
  keep.scale = "all"
)

p + theme(
  legend.position = "right",
  legend.direction = "vertical",
)

p <- FeaturePlot(
  ASFLEX_combined_downsampled_VIC,
  features = "RAET1E",
  split.by = "ASAR",
  min.cutoff = 0.1 ,
  max.cutoff = "q90",
  keep.scale = "all"
)

p + theme(
  legend.position = "right",
  legend.direction = "vertical",
)

# volcano VIC------
ASFLEX_combined_downsampled$celltype <- Idents(ASFLEX_combined_downsampled)
ASFLEX_combined_downsampled$celltype.stim <- paste(Idents(ASFLEX_combined_downsampled), ASFLEX_combined_downsampled$ASAR, sep="_")
Idents(ASFLEX_combined_downsampled) <- "celltype.stim" 
levels(ASFLEX_combined_downsampled) 

ASFLEX_combined_downsampled.table <- FindMarkers(ASFLEX_combined_downsampled, ident.1 = c("0_AS","1_AS","2_AS"), 
                                     ident.2 =c("0_AR","1_AR","2_AR"), verbose = FALSE, logfc.threshold = 0)
ASFLEX_combined_downsampled.table$logp <- -log10(ASFLEX_combined_downsampled.table$p_val)
ASFLEX_combined_downsampled_filtered_left = subset(ASFLEX_combined_downsampled.table, logp>=2 & avg_log2FC <= -0.8)
ASFLEX_combined_downsampled_filtered_right = subset(ASFLEX_combined_downsampled.table, logp>=2 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(ASFLEX_combined_downsampled_filtered_left)
genes.to.label.right <- rownames(ASFLEX_combined_downsampled_filtered_right)

p1 <- ggplot(ASFLEX_combined_downsampled.table, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

###Split into ASAR
ASFLEX_combined_downsampled_VIC_AS <- subset(ASFLEX_combined_downsampled_VIC, subset = ASAR =="AS")
ASFLEX_combined_downsampled_VIC_AR <- subset(ASFLEX_combined_downsampled_VIC, subset = ASAR =="AR")

###Count cell number of the subcluster
ASFLEX_combined_downsampled_VIC_AS_0 <- subset(ASFLEX_combined_downsampled_VIC_AS, idents = 0)
ASFLEX_combined_downsampled_VIC_AS_1 <- subset(ASFLEX_combined_downsampled_VIC_AS, idents = 1)
ASFLEX_combined_downsampled_VIC_AS_2 <- subset(ASFLEX_combined_downsampled_VIC_AS, idents = 2)
ASFLEX_combined_downsampled_VIC_AS_3 <- subset(ASFLEX_combined_downsampled_VIC_AS, idents = 3)
ASFLEX_combined_downsampled_VIC_AS_4 <- subset(ASFLEX_combined_downsampled_VIC_AS, idents = 4)
ASFLEX_combined_downsampled_VIC_AS_5 <- subset(ASFLEX_combined_downsampled_VIC_AS, idents = 5)

###Count cell number of the subcluster
ASFLEX_combined_downsampled_VIC_AR_0 <- subset(ASFLEX_combined_downsampled_VIC_AR, idents = 0)
ASFLEX_combined_downsampled_VIC_AR_1 <- subset(ASFLEX_combined_downsampled_VIC_AR, idents = 1)
ASFLEX_combined_downsampled_VIC_AR_2 <- subset(ASFLEX_combined_downsampled_VIC_AR, idents = 2)
ASFLEX_combined_downsampled_VIC_AR_3 <- subset(ASFLEX_combined_downsampled_VIC_AR, idents = 3)
ASFLEX_combined_downsampled_VIC_AR_4 <- subset(ASFLEX_combined_downsampled_VIC_AR, idents = 4)
ASFLEX_combined_downsampled_VIC_AR_5 <- subset(ASFLEX_combined_downsampled_VIC_AR, idents = 5)

# Save as RDS file (VIC)
saveRDS(ASFLEX_combined_downsampled_VIC, "~/Desktop/Human_AS/RDSfiles/ASFLEX_vic_downsampling_240804.rds")
ASFLEX_combined_downsampled_VIC<- readRDS ("~/Desktop/Human_AS/RDSfiles/ASFLEX_vic_downsampling_240804.rds")

#EC-------
ASFLEX_combined_downsampled_EC <- subset(ASFLEX_combined_downsampled, idents = c(6))
ASFLEX_combined_downsampled_EC <- ScaleData(ASFLEX_combined_downsampled_EC, vedsbose = FALSE)
ASFLEX_combined_downsampled_EC<- RunPCA(ASFLEX_combined_downsampled_EC, npcs = 30, verbose = FALSE)
ASFLEX_combined_downsampled_EC <- RunUMAP(ASFLEX_combined_downsampled_EC , reduction = "pca", dims = 1:15)
ASFLEX_combined_downsampled_EC <- FindNeighbors(ASFLEX_combined_downsampled_EC, reduction = "pca", dims = 1:15)
ASFLEX_combined_downsampled_EC <- FindClusters(ASFLEX_combined_downsampled_EC, resolution = 0.15)
DimPlot(ASFLEX_combined_downsampled_EC, reduction = "umap", split.by = "ASAR")

VlnPlot(ASFLEX_combined_downsampled_EC, pt.size = 0, features = c("ICAM1","VCAM1","TNFRSF9","BMP2","IL1R1","TNFRSF1A","CD68","COL1A1","RUNX2","IBSP","BGLAP","BMP4"))
VlnPlot(ASFLEX_combined_downsampled_EC, assay ="RNA", pt.size = 0, features = c("STAT1","IRF1","CXCL9","CXCL10","CXCL11","HLA-DRA","HLA-DRB1","B2M","TAP1"), split.by = "ASAR")


ASFLEX_downsampled_combined.table <- FindMarkers(ASFLEX_downsampled_combined, ident.1 = c("6_AS"), 
                                     ident.2 =c("6_AR"), verbose = FALSE, logfc.threshold = 0)
ASFLEX_downsampled_combined.table$logp <- -log10(ASFLEX_downsampled_combined.table$p_val)
ASFLEX_downsampled_combined_filtered_left = subset(ASFLEX_downsampled_combined.table, logp>=5 & avg_log2FC <= -1)
ASFLEX_downsampled_combined_filtered_right = subset(ASFLEX_downsampled_combined.table, logp>=5 & avg_log2FC >= 1)

genes.to.label.left <- rownames(ASFLEX_downsampled_combined_filtered_left)
genes.to.label.right <- rownames(ASFLEX_downsampled_combined_filtered_right)

p1 <- ggplot(ASFLEX_downsampled_combined.table, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1
