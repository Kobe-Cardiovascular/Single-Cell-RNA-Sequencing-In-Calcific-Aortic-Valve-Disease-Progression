#Import library------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#T cell RDS file 
saveRDS(AS.combined_2, "~/Desktop/Human_AS/RDSfiles/AS_combined_Tcell_240125.rds")
AS.combined_2<- readRDS ("~/Desktop/Human_AS/RDSfiles/AS_combined_Tcell_240125.rds")

DimPlot(AS.combined_2, reduction = "umap", label = TRUE,  repel = TRUE)


##Split into AS and AR
AS.combined_2_AS <- subset(AS.combined_2, subset = ASAR =="AS")
AS.combined_2_AR <- subset(AS.combined_2, subset = ASAR =="AR")

###Count the cell number of subcluster
AS.combined_2_AS_0 <- subset(AS.combined_2_AS, idents = 0)
AS.combined_2_AS_1 <- subset(AS.combined_2_AS, idents = 1)
AS.combined_2_AS_2 <- subset(AS.combined_2_AS, idents = 2)
AS.combined_2_AS_3 <- subset(AS.combined_2_AS, idents = 3)
AS.combined_2_AS_4 <- subset(AS.combined_2_AS, idents = 4)
AS.combined_2_AS_5 <- subset(AS.combined_2_AS, idents = 5)

AS.combined_2_AR_0 <- subset(AS.combined_2_AR, idents = 0)
AS.combined_2_AR_1 <- subset(AS.combined_2_AR, idents = 1)
AS.combined_2_AR_2 <- subset(AS.combined_2_AR, idents = 2)
AS.combined_2_AR_3 <- subset(AS.combined_2_AR, idents = 3)
AS.combined_2_AR_4 <- subset(AS.combined_2_AR, idents = 4)
AS.combined_2_AR_5 <- subset(AS.combined_2_AR, idents = 5)

#Proportion
AS.combined_2_prop. <- prop.table(table(Idents(AS.combined_2),AS.combined_2$ASAR))
write.table(AS.combined_2_prop., "~/Desktop/AS.combined_2_prop_20250226.txt", quote=F, col.names=F, append=T)
write.csv(AS.combined_2_prop., "AS.combined_2_20250226_prop.csv")

#TNF NFkB1 AS vs. AR
VlnPlot(AS.combined_2,pt.size = 0.1,assay ="RNA", features = c("TNF"),split.by = "ASAR")
VlnPlot(AS.combined_2,pt.size = 0.1,assay ="RNA", features = c("NFKB1"),split.by = "ASAR")

#Heatmap
AS.combined_2.markers <- FindAllMarkers(AS.combined_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AS.combined_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AS.combined_2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AS.combined_2, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


AS.combined_2_trans <- SCTransform (AS.combined_2)

cd_genes <- c("CD3E","CD4","CD8A","IL7R","MKI67","NKG7","GNLY","GZMB","TYROBP","FCGR3A","XCL1","XCL2","TRDV2","TRGV9","GATA3")
DotPlot(AS.combined_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()


#CD4----
AS.combined_2_CD4<- subset(AS.combined_2, idents = c(0,1))
AS.combined_2_CD4<- ScaleData(AS.combined_2_CD4, verbose = FALSE)
AS.combined_2_CD4<- RunPCA(AS.combined_2_CD4, npcs = 30, verbose = FALSE)
AS.combined_2_CD4<- RunUMAP(AS.combined_2_CD4, reduction = "pca", dims = 1:10)
AS.combined_2_CD4<- FindNeighbors(AS.combined_2_CD4, reduction = "pca", dims = 1:10)
AS.combined_2_CD4<- FindClusters(AS.combined_2_CD4, resolution = 0.32)

DimPlot(AS.combined_2_CD4, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD4, reduction = "umap", split.by = "ASAR", label = TRUE)
DimPlot(AS.combined_2_CD4, reduction = "umap", split.by = "ASAR", label = FALSE)

levels(AS.combined_2_CD4) <- c("0","1","2","3","4","5","7","6")
new.cluster.ids <- c("0","1","2","3","4","5","6","7")
names(new.cluster.ids) <- levels(AS.combined_2_CD4)
AS.combined_2_CD4<- RenameIdents(AS.combined_2_CD4, new.cluster.ids)


VlnPlot(AS.combined_2_CD4,pt.size = 0, assay = 'RNA', features = c("GZMB","GZMK","CD69","ITGB2","CD8A","PLAC8"))
VlnPlot(AS.combined_2_CD4,pt.size = 0, features = c("IL7R","CD28","GZMB","PRF1","CCR7","SELL"))
VlnPlot(AS.combined_2_CD4,pt.size = 0, features = c("LEF1","FOXP3","CTLA4","IFNG","GATA3","SELL"))
VlnPlot(AS.combined_2_CD4,pt.size = 0, assay = 'RNA', features = c("PDCD1","CXCR5","CD28","CCR7"))


FeaturePlot(AS.combined_2_CD4, features=c('CD8A'))
FeaturePlot(AS.combined_2_CD4, features=c('GZMK'))
FeaturePlot(AS.combined_2_CD4, features=c('KLRD1'))
FeaturePlot(AS.combined_2_CD4, features=c('CD8A'), split.by = "ASAR")

FeaturePlot(AS.combined_2_CD4, features=c('KLRK1'), split.by = "ASAR")
FeaturePlot(AS.combined_2_CD4, features=c('FASLG'), split.by = "ASAR")


Naive
FeaturePlot(AS.combined_2_CD4, features=c('CCR7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('LEF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('SELL'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('TCF7'), min.cutoff=0.5, max.cutoff='q90')
Cytotoxicity
FeaturePlot(AS.combined_2_CD4, features=c('GZMA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('GZMB'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('GZMH'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('GZMK'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('GNLY'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('PRF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('NKG7'), min.cutoff=0.5, max.cutoff='q90')
TRM
FeaturePlot(AS.combined_2_CD4, features=c('CD69'), min.cutoff=3.0, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CXCR6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('ITGA1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('ZNF683'), min.cutoff=0.5, max.cutoff='q90')
Inhibitors
FeaturePlot(AS.combined_2_CD4, features=c('CTLA4'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('HAVCR2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('FOXP3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('LAG3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('PDCD1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('TIGIT'), min.cutoff=0.5, max.cutoff='q90')
MAIT
FeaturePlot(AS.combined_2_CD4, features=c('TRAV1-2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('SLC4A10'), min.cutoff=0.5, max.cutoff='q90')
LncRNA
FeaturePlot(AS.combined_2_CD4, features=c('MALAT1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('NEAT1'), min.cutoff=0.5, max.cutoff='q90')

#stress
#FOS, JUN, HSP
FeaturePlot(AS.combined_2_CD4, features=c('FOS'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('JUN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('HSPA1B'), min.cutoff=0.5, max.cutoff='q90')
remodeling
#DCN, COL1A2, CLU
FeaturePlot(AS.combined_2_CD4, features=c('DCN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('COL1A2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CLU'), min.cutoff=0.5, max.cutoff='q90')
#regulatory
#NFKB1, NFKBID
FeaturePlot(AS.combined_2_CD4, features=c('NFKB1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('NFKBID'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AS.combined_2_CD4, features=c('KLRG1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('ZEB2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('TNF'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CCL5'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CD69'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CD74'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CD55'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('SELL'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('IKZF2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CCL5'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('ITGB2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('GIMAP7'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AS.combined_2_CD4, features=c('SPON2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('FCGR3A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CX3CR1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD4, features=c('CD28'), min.cutoff=0.5, max.cutoff='q90')


# Save RDS file (CD4)
saveRDS(AS.combined_2_CD4, "~/Desktop/Human_AS/RDSfiles/AS_combined_CD4Tcell_251218_male.rds")
AS.combined_2_CD4<- readRDS ("~/Desktop/Human_AS/RDSfiles/AS_combined_CD4Tcell_251218_male.rds")


# Dotplot_CD4
AS.combined_2_CD4_trans <- SCTransform (AS.combined_2_CD4)

cd_genes <- c("GZMA","CD69","TNFAIP3","HSPH1","CYBA","MALAT1","CCR7","SELL","FOXP3","CTLA4","LEF1","PRF1","KLRB1","CCL4","IFNG")
DotPlot(AS.combined_2_CD4_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD69","CD74","CD28","CCL5","GZMA","GZMK","CXCR6","PDCD1","FOXP3","CTLA4","TNF","FOS","TNFAIP3","HSPA1A","IL7R","SELL","CCR7","LEF1","CXCR3","ZEB2","KLRG1","GZMB","PRF1")
DotPlot(AS.combined_2_CD4_trans,features = cd_genes)+RotatedAxis()+coord_flip()

# TNF, NFkB AS vs. AR
VlnPlot(AS.combined_2_CD4,pt.size = 0.1,assay ="RNA", features = c("TNF"),split.by = "ASAR")
VlnPlot(AS.combined_2_CD4,pt.size = 0.1,assay ="RNA", features = c("NFKB1"),split.by = "ASAR")

# Merge clusters based on the previous clustering results
new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")
new.cluster.ids <- c("0","0","1","2","3","4","0","0","5","6","7")
names(new.cluster.ids) <- levels(AS.combined_2_CD4)
AS.combined_2_CD4<- RenameIdents(AS.combined_2_CD4, new.cluster.ids)

DimPlot(AS.combined_2_CD4, reduction = "umap", label = TRUE,  repel = TRUE)

###Heatmap
AS.combined_2_CD4.markers <- FindAllMarkers(AS.combined_2_CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- AS.combined_2_CD4.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AS.combined_2_CD4, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#Reorder
levels(AS.combined_2_CD4) # [1] "0""1""2""3""4""5""6""7"
levels(AS.combined_2_CD4) <- c("2","5","1","0","4","3","6","7")
new.cluster.ids <- c("0","1","2","3","4","5","6","7")
names(new.cluster.ids) <- levels(AS.combined_2_CD4)
AS.combined_2_CD4<- RenameIdents(AS.combined_2_CD4, new.cluster.ids)

#After reordering_AS.combined_2_CD4
saveRDS(AS.combined_2_CD4, "~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD4_並び替え後.rds")
AS.combined_2_CD4<- readRDS ("~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD4_並び替え後.rds")

###Split into AS and AR
AS.combined_2_CD4_AS <- subset(AS.combined_2_CD4, subset = ASAR =="AS")
AS.combined_2_CD4_AR <- subset(AS.combined_2_CD4, subset = ASAR =="AR")

DimPlot(AS.combined_2_CD4, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD4_AS, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD4_AR, reduction = "umap", label = TRUE,  repel = TRUE)


#Exclude 6, 7
AS.combined_2_CD4_2<- subset(AS.combined_2_CD4, idents = c(0,1,2,3,4,5))
AS.combined_2_CD4_2<- ScaleData(AS.combined_2_CD4_2, verbose = FALSE)
AS.combined_2_CD4_2<- RunPCA(AS.combined_2_CD4_2, npcs = 30, verbose = FALSE)
AS.combined_2_CD4_2<- RunUMAP(AS.combined_2_CD4_2, reduction = "pca", dims = 1:8)
AS.combined_2_CD4_2<- FindNeighbors(AS.combined_2_CD4_2, reduction = "pca", dims = 1:8)
AS.combined_2_CD4_2<- FindClusters(AS.combined_2_CD4_2, resolution = 0.6)

DimPlot(AS.combined_2_CD4_2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD4_2, reduction = "umap", split.by = "ASAR", label = TRUE)

# Merge clusters based on the previous clustering results
new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")
new.cluster.ids <- c("0","1","2","2","3","4","3","3","3","5","2")
names(new.cluster.ids) <- levels(AS.combined_2_CD4_2)
AS.combined_2_CD4_2<- RenameIdents(AS.combined_2_CD4_2, new.cluster.ids)

#Afer merging Save as RDS 
saveRDS(AS.combined_2_CD4_2, "~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD4_2まとめ後.rds")
AS.combined_2_CD4_2<- readRDS ("~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD4_2まとめ後.rds")

#Reorder
levels(AS.combined_2_CD4_2) # [1] "0""1""2""3""4""5"
levels(AS.combined_2_CD4_2) <- c("2","1","3","0","5","4")
new.cluster.ids <- c("0","1","2","3","4","5")
names(new.cluster.ids) <- levels(AS.combined_2_CD4_2)
AS.combined_2_CD4_2<- RenameIdents(AS.combined_2_CD4_2, new.cluster.ids)

###After reordering Save as RDS file (AS.combined_2_CD4_2) This is the final version!!###
saveRDS(AS.combined_2_CD4_2, "~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD4_2並び替え後.rds")
AS.combined_2_CD4_2<- readRDS ("~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD4_2並び替え後.rds")

# Dotplot_CD4_2
AS.combined_2_CD4_2_trans <- SCTransform (AS.combined_2_CD4_2)

cd_genes <- c("CD3E","CD69","CD28","SELL","LEF1","CXCR6","PDCD1","NFKB1","PRF1","NKG7","GZMK","GZMA","TNF","FOS","HSPA1B","FOXP3","CTLA4")
DotPlot(AS.combined_2_CD4_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD69","CD74","CD28","SELL","CCR7","LEF1","FOXP3","IKZF2","CTLA4","CXCR6","PDCD1","NFKB1","TNF","FOS","HSPA1B","PRF1","NKG7","GZMK","GZMA")
DotPlot(AS.combined_2_CD4_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

# TNF, NFkB, IFNg AS vs. AR
VlnPlot(AS.combined_2_CD4_2,pt.size = 0.05,assay ="RNA", features = c("TNF"),split.by = "ASAR", cols = c("AR" = "blue", "AS" = "red")) 
VlnPlot(AS.combined_2_CD4_2,pt.size = 0.05,assay ="RNA", features = c("NFKB1"),split.by = "ASAR",cols = c("AR" = "blue", "AS" = "red"))
VlnPlot(AS.combined_2_CD4_2,pt.size = 0.05,assay ="RNA", features = c("IFNG"),split.by = "ASAR",cols = c("AR" = "blue", "AS" = "red"))

###Split into AS and AR
AS.combined_2_CD4_2_AS <- subset(AS.combined_2_CD4_2, subset = ASAR =="AS")
AS.combined_2_CD4_2_AR <- subset(AS.combined_2_CD4_2, subset = ASAR =="AR")

DimPlot(AS.combined_2_CD4_2_AS, reduction = "umap", label = FALSE,  repel = TRUE)
DimPlot(AS.combined_2_CD4_2_AR, reduction = "umap", label = FALSE,  repel = TRUE)
DimPlot(AS.combined_2_CD4_2_AS, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD4_2_AR, reduction = "umap", label = TRUE,  repel = TRUE)

###Count cell number of the subcluster (AS)
AS.combined_2_CD4_2_AS_0 <- subset(AS.combined_2_CD4_2_AS, idents = 0)
AS.combined_2_CD4_2_AS_1 <- subset(AS.combined_2_CD4_2_AS, idents = 1)
AS.combined_2_CD4_2_AS_2 <- subset(AS.combined_2_CD4_2_AS, idents = 2)
AS.combined_2_CD4_2_AS_3 <- subset(AS.combined_2_CD4_2_AS, idents = 3)
AS.combined_2_CD4_2_AS_4 <- subset(AS.combined_2_CD4_2_AS, idents = 4)
AS.combined_2_CD4_2_AS_5 <- subset(AS.combined_2_CD4_2_AS, idents = 5)

###Count cell number of the subcluster (AR)
AS.combined_2_CD4_2_AR_0 <- subset(AS.combined_2_CD4_2_AR, idents = 0)
AS.combined_2_CD4_2_AR_1 <- subset(AS.combined_2_CD4_2_AR, idents = 1)
AS.combined_2_CD4_2_AR_2 <- subset(AS.combined_2_CD4_2_AR, idents = 2)
AS.combined_2_CD4_2_AR_3 <- subset(AS.combined_2_CD4_2_AR, idents = 3)
AS.combined_2_CD4_2_AR_4 <- subset(AS.combined_2_CD4_2_AR, idents = 4)
AS.combined_2_CD4_2_AR_5 <- subset(AS.combined_2_CD4_2_AR, idents = 5)


# volcano  CD4------
AS.combined_2_CD4_2$celltype <- Idents(AS.combined_2_CD4_2)
AS.combined_2_CD4_2$celltype.stim <- paste(Idents(AS.combined_2_CD4_2), AS.combined_2_CD4_2$ASAR, sep="_")
Idents(AS.combined_2_CD4_2) <- "celltype.stim" # Update the identities using the Idents() function
levels(AS.combined_2_CD4_2) # Check that the identities have been updated


AS.combined_2_CD4_2.table <- FindMarkers(AS.combined_2_CD4_2, ident.1 = c("0_AS","5_AS", "1_AS", "2_AS", "3_AS", "4_AS"), 
                                     ident.2 =c("1_AR", "5_AR", "3_AR", "2_AR", "0_AR", "4_AR"), verbose = FALSE, logfc.threshold = 0)
AS.combined_2_CD4_2.table$logp <- -log10(AS.combined_2_CD4_2.table$p_val)
AS.combined_2_CD4_2_filtered_left = subset(AS.combined_2_CD4_2.table, logp>=1.0 & avg_log2FC <= -0.5)
AS.combined_2_CD4_2_filtered_right = subset(AS.combined_2_CD4_2.table, logp>=1.0 & avg_log2FC >= 0.5)

genes.to.label.left <- rownames(AS.combined_2_CD4_2_filtered_left)
genes.to.label.right <- rownames(AS.combined_2_CD4_2_filtered_right)

p1 <- ggplot(AS.combined_2_CD4_2.table, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1


#CD8----
AS.combined_2_CD8<- subset(AS.combined_2, idents = c(3,4,5,6,7))
AS.combined_2_CD8<- ScaleData(AS.combined_2_CD8, verbose = FALSE)
AS.combined_2_CD8<- RunPCA(AS.combined_2_CD8, npcs = 30, verbose = FALSE)
AS.combined_2_CD8<- RunUMAP(AS.combined_2_CD8, reduction = "pca", dims = 1:8)
AS.combined_2_CD8<- FindNeighbors(AS.combined_2_CD8, reduction = "pca", dims = 1:8)
AS.combined_2_CD8<- FindClusters(AS.combined_2_CD8, resolution = 2.0)

DimPlot(AS.combined_2_CD8, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD8, reduction = "umap", split.by = "ASAR", label = TRUE)

FeaturePlot(AS.combined_2_CD8, features=c('CD8A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('F2R'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('KLRG1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('CD55'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('CORO1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('OAS3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('LAG3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('CCR6'), min.cutoff=0.5, max.cutoff='q90')

#Naive
FeaturePlot(AS.combined_2_CD8, features=c('CCR7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('LEF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('SELL'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('TCF7'), min.cutoff=0.5, max.cutoff='q90')
#Cytotoxicity
FeaturePlot(AS.combined_2_CD8, features=c('GZMA'), min.cutoff=1.0, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('GZMB'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('GZMH'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('GZMK'), min.cutoff=1.0, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('GNLY'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('PRF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('NKG7'), min.cutoff=0.5, max.cutoff='q90')
#MAIT
FeaturePlot(AS.combined_2_CD8, features=c('TRAV1-2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('SLC4A10'), min.cutoff=0.5, max.cutoff='q90')
#LncRNA
FeaturePlot(AS.combined_2_CD8, features=c('MALAT1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('NEAT1'), min.cutoff=0.5, max.cutoff='q90')
#Proliferation
FeaturePlot(AS.combined_2_CD8, features=c('MKI67'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('PCNA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('STMN1'), min.cutoff=0.5, max.cutoff='q90')
#TRM
FeaturePlot(AS.combined_2_CD8, features=c('CD69'), min.cutoff=3.0, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('CXCR6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('ITGA1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('ZNF683'), min.cutoff=0.5, max.cutoff='q90')
Inhibitors
FeaturePlot(AS.combined_2_CD8, features=c('CTLA4'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('HAVCR2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('LAG3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('PDCD1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('TIGIT'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AS.combined_2_CD8, features=c('PDCD1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('JAML'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('IL7R'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('NEAT1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('MALAT1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('MKI67'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('HSPA1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('HSPA1A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('DNAJB1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('DUSP1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('CTSW'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('JUN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('FOS'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('KLRG1'), min.cutoff=0.5, max.cutoff='q90')


#Naive
#LEF1, CCR7, AQP3, SELL
FeaturePlot(AS.combined_2_CD8, features=c('AQP3'), min.cutoff=0.5, max.cutoff='q90')
#Exhausion
#HOPX, LAG3, TOX, TIGHT, PDCD1
FeaturePlot(AS.combined_2_CD8, features=c('HOPX'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('TOX'), min.cutoff=0.5, max.cutoff='q90')
#Cytotoxic
#GZMA, GZMH, GNLY, PRF1, NKG7
#Activated
#CD247, CD69, IL2RA
FeaturePlot(AS.combined_2_CD8, features=c('CD247'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('IL2RA'), min.cutoff=0.5, max.cutoff='q90')
#MAIT
#KLRB1, CCR6, CXCR6
FeaturePlot(AS.combined_2_CD8, features=c('KLRB1'), min.cutoff=0.5, max.cutoff='q90')

#stress
#FOS, JUN, HSP
FeaturePlot(AS.combined_2_CD8, features=c('FOS'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('JUN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('HSPA1B'), min.cutoff=0.5, max.cutoff='q90')
#remodeling
#DCN, COL1A2, CLU
FeaturePlot(AS.combined_2_CD8, features=c('DCN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('COL1A2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('CLU'), min.cutoff=0.5, max.cutoff='q90')
#regulatory
#NFKB1, NFKBID
FeaturePlot(AS.combined_2_CD8, features=c('NFKB1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8, features=c('NFKBID'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(AS.combined_2_CD8, pt.size = 0, assay = 'RNA', features = c("GZMB","HSPA1B","FOS","GNLY") )
VlnPlot(AS.combined_2_CD8,pt.size = 0, features = c("CD8A","GZMA","GZMK","CXCR6","CD27"))
VlnPlot(AS.combined_2_CD8,pt.size = 0, features = c("GZMA","GZMK","CXCR6","CD27","TNF","CD74","IFNG","PRF1","GZMB","TBX21","CX3CR1","NKG7","SELL","CCR7","CD69"))
#VlnPlot(AS.combined_2_CD8,pt.size = 0, features = c("CELL","IL7R","LEF1","C1QA","TNF","CXCL8","LYVE1","TET2","CD163","FOXP3","CTLA4") )
#VlnPlot(AS.combined_2_CD8,pt.size = 0, features = c("GZMA","GZMK","CD69","CCL4","TNFAIP3","CCL5","IL7R","CD28","GZMB","PRF1","CCR7","SELL","LEF1","FOXP3","CTLA4","IFNG","IL4","GATA3"))

# Myeloid----
VlnPlot(AS.combined_2_CD8, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(AS.combined_2_CD8, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","TNF","CXCL8") )
VlnPlot(AS.combined_2_CD8, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AS.combined_2_CD8, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
# Bcells-----
VlnPlot(AS.combined_2_CD8, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(AS.combined_2_CD8, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )
cDC1
VlnPlot(AS.combined_2_CD8, pt.size = 0 , features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(AS.combined_2_CD8, pt.size = 0 , features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(AS.combined_2_CD8, pt.size = 0 , features = c("IL3RA","CLEC4C") )

#dotplot
AS.combined_2_CD8_trans <- SCTransform (AS.combined_2_CD8)

cd_genes <- c("CD3E","CD8A","CD28","CD69","CD74","CCL5","GZMA","GZMK","PDCD1","CXCR6","JAML","IL7R","NEAT1","MALAT1","CCR7","TNF","GZMB","NKG7","PRF1","CX3CR1","FCGR3A","EOMES","KLRG1","MKI67")
DotPlot(AS.combined_2_CD8_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD8A","CD28","CD69","CD74","CCL5","GZMA","GZMK","PDCD1","CXCR6","IL7R","CCR7","TNF","GZMB","PRF1","EOMES","KLRG1","MKI67")
DotPlot(AS.combined_2_CD8_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD8A","CD28","CD69","GZMB","GNLY","CTSW","FOS","JUN","HSPA1B","IL7R","AREG","PDCD1","CREM","KLRB1","MKI67")
DotPlot(AS.combined_2_CD8_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CTLA4","LAG3","PDCD1","CD69","CXCR6","MKI67","STMN1","GZMK","GZMB","GNLY","PRF1","NKG7","CCR7","LEF1","TCF7")
DotPlot(AS.combined_2_CD8_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD8A","CD28","CD69","PDCD1","CXCR6","CD55","KLRB1","NFKB1","IL7R","AREG","FOS","JUN","HSPA1B","DUSP1","GIMAP7","GZMA","GZMB","PRF1","MKI67","STMN1")
DotPlot(AS.combined_2_CD8_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD8A","CD28","CD69","PDCD1","CXCR6","CD55","NFKB1","IL7R","AREG","GZMB","PRF1","FOS","JUN","HSPA1B","DUSP1","GZMA","GZMK","KLRG1","MKI67","STMN1")
DotPlot(AS.combined_2_CD8_trans,features = cd_genes)+RotatedAxis()+coord_flip()


###Heatmap
AS.combined_2_CD8.markers <- FindAllMarkers(AS.combined_2_CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- AS.combined_2_CD8.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AS.combined_2_CD8, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

FeaturePlot(AS.combined_2_CD8, split.by = "ASAR", features=c('CXCR6'), min.cutoff=0.5, max.cutoff='q90')

#Merge clusters based on the previous clustering results
new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
new.cluster.ids <- c("0","1","1","1","2","1","3","3","1","3","1","0","1","1","1","2","2","3","2","4","5","0","6","7","8","1","3","5")
names(new.cluster.ids) <- levels(AS.combined_2_CD8)
AS.combined_2_CD8<- RenameIdents(AS.combined_2_CD8, new.cluster.ids)

#Reordering
levels(AS.combined_2_CD8) # [1] "0""1""2""3""4""5""6""7""8"
levels(AS.combined_2_CD8) <- c("0","3","1","2","5","4","6","7","8")
new.cluster.ids <- c("0","1","2","3","4","5","6","7","8")
names(new.cluster.ids) <- levels(AS.combined_2_CD8)
AS.combined_2_CD8<- RenameIdents(AS.combined_2_CD8, new.cluster.ids)

DimPlot(AS.combined_2_CD8, reduction = "umap", label = TRUE,  repel = TRUE)
FeaturePlot(AS.combined_2_CD8, features=c('STMN1'), min.cutoff=0.5, max.cutoff='q90')

#After reordering
saveRDS(AS.combined_2_CD8, "~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD8_並び替え後.rds")
AS.combined_2_CD8<- readRDS ("~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD8_並び替え後.rds")

#CD8_2----Exclude CD8 5,6,7,8
AS.combined_2_CD8_2<- subset(AS.combined_2_CD8, idents = c(0,1,2,3,4))
AS.combined_2_CD8_2<- ScaleData(AS.combined_2_CD8_2, verbose = FALSE)
AS.combined_2_CD8_2<- RunPCA(AS.combined_2_CD8_2, npcs = 30, verbose = FALSE)
AS.combined_2_CD8_2<- RunUMAP(AS.combined_2_CD8_2, reduction = "pca", dims = 1:7)
AS.combined_2_CD8_2<- FindNeighbors(AS.combined_2_CD8_2, reduction = "pca", dims = 1:7)
AS.combined_2_CD8_2<- FindClusters(AS.combined_2_CD8_2, resolution = 0.1)

DimPlot(AS.combined_2_CD8_2, reduction = "umap", label = TRUE,  repel = TRUE)

Naive
FeaturePlot(AS.combined_2_CD8_2, features=c('CCR7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('LEF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('SELL'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('TCF7'), min.cutoff=0.5, max.cutoff='q90')
Cytotoxicity
FeaturePlot(AS.combined_2_CD8_2, features=c('GZMA'), min.cutoff=1.0, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('GZMB'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('GZMH'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('GZMK'), min.cutoff=1.0, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('GNLY'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('PRF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('NKG7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('MKI67'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_2, features=c('HSPA1B'), min.cutoff=0.5, max.cutoff='q90')

#CD8_3----Exclude CD8_2(pca7, reso0.1) cluster4
AS.combined_2_CD8_3<- subset(AS.combined_2_CD8_2, idents = c(0,1,2,3))
AS.combined_2_CD8_3<- ScaleData(AS.combined_2_CD8_3, verbose = FALSE)
AS.combined_2_CD8_3<- RunPCA(AS.combined_2_CD8_3, npcs = 30, verbose = FALSE)
AS.combined_2_CD8_3<- RunUMAP(AS.combined_2_CD8_3, reduction = "pca", dims = 1:6)
AS.combined_2_CD8_3<- FindNeighbors(AS.combined_2_CD8_3, reduction = "pca", dims = 1:6)
AS.combined_2_CD8_3<- FindClusters(AS.combined_2_CD8_3, resolution = 1.4)

DimPlot(AS.combined_2_CD8_3, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD8_3, reduction = "umap", split.by = "ASAR", label = TRUE)

#Merge clusters based on the previous clustering results
new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")
new.cluster.ids <- c("0","0","1","2","3","0","1","3","1","0","0","2","2","3","3","4","0","3","4")
names(new.cluster.ids) <- levels(AS.combined_2_CD8_3)
AS.combined_2_CD8_3<- RenameIdents(AS.combined_2_CD8_3, new.cluster.ids)

#Reordering
levels(AS.combined_2_CD8_3) # [1] "0""1""2""3""4"
levels(AS.combined_2_CD8_3) <- c("3","1","0","2","4")
new.cluster.ids <- c("0","1","2","3","4")
names(new.cluster.ids) <- levels(AS.combined_2_CD8_3)
AS.combined_2_CD8_3<- RenameIdents(AS.combined_2_CD8_3, new.cluster.ids)

DimPlot(AS.combined_2_CD8_3, reduction = "umap", split.by = "ASAR", label = TRUE)

FeaturePlot(AS.combined_2_CD8_3, features=c('CD8A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('STMN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('MKI67'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('GZMB'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('PRF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('NKG7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('FOS'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('JUN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('HSPA1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('KLRG1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('NFKB1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('CXCR6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('GZMK'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('GZMA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('LAG3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('TOX'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('TIGHT'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('TNF'), min.cutoff=0.5, max.cutoff='q90')
TRM
FeaturePlot(AS.combined_2_CD8_3, features=c('CD69'), min.cutoff=3.0, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('CXCR6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('ITGA1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('ZNF683'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('PDCD1'), min.cutoff=0.5, max.cutoff='q90')
Naive
FeaturePlot(AS.combined_2_CD8_3, features=c('CCR7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('LEF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('SELL'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AS.combined_2_CD8_3, features=c('TCF7'), min.cutoff=0.5, max.cutoff='q90')

#dotplot
AS.combined_2_CD8_3_trans <- SCTransform (AS.combined_2_CD8_3)

cd_genes <- c("CD3E","CD28","CD69","PDCD1","CXCR6","IL7R","GZMB","PRF1","FOS","JUN","HSPA1B","GZMK","GZMA","NKG7","GIMAP7","KLRG1","MKI67","STMN1")
DotPlot(AS.combined_2_CD8_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD28","CD69","TNF","JUN","HSPA1B","PDCD1","CXCR6","IL7R","GZMB","PRF1","GZMK","GZMA","GIMAP7","KLRG1","MKI67","STMN1")
DotPlot(AS.combined_2_CD8_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD8A","CD28","CD69","FOS","JUN","HSPA1B","PDCD1","CXCR6","CD55","NFKB1","IL7R","GZMB","PRF1","GZMA","GZMK","GIMAP7","KLRG1","MKI67","STMN1")
DotPlot(AS.combined_2_CD8_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()


#After reordering, save as RDS file (!!This is the final version!!)
saveRDS(AS.combined_2_CD8_3, "~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD8_3.rds")
AS.combined_2_CD8_3<- readRDS ("~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD8_3.rds")

DimPlot(AS.combined_2_CD8_3, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD8_3, reduction = "umap", label = FALSE,  repel = TRUE)

###Split into AS and AR
AS.combined_2_CD8_3_AS <- subset(AS.combined_2_CD8_3, subset = ASAR =="AS")
AS.combined_2_CD8_3_AR <- subset(AS.combined_2_CD8_3, subset = ASAR =="AR")

DimPlot(AS.combined_2_CD8_3_AS, reduction = "umap", label = FALSE,  repel = TRUE)
DimPlot(AS.combined_2_CD8_3_AR, reduction = "umap", label = FALSE,  repel = TRUE)
DimPlot(AS.combined_2_CD8_3_AS, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AS.combined_2_CD8_3_AR, reduction = "umap", label = TRUE,  repel = TRUE)

###Count cell number of the subcluster (AS)
AS.combined_2_CD8_3_AS_0 <- subset(AS.combined_2_CD8_3_AS, idents = 0)
AS.combined_2_CD8_3_AS_1 <- subset(AS.combined_2_CD8_3_AS, idents = 1)
AS.combined_2_CD8_3_AS_2 <- subset(AS.combined_2_CD8_3_AS, idents = 2)
AS.combined_2_CD8_3_AS_3 <- subset(AS.combined_2_CD8_3_AS, idents = 3)
AS.combined_2_CD8_3_AS_4 <- subset(AS.combined_2_CD8_3_AS, idents = 4)

###Count cell number of the subcluster (AR)
AS.combined_2_CD8_3_AR_0 <- subset(AS.combined_2_CD8_3_AR, idents = 0)
AS.combined_2_CD8_3_AR_1 <- subset(AS.combined_2_CD8_3_AR, idents = 1)
AS.combined_2_CD8_3_AR_2 <- subset(AS.combined_2_CD8_3_AR, idents = 2)
AS.combined_2_CD8_3_AR_3 <- subset(AS.combined_2_CD8_3_AR, idents = 3)
AS.combined_2_CD8_3_AR_4 <- subset(AS.combined_2_CD8_3_AR, idents = 4)

# TNF NFkB1 AS vs AR
VlnPlot(AS.combined_2_CD8_3,pt.size = 0.05,assay ="RNA", features = c("TNF"),split.by = "ASAR", cols = c("AR" = "blue", "AS" = "red")) 
VlnPlot(AS.combined_2_CD8_3,pt.size = 0.05,assay ="RNA", features = c("NFKB1"),split.by = "ASAR",cols = c("AR" = "blue", "AS" = "red"))
VlnPlot(AS.combined_2_CD8_3,pt.size = 0.05,assay ="RNA", features = c("IFNG"),split.by = "ASAR",cols = c("AR" = "blue", "AS" = "red"))

# volcano  CD8------
AS.combined_2_CD8_3$celltype <- Idents(AS.combined_2_CD8_3)
AS.combined_2_CD8_3$celltype.stim <- paste(Idents(AS.combined_2_CD8_3), AS.combined_2_CD8_3$ASAR, sep="_")
Idents(AS.combined_2_CD8_3) <- "celltype.stim" # Update the identities using the Idents() function
levels(AS.combined_2_CD8_3) # Check that the identities have been updated

AS.combined_2_CD8_3.table <- FindMarkers(AS.combined_2_CD8_3, ident.1 = c("0_AS", "1_AS", "2_AS", "3_AS", "4_AS"), 
                                         ident.2 =c("1_AR","3_AR", "2_AR", "0_AR", "4_AR"), verbose = FALSE, logfc.threshold = 0)
AS.combined_2_CD8_3.table$logp <- -log10(AS.combined_2_CD8_3.table$p_val)
AS.combined_2_CD8_3_filtered_left = subset(AS.combined_2_CD8_3.table, logp>=1.0 & avg_log2FC <= -0.5)
AS.combined_2_CD8_3_filtered_right = subset(AS.combined_2_CD8_3.table, logp>=1.0 & avg_log2FC >= 0.5)

genes.to.label.left <- rownames(AS.combined_2_CD8_3_filtered_left)
genes.to.label.right <- rownames(AS.combined_2_CD8_3_filtered_right)

p1 <- ggplot(AS.combined_2_CD8_3.table, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1