#library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scMetabolism)
library(rsvd)

AS_Mye <- readRDS ("~/Desktop/Human_AS/RDSfiles/AS_combined_myeloid_240125.rds")

levels(AS_Mye)

DimPlot(AS_Mye, reduction = "umap", label = TRUE, repel = TRUE)

#Change Cluster name
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9")
names(new.cluster.ids) <- levels(CARO.combined_Mye)
CARO.combined_Mye<- RenameIdents(CARO.combined_Mye, new.cluster.ids)

AS_Mye_AR <- subset(AS_Mye, ASAR == "AR")
AS_Mye_AS <- subset(AS_Mye, ASAR == "AS")

DimPlot(AS_Mye_AS, reduction = "umap", label = TRUE, repel = TRUE)

#scMeta---
DefaultAssay(AS_Mye_AR) <- "RNA"
countexp.Seurat<-sc.metabolism.Seurat(obj = AS_Mye_AR, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
saveRDS(countexp.Seurat,"~/Desktop/Human_AS/RDSfiles/ARsubset_scMetabolism_240520.rds.rds")

DefaultAssay(AS_Mye_AS) <- "RNA"
countexp.Seurat.1<-sc.metabolism.Seurat(obj = AS_Mye_AS, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
saveRDS(countexp.Seurat.1,"~/Desktop/Human_AS/RDSfiles/ASsubset_scMetabolism_240520.rds.rds")

# countexp.Seurat <- readRDS("")

input.pathway<-c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)","Folate biosynthesis",
                 "Fatty acid biosynthesis","Fatty acid elongation","Fatty acid degradation","Oxidative phosphorylation",
                 "Fructose and mannose metabolism", "Tyrosine metabolism", "Riboflavin metabolism", "Phosphonate and phosphinate metabolism",
                 "Porphyrin and chlorophyll metabolism","Folate biosynthesis", "Glycosaminoglycan degradation",
                 "Amino sugar and nucleotide sugar metabolism", "Ascorbate and aldarate metabolism",
                 "Other glycan degradation","Pentose and glucuronate interconversions","Phenylalanine metabolism")

input.pathway<-c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)",
                 "Fatty acid biosynthesis",
                 "Glycosaminoglycan degradation","Tyrosine metabolism",
                 "Amino sugar and nucleotide sugar metabolism", 
                "Phenylalanine metabolism")


DimPlot(countexp.Seurat, reduction = "umap", label = TRUE, repel = TRUE,raster = FALSE)
DimPlot(countexp.Seurat, reduction = "umap", split.by = "sample")
DimPlot(countexp.Seurat, reduction = "umap", split.by = "stim")

DimPlot.metabolism(obj = countexp.Seurat, pathway = "Citrate cycle (TCA cycle)", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)
DimPlot.metabolism(obj = countexp.Seurat.1, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)
DimPlot.metabolism(obj = countexp.Seurat.1, pathway = "Citrate cycle (TCA cycle)", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)


DimPlot.metabolism(obj = countexp.Seurat, pathway = "Folate biosynthesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Fatty acid biosynthesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Fatty acid elongation", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Fatty acid degradation", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Oxidative phosphorylation", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)



countexp.Seurat@active.ident

countexp.Seurat$ASAR

paste(countexp.Seurat$ASAR, countexp.Seurat@active.ident, sep="_")

countexp.Seurat<- AddMetaData(countexp.Seurat , countexp.Seurat@active.ident, col.name = "ident")

countexp.Seurat<- AddMetaData(countexp.Seurat , paste(countexp.Seurat$ASAR, countexp.Seurat@active.ident, sep="_"), col.name = "ASAR_ident")

DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "ASAR_ident", norm = "y")
