library(Rcpp)
library(Seurat)
library(dplyr)

#Import RDS file
data <- readRDS (file =  "~/Desktop/Human_AS/RDSfiles/AS_combined_myeloid_240125.rds")

# The following code extracts upregulated genes with avg_log2FC > 0.
# To extract downregulated genes, use avg_log2FC < 0 instead.

# Identify differentially expressed genes for Type 0
type0.markers <- FindMarkers(data, ident.1 = "0", ident.2 = c("1","2","3","4"))
type0.markers <- type0.markers[type0.markers$p_val_adj < 0.05 & type0.markers$avg_log2FC > 0,]

# Identify differentially expressed genes for Type 1
type1.markers <- FindMarkers(data, ident.1 = "1", ident.2 = c("0","2","3","4"))
type1.markers <- type1.markers[type1.markers$p_val_adj < 0.05 & type1.markers$avg_log2FC > 0,]

# Identify differentially expressed genes for Type 2
type2.markers <- FindMarkers(data, ident.1 = "2", ident.2 = c("0","1","3","4"))
type2.markers <- type2.markers[type2.markers$p_val_adj < 0.05 & type2.markers$avg_log2FC > 0,]

# Identify differentially expressed genes for Type 3
type3.markers <- FindMarkers(data, ident.1 = "3", ident.2 = c("0","1","2","4"))
type3.markers <- type3.markers[type3.markers$p_val_adj < 0.05 & type3.markers$avg_log2FC > 0,]

# Identify differentially expressed genes for Type 4
type4.markers <- FindMarkers(data, ident.1 = "4", ident.2 = c("0","1","2","3"))
type4.markers <- type4.markers[type4.markers$p_val_adj < 0.05 & type4.markers$avg_log2FC > 0,]

# Obtain a mapping table from gene symbols to gene IDs
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

library(clusterProfiler)

# Convert gene symbols of differentially expressed genes in Type 0 to Entrez Gene IDs
type0.gene_SYMBOLs <- rownames(type0.markers)
type0.gene_IDs <- AnnotationDbi::select(hs, keys = type0.gene_SYMBOLs,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")$ENTREZID

# Convert gene symbols of differentially expressed genes in Type 1 to Entrez Gene IDs
type1.gene_SYMBOLs <- rownames(type1.markers)
type1.gene_IDs <- AnnotationDbi::select(hs, keys = type1.gene_SYMBOLs,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")$ENTREZID

# Convert gene symbols of differentially expressed genes in Type 2 to Entrez Gene IDs
type2.gene_SYMBOLs <- rownames(type2.markers)
type2.gene_IDs <- AnnotationDbi::select(hs, keys = type2.gene_SYMBOLs,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")$ENTREZID

# Convert gene symbols of differentially expressed genes in Type 3 to Entrez Gene IDs
type3.gene_SYMBOLs <- rownames(type3.markers)
type3.gene_IDs <- AnnotationDbi::select(hs, keys = type3.gene_SYMBOLs,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")$ENTREZID

# Convert gene symbols of differentially expressed genes in Type 4 to Entrez Gene IDs
type4.gene_SYMBOLs <- rownames(type4.markers)
type4.gene_IDs <- AnnotationDbi::select(hs, keys = type4.gene_SYMBOLs,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL")$ENTREZID

# compare cluster 
type0.gene_SYMBOLs <- rownames(type0.markers)
type0.gene_IDs <- AnnotationDbi::select(hs, keys=type0.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type1.gene_SYMBOLs <- rownames(type1.markers)
type1.gene_IDs <- AnnotationDbi::select(hs, keys=type1.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type2.gene_SYMBOLs <- rownames(type2.markers)
type2.gene_IDs <- AnnotationDbi::select(hs, keys=type2.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type3.gene_SYMBOLs <- rownames(type3.markers)
type3.gene_IDs <- AnnotationDbi::select(hs, keys=type3.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type4.gene_SYMBOLs <- rownames(type4.markers)
type4.gene_IDs <- AnnotationDbi::select(hs, keys=type4.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID


genelist <- list(type0.gene_IDs, type1.gene_IDs, type2.gene_IDs, type3.gene_IDs,type4.gene_IDs)
names(genelist) <- c("type0", "type1", "type2","type3","type4")

# GO:BP enrichment analysis
cgBP <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="BP",OrgDb='org.Hs.eg.db')
dotplot(cgBP)

# GO:MF enrichment analysis
cgMF <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="MF",OrgDb='org.Hs.eg.db')
dotplot(cgMF)

# GO:CC enrichment analysis
cgCC <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="CC",OrgDb='org.Hs.eg.db')
dotplot(cgCC)


