library(edgeR)
library(openxlsx)

#### Count data import 
data <- read.xlsx("~/Desktop/Human_AS/Bulkdata/PR4211/expression/PR4211_gene_expression.xlsx")


View(data)
gene_row_name = "Gene.Name"
sample_use = c("PR4211_01_a.TPM",
               "PR4211_02_a.TPM",
               "PR4211_03_a.TPM",
               "PR4211_04_a.TPM",
               "PR4211_05_a.TPM",
               "PR4211_06_a.TPM",
               "PR4211_07_a.TPM",
               "PR4211_08_a.TPM",
               "PR4211_09_a.TPM",
               "PR4211_10_a.TPM",
               "PR4211_12_a.TPM",
               "PR4211_13_a.TPM",
               "PR4211_14_a.TPM",
               "Gene.ID")

sample_group = c("1",
                 "1",
                 "1",
                 "1",
                 "2",
                 "2",
                 "2",
                 "2",
                 "2",
                 "2",
                 "2",
                 "2",
                 "2")

pair = c("1", "2")


#### Extract data ####
columns_to_use <- c(gene_row_name,
                    sample_use)

data.filtered <- data[columns_to_use]

count.data<-as.matrix(data.filtered[2:(ncol(data.filtered)-1)])
gene_name <-as.character(as.vector(data.filtered[,gene_row_name]))
gene_id <-as.character(as.vector(data.filtered[,"Gene.ID"])) 
rownames(count.data) <- gene_name
count.data <- as.data.frame(count.data)
      
####  Preprocessing ####
count.keep <- rowSums(count.data) > 0
count.data["Gene.ID"] <- gene_id
count.data <- count.data[count.keep,]

gene_id <- count.data["Gene.ID"]

#### DEG comparison ####
group <- factor(sample_group)
count.dge <- DGEList(counts=count.data[1:13], group=group)
count.dge <- calcNormFactors(count.dge)
count.dge <- estimateCommonDisp(count.dge)
count.dge <- estimateTagwiseDisp(count.dge)

result<- exactTest(count.dge, pair=pair)
result$table["FDR"] <- p.adjust(result$table$PValue, "fdr")

# Output
#dir.create("output")
write.csv(result$table, "~/Desktop/Human_AS/Bulkdata/PR4211/allgene_result.AR.vs.AS.240719.csv")

# DEG file
result$table["Gene.ID"] <- gene_id
write.csv(result$table, "~/Desktop/Human_AS/Results/DEG.AR_vs.AS_3_240809.csv")

## Rename
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")

library("org.Mm.eg.db") 
library("org.Hs.eg.db")
library("clusterProfiler")
library(enrichplot)

DEG <- read.csv("~/Desktop/Human_AS/Results/DEG.AS_vs.AR.csv",row.names = 1, sep=",")

Names.all <- rownames(DEG)

my.ensmusg <- Names.all
ensmusg <- unlist(as.list(org.Hs.egALIAS2EG))
my.ensmusg <- sapply(strsplit(my.ensmusg, "\\."), function(x) x[[1]])
my.eg <- ensmusg[my.ensmusg] 
Names.all <- my.eg


# GO analysis
# Figure S10B
ego <- enrichGO(gene = Names.all, 
                ont ="CC", 
                OrgDb='org.Hs.eg.db')
res <- summary(ego)

#Barplot and dotplot
barplot(ego, showCategory=20)
dotplot(ego, showCategory=20)


##GO_stat csv output
write.csv(res, "output/GO_stat.csv")

#cnetplot
# Figure S10C
p1 <- cnetplot(egox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(egox, node_label="gene", 
               cex_label_gene = 0.8, showCategory = 20) 
p3 <- cnetplot(egox, node_label="all") 
p4 <- cnetplot(egox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
