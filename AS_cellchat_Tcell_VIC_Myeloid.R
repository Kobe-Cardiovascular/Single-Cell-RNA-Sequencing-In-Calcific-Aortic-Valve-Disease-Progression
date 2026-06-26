library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

AS_FLEX <- readRDS( "~/Desktop/Human_AS/RDSfiles/ASFLEX_vic_downsampling_240804.rds" )
AS_CD4 <- readRDS("~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD4_2並び替え後.rds")
AS_CD8 <- readRDS("~/Desktop/Human_AS/RDSfiles/shoda_AS.combined_2_CD8_3.rds")
AS_Myeloid <- readRDS("~/Desktop/Human_AS/RDSfiles/AS_combined_myeloid_240125.rds")


#Cluster name change
new.cluster.ids <- c("V.0","V.1","V.2","V.3","V.4","V.5")
names(new.cluster.ids) <- levels(AS_FLEX)
AS_FLEX<- RenameIdents(AS_FLEX, new.cluster.ids)

new.cluster.ids <- c("CD4.0","CD4.1","CD4.2","CD4.3","CD4.4","CD4.5")
names(new.cluster.ids) <- levels(AS_CD4)
AS_CD4<- RenameIdents(AS_CD4, new.cluster.ids)

new.cluster.ids <- c("CD8.0","CD8.1","CD8.2","CD8.3","CD8.4")
names(new.cluster.ids) <- levels(AS_CD8)
AS_CD8<- RenameIdents(AS_CD8, new.cluster.ids)

new.cluster.ids <- c("My.0","My.1","My.2","My.3","My.4","My.5","My.6")
names(new.cluster.ids) <- levels(AS_Myeloid)
AS_Myeloid<- RenameIdents(AS_Myeloid, new.cluster.ids)

#VIC + Tcell(CD4, CD8)+ Myeloid

VIC_Myeloid <- merge(AS_FLEX, y = AS_Myeloid)
VIC_CD4_T_AS <- merge(VIC_Myeloid, y = AS_CD4)
VIC_T_Mye_AS <- merge(VIC_CD4_T_AS, y = AS_CD8)

levels(VIC_T_Mye_AS)

VIC_T_Mye_AS <- NormalizeData(VIC_T_Mye_AS)

DefaultAssay(VIC_T_Mye_AS) <- "RNA"
cellchat <- createCellChat(object = VIC_T_Mye_AS, group.by = "ident")


CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB[["interaction"]]$receptor[763] <- "IL1R1"
CellChatDB[["interaction"]]$receptor[763] 

CellChatDB[["interaction"]]$receptor[762] <- "IL1R1"
CellChatDB[["interaction"]]$receptor[762] 

cellchat@DB <- CellChatDB

#### Preprocessing gene expression data for cell–cell communication analysis ####
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

# Store signaling genes overexpressed in each cluster in object@var.features
cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.1)
# View(cellchat@var.features$features.info)
# Identify overexpressed ligand–receptor interactions and save them to object@LR
cellchat <- identifyOverExpressedInteractions(cellchat)
# View(cellchat@LR$LRsig)

#### Compute communication probabilities and infer the cell–cell communication network ####
# Estimate communication probabilities and store them in object@net
cellchat <- computeCommunProb(cellchat, type="truncatedMean")

# Filter out communications involving only a small number of cells
cellchat <- filterCommunication(cellchat, min.cells = 10)

#### Infer cell–cell communication at the signaling pathway level ####
# Infer signaling pathway–level communication and store the results in object@netP
cellchat <- computeCommunProbPathway(cellchat)
#write.csv(cellchat@netP$prob, "~/Desktop/JRAS/R-CellChat/output/netP.csv")

#### Aggregate the cell–cell communication network ####
# Store the aggregated cell–cell communication network in object@net$count and object@net$weight
cellchat <- aggregateNet(cellchat)
# View(cellchat@net$count)
# View(cellchat@net$weight)

#AS_Tcell+VIC
saveRDS(cellchat, "~/Desktop/Human_AS/RDSfiles/AS_cellchat_20250805.rds")
cellchat <- readRDS("~/Desktop/Human_AS/RDSfiles/AS_cellchat_20250805.rds")

#---------------------------------
groupSize <- as.numeric(table(cellchat@idents)) # Store the number of cells in each cluster
par(mfrow = c(1,2), xpd = TRUE)  # Set plotting parameters to display the next two plots side by side
# Plot the number and strength (weight) of cell–cell communications
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Visualize cell–cell communications for each cluster separately
mat <- cellchat@net$weight
#dev.off()
par(mfrow = c(1,2), xpd=TRUE) 
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], arrow.size = 0.0001)
}

#### Visualize individual signaling pathways (Hierarchy plot, Circle plot, and Chord diagram) ####
# Select the signaling pathway(s) to visualize from the list below
cellchat@netP$pathways
pathways.show <- c("NKG2D")  # Specify the signaling pathway(s) selected above

#print(cellchat@netP$pathways)
#print(cellchat@netP$prob[,,15])

# Visualize cell–cell communication
# Check the order of clusters
levels(cellchat@idents)

# Specify the cells to display on the left side using their indices in cellchat@idents
vertex.receiver = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18) # Numeric vector; myeloid clusters

# Hierarchy plot
# Circle plot
par(mfrow = c(1,1))
library(scales) # Library used to obtain color information
n_cluster <- 24 # Specify the number of clusters
palette <- hue_pal()(n_cluster) # Generate color palette

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",
                    color.use = palette, show.legend = TRUE, arrow.size = 1)

netVisual_aggregate(cellchat, signaling = pathways.show,
                    layout = "hierarchy", vertex.receiver = vertex.receiver)

# Circle plot
# Visualize signaling interactions as arrows between clusters arranged in a circle
par(mfrow = c(1,3))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")