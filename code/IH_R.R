library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2) 
library(Matrix)

IH = readRDS('E:/SpyderCode/project_9/sample1-3/IH_pro.RDS')
IH$cell_type <- NA
IH_celltype = read.csv('E:/SpyderCode/project_9/sample1-3/IH_pro_celltype.csv')
IH$cell_type = IH_celltype$type
features <- c("CD79A", "PECAM1", "DSP", "OGN","NOTCH3",
              "C1QB", "TPSAB1", "CD1C","S100B","KLRB1")
colors <- list(rgb(0.2118, 0.6235, 0.1765), rgb(0.9882, 0.5019, 0.0078), rgb(0.5725, 0.7608, 0.8667), rgb(0.0863, 0.3882, 0.6627), rgb(0.7765,0.1922,0.2078),rgb(0.3804, 0.2509, 0.6000)
              ,rgb(0.9804, 0.7804, 0.7019),rgb(0.6784, 0.8588, 0.5333),rgb(0.5059, 0.5059, 0.5353),rgb(0.7059, 0.7059, 0.8353))
# pdf(file = "1-c.pdf", width = 10, height = 10)
p2 <- VlnPlot(IH, features, stack = TRUE,
              sort = FALSE, flip = TRUE, cols = colors, group.by="cell_type") +
  theme(legend.position = "none",
        text = element_text(size = 20))
p2
# dev.off()
# print('ok')


# RNA velocity

library(circlize)
library(gatepoints)
library(stringr)
library(igraph)
library(gmodels)
library(Seurat)

IH.combined = IH
IH.combined$cell_type <- NA
IH_celltype = read.csv('IH_pro_celltype.csv')
IH.combined$cell_type = IH_celltype$type
DefaultAssay(IH.combined) <- "integrated"
IH.combined <- ScaleData(IH.combined, verbose = FALSE)
IH.combined <- RunPCA(IH.combined, npcs = 30, verbose = FALSE)
IH.combined <- RunUMAP(IH.combined, reduction = "pca", dims = 1:30)
DimPlot(IH.combined, reduction = "umap")
IH.sub = IH.combined # subset(IH.combined, subset = cell_type == "Endothelial cell") #


VEC = IH.sub@reductions$umap@cell.embeddings
rownames(VEC) = colnames(IH.sub)
PCA = IH.sub@reductions$pca@cell.embeddings
source('https://raw.githubusercontent.com/jumphone/Vector/master/Vector.R')
PCA=vector.rankPCA(PCA)
OUT=vector.buildGrid(VEC, N=35,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=25, SHOW=TRUE) #35/25-EC
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.99,SHOW=TRUE) #0.99-EC，0.9
OUT=vector.drawArrow(OUT,P=0.99,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=FALSE)

# CellChat
library(CellChat)
library(scRNAseq)
library(SummarizedExperiment)
options(stringsAsFactors = FALSE)

IH.combined = readRDS('E:/SpyderCode/project_9/sample1-3/IH_pro.RDS')
IH_celltype = read.csv('E:/SpyderCode/project_9/sample1-3/IH_pro_celltype.csv')
data = IH.combined@assays[["RNA"]]@data
cellchat <- createCellChat(object = data)
cellchat <- addMeta(cellchat, meta = IH_celltype)
celltype_ident <- as.factor(IH_celltype$type)
cellchat@idents <- celltype_ident
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
DB <- dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") #Secreted Signaling #Cell-Cell Contact #ECM-Receptor
cellchat@DB <- CellChatDB.use
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)#取出表达数据
cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因
cellchat <- identifyOverExpressedInteractions(cellchat)#寻找高表达的通路
cellchat <- computeCommunProb(cellchat,population.size=TRUE) #数量分数
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# # write.csv(CellChatDB[["interaction"]],'CellChatDB.csv')

colors <- c(rgb(0.9804, 0.7804, 0.7019), rgb(0.9882, 0.5019, 0.0078), rgb(0.597, 0.726, 0.871), rgb(0.0863, 0.3882, 0.662), rgb(0.7765,0.1922,0.2078),rgb(0.3804, 0.2509, 0.6000)
            ,rgb(0.675, 0.644, 0.107),rgb(0.6784, 0.8588, 0.5333),rgb(0.5059, 0.5059, 0.5353), rgb(0.2118, 0.6235, 0.1765))

groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, color.use = colors,vertex.weight = groupSize,  vertex.size.max = 15, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, color.use = colors,vertex.weight = groupSize,  vertex.size.max = 15, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat@netP$pathways
pathways.show <- c("VEGF")
netVisual_aggregate(cellchat, color.use = colors,signaling = pathways.show, vertex.size.max = 15,layout = "circle",vertex.label.cex = 2)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, color.use = colors,signaling = pathways.show, width = 8, height = 2.5, font.size = 12, font.size.title = 14)

# Sub_cluster
IH.combined = readRDS('IH.proliferation.RDS')
IH_celltype = read.csv('IH_pro_celltype.csv')
IH.combined$cell_type = IH_celltype$sub
IH.sub = subset(IH.combined, subset = cell_type %in% c("EC1","EC2","EC3","Hem-MC1","Hem-MC2","Hem-MC3"))
IH_celltype_sub = IH_celltype[IH_celltype$sub %in% c("EC1","EC2","EC3","Hem-MC1","Hem-MC2","Hem-MC3"),]

data = IH.sub@assays[["RNA"]]@data
cellchat <- createCellChat(object = data)
cellchat <- addMeta(cellchat, meta = IH_celltype_sub)
celltype_ident <- as.factor(IH_celltype_sub$sub)
cellchat@idents <- celltype_ident
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
DB <- dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") #Secreted Signaling #Cell-Cell Contact #ECM-Receptor
cellchat@DB <- CellChatDB.use
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)#取出表达数据
cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因
cellchat <- identifyOverExpressedInteractions(cellchat)#寻找高表达的通路
cellchat <- computeCommunProb(cellchat,population.size=TRUE) #数量分数
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

unique(IH_celltype_sub$sub)
colors <- c(rgb(0.0863, 0.3882, 0.6627), rgb(0.0000, 0.4980, 0.3294), rgb(0.484, 0.535, 1),
            rgb(0.8706, 0.3471, 0.1647), rgb(0.047, 0.718, 0.876), rgb(0.992, 0.847, 0.271))
cellchat@netP$pathways

pathways.show <- c("VEGF")
netVisual_aggregate(cellchat, color.use = colors,signaling = pathways.show, vertex.size.max = 15,layout = "circle",vertex.label.cex = 2)

pathways.show <- c("PDGF")
netVisual_aggregate(cellchat, color.use = colors,signaling = pathways.show, vertex.size.max = 15,layout = "circle",vertex.label.cex = 2)

pathways.show <- c("TGFb")
netVisual_aggregate(cellchat, color.use = colors,signaling = pathways.show, vertex.size.max = 15,layout = "circle",vertex.label.cex = 2)

pathways.show <- c("NOTCH")
netVisual_aggregate(cellchat, color.use = colors,signaling = pathways.show, vertex.size.max = 15,layout = "circle",vertex.label.cex = 2)


pathways.show <- c("ANGPT")
netVisual_aggregate(cellchat, color.use = colors,signaling = pathways.show, vertex.size.max = 15,layout = "circle",vertex.label.cex = 2)

