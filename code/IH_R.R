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