library(Matrix, lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(farver,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(ggplot2,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(vctrs,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(plyr, lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(dplyr,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(Seurat,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(patchwork,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(plotly,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(reticulate,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(future,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(future.apply,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(cowplot,lib="/oak/stanford/groups/icobos/inma/software_packages/")
library(gridExtra, lib="/oak/stanford/groups/icobos/inma/software_packages/")

plan("multiprocess", workers = 2 )
options(future.globals.maxSize = 256000 * 1024^4)

wdir <- ("/oak/stanford/groups/icobos/inma/data_output/")
setwd(wdir)
outdir <- ("/oak/stanford/groups/icobos/inma/data_output/")
setwd(outdir)


#REFERENCE BASED
Ex10xDq2 <- readRDS('/oak/stanford/groups/icobos/inma/data_output/EI_10xDq_AfterQC_WithMeta_IC04282020.rds')
Ex10xDq2.list <- SplitObject(Ex10xDq2, split.by = "Assay.Region")
for (i in names(Ex10xDq2.list)) {
  Ex10xDq2.list[[i]] <- SCTransform(Ex10xDq2.list[[i]], vars.to.regress = c("percent.mt", "nCount_RNA"), variable.features.n = 3000, verbose = TRUE)
}
Ex10xDq2.features <- SelectIntegrationFeatures(object.list = Ex10xDq2.list, nfeatures = 3000)
Ex10xDq2.list <- PrepSCTIntegration(object.list = Ex10xDq2.list, anchor.features = Ex10xDq2.features)

reference_dataset <- which(names(Ex10xDq2.list) == "v3.F")

Ex10xDq2.anchors <- FindIntegrationAnchors(object.list = Ex10xDq2.list, normalization.method = "SCT", 
                                          anchor.features = Ex10xDq2.features, reference = reference_dataset)
Ex10xDq2.integrated <- IntegrateData(anchorset = Ex10xDq2.anchors, normalization.method = "SCT")

Ex10xDq2.integrated <- RunPCA(object = Ex10xDq2.integrated, verbose = TRUE)
Ex10xDq2.integrated <- RunUMAP(object = Ex10xDq2.integrated, dims = 1:20)

saveRDS(Ex10xDq2.integrated, file = "/oak/stanford/groups/icobos/inma/data_output/EI_10xDq_AfterQSCTIntegrIC04282020.rds")


Ex10xDq2.integrated <- FindNeighbors(Ex10xDq2.integrated, dims = 1:20, verbose = FALSE)
Ex10xDq2.integrated <- FindClusters(Ex10xDq2.integrated, verbose = FALSE, resolution = 0.8)


png(filename="plot1d.png", width=3000, height=3000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "SORT.y"))
dev.off()

png(filename="plot2d.png", width=3000, height=3000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "Sequencing.batch.y"))
dev.off()

png(filename="plot3d.png", width=3000, height=3000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, label = TRUE))
dev.off()

png(filename="plot4d.png", width=3000, height=3000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "Brain.Region.y"))
dev.off()

png(filename="plot5d.png", width=3000, height=3000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "Assay.y"))
dev.off()

png(filename="plot6d.png", width=3000, height=3000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "Braak"))
dev.off()

png(filename="plot7d.png", width=3000, height=3000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "Assay.Region"))
dev.off()

png(filename="plot8d.png", width=6000, height=2000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "Assay.Region", split.by = "Brain.Region.y"))
dev.off()

png(filename="plot9d.png", width=6000, height=2000, pointsize=3)
print(DimPlot(Ex10xDq2.integrated, pt.size = 0.1, group.by = "Assay.Region", split.by = "Assay.y"))
dev.off()

png(filename="plot10d.png", width=4000, height=2000, pointsize=3)
print(FeaturePlot(Ex10xDq2.integrated, pt.size = 0.1, features = c("nCount_RNA", "percent.mt")))
dev.off()

