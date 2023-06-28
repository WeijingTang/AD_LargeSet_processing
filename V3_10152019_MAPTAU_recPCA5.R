library(future)
library(future.apply)
library(Matrix)
library(dplyr)
library(Seurat, lib.loc="/home/users/ogarmar3/R/Seurat_v3.1/")
library(ggplot2)
library(cowplot)
library(gtools)
library(lattice)
#library(loomR)
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 128000 * 1024^2)


wdir <- ("/scratch/users/ogarmar3/10x_MAPTAU_v3_10152019_recPCA4/")
dir.create(wdir)
outdir <- paste(wdir, "out/", sep="")
dir.create(outdir)
datadir <- ("/scratch/users/ogarmar3/_data_10x_MAPTAUCont/")


# Metadata table
mdtable<-read.csv("/oak/stanford/groups/icobos/shared/Single_Cell/Processed_DGE/10x_Processed_cellranger-v3_summarized/metadatatable_10x_10082019.csv")

setwd(datadir)
file=list.files()
fileX<-file[1:16]
names <- fileX 
#assign(paste(gsub('.dge.txt.gz', '', file2[i]), sep=""), read.table(file = file2[[i]], sep = "\t", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 3000))))
#data<-read.table(file = paste(names[i],".dge.txt.gz", sep=""), sep = "\t", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 3000)))

for (i in 1:length(fileX)){
  assign(paste(gsub('', '', fileX[i]), sep=""),Read10X(paste(datadir, fileX[i],"/outs/filtered_feature_bc_matrix/",sep="")))
}

mdtableX<-mdtable[(mdtable$Unique.identifier %in% fileX),]

setwd(outdir)

list <- mget(names)

for  (i in 1:length(unique(names))){
  rm(list=names[i])
}
gc()

for (i in 1:length(unique(list))){
  dataf <- list[[i]]
  # Filter out ucharacterized RP* genes
  #dataf <- dataf[!grepl("RP1-|RP11-|RP13-|RP3-|RP4-|RP5-|RP6-", row.names(dataf)),]
  #dataf <- dataf[!grepl("MT-", row.names(dataf)),]
  #dataf <- dataf[!grepl("MALAT", row.names(dataf)),]
  #dataf <- dataf[!grepl("MIAT", row.names(dataf)),]
  #dataf <- dataf[!grepl("MEG3", row.names(dataf)),]
  #dataf <- dataf[!grepl("RP[SL][[:digit:]]", row.names(dataf)),]
  # Assign unique names to cells
  colnames(dataf) <- paste(names[i], colnames(dataf), sep = "-")
  # Create seurat object
  dataf <- CreateSeuratObject(dataf, min.cells = 3, min.features = 250, project = names[i])
  
  # Add metadata
  mdtablei<-mdtableX[i,]
  mdsample<-dataf@meta.data
  comb<-mdsample
  comb[colnames(mdtable)]<-mdtablei
  #comb<-rbind.fill(comb,mdtablei)
  dataf@meta.data<-comb
  
  # Get percent.mt and percent.ribo
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(object = dataf)), value = TRUE)
  percent.mt <- Matrix::colSums(GetAssayData(object = dataf, slot = "counts")[mito.genes, ])/Matrix::colSums(GetAssayData(object = dataf, slot = "counts"))
  dataf$percent.mt <- percent.mt
  ribo.genes <- grep(pattern = "RP[SL][[:digit:]]", x = rownames(x = GetAssayData(object = dataf)), value = TRUE)
  percent.ribo <- Matrix::colSums(GetAssayData(object = dataf, slot = "counts")[ribo.genes, ])/Matrix::colSums(GetAssayData(object = dataf, slot = "counts"))
  dataf$percent.ribo <- percent.ribo
  
  # Assign object
  #assign(paste(gsub('-MAP', '', names[i]), sep=""), data)
  assign(paste(names[i], sep=""), dataf)
  #rm(list=names[i])
  rm(dataf)
}

list <- mget(names)
for  (i in 1:length(unique(names))){
  rm(list=names[i])
}
gc()

for (i in 1:length(names)){
  x <- list[[i]]
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
  assign(paste(names[i], sep=""), x)
  rm(x)
}

list <- mget(names)
for  (i in 1:length(unique(names))){
  rm(list=names[i])
}
gc()

features <- SelectIntegrationFeatures(object.list = list)

for (i in 1:length(list)){
  x <- list[[i]]
  x <- ScaleData(x, features = features, verbose = T)#, vars.to.regress = c("nCounts_RNA", "percent.mt"))
  x <- RunPCA(x, features = features, verbose = T)
  assign(paste(names(list)[i], sep=""), x)
  rm(x)
}

list <- mget(names)
for  (i in 1:length(unique(names))){
  rm(list=names[i])
}
gc()


anchors <- FindIntegrationAnchors(object.list = list, reference = c(5,9,11,15), reduction = "rpca", 
                                dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
integrated <- ScaleData(integrated, verbose = FALSE, vars.to.regress = c("percent.mt", "nCount_RNA"))
integrated <- RunPCA(integrated, verbose = FALSE)
png(filename="ElbowPlot.png", width=800, height=600, pointsize=3)
ElbowPlot(integrated, ndims=40)
dev.off()

##########
pcs <- 25
res <- 0.7
##########

integrated <- RunUMAP(integrated, dims = 1:pcs)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:pcs)
integrated <- FindClusters(integrated, resolution = res)

# Visualization
p1 <- DimPlot(integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE)
p3 <- DimPlot(integrated, reduction = "umap", group.by = "SORT")
p3.2 <- DimPlot(integrated, reduction = "umap", group.by = "NextSeq_01172018")
p4 <- DimPlot(integrated, reduction = "umap", group.by = "NP.Diagnosis")
p4.2 <- DimPlot(integrated, reduction = "umap", group.by = "Braak.stage")
p5 <- DimPlot(integrated, reduction = "umap", group.by = "Amyloid")
p6 <- DimPlot(integrated, reduction = "umap", group.by = "Age")
p7 <- DimPlot(integrated, reduction = "umap", group.by = "Gender")
p8 <- DimPlot(integrated, reduction = "umap", group.by = "RIN")
p9 <- DimPlot(integrated, reduction = "umap", group.by = "Assay")
p14 <- FeaturePlot(object=integrated, features = "nCount_RNA")
png(filename=paste(pcs,"c",res,"r_UMAP.png",sep=""), width=1300, height=1200, pointsize=3)
plot_grid(p1, p2, p3, p4.2)
dev.off()
png(filename=paste(pcs,"c",res,"r_metadata.png",sep=""), width=3300, height=3200, pointsize=3)
plot_grid(p1, p2, p3, p3.2, p4.2, p4, p5, p6, p7, p8, p9, p14)
dev.off()

setwd(outdir)
saveRDS(integrated, paste(pcs,"c",res,"r_Integrated.rds",sep=""))
#integrated <- readRDS("25c0.7r_Integrated.rds")

DefaultAssay(integrated) <- "RNA"

png(filename=paste(pcs,"c",res,"r_UMAP-Glia.png",sep=""), width=1700, height=1400, pointsize=5)
print(FeaturePlot(object = integrated, features = c("PLP1", "MBP", "CD74", "PDGFRA", "APOE", "GFAP", "SOX10", "OLIG1"), max.cutoff = "q90"))
dev.off()
png(filename=paste(pcs,"c",res,"r_UMAP-EN.png",sep=""), width=1700, height=1600, pointsize=5)
print(FeaturePlot(object = integrated, features = c("SLC17A7", "CUX2", "LAMP5", "COL5A2", "CARTPT", "RORB", "GAL", "GABRG1", "RPRM", "PCP4", "NFIA", "THEMIS", "ROBO3", "NR4A2", "NTNG2", "CTGF"), max.cutoff = "q90"))
dev.off()
png(filename=paste(pcs,"c",res,"r_UMAP-IN.png",sep=""), width=1700, height=1600, pointsize=5)
print(FeaturePlot(object = integrated, features = c("GAD1", "GAD2", "LHX6", "PVALB", "SCUBE3", "SST", "NPY", "ADARB2", "LAMP5", "CXCL14", "CALB2", "VIP"), max.cutoff = "q90"))
dev.off()
png(filename=paste(pcs,"c",res,"r_UMAP-nCounts_RNA.png",sep=""), width=700, height=600, pointsize=5)
print(FeaturePlot(object=integrated, features = "nCount_RNA"), max.cutoff = "q90")
dev.off()

cluster.ident <- table(integrated@active.ident)
cluster.ident <- as.data.frame(cluster.ident, sep = ",", col.names=c("A","B"))
colnames(cluster.ident) <- c("CLUSTER", "TOTAL CELLS")
cluster.ident <- t(cluster.ident)
case.ident <- table(integrated@meta.data$orig.ident, integrated@active.ident)
combined <- rbind(cluster.ident, case.ident)
write.table(combined, paste(pcs,"c",res,"r_CellsPerCluster.csv",sep=""), sep = ",", col.names = NA)

integrated.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#write.csv(integrated.markers, "ClusterMarkers.csv")
#write.table(integrated.markers, "ClusterMarkers_t.csv", sep = ",", col.names = NA)
write.csv(integrated.markers, paste(pcs,"c",res,"r_ClusterMarkers.csv",sep=""))




quit(save="no")











###########
# Subsets # 25c07r
###########
datadir <- ("/scratch/users/ogarmar3/_data_10x_MAPTAUCont/")
# Metadata table
mdtable<-read.csv("/oak/stanford/groups/icobos/shared/Single_Cell/Processed_DGE/10x_Processed_cellranger-v3_summarized/metadatatable_10x_10082019.csv")
setwd(datadir)
file=list.files()
fileX<-file#[1:2]
names <- fileX 


setwd(outdir)
integrated <- readRDS("25c0.7r_Integrated.rds")


# 25c07r
excit <- subset(x = integrated, idents = c("4","8","6","3","19","7","2", "24", "0"), invert = F)
rm(integrated)
all.genes <- rownames(excit)

gc()

#Excitatory
########################################################
outdirEx <- paste(outdir, paste("rerun_Ex_25c07r/", sep=""), sep="")
dir.create(outdirEx)
setwd(outdirEx)

p1 <- DimPlot(excit, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(excit, reduction = "umap", label = TRUE)
p3 <- DimPlot(excit, reduction = "umap", group.by = "SORT")
p3.2 <- DimPlot(excit, reduction = "umap", group.by = "NextSeq_01172018")
p4.2 <- DimPlot(excit, reduction = "umap", group.by = "Braak.stage")
png(filename="PRERUN_UMAP.png", width=1300, height=1200, pointsize=3)
plot_grid(p1, p2, p3, p4.2)
dev.off()
###############################

list<- SplitObject(excit, split.by = "orig.ident")

for (i in 1:length(names)){
  x <- list[[i]]
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
  assign(paste(names[i], sep=""), x)
  rm(x)
}

list <- mget(names)
for  (i in 1:length(unique(names))){
  rm(list=names[i])
}
gc()

features <- SelectIntegrationFeatures(object.list = list)

for (i in 1:length(names)){
  x <- list[[i]]
  x <- ScaleData(x, features = features, verbose = T)#, vars.to.regress = c("nCounts_RNA", "percent.mt"))
  x <- RunPCA(x, features = features, verbose = T)
  assign(paste(names[i], sep=""), x)
  rm(x)
}

list <- mget(names)
for  (i in 1:length(unique(names))){
  rm(list=names[i])
}
gc()


anchors <- FindIntegrationAnchors(object.list = list, reference = c(5,9,11,15,17,19,21,23), reduction = "rpca", 
                                  dims = 1:30)
saveRDS(anchors, "anchors_excit.rds")
#anchors<-readRDS("anchors_excit.rds")
rm(list)
rm(excit)
gc()

integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
integrated <- ScaleData(integrated, verbose = FALSE, vars.to.regress = c("percent.mt", "nCount_RNA"))
integrated <- RunPCA(integrated, verbose = FALSE)
png(filename="ElbowPlot.png", width=800, height=600, pointsize=3)
ElbowPlot(integrated, ndims=30)
dev.off()

##########
pcs<-20
res<-0.5
##########

excit <- RunUMAP(excit, dims = 1:pcs)
excit <- FindNeighbors(excit, reduction = "pca", dims = 1:pcs)
excit <- FindClusters(excit, resolution = res)

# Visualization
p1 <- DimPlot(excit, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(excit, reduction = "umap", label = TRUE)
p3 <- DimPlot(excit, reduction = "umap", group.by = "SORT")
p3.2 <- DimPlot(excit, reduction = "umap", group.by = "Sequencing.batch")
p4 <- DimPlot(excit, reduction = "umap", group.by = "NP.Diagnosis")
p4.2 <- DimPlot(excit, reduction = "umap", group.by = "Braak.stage")
p5 <- DimPlot(excit, reduction = "umap", group.by = "Amyloid")
p6 <- DimPlot(excit, reduction = "umap", group.by = "Age")
p7 <- DimPlot(excit, reduction = "umap", group.by = "Gender")
p8 <- DimPlot(excit, reduction = "umap", group.by = "RIN")
p9 <- DimPlot(excit, reduction = "umap", group.by = "Assay")
p14 <- FeaturePlot(object=excit, features = "nCount_RNA")
png(filename=paste(pcs,"c",res,"r_UMAP.png",sep=""), width=1300, height=600, pointsize=3)
plot_grid(p1, p2, p3, p4.2)
dev.off()
png(filename=paste(pcs,"c",res,"r_metadata.png",sep=""), width=3300, height=3200, pointsize=3)
plot_grid(p1, p2, p3, p3.2, p4.2, p4, p5, p6, p7, p8, p9, p14)
dev.off()

setwd(outdir)
saveRDS(excit, paste(pcs,"c",res,"r_excit.rds",sep=""))
#excit <- readRDS("30c07r_excit.rds")

DefaultAssay(excit) <- "RNA"

png(filename=paste(pcs,"c",res,"r_UMAP-Glia.png",sep=""), width=1700, height=1400, pointsize=5)
print(FeaturePlot(object = excit, features = c("PLP1", "MBP", "CD74", "PDGFRA", "APOE", "GFAP", "SOX10", "OLIG1")))#, cols.use = c("grey", "blue"), reduction.use = "UMAP")
dev.off()
png(filename=paste(pcs,"c",res,"r_UMAP-EN.png",sep=""), width=1700, height=1600, pointsize=5)
print(FeaturePlot(object = excit, features = c("SLC17A7", "CUX2", "LAMP5", "COL5A2", "CARTPT", "RORB", "GAL", "GABRG1", "RPRM", "PCP4", "NFIA", "THEMIS", "ROBO3", "NR4A2", "NTNG2", "CTGF")))#)#
dev.off()
png(filename=paste(pcs,"c",res,"r_UMAP-IN.png",sep=""), width=1700, height=1600, pointsize=5)
print(FeaturePlot(object = excit, features = c("GAD1", "GAD2", "LHX6", "PVALB", "SCUBE3", "SST", "NPY", "ADARB2", "LAMP5", "CXCL14", "CALB2", "VIP")))#
dev.off()
png(filename=paste(pcs,"c",res,"r_UMAP-nCounts_RNA.png",sep=""), width=700, height=600, pointsize=5)
print(FeaturePlot(object=excit, features = "nCount_RNA"))#, cols.use = c("grey","blue"), pt.size=5, nCol=1)
dev.off()

cluster.ident <- table(excit@active.ident)
cluster.ident <- as.data.frame(cluster.ident, sep = ",", col.names=c("A","B"))
colnames(cluster.ident) <- c("CLUSTER", "TOTAL CELLS")
cluster.ident <- t(cluster.ident)
case.ident <- table(excit@meta.data$orig.ident, excit@active.ident)
combined <- rbind(cluster.ident, case.ident)
write.table(combined, paste(pcs,"c",res,"r_CellsPerCluster.csv",sep=""), sep = ",", col.names = NA)

excit.markers <- FindAllMarkers(excit, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
excit.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#write.csv(excit.markers, "ClusterMarkers.csv")
#write.table(excit.markers, "ClusterMarkers_t.csv", sep = ",", col.names = NA)
write.csv(excit.markers, paste(pcs,"c",res,"r_ClusterMarkers.csv",sep=""))





quit(save = "no")

