library(Matrix)
library(dplyr)
library(Seurat)#, lib.loc="C:/Users/ogarmar3/Box/_Cobos Lab/Marcos/Seurat_BOX/Seurat-v2.3/") 
library(DropletUtils)
library(DoubletFinder)
library(lattice)

datadir<-("C:/Users/ogarmar3/Box/_Cobos Lab/Marcos/Seurat_BOX/_DGEs_pooled_09192019/")
setwd(datadir)

file=list.files()
file2 <- file#[c(1:63)]

wdir<-("C:/Users/ogarmar3/Box/_Cobos Lab/Marcos/Seurat_BOX/_DGEs_09192019_emptyDrops10012019/")
dir.create(wdir)
setwd(wdir)


for (i in 191:length(unique(file2))){
  setwd(datadir)
  #assign(paste(gsub('.dge.txt.gz', '', file2[i]), sep=""), read.table(file = file2[[i]], sep = "\t", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 3000))))
  assign("mt", read.table(file = file2[[i]], sep = "\t", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 3000))))
  
  outdir<-(paste(wdir,gsub('.dge.txt.gz', '', file2[i]), "_10012019", sep=""))
  dir.create(outdir)
  setwd(outdir)
  
  tag<-paste(gsub('.dge.txt.gz', '', file2[i]), sep="")
  ms<-CreateSeuratObject(mt)
  
  ## barcoderanks ##
  ms.out<-barcodeRanks(ms, lower=500)
  names(ms.out)
  head(ms.out)
  
  png(filename="kneeplot.png", width=1200, height=1000, pointsize=15)
  print(plot(ms.out$rank, ms.out$total, log="xy", xlab="Rank", ylab="Total"))
  o <- order(ms.out$rank)
  lines(ms.out$rank[o], ms.out$fitted[o], col="red")
  abline(h=metadata(ms.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(ms.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  dev.off()
  
  ms.out
  metadata(ms.out)$knee
  metadata(ms.out)$inflection
  
  threshod.inf<-(metadata(ms.out)$inflection)
  
  ## emptydrops ##
  me.out<-emptyDrops(ms@assays$RNA@counts, lower=0.5*threshod.inf, ignore=0.05*threshod.inf)
  is.cell <- me.out$FDR <= 0.01
  sum(is.cell, na.rm=TRUE)
  
  is.cell2<-is.cell
  is.cell2[is.na(is.cell2)]<-FALSE
  
  mf<-mt[,is.cell2]
  
  dataf<-CreateSeuratObject(mf)
  datan<-NormalizeData(dataf)
  
  # Seurat analysis
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(object = dataf)), value = TRUE)
  percent.mt <- Matrix::colSums(GetAssayData(object = dataf,  slot = "counts")[mito.genes, ])/Matrix::colSums(GetAssayData(object = dataf, slot = "counts"))
  dataf$percent.mt <- percent.mt
  ribo.genes <- grep(pattern = "RP[SL][[:digit:]]", x = rownames(x = GetAssayData(object =dataf)), value = TRUE)
  percent.ribo <- Matrix::colSums(GetAssayData(object =dataf, slot = "counts")[ribo.genes, ])/Matrix::colSums(GetAssayData(object =dataf, slot = "counts"))
  dataf$percent.ribo <- percent.ribo
  
  png(filename="Features_vln.png", width=800, height=600, pointsize=3)
  print(VlnPlot(dataf, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  plot1 <- FeatureScatter(dataf, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(dataf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  png(filename="Features_scatter.png", width=800, height=500, pointsize=3)
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
  
  # Metrics
  Totalcells<-length(dataf@meta.data$nFeature_RNA)
  cells250g<-length(dataf@meta.data$nFeature_RNA[dataf@meta.data$nFeature_RNA>250])
  cells100g<-length(dataf@meta.data$nFeature_RNA[dataf@meta.data$nFeature_RNA>100])
  pct250.100g<-(cells250g/cells100g*100)
  
  cells12kU<-length(dataf@meta.data$nCount_RNA[dataf@meta.data$nCount_RNA>12000])
  cells5kU<-length(dataf@meta.data$nCount_RNA[dataf@meta.data$nCount_RNA>5000])
  pct12kU<-(cells12kU/Totalcells*100)
  
  cells5mt<-length(dataf@meta.data$percent.mt[dataf@meta.data$percent.mt>0.05])
  pct5mt<-(cells5mt/Totalcells*100)
  
  f.med.g<-median(dataf@meta.data$nFeature_RNA)#[dataf@meta.data$nFeature_RNA>250])
  f.mean.g<-mean(dataf@meta.data$nFeature_RNA)#[dataf@meta.data$nFeature_RNA>250])
  
  f.med.U<-median(dataf@meta.data$nCount_RNA)#[dataf@meta.data$nFeature_RNA>250])
  f.mean.U<-mean(dataf@meta.data$nCount_RNA)#[dataf@meta.data$nFeature_RNA>250])
  
  combinedf<-data.frame(Totalcells,cells250g,cells100g,pct250.100g,cells12kU,cells5kU,pct12kU,cells5mt,pct5mt,
                        f.med.g,f.mean.g,f.med.U,f.mean.U)
  row.names(combinedf)<-gsub(".dge.txt.gz","",file2[i])
  write.csv(combinedf, "gene_metrics.csv")

  #dataf<-mf# subset(mf, subset = nFeature_RNA > 250 & nCount_RNA < 12000 & percent.mt < 5)
  #Filtered250g12kU5mt.cells<-length(dataf@active.ident)
  
  dataf <- NormalizeData(dataf)
  dataf <- FindVariableFeatures(dataf, selection.method = "vst", nfeatures = 4000)
  
  top20 <- head(VariableFeatures(dataf), 20)
  plot1<-VariableFeaturePlot(dataf)
  png(filename="VariableGenes.png", width=400, height=500, pointsize=3)
  print(LabelPoints(plot = plot1, points = top20, repel = TRUE))
  dev.off()
  
  dataf <- ScaleData(dataf, verbose = FALSE, vars.to.regress = c("percent.mt", "nCount_RNA"))
  dataf <- RunPCA(dataf, npcs = 40, verbose = FALSE)
  dataf <- RunUMAP(dataf, reduction = "pca", dims = 1:40)
  dataf <- FindNeighbors(dataf, reduction = "pca", dims = 1:40)
  dataf <- FindClusters(dataf, resolution = 0.8)
  
  png(filename="UMAP_plot.png", width=700, height=600, pointsize=3)
  print(DimPlot(dataf, reduction = "umap", label = TRUE))
  dev.off()
  
  png(filename="MarkersGlia_UMAP.png", width=1700, height=1600, pointsize=5)
  print(FeaturePlot(object = dataf, features = c("PLP1", "MBP", "CD74", "PDGFRA", "APOE", "GFAP", "SOX10", "OLIG1")))#, cols.use = c("grey", "blue"), reduction.use = "UMAP")
  dev.off()
  png(filename="Markers_EN_UMAP.png", width=1700, height=1600, pointsize=5)
  print(FeaturePlot(object = dataf, features = c("SLC17A7", "CUX2", "LAMP5", "COL5A2", "RORB", "GAL", "GABRG1", "RPRM", "PCP4", "NFIA", "THEMIS", "ROBO3", "FEZF2", "NR4A2", "NTNG2", "CTGF")))#
  dev.off()
  png(filename="Markers_IN_UMAP.png", width=1700, height=1600, pointsize=5)
  print(FeaturePlot(object = dataf, features = c("GAD1", "GAD2", "LHX6", "PVALB", "SCUBE3", "SST", "ADARB2", "LAMP5", "CXCL14", "CALB2", "VIP")))#
  dev.off()
  
  png(filename="Vln_Glia.png", width=2000, height=800, pointsize=1)
  print(VlnPlot(object = dataf, features = c("PLP1", "MBP", "CD74", "PDGFRA", "APOE", "GFAP", "SOX10", "OLIG1")))
  dev.off()
  png(filename="Vln_EN.png", width=2000, height=800, pointsize=1)
  print(VlnPlot(object = dataf, features = c("SLC17A7", "CUX2", "LAMP5", "COL5A2", "RORB", "GAL", "GABRG1", "RPRM", "PCP4", "NFIA", "THEMIS", "ROBO3", "FEZF2", "NR4A2", "NTNG2", "CTGF")))
  dev.off()
  png(filename="Vln_IN.png", width=2000, height=800, pointsize=1)
  print(VlnPlot(object = dataf, features = c("GAD1", "GAD2", "LHX6", "PVALB", "SCUBE3", "SST", "ADARB2", "LAMP5", "CXCL14", "CALB2", "VIP")))
  dev.off()
  
  png(filename="nUMI_UMAP.png", width=700, height=600, pointsize=5)
  print(FeaturePlot(object=dataf, features = "nCount_RNA"))#, cols.use = c("grey","blue"), pt.size=5, nCol=1))
  dev.off()
  
  cluster.ident <- table(dataf@active.ident)
  cluster.ident <- as.data.frame(cluster.ident, sep = ",", col.names=c("A","B"))
  colnames(cluster.ident) <- c("CLUSTER", "TOTAL CELLS")
  cluster.ident <- t(cluster.ident)
  
  avg.exp<-as.data.frame(AverageExpression(dataf, features = c("MBP","SLC17A7","GAD1","GFAP")))
  clusternum<-as.data.frame(cluster.ident)
  colnames(clusternum)<-colnames(avg.exp)
  
  binded<-rbind(avg.exp,clusternum)
  binded
  write.table(binded, "BindedCellsPerCluser.csv", sep = ",", col.names = NA)
  
  MBP<-as.data.frame(AverageExpression(datan, features = c("MBP")))
  SLC17A7<-as.data.frame(AverageExpression(datan, features = c("SLC17A7")))
  GAD1<-as.data.frame(AverageExpression(datan, features = c("GAD1")))
  binded2<-rbind(MBP,SLC17A7,GAD1)
  binded2<-t(binded2)
  row.names(binded2)<-paste(gsub('.dge.txt.gz', '', file2[i]), sep="")
  write.table(binded2, "Markers_average_exp.csv", sep = ",", col.names = NA)
  
  dataf.markers <- FindAllMarkers(dataf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  dataf.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  write.csv(dataf.markers, "ClusterMarkers.csv")
}
  



## Metrics only ##

for (i in 1:length(unique(file2))){
  setwd(datadir)
  #assign(paste(gsub('.dge.txt.gz', '', file2[i]), sep=""), read.table(file = file2[[i]], sep = "\t", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 3000))))
  assign("mt", read.table(file = file2[[i]], sep = "\t", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric", 3000))))
  
  outdir<-(paste(wdir,gsub('.dge.txt.gz', '', file2[i]), "_10012019", sep=""))
  dir.create(outdir)
  setwd(outdir)
  
  tag<-paste(gsub('.dge.txt.gz', '', file2[i]), sep="")
  ms<-CreateSeuratObject(mt)
  
  ## barcoderanks ##
  ms.out<-barcodeRanks(ms, lower=500)
  names(ms.out)
  head(ms.out)
  
  ms.out
  metadata(ms.out)$knee
  metadata(ms.out)$inflection
  
  threshod.inf<-(metadata(ms.out)$inflection)
  
  ## emptydrops ##
  me.out<-emptyDrops(ms@assays$RNA@counts, lower=0.5*threshod.inf, ignore=0.05*threshod.inf)
  is.cell <- me.out$FDR <= 0.01
  sum(is.cell, na.rm=TRUE)
  
  is.cell2<-is.cell
  is.cell2[is.na(is.cell2)]<-FALSE
  
  mf<-mt[,is.cell2]
  
  dataf<-CreateSeuratObject(mf)
  
  # Seurat analysis
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(object = dataf)), value = TRUE)
  percent.mt <- Matrix::colSums(GetAssayData(object = dataf,  slot = "counts")[mito.genes, ])/Matrix::colSums(GetAssayData(object = dataf, slot = "counts"))
  dataf$percent.mt <- percent.mt
  ribo.genes <- grep(pattern = "RP[SL][[:digit:]]", x = rownames(x = GetAssayData(object =dataf)), value = TRUE)
  percent.ribo <- Matrix::colSums(GetAssayData(object =dataf, slot = "counts")[ribo.genes, ])/Matrix::colSums(GetAssayData(object =dataf, slot = "counts"))
  dataf$percent.ribo <- percent.ribo
  
  # Metrics
  Totalcells<-length(dataf@meta.data$nFeature_RNA)
  cells250g<-length(dataf@meta.data$nFeature_RNA[dataf@meta.data$nFeature_RNA>250])
  cells100g<-length(dataf@meta.data$nFeature_RNA[dataf@meta.data$nFeature_RNA>100])
  pct250.100g<-(cells250g/cells100g*100)
  
  cells12kU<-length(dataf@meta.data$nCount_RNA[dataf@meta.data$nCount_RNA>12000])
  cells5kU<-length(dataf@meta.data$nCount_RNA[dataf@meta.data$nCount_RNA>5000])
  pct12kU<-(cells12kU/Totalcells*100)
  
  cells5mt<-length(dataf@meta.data$percent.mt[dataf@meta.data$percent.mt>0.05])
  pct5mt<-(cells5mt/Totalcells*100)
  
  f.med.g<-median(dataf@meta.data$nFeature_RNA)#[dataf@meta.data$nFeature_RNA>250])
  f.mean.g<-mean(dataf@meta.data$nFeature_RNA)#[dataf@meta.data$nFeature_RNA>250])
  
  f.med.U<-median(dataf@meta.data$nCount_RNA)#[dataf@meta.data$nFeature_RNA>250])
  f.mean.U<-mean(dataf@meta.data$nCount_RNA)#[dataf@meta.data$nFeature_RNA>250])
  
  combinedf<-data.frame(Totalcells,cells250g,cells100g,pct250.100g,cells12kU,cells5kU,pct12kU,cells5mt,pct5mt,
                        f.med.g,f.mean.g,f.med.U,f.mean.U)
  row.names(combinedf)<-gsub(".dge.txt.gz","",file2[i])
  write.csv(combinedf, "gene_metrics.csv")
}










###################
## FIND DOUBLETS ##   Not accurate
###################

annotations <- dataf@active.ident
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- dataf@meta.data$ClusteringResults
nExp_poi <- round(0.05*length(colnames(dataf)))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

dataf <- doubletFinder_v3(dataf, PCs = 1:40, pN = 0.25, pK = 0.5, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
dataf@meta.data[,"CellTypes_DF"] <- dataf@meta.data$DF.classifications_0.25_0.5_35

png(filename="Doublets.png", width=700, height=600, pointsize=5)
DimPlot(dataf, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("grey", "blue"))
dev.off()



dataf@meta.data$CellTypes_DF[which(dataf@meta.data$DF.classifications_0.25_0.5_52 == "Doublet")] <- 1
dataf@meta.data$CellTypes_DF[which(dataf@meta.data$DF.classifications_0.25_0.5_52 != "Doublet")] <- 2

data.names<-colnames(dataf)[which(dataf@meta.data$DF.classifications_0.25_0.5_35 != "Doublet")]

#dataf@active.ident <- as.factor(dataf@meta.data$CellTypes_DF)
#data.names2 <- WhichCells(dataf, idents = "Singlet")

datai <- SubsetData(dataf, cells=data.names)



