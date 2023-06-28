library(Matrix)
library(dplyr)
library(Seurat)#, lib.loc="C:/Users/ogarmar3/Box/_Cobos Lab/Marcos/Seurat_BOX/Seurat-v2.3/") 
library(DropletUtils)
library(DoubletFinder)


ms<-CreateSeuratObject(D0202)


## barcoderanks ##

ms.out<-barcodeRanks(ms, lower=300)
names(ms.out)
head(ms.out)


png(filename="kneeplot.png", width=1200, height=1000, pointsize=15)
plot(ms.out$rank, ms.out$total, log="xy", xlab="Rank", ylab="Total")
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

D0202f<-D0202[,is.cell2]

mf<-CreateSeuratObject(D0202f)

##########################
## Seurat analysis here ##
##########################

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



