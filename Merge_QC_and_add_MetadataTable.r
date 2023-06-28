## merge 


## merge more than two Seurat objects
for (i in 1:number_Of_samples){
  a <- Read10X(paste(datadir, sample_name,sep=""))
  Oa <- CreateSeuratObject(counts = a, project = sample_name)
  Omerged10x <- merge(x = Omerged10x, y = Oa)
}

## add metadata table
## this is the current version, will change to a better one later

metadata.table <- read.csv("metadata_table.csv", header = TRUE)
metadata.table[is.na(metadata.table)]<-"NA"
mdtable <- metadata.table
mdsample<-seurat_object@meta.data
mdsample <- cbind(rownames(mdsample), mdsample)
rownames(mdsample) <- c()
sampleID <- metadata.table[,1]
mdsample$wpt <- as.character(mdsample$orig.ident)
mdtable$wpt <- as.character(mdtable$Unique.identifier)
A <- left_join(data.frame(mdsample, row.names=NULL), data.frame(mdtable, row.names=NULL), by = "wpt", all = TRUE)[-1]
B <- cbind(mdsample[,1], A)
cellcodes <- as.data.frame(seurat_object@assays$RNA@data@Dimnames[[2]])
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes
cellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", cellcodes$barcodes))
celllib <- as.data.frame(cellcodes$libcodes, row.names = row.names(cellcodes))
seurat_object <- AddMetaData(object = seurat_object, metadata = celllib, col.name = "libcodes")
for (i in 5:length(seurat_object@meta.data)){
  W <- names(seurat_object@meta.data[5])
  seurat_object@meta.data <- seurat_object@meta.data[, -which(colnames(seurat_object@meta.data) %in% W)]
}

for (i in 1:length(B[0,])){
  V <- B[0,][i]
  #V <- gsub(".y","",names(V))
  cellcodes[,names(V)] <- as.vector(B[,i])
  seurat_object <- AddMetaData(object = seurat_object, metadata = cellcodes[,names(V)], col.name = names(V))#col.name = gsub(".y","",names(V))
}
 
 
## QC
q.nFeat <- quantile(Eonly_10x@meta.data$nFeature_RNA, .98)
q.nCount <- quantile(Eonly_10x@meta.data$nCount_RNA, .98)
q.mt <- quantile(Eonly_10x@meta.data$percent.mt, .98)
q.nFeat; q.nCount; q.mt
Eonly_10x <- subset(Eonly_10x, subset = nFeature_RNA > 450 & nFeature_RNA < q.nFeat & percent.mt < q.mt & nCount_RNA > 675 & nCount_RNA < q.nCount)
