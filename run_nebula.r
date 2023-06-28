library(ggplot2)
library(Seurat)
#library(GGally)
library(GSEABase)
library(limma)
library(reshape2)
library(data.table)
library(knitr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(NMF)
library(rsvd)
library(RColorBrewer)
library(MAST)
library(nebula)
library(Seurat)
library(SeuratData)
library(Matrix)

obj <- readRDS("/oak/stanford/groups/icobos/TAU_Manuscript_Oct2020/Inhibitory_Subset_Annotated/3wayInClusters12dim03res_Annotated.rds")
objSub <- obj[, sample(colnames(obj), size = 3000)]
#data_g = group_cell(count=sample_data$count,id=sample_data$sid,pred=df)

#"percent.mt"
#"PMI.hr."               
#"Gender"
#"nCount_RNA"
#"ApoE"
#"RIN"
#"Age"
#"Sequencing.batch"
dataf <- data.frame(SORT = objSub$SORT,  # Create example data
                   Age = objSub$Age,
                   Gender = objSub$Gender,
                   PMI = objSub$PMI.hr.,
                   RIN = objSub$RIN,
                   APOE = objSub$ApoE,
                   CCount = objSub$nCount_RNA,
                   Pmt = objSub$percent.mt,
                   Sbatch = objSub$Sequencing.batch)

head(dataf)
df = model.matrix(~SORT+Age+Gender+PMI+RIN+APOE+CCount+Pmt+Sbatch, data=dataf)
data_g = group_cell(count=objSub@assays$RNA@counts,id=objSub$Unique.identifier,pred=df)
re = nebula(data_g$count, data_g$id, pred=df)
print(re)
