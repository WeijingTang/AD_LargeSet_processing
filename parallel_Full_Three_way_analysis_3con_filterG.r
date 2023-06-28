library(Matrix)
library('vctrs',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library('Seurat',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library('dplyr',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
require('plyr', lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
#library('devtools',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library('ggplot2',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
#library('Matrix',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
#library('biomaRt',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library('plotly',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
#library("reshape",lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library(reticulate,lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
#library(umap,lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library(future,lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library(future.apply,lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library(cowplot,lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library(gridExtra,lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library("farver",lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat)
#install.packages('glmmTMB',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
library('glmmTMB',lib="/oak/stanford/groups/icobos/Single_Cell_data/Analyzed_04122020/test/")
require(gdata)
plan("multiprocess", workers = 16)
options(future.globals.maxSize = Inf)
date <- format(Sys.Date(),format='%m.%d.%y')

#TclusterSet <- readRDS("/home/groups/icobos/Weijing/Differential_Gene_Expression/testdata/C6_Three_conditions_2020_CLUSTERSET.rds")
obj <- readRDS("/oak/stanford/groups/icobos/Sam/Data/TAU_Manuscript_Final/Debia_3way/Inma_Clustering/MultiAssay10dim05res_Excitatory/Export_Import_3way_Integration/3wayClusters10dim05res_imported_Integrated_Annotated.rds")

DefaultAssay(obj) <- "Exon"
#Idents(obj) <- obj$active.ident
sepobj <- obj
DefaultAssay(sepobj) <- "Exon"
sepobj <- NormalizeData(sepobj, verbose = TRUE)

TclusterSet = sepobj@meta.data$seurat_clusters
TclusterSet = as.numeric(as.character(TclusterSet))
names(TclusterSet) = colnames(sepobj)

Idents(sepobj) <- sepobj$SORT
objMAP2 <- subset(sepobj, idents = c("MAP2"), invert = FALSE)
objAT8 <- subset(sepobj, idents = c("AT8"), invert = FALSE)
objctr <- subset(sepobj, idents = c("MAP2control"), invert = FALSE)

for(i in 1:length(TclusterSet))
{
  if(names(TclusterSet)[i] %in% colnames(objMAP2))
  {
    TclusterSet[i] = 0
  }
  if(names(TclusterSet)[i] %in% colnames(objAT8))
  {
    TclusterSet[i] = 1
  }
  if(names(TclusterSet)[i] %in% colnames(objctr))
  {
    TclusterSet[i] = 2
  }
}
saveRDS(TclusterSet, file = '/home/groups/icobos/Weijing/Differential_Gene_Expression/testdata/DGE_C6_test/Full_TclusterSet_Three_conditions_2020_filtered.rds')
#Idents(obj) <- obj$seurat_clusters
#sepobj <- subset(obj, idents = c("6"), invert = FALSE)
sepobj <- obj
DefaultAssay(sepobj) <- "Exon"
sepobj <- NormalizeData(sepobj, verbose = TRUE)

mx <- GetAssayData(sepobj, slot='counts')
dim(mx)
mx1 <- mx[1:16000,]
mx1 <- as.matrix(mx1)
Trowfilter1 <- apply(mx1,1,function(x){ if(sum(x == 0) < (length(x)-0.2*length(x))){ return(T)} else{ return(F) }}) 
mx1 <- mx1[Trowfilter1, ]
mx2 <- mx[16001:33538,]
mx2 <- as.matrix(mx2)

Trowfilter2 <- apply(mx2,1,function(x){ if(sum(x == 0) < (length(x)-0.2*length(x))){ return(T)} else{ return(F) }})
mx2 <- mx2[Trowfilter2, ]
Tdata <- rbind(mx1,mx2)

#Tdata <- as.matrix(GetAssayData(sepobj, slot='counts'))
print("check_where")

save(Tdata, file = '/home/groups/icobos/Weijing/Differential_Gene_Expression/testdata/Full_Three_conditions_2020_filtered.Rdata')
TclusterSet = TclusterSet[colnames(Tdata)]
Tsubjs <- unlist(lapply(colnames(Tdata),function(x){ return(unlist(strsplit(x,'_'))[1]) }))

mttable <- read.csv("/home/groups/icobos/Weijing/Differential_Gene_Expression/Metadatatable_PremRNA_3wayIC08172020.csv")
key <- mttable
colnames(key)[1] <- 'Sample'
TDAssay <- key[,c('Sample','Assay')]
TADstatus <- key[,c('Sample','NP.Diagnosis')]
TADBraakstatus <- key[,c('Sample','Braak.stage')]
#OTstatus <- key[,c('Sample','Final_Status')]
TAODstatus <- key[,c('Sample','Age')]
TSEXstatus <- key[,c('Sample','Gender')]
TMS4Astatus<- key[,c('Sample','ApoE')]
TPMIstatus <- key[,c('Sample','PMI.hr.')]
#cellCounts <- key[,c('Sample','count_clusterMicro')]
TAmyloidstatus <- key[,c('Sample','Amyloid')]
TSbatchstatus <- key[,c('Sample','SORT.BATCH')]
#TADstatus <- key[,c('Sample','SORT')]
Tnm <- TADstatus[,1]

cellCount <- as.data.frame(table(Tsubjs))
cc <- as.character(cellCount[,-1])
cc <- as.numeric(cc) - (mean(as.numeric(cc)))#center around the mean
cc <- scale(cc)
names(cc) <- Tnm
#TAODstatus <- as.character(scale(as.numeric(TAODstatus)))
cCount <- unlist(lapply(Tsubjs,function(x){ return(cc[x]) }))
cCount <- as.numeric(cCount)

#apply subject SEX status to the cells

TSEXstatus <- as.character(TSEXstatus[,-1])
names(TSEXstatus) <- Tnm
TSEXstatus[TSEXstatus == 'M'] <- 0
TSEXstatus[TSEXstatus == 'F'] <- 1
SEX <- unlist(lapply(Tsubjs,function(x){ return(TSEXstatus[x]) }))
SEX<- as.numeric(SEX)
#apply subject AOD status to the cells

TAODstatus <- as.character(TAODstatus[,-1])

TAODstatus <- as.numeric(TAODstatus) - (mean(as.numeric(TAODstatus)))#center around the mean
TAODstatus <- scale(TAODstatus)
names(TAODstatus) <- Tnm
#TAODstatus <- as.character(scale(as.numeric(TAODstatus)))
AOD <- unlist(lapply(Tsubjs,function(x){ return(TAODstatus[x]) }))
AOD <- as.numeric(AOD)



Tsubjs <- as.numeric(factor(Tsubjs))
TclusterSet = as.numeric(TclusterSet)

sub1 = '0'
sub2 = '1'
Tmainfilename <- paste('/home/groups/icobos/Weijing/Differential_Gene_Expression/testdata/DGE_C6_test/parallel_Full_ThreeCondition_Gender_Age_Ccount_gene_1018filter',sub1,sub2, sep='_')
outfile <- paste(Tmainfilename,date,sep='_')

write('Gene\tEstimate\tStd.Error\tz.Value\tPr(>|z|)\n',file=paste0(outfile,'.out')) #write the header line
results <- future_apply(cbind(seq_len(nrow(Tdata)),Tdata),1,function(x){
  y <- x[-1]#the gene row position is passed so it needs to be removed
  re = NA
  re <- tryCatch({
    re <- glmmTMB(y~ TclusterSet + SEX + AOD + cCount + (1|Tsubjs),ziformula=~1,family=nbinom2)
  },warning = function(warn){return(NA)},error = function(Err){return(NA)},finally={
    ##write the effects and pvalue from the carrier variable to a file as it goes
    if (!is.na(re[1])){
      carrierResult <- summary(re)$coefficients$cond[2,]
      #print(paste0(rownames(Tdata),'\t',paste(carrierResult,collapse='\t')))
      write(paste0(rownames(Tdata)[x[1]],'\t',paste(carrierResult,collapse='\t')), file=paste0(outfile,'.out'), append=T)
    }else{
      #if there is a convergence error keep it in the output to show we tested the gene
      write(paste0(rownames(Tdata)[x[1]],'_ConvergenceERROR','\t',paste(rep(1.00,4),collapse='\t')), file=paste0(outfile,'.out'), append=T)}})
  return(NA)
})


