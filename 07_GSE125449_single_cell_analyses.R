
library(Seurat)
library(Matrix)

#scan the set1 and create first Seurat object
data1 <- Read10X(data.dir = "set1")
set1  <- CreateSeuratObject(counts = data1, min.cells = 3, project = "set1")
meta1<-read.table("set1/samples.txt",h=T,sep="\t")
row.names(meta1)<-meta1$Cell.Barcode
set1 <- AddMetaData(object = set1,  metadata = meta1)
head(set1[[]])



#scan the set2 and create second Seurat object
data2 <- Read10X(data.dir = "set2")
set2  <- CreateSeuratObject(counts = data2, min.cells = 3, project = "set2")
meta2<-read.table("set2/samples.txt",h=T,sep="\t")
row.names(meta2)<-meta2$Cell.Barcode
set2 <- AddMetaData(object = set2,  metadata = meta2)
head(set2[[]])



# merge seurat objets
all <- merge(set1, y = set2, project = "liver")

all
An object of class Seurat 
21287 features across 9946 samples within 1 assay 
Active assay: RNA (21287 features, 0 variable features)

#description of the object by their cell origin
list <- SplitObject(all, split.by = "orig.ident")



#normalization
list <- lapply(X = list, FUN = function(x) {
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", )
})
nfeatures = 2000

#find common anchors 
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) <- "integrated"




#scale data
combined <- ScaleData(combined, verbose = FALSE)
#dimensionnal reductions
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
ElbowPlot(combined,n=30)

combined <- RunUMAP(combined, reduction = "pca", dims = 1:15)



combined <- RunTSNE(combined, npcs = 15, verbose = FALSE)

DimPlot(combined, reduction = "umap",group.by ="orig.ident")

head(combined[[]])

DimPlot(combined, reduction = "umap",group.by ="Type")
ica<-subset(combined,subset = Diagnosis == "iCCA")

DimPlot(ica, reduction = "tsne",group.by ="ID", cols=cols25(19),pt.size=1)
DimPlot(ica, reduction = "tsne",group.by ="conca", cols=cols25(),pt.size=1)

library(viridis)
DimPlot(ica, reduction = "tsne",group.by ="conca", cols=inferno(16),pt.size=1.5)
library(pals)
DimPlot(ica, reduction = "tsne",group.by ="ident", cols=cols25(16),pt.size=1)
DimPlot(ica, reduction = "tsne",group.by ="conca", cols=cols25(19),pt.size=1.5)

features<-c("PLOD2","PLOD1","FH","MAT2B","NT5DC3","PDE6D","ALDOC")

DotPlot(ica, features = features,cols=inferno(8),split.by = "conca",assay = "RNA") + RotatedAxis()

x<-ica[[]]



library(Seurat)
library(viridis)
library(RColorBrewer)
library(ggplot2)



##compute metabolic score in seurat 
meta.score <- list(c("PLOD2","PLOD1","FH","MAT2B","NT5DC3","PDE6D","ALDOC"))

ica<- AddModuleScore(object = ica,
  features = meta.score, name = 'metabolism.score')

FeaturePlot(ica,"metabolism.score1",reduction = "tsne",pt.size=0.75,cols=c("grey90","darkblue"),min.cut="q9") 


VlnPlot(ica, features = c("metabolism.score1"), slot = "data", log = TRUE,pt.size=1,split.by="conca",col=cols25())

save(ica,file="ica.rda")





