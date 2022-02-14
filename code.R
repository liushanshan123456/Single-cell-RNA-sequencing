install.packages("Seurat")
install.packages("outliers")
install.packages("pbmcapply")
install.packages("doFuture")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("singscore")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("limma")

install.packages("devtools")
library(devtools)
devtools::install_github('dviraran/SingleR')



#2.1 Data download, filtering and standardization
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)

setwd("C:\\Users\\LSS\\Desktop\\µ¥Ï¸°û")          

rt=read.table("geneMatrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)


pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 2, min.features = 10, names.delim = "_",)
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="Fig.1 abc.pdf",width=10,height=6)         
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 8000 & percent.mt < 15)    


pdf(file="Fig.1 de.pdf",width=11,height=6)             
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.1)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.1)
CombinePlots(plots = list(plot1, plot2))
dev.off()


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2500)
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="Fig.1 fg.pdf",width=12,height=6)            
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()





#2.2 PCA (principal component analysis)
pbmc=ScaleData(pbmc)                    
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))   


pdf(file="Fig.2 a.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()


pdf(file="Fig.2 b.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()


pdf(file="Fig.2 c.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()


pbmc <- JackStraw(object = pbmc, num.replicate = 45)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:11)
pdf(file="Fig.2 d.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:11)
dev.off()





#2.3 TSNE cluster analysis and marker gene
pcSelect=6
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)               
pbmc <- FindClusters(object = pbmc, resolution = 1)                  
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)                      
pdf(file="Fig.3 a.pdf",width=6,height=6)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 2, label = TRUE)   
dev.off()
write.table(pbmc$seurat_clusters,file="tsneCluster.txt",quote=F,sep="\t",col.names=F)

logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="markers.xls",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf(file="Fig.3 b.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()


pdf(file="Fig.3 c.pdf",width=12,height=3)
VlnPlot(object = pbmc, features = c("SPP1", "CD36","CD44", "AREG"),ncol = 4)
dev.off()


pdf(file="Fig.3 d.pdf",width=12,height=3)
FeaturePlot(object = pbmc, features = c("SPP1", "CD36","CD44", "AREG"),cols = c("yellow", "red"),ncol = 4)
dev.off()


pdf(file="Fig.3 e.pdf",width=12,height=6)
cluster10Marker=c("GNAT3", "PAGE2", "PRAME", "HGF", "SSX1", "MAGEA6", "MAEL", "DCBLD2", "EREG", "MMP1")
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()





#2.4 Cell type annotation and trajectory analysis
library(SingleR)
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
singler = CreateSinglerObject(counts, annot = ann, "pbmc", min.genes = 0,
  species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = F, do.signatures = T, clusters = clusters, do.main.types = T,
  reduce.file.size = T, numCores = 1)
singler$seurat = pbmc
singler$meta.data$xy = pbmc@reductions$tsne@cell.embeddings
clusterAnn=singler$singler[[2]]$SingleR.clusters.main$labels
write.table(clusterAnn,file="clusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(singler$other,file="cellAnn.txt",quote=F,sep="\t",col.names=F)


monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(singler$other,file="monocleClusterAnn.txt",quote=F,sep="\t",col.names=F)
write.table(sig.markers,file="0monocleMarkers.txt",sep="\t",row.names=F,quote=F)

