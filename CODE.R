####1. single cell analysis of all sample ####
####merge the scRNA-data####
rm(list=ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
#BiocManager::install("clusterProfiler")
#BiocManager::install("msigdbr")

"plan(""multicore"", workers = 12) ###set the compute core"
options(future.globals.maxSize = 10000 * 1024^2)
getwd()

setwd("F:/IPF1/Single Cell")
###########
gc()
"sce.mergeTEN1<- SplitObject(sce.mergeTEN, split.by = ""orig.ident"") "
"sce.mergeTEN2<- merge(x=sce.mergeTEN1[[1]],y=c(sce.mergeTEN1[2],sce.mergeTEN1[3],sce.mergeTEN1[4],"
"                      sce.mergeTEN1[5],sce.mergeTEN1[18],sce.mergeTEN1[19],sce.mergeTEN1[20],"
"                      sce.mergeTEN1[21],sce.mergeTEN1[22],sce.mergeTEN1[23],sce.mergeTEN1[24]))"
meta.data1 <- sce.mergeTEN2@meta.data
"write.csv(meta.data1,""./control.csv"")"
"saveRDS(sce.mergeTEN2,file=""./sce.mergeTEN_control.rds"")"

"sce.mergeTEN3 <- merge(x=sce.mergeTEN1[[6]],y=c(sce.mergeTEN1[10],sce.mergeTEN1[11],sce.mergeTEN1[7],"
"                                               sce.mergeTEN1[8],sce.mergeTEN1[9],sce.mergeTEN1[14],sce.mergeTEN1[15],"
"                                               sce.mergeTEN1[26],sce.mergeTEN1[32],sce.mergeTEN1[33],sce.mergeTEN1[34],"
"                                               sce.mergeTEN1[35],sce.mergeTEN1[36],sce.mergeTEN1[37],sce.mergeTEN1[40],"
"                                               sce.mergeTEN1[41],sce.mergeTEN1[42],sce.mergeTEN1[43]))"
meta.data2 <- sce.mergeTEN3@meta.data
"write.csv(meta.data2,""./ipf.csv"")"
"saveRDS(sce.mergeTEN3,file=""./sce.mergeTEN_ipf.rds"")"

"rm(sce.mergeTEN1,meta.data1,meta.data2)"
sce.mergeTEN2$tissue_type <- "Control"
sce.mergeTEN3$tissue_type <- "IPF"
sce.mergeTEN2@meta.data

"sce.mergeTEN1<- merge(x=sce.mergeTEN2,y=sce.mergeTEN3,project = ""scTEN"")"
sce.mergeTEN1@meta.data
"rm(sce.mergeTEN,sce.mergeTEN2,sce.mergeTEN3)"
"saveRDS(sce.mergeTEN1,file=""./sce.mergeTEN.rds"")"
sce.mergeTEN <- sce.mergeTEN1
rm(sce.mergeTEN1)

sce.mergeTEN<-readRDS(file="./sce.mergeTEN.rds")
####integration with harmony####
library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)
gc()
"sce.mergeTEN@meta.data[1:5,]"


"sce.mergeTEN[[""percent.mt""]] <- PercentageFeatureSet(sce.mergeTEN, pattern = ""^MT-"")"
hist(sce.mergeTEN[["percent.mt"]]$percent.mt)

"VlnPlot(sce.mergeTEN, features = c(""nFeature_RNA"", ""nCount_RNA"", ""percent.mt""), ncol = 3)"
sce.mergeTEN
"plot1 <- FeatureScatter(sce.mergeTEN, feature1 = ""nCount_RNA"", feature2 = ""percent.mt"")"
"plot2 <- FeatureScatter(sce.mergeTEN, feature1 = ""nCount_RNA"", feature2 = ""nFeature_RNA"")"
"CombinePlots(plots = list(plot1,plot2))"
sce.mergeTEN

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
"sce.mergeTEN <- CellCycleScoring(sce.mergeTEN, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)"
"sce.mergeTEN@meta.data[1:5,]"
"sce.mergeTEN<-NormalizeData(sce.mergeTEN,verbose = T) "
"sce.mergeTEN<-FindVariableFeatures(sce.mergeTEN,selection.method = ""vst"", nfeatures = 2000)"
"sce.mergeTEN<-ScaleData(sce.mergeTEN,vars.to.regress = c(""percent.mt"",""S.Score"",""G2M.Score""),verbose = FALSE)"
"sce.mergeTEN<-RunPCA(sce.mergeTEN,verbose = T,npcs = 50)"
"ElbowPlot(sce.mergeTEN,ndims = 50)"

"p1 <- DimPlot(object = sce.mergeTEN, reduction = ""pca"", pt.size = .1, group.by = ""orig.ident"",raster=FALSE)"
"p2 <- VlnPlot(object = sce.mergeTEN, features = ""PC_1"", group.by = ""orig.ident"", pt.size = .1,raster=FALSE)"
"CombinePlots(plots=list(p1,p2))"

"sce.mergeTEN<-RunHarmony(sce.mergeTEN,group.by.vars = c(""orig.ident""), plot_convergence = TRUE)"
"harmony_embeddings <- Embeddings(sce.mergeTEN, 'harmony')"
dim(harmony_embeddings)

"p3 <- DimPlot(object = sce.mergeTEN, reduction = ""harmony"", pt.size = .1, group.by = ""orig.ident"",raster=FALSE)"
"p4 <- VlnPlot(object = sce.mergeTEN, features = ""harmony_1"", group.by = ""orig.ident"", pt.size = .1,raster=FALSE)"
"CombinePlots(plots=list(p3,p4))"
p1+p2+p3+p4
#跑ump和tsne
sce.mergeTEN <- sce.mergeTEN %>% 
"  RunUMAP(reduction = ""harmony"", dims = 1:50) %>% "
"  RunTSNE(reduction = ""harmony"", dims = 1:50) %>%"
"  FindNeighbors(reduction = ""harmony"", dims = 1:50)"

"sce.mergeTEN<-FindClusters(sce.mergeTEN,resolution = 0.5)#大样本想要分得更细，resolution可以调的大一点"
table(Idents(sce.mergeTEN))
#Idents(sce_inter)<-sce_inter$seurat_clusters
"sce.mergeTEN@meta.data[1:5,]"
metadata <- sce.mergeTEN@meta.data

"p1<-DimPlot(sce.mergeTEN,reduction = ""tsne"",label = T,raster=FALSE)"
p1
"DimPlot(sce.mergeTEN,reduction = ""umap"",label = T,raster=FALSE)"
"DimPlot(sce.mergeTEN,reduction = ""tsne"",label = T,split.by = ""tissue_type"",raster=FALSE)"
"DimPlot(sce.mergeTEN,reduction = ""umap"",label = T,split.by = ""tissue_type"")"
"sce.mergeTEN@meta.data[1:5,]"

####marker####
#epithelial cells
"FeaturePlot(sce.mergeTEN,features = c(""EPCAM"",""SCGB3A2"",""KRT5"",""FOXJ1"",""SFTPB""),reduction = ""tsne"",label = T)"
#immune cells
"FeaturePlot(sce.mergeTEN,features = c(""PTPRC"",""CD79A"",""IGHG1"",""CPA3""),reduction = ""tsne"",label = T)"
"FeaturePlot(sce.mergeTEN,features = c(""LYZ"",""S100A9"",""CCL5"",""C1QA""),reduction = ""tsne"",label = T)"
"FeaturePlot(sce.mergeTEN,features = c(""APOC1"",""NKG7"",""CD3E"",""FCER1A""),reduction = ""tsne"",label = T)"
#immune cells
"FeaturePlot(sce.mergeTEN,features = c(""PTPRC"",""LYZ"",""S100A9"",""APOC1"",""NKG7""),reduction = ""tsne"",label = T)"

#fibroblasts markers
"FeaturePlot(sce.mergeTEN,features = c(""LUM"", ""DCN"", ""TAGLN"",""ACTA2"",""THY1""),reduction = ""tsne"",label = T)"
#Endothelial markers
"FeaturePlot(sce.mergeTEN,features = c(""PECAM1"",""VWF"", ""VIPR1"",""EMCN"",""ACKR1""),reduction = ""tsne"",label = T)"

"sce.mergeTEN<-RenameIdents(sce.mergeTEN,""2""=""Epithelial cells"",""4""=""Epithelial cells"",""7""=""Epithelial cells"",""10""=""Epithelial cells"","
"                                             ""13""=""Epithelial cells"",""14""=""Epithelial cells"",""22""=""Epithelial cells"",""23""=""Epithelial cells"","
"                                             ""27""=""Epithelial cells"",""28""=""Epithelial cells"","
"                                             ""0""=""Immune cells"",""1""=""Immune cells"",""3""=""Immune cells"",""6""=""Immune cells"",""8""=""Immune cells"","
"                                             ""9""=""Immune cells"",""12""=""Immune cells"",""17""=""Immune cells"",""18""=""Immune cells"",""19""=""Immune cells"","
"                                             ""20""=""Immune cells"",""24""=""Immune cells"",""26""=""Immune cells"",""29""=""Immune cells"","
"                                             ""11""=""fibroblasts"",""15""=""fibroblasts"","
"                                             ""16""=""Endothelial cells"",""5""=""Endothelial cells"",""21""=""Endothelial cells"",""25""=""Unknown"")"

"DimPlot(sce.mergeTEN,reduction = ""tsne"",label = F)"
"DimPlot(sce.mergeTEN,reduction = ""umap"",label = T)"

sce.mergeTEN$Seurat_harmony<-Idents(sce.mergeTEN)
"sce.mergeTEN@meta.data[1:5,]"

"markers<-FindAllMarkers(sce.mergeTEN,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)"
"write.csv(markers,file=""./markers.csv"")#"

"saveRDS(sce.mergeTEN,file=""./sce.mergeTEN_markered.rds"")"

#####draw figure####
#BiocManager::install("ggsci")
library(ggsci)
##define the color
??ggsci

#pal_npg()(10)
cors<-pal_igv()(11)
# <- pal_aaas()(10)
"#cors <- c(cors,""grey"")"
##Reduction plot##
"DimPlot(sce.mergeTEN,reduction=""tsne"",cols=cors)"
"DimPlot(sce.mergeTEN,reduction=""umap"",cols=cors)"

####draw figure####
#sample
"DimPlot(sce.mergeTEN,reduction=""tsne"",cols = cors,split.by = ""orig.ident"",ncol = 3)"
"DimPlot(sce.mergeTEN,reduction=""umap"",cols = cors,split.by = ""orig.ident"",ncol = 3)"
#type
sce.mergeTEN$tissue_type
"DimPlot(sce.mergeTEN,reduction=""tsne"",cols = cors,split.by = ""tissue_type"")"
"#DimPlot(sce_inter,reduction=""umap"",cols = cors,label = F)"
"#DimPlot(sce_inter,reduction=""umap"",cols = cors,split.by = ""orig.ident"",ncol = 3)"

"#markers <- FindAllMarkers(sce.mergeTEN, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)"
"write.csv(markers,file=""./markers.csv"")"

"features = c(""EPCAM"",""SCGB3A2"",""KRT5"",""FOXJ1"",""SFTPB"","
"             ""PTPRC"",""LYZ"",""S100A9"",""APOC1"",""NKG7"","
"             ""LUM"", ""DCN"", ""TAGLN"",""ACTA2"",""THY1"","
"             ""PECAM1"",""VWF"", ""VIPR1"",""EMCN"",""ACKR1"")"
"features = c(""CXCL14"",""TM4SF1"",""CYTL1"",""SOD3"",""MMP10"")"

length(features)
"VlnPlot(sce.mergeTEN,features = features[1:5],cols = cors,pt.size = 0)"
"VlnPlot(sce.mergeTEN,features = features[6:10],cols = cors,pt.size = 0)"
"VlnPlot(sce.mergeTEN,features = features[11:15],cols = cors,pt.size = 0)"
"VlnPlot(sce.mergeTEN,features = features[16:20],cols = cors,pt.size = 0)"

"VlnPlot(sce.mergeTEN,features = c(""CD14"",""VEGFA""),cols = cors,pt.size = 0)"

"DotPlot(sce.mergeTEN,features = features,cols = c(""grey"",""darkblue""))"

####different cells in all samples####
library(ggplot2)
"table(sce.mergeTEN$tissue_type,sce.mergeTEN$orig.ident)"

"tab<-table(Idents(sce.mergeTEN),sce.mergeTEN$orig.ident)"
tab

tab<-as.data.frame(tab)
tab
"ggplot(tab,aes(x=Var2,y=Freq,fill=Var1))+ "
"  geom_bar(stat='identity',position='stack',alpha=.5)+ "
"  labs(title='',x='group',y='Cell Number')+ "
"  theme(legend.justification = 'right', "
"        legend.position = 'right', "
"        legend.key.height = unit(0.1,'cm'),"
"        panel.background = element_blank(),"
"        axis.line=element_line(size=0.5,colour=""black"")"
  )+ scale_fill_manual(values=cors)

"tab<-table(Idents(sce.mergeTEN),sce.mergeTEN$orig.ident)"
tab
"tab<-prop.table(tab,2)*100"
tab
tab<-as.data.frame(tab)
"ggplot(tab,aes(x=Var2,y=Freq,fill=Var1))+ "
"  geom_bar(stat='identity',position='stack',alpha=.5)+ "
"  labs(title='',x='group',y='Cell Proportion (%)')+ "
"  theme(legend.justification = 'right', "
"        legend.position = 'right', "
"        legend.key.height = unit(0.1,'cm'),"
"        panel.background = element_blank(),"
"        axis.line=element_line(size=0.5,colour=""black"")"
  )+ scale_fill_manual(values=cors)


####figure####
#####Doheatmap#####
"features<-unique(markers$gene[markers$p_val_adj<0.0001&markers$avg_log2FC>1.5]) #通常p_val_adj<0.05,avg_log2FC>0.25"
features
getwd()

"#sce.mergeTEN<-ScaleData(sce.mergeTEN, vars.to.regress = c(""percent.mt"",""S.Score"",""G2M.Score""))"
"sce.mergeTEN_DR<- ScaleData(sce.mergeTEN, vars.to.regress = c(""percent.mt"",""S.Score"",""G2M.Score""),features = features)"

"pdf(file=""DO_heatmap.pdf"",width = 8,height = 10)"
"DoHeatmap(sce.mergeTEN_DR,features = features,group.colors = cors)#It take several mins"
dev.off()

#####heatmap of marker gene#####
#cors
####matrix
#genemeanMatrix<-AverageExpression(sce.mergeTEN)

"#features_to_anno = c(""PSCA"",""CD24"",""UPK2"",""CD3D"",""CD8A"",""CD79A"",""JCHAIN"",""CD74"",""CD14"",""FCGR3A"",""VCAN"",""CPA3"",""COL1A1"",""LUM"",""VWF"",""PECAM1"",""RGS5"",""PDGFRB"",""S100B"",""PLP1"")"
#length(features_to_anno)
"#geneMatrix<-genemeanMatrix[features_to_anno,]"
#rownames
#cc<-rownames(geneMatrix)
#cc
#cc[which(!(cc%in%features_to_anno))]<-""
#cc
#table(Idents(sce.mergeTEN))
#levels(Idents(sce.mergeTEN))

#annotation_col = data.frame(
"#  CellType = c(""Epithelial"", ""T/NK"",""B"",""Plasma_B"",""Myeloid"",""MAST"",""Fibroblast"",""Endothelial"", ""Pericytes"", ""Neuron"",""Doublets"")"
#)
#rownames(annotation_col) = colnames(geneMatrix)
#define the seq
"#annotation_col$CellType<-factor(annotation_col$CellType,levels= c(""Epithelial"", ""T/NK"",""B"",""Plasma_B"",""Myeloid"",""MAST"",""Fibroblast"",""Endothelial"", ""Pericytes"", ""Neuron"",""Doublets""))"
#levels(Idents(sce.mergeTEN))
#cors
"#paste(c(""Epithelial"", ""T/NK"",""B"",""Plasma_B"",""Myeloid"",""MAST"",""Fibroblast"",""Endothelial"", ""Pericytes"", ""Neuron"",""Doublets""),cors,sep='""=""')"
#color to use in heatmap
#ann_colors = list(
"#  CellType = c(""Epithelial""=""#5050FFFF"", ""T/NK""=""#CE3D32FF"",        ""B""=""#749B58FF"",           ""Plasma_B""=""#F0E685FF"","
"#               ""Myeloid""=""#466983FF"",     ""MAST""=""#BA6338FF"",        ""Fibroblast""=""#5DB1DDFF"",  ""Endothelial""=""#802268FF"","
"#               ""Pericytes""=""#6BD76BFF"",   ""Neuron""=""#D595A7FF"",      ""Doublets""=""#924822FF"" )"
"  #GeneClass = c(""CD4+"" = ""#7570B3"", ""CD8+"" = ""#E7298A"",""NK"" = ""#66A61E"")"
#)
#?pheatmap
#library(RColorBrewer)
#dev.off()
#dev.new()
"#ComplexHeatmap::pheatmap(as.matrix(geneMatrix), annotation_col = annotation_col,cluster_rows = F,labels_row  = cc,"
"#                         annotation_colors = ann_colors,cluster_cols = F,scale = ""row"",color = colorRampPalette(rev(brewer.pal(n = 20, name =""RdYlBu"")))(50))"
"#gaps_row = c(8, 14) )"



####2. single cell analysis of fibroblasts####
library(Seurat)
library(dyno)
library(monocle3)
library(ggplot2)
library(dplyr)
library(CellChat)
#devtools::install_github("sqjin/CellChat")
"#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils'))"
#install.packages("devtools")
#library(speedglm)
#setwd("C:/Program Files/R/R-4.2.3/library/speedglm_0.3-4")
#devtools::load_all("C:/Program Files/R/R-4.2.3/library/speedglm")
"#devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')"
"#devtools::install_github('cole-trapnell-lab/leidenbase',force = TRUE)"
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github("dynverse/dyno")
seurat.obj <- readRDS("sce.mergeTEN.rds")
seurat.obj <- sce.mergeTEN
"rm(plot1,plot2)"
setwd("F:/IPF1/Single Cell")
#subcluster
"fibroblast <- subset(seurat.obj,Seurat_harmony==""fibroblasts"")"
"fibroblast <- FindVariableFeatures(fibroblast, selection.method = ""vst"", nfeatures = 2000)"
"fibroblast <- RunPCA(fibroblast, features = VariableFeatures(object = fibroblast))"
"ElbowPlot(fibroblast,ndims = 20)#有批次效应时，批次效应选的尽量大"
"fibroblast <- FindNeighbors(fibroblast, dims = 1:20)"
"fibroblast <- FindClusters(fibroblast, resolution = 0.5)"
"fibroblast <- RunUMAP(fibroblast, dims = 1:20)"
"fibroblast <- RunTSNE(fibroblast, dims = 1:20)"

"DimPlot(fibroblast, reduction = ""tsne"")"
"DimPlot(fibroblast, reduction = ""tsne"",split.by = ""tissue_type"",label = F,raster=FALSE)"

"#0  ""MFAP5"",""PCOLCE2"",""PI16"""
"# Smooth Muscle Cells (SMCs) 6,7,8,0"
"FeaturePlot(fibroblast,features = c(""ACTA2"", ""PDGFRB"", ""MYH11"", ""TAGLN""),reduction = ""tsne"",label = T)"
"#9 Mesothelial ,9"
"#FeaturePlot(fibroblast,features = c(""MSLN"", ""UPK3B"", ""HP"", ""WT1""),reduction = ""tsne"",label = T)"
#1 "F3_related " 
"FeaturePlot(fibroblast,features = c(""ERRFI1"", ""SPSB1"",""PDPN"", ""PLA2G2A"",""F3""),reduction = ""tsne"",label = T)"
#3.4 "ROBO2_related"
"FeaturePlot(fibroblast,features = c(""ROBO2"" ,""BMP5"",""ASPN"",""POSTN""),reduction = ""tsne"",label = T)"
#2.5 undetermined
"FeaturePlot(fibroblast,features = c(""MALAT1"" , ""SPARC"",""SFRP2"",""CFD"" ),reduction = ""tsne"",label = T)"
"#1,6  ""BMP5"",""LTBP2"",""TCF21"
"#FeaturePlot(fibroblast,features = c(""BMP5"",""LBH"",""NPNT"",""TCF21""),reduction = ""tsne"",label = T)"
"#2,5  ""MYH11"",""PCSK1N"",""WIF1"" COLLAGEN"
"#FeaturePlot(fibroblast,features = c(""COL11A1"",""COL10A1"",""COL1A2"",""FN1""),reduction = ""tsne"",label = T)"
"#FeaturePlot(fibroblast,features = c(""FN1"",""COL1A2"",""TPM1"",""SPARC""),reduction = ""umap"",label = T)"
"#3  ""IGHG4"",""IGLC3"",""IGHG1"",""IGHG3"" ICAF"
"#FeaturePlot(fibroblast,features = c(""IGHG4"",""IGHG1"",""IGHG3"",""IGLC3""),reduction = ""tsne"",label = F)"
"#4,8  ""MYH11"",""ACTA2"",""TINAGL1"",""ACTG2"",""DES"" MYCAF"
"#FeaturePlot(fibroblast,features = c(""TINAGL1"",""MYH11"",""ACTG2"",""TAGLN""),reduction = ""tsne"",label = T)"
"#5   ""POSTN"",""COL11A1"",""MMP11"""
"#FeaturePlot(fibroblast,features = c(""POSTN"",""MMP11"",""CTHRC1""),reduction = ""tsne"",label = F)"
"#6  ""FGFR4"",""FIGF"",""TCF21"" "
"#FeaturePlot(fibroblast,features = c(""FGFR4"",""FIGF"",""RGCC"",""LIMCH1"",""TCF21""),reduction = ""umap"",label = T)"
"#7  ""MYC"",""ADAMTS1"",""IL6"" INCAT"
"#FeaturePlot(fibroblast,features = c(""MYC"",""IL6"",""CXCL2"",""NAMPT""),reduction = ""tsne"",label = T)"
"#8  ""MYH11"",""PCSK1N"",""WIF1"""
"#FeaturePlot(fibroblast,features = c(""WIF1"",""FAM150A"",""PCSK1N"",""SCX"",""MYH11""),reduction = ""umap"",label = T)"
"#2,5  ""MYH11"",""PCSK1N"",""WIF1"""
"#FeaturePlot(fibroblast,features = c(""COL11A1"",""COL8A1"",""COL10A1"",""CTHRC1""),reduction = ""umap"",label = T)"
#markers
n_clust <- names(table(fibroblast@active.ident)[table(fibroblast@active.ident)!=0])
marker.list <- list()
for(i in n_clust){
  ident1 <- i
"  markers <- FindMarkers(fibroblast,"
"                         ident.1 = ident1, "
"                         only.pos = FALSE, "
"                         min.pct = 0.25, "
                         logfc.threshold = 0.25)
"  marker.list[[i]] <- rownames(markers[1:20,])  "
}

"new_cluster <- c(""Smooth Muscle Cells"",""F3+ Fibroblast"",""Undetermined"","
"                 ""ROBO2+ Fibroblast"",""ROBO2+ Fibroblast"",""Undetermined"","
"                 ""Smooth Muscle Cells"",""Smooth Muscle Cells"",""Smooth Muscle Cells"",""Mesothelial Cells"")"

names(new_cluster) <- levels(fibroblast)
"fibroblast <- RenameIdents(fibroblast,new_cluster)"
"DimPlot(fibroblast, "
"        repel = TRUE,"
"        reduction = ""tsne"", "
"        label=TRUE, "
        pt.size = .1) + NoLegend()

"features = c(""F3"",""ROBO2"")"
length(features)
"VlnPlot(fibroblast,features = features[1:2],cols = cors,pt.size = 0)"

"saveRDS(fibroblast, file=""fb_seurat.rds"")"

####3. bar plot of F3+ fibroblasts and ROBO2+ fibroblasts in normal and IPF####
library(Seurat)
library(CellChat)
library(reshape2)

fibroblast <- readRDS("/fb_seurat.rds")
fibroblast$celltype <- fibroblast@active.ident

"count_df = table(fibroblast@active.ident,fibroblast$tissue_type)"
"count_df <- melt(count_df,value.name = ""Count"",varnames = c(""Cell Type"",""sample""))"
count_df$sample <- as.character(count_df$sample)
"count_df <- count_df[order(count_df$sample,count_df$`Cell Type`),]"
count_df <- count_df %>% filter(`Cell Type`=='F3+ Fibroblast' | `Cell Type`=='ROBO2+ Fibroblast')
percent_df = count_df %>% group_by(`sample`) %>% mutate(Percentage=Count/sum(Count))
"ggplot(percent_df, aes(fill=`Cell Type`, y=Percentage, x=sample)) + "
"  geom_bar(position=""fill"", stat=""identity"")"


####4. monocle and GO/KEGG####
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(CellChat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
"rm(count_df,marker.list,markers,percent_df)"
setwd("F:/IPF1/Single Cell")
fibroblast <- readRDS("./fb_seurat.rds")
f1 <- fibroblast@meta.data
fibroblast$celltype <- fibroblast@active.ident
"seuratObj <- subset(fibroblast, subset = celltype %in% c(""F3+ Fibroblast"",""ROBO2+ Fibroblast""))"
"data <- GetAssayData(seuratObj, assay = 'RNA', slot = 'counts')"
cell_metadata <- seuratObj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
"cds <- new_cell_data_set(data,"
"                         cell_metadata = cell_metadata,"
                         gene_metadata = gene_annotation)
"cds <- preprocess_cds(cds, num_dim = 50)"

#umap dimension reduction
"cds <- reduce_dimension(cds, preprocess_method = ""PCA"")"

## import umap embeddings from seurat
cds.embed <- cds@int_colData$reducedDims$UMAP
"int.embed <- Embeddings(seuratObj, reduction = ""umap"")"
"int.embed <- int.embed[rownames(cds.embed),]"
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
"cds <- learn_graph(cds,use_partition = FALSE)"
cds <- order_cells(cds)
"p2 <- plot_cells(cds,"
"                 color_cells_by = ""pseudotime"","
"                 label_cell_groups=FALSE,"
"                 label_leaves=FALSE,"
"                 label_branch_points=FALSE,"
                 graph_label_size=1.5)
"p1 <- plot_cells(cds, reduction_method=""UMAP"", color_cells_by=""celltype"") + ggtitle('cds.umap')"
"Track_genes <- graph_test(cds, neighbor_graph=""principal_graph"", cores=8, alternative = ""two.sided"")"
"write.table(Track_genes,""./monocle and cellchat/time.txt"",quote = FALSE,row.names=FALSE)"
"write.csv(Track_genes,""./monocle and cellchat/time.csv"")"

"#Track_genes <- graph_test(cds, neighbor_graph=""principal_graph"", cores=8, alternative = ""two.sided"")"
"Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%"
  pull(gene_short_name) %>% as.character()
"p3 = plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by=""celltype"", "
"                              min_expr=0.5, ncol = 2)"

# gene expression trend
"p4 = plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by=""pseudotime"", "
"                              min_expr=0.5, ncol = 2)"

p3+p4
"p5 = plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,"
"                label_cell_groups=FALSE,  label_leaves=FALSE)"

# violin plot
"p6 = plot_genes_violin(cds[Track_genes_sig], group_cells_by=""celltype"", ncol=2) +"
"  theme(axis.text.x=element_text(angle=45, hjust=1))"
p5+p6
#de analysis
"fb.de.markers <- FindMarkers(fibroblast, ident.1 = ""F3+ Fibroblast"",ident.2 = ""ROBO2+ Fibroblast"")"

mcaf.marker <- fb.de.markers %>% filter(p_val_adj<0.01 & avg_log2FC>1) %>% rownames()
immune.marker <- fb.de.markers %>% filter(p_val_adj<0.01 & avg_log2FC< -1) %>% rownames()

"write.csv(mcaf.marker,""./F3.csv"")"
"write.csv(immune.marker,""./ROBO2.csv"")"


"mcaf.marker <- mapIds(org.Hs.eg.db, mcaf.marker, 'ENTREZID', 'SYMBOL')"
"immune.marker <- mapIds(org.Hs.eg.db, immune.marker, 'ENTREZID', 'SYMBOL')"

"mcaf.go <- enrichGO(gene = mcaf.marker,"
"                    OrgDb=org.Hs.eg.db,"
"                    pAdjustMethod =""BH"", minGSSize=1,"
"                    pvalueCutoff=0.1, qvalueCutoff=0.05, readable=TRUE)"
"mcaf.kk<-enrichKEGG(gene=mcaf.marker, organism = 'hsa', "
"                    pvalueCutoff =0.9,use_internal_data =F)"

"immune.go <- enrichGO(gene = immune.marker,"
"                      OrgDb=org.Hs.eg.db,"
"                      pAdjustMethod =""BH"", minGSSize=1,"
"                      pvalueCutoff=0.1, qvalueCutoff=0.05, readable=TRUE)"
"immune.kk<-enrichKEGG(gene=immune.marker, organism = 'hsa', "
"                      pvalueCutoff = 0.9,use_internal_data =F)"
"dotplot(mcaf.go,title=""F3+ Fibroblasts markers Enrichment GO dot"")"
"dotplot(immune.go,title=""ROBO2+ Fibroblasts markers Enrichment GO dot"")"
"dotplot(mcaf.kk,title=""F3+ Fibroblasts markers Enrichment KEGG dot"")"
"dotplot(immune.kk,title=""ROBO2+ Fibroblasts markers Enrichment KEGG dot"")"

## find co-expression modules
"genelist <- pull(Track_genes, gene_short_name) %>% as.character()"
"gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)"


####5. analysis of bulk RNA_seq data####
data <- normalize
"Type=gsub(""(.*)\\_(.*)"", ""\\2"", colnames(data))"
####deal with the raw data####
setwd("F:/IPF1/Training Set")
library(tidyverse)
BiocManager::install('GEOquery')
library(GEOquery)

chooseBioCmirror()
"gset = getGEO('GSE70866', destdir=""."", AnnotGPL = F, getGPL = F)"

Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
class(gset)

#install.packages('ggfortify')
library(ggfortify)
gset[[1]]
gset[1]

pdata <- pData(gset[[1]])
table(pdata$title)
library(stringr)

"group_list <- ifelse(str_detect(pdata$title, ""BAL_IPF""), ""IPF"","
                     "normal")

"group_list = factor(group_list,"
"                    levels = c(""normal"",""IPF""))"

exp <- exprs(gset[[1]])
"boxplot(exp,outline=FALSE, notch=T,las=2)"
dev.off()

dim(exp)
library(limma) 
exp=normalizeBetweenArrays(exp)
"boxplot(exp,outline=FALSE, notch=T, las=2)"
range(exp)
exp <- log2(exp+1)
range(exp)
dev.off()

"#write.table(exp,file=""./exp.txt"",sep = ""\t"",row.names = F,col.names = T,quote = T)"
"write.table(pdata,file=""./cli.txt"",sep = ""\t"",row.names = F,col.names = T,quote = T)"

#index = gset[[1]]@annotation
#if(!require("hgu133plus2.db"))
#BiocManager::install("hgu133plus2.db")
#library(hgu133plus2.db)
#ls("package:hgu133plus2.db")
#ids <- toTable(hgu133plus2SYMBOL)
#head(ids)
#length(unique(ids$symbol))
#table(sort(table(ids$symbol)))

#library(tidyverse)
#exp <- as.data.frame(exp)
#exp <- exp %>% mutate(probe_id=rownames(exp))
"#exp <- exp %>% inner_join(ids,by=""probe_id"") "
"#exp <- exp[!duplicated(exp$symbol),]"
#rownames(exp) <- exp$symbol
"#exp <- exp[,-(29:30)]"
"#write.table(exp, file = ""exp.txt"",sep = ""\t"",row.names = T,col.names = NA,quote = F)"

exp <- as.data.frame(exp)

"comname <- intersect(rownames(exp),rownames(GPL))"
"exp <- exp[comname,]"
"GPL <- GPL[comname,]"
"exp1 <- cbind(GPL,exp)"
"exp1 <- exp1[!duplicated(exp1$GENE_SYMBOL),]"
rownames(exp1) <- exp1$GENE_SYMBOL
"exp1 <- exp1[,-(1:16)]"
"write.table(exp1, file = ""./geneMatrix"",sep = ""\t"",row.names = T,col.names = NA,quote = F)"


####quality control of the raw data####

"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)               
expFile="geneMatrix.txt"     
conFile="s1.txt"             
treatFile="s2.txt"          
setwd("C:\Training Set")      

"rt=read.table(expFile, header=T, sep=""\t"", check.names=F)"
rt=as.matrix(rt)
"rownames(rt)=rt[,1]"
"exp=rt[,2:ncol(rt)]"
"dimnames=list(rownames(exp),colnames(exp))"
"data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)"
data=avereps(data)

"s1=read.table(conFile, header=F, sep=""\t"", check.names=F)"
"sampleName1=as.vector(s1[,1])"
"conData=data[,sampleName1]"

"s2=read.table(treatFile, header=F, sep=""\t"", check.names=F)"
"sampleName2=as.vector(s2[,1])"
"treatData=data[,sampleName2]"

"rt=cbind(conData, treatData)"
rt <- exp1
"qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))"
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)
conNum=ncol(conData)
treatNum=ncol(treatData)
"Type=c(rep(""Control"",conNum),rep(""Treat"",treatNum))"
"outData=rbind(id=paste0(colnames(data),""_"",Type),data)"
"write.table(outData, file=""normalize.txt"", sep=""\t"", quote=F, col.names=F)"
"boxplot(normalize,outline=FALSE, notch=T, las=2)"
####got the DEGs####
exp <- normalize
library(limma)
design=model.matrix(~group_list)
"fit=lmFit(exp,design)"
fit=eBayes(fit)
"deg=topTable(fit,coef=2,number = Inf)"
"write.table(deg, file = ""deg_all.txt"",sep = ""\t"",row.names = T,col.names = NA,quote = F)"

deg$diffexpressed <- 'Not significant'
genes <- deg %>% filter(adj.P.Val<0.05&logFC< -1.5) %>% rownames()
"deg[genes,""diffexpressed""] <- 'Downregulated'"
genes <- deg %>% filter(adj.P.Val<0.05&logFC> 1.5) %>% rownames()
"deg[genes,""diffexpressed""] <- 'Upregulated'"
"ggplot(data = deg, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed)) +"
"  geom_vline(xintercept = c(-1.5, 1.5), col = ""black"", linetype = 'dashed') +"
"  geom_hline(yintercept = -log10(0.05), col = ""black"", linetype = 'dashed') +"
  geom_point(size = 2) +
"  scale_color_manual(values = c(""#00AFBB"", ""grey"", ""#bb0c00""), "
"                     labels = c(""Downregulated"", ""Not significant"", ""Upregulated"")) + "
"  coord_cartesian(ylim = c(0, 10), xlim = c(-5, 5)) + "
"  labs(color = 'Severe', "
"       x = expression(""log""[2]*""FC""), y = expression(""-log""[10]*""q-value"")) +"
"  scale_x_continuous(breaks = seq(-10, 10, 2)) +"
  ggtitle('Cancer tissue VS Normal lung tissue')
diff.gene <- deg %>% filter(abs(logFC)>1.5&adj.P.Val<0.05) %>% rownames()
"DEGs <- deg[diff.gene,]"
"write.table(DEGs, file = ""DEGs.txt"",sep = ""\t"",row.names = T,col.names = NA,quote = F)"
"write.csv(DEGs,file = ""DEGs.csv"")"
####consesus clustering####

"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")
geneFile="gene.txt"
setwd("F:/IPF1/Consesus Clustering")    
"geneRT=read.table(geneFile, header=F, sep=""\t"", check.names=F)"
"data=normalize[as.vector(geneRT[,1]),]"
"row.names(data)=gsub(""-"", ""_"", row.names(data))"

library(ConsensusClusterPlus)      
workDir="F:/IPF1/Consesus Clustering"     
setwd(workDir)     


"#data=read.table(expFile, header=T, sep=""\t"", check.names=F, row.names=1)"
data=as.matrix(data)


"group=sapply(strsplit(colnames(data),""\\_""), ""["", 2)"
"data=data[,group==""Treat""]"


maxK=9    
"results=ConsensusClusterPlus(data,"
"                             maxK=maxK,"
"                             reps=50,"
"                             pItem=0.8,"
"                             pFeature=1,"
"                             title=workDir,"
"                             clusterAlg=""km"","
"                             distance=""euclidean"","
"                             seed=123456,"
                             plot="png")


"calcICL(results, title=""consensusScore"", plot=""png"")"


clusterNum=2       
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("Cluster")
"cluster$Cluster=paste0(""C"", cluster$Cluster)"
"outTab=cbind(t(data), cluster)"
"outTab=rbind(ID=colnames(outTab), outTab)"
"write.table(outTab, file=""cluster.txt"", sep=""\t"", quote=F, col.names=F)"


####The marker genes in the 2 clusters####
"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")

library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

clusterFile="cluster.txt"      
setwd("F:/IPF1/cluster")     


"rt=read.table(clusterFile, header=T, sep=""\t"", check.names=F, row.names=1)"
"rt=rt[order(rt$Cluster),]"


"data=t(rt[,1:(ncol(rt)-1),drop=F])"
"Type=rt[,ncol(rt),drop=F]"


"bioCol=c(""#0066FF"",""#FF0000"",""#FF9900"",""#6E568C"",""#7CC767"",""#223D6C"",""#D20A13"",""#FFD121"",""#088247"",""#11AA4D"")"
ann_colors=list()
crgCluCol=bioCol[1:length(levels(factor(Type$Cluster)))]
names(crgCluCol)=levels(factor(Type$Cluster))
ann_colors[["Cluster"]]=crgCluCol


"pdf(""heatmap.pdf"", width=7, height=4.5)"
"pheatmap(data,"
"         annotation=Type,"
"         annotation_colors = ann_colors,"
"         color = colorRampPalette(c(rep(""blue"",2), ""white"", rep(""red"",2)))(100),"
"         cluster_cols =F,"
"         cluster_rows =T,"
"         scale=""row"","
"         show_colnames=F,"
"         show_rownames=T,"
"         fontsize=7,"
"         fontsize_row=7,"
         fontsize_col=7)
dev.off()


"data=melt(rt, id.vars=c(""Cluster""))"
"colnames(data)=c(""Cluster"", ""Gene"", ""Expression"")"


"p=ggboxplot(data, x=""Gene"", y=""Expression"", color = ""Cluster"","
"            xlab="""","
"            ylab=""Gene expression"","
"            legend.title=""Cluster"","
"            palette = crgCluCol,"
"            width=0.8,"
            add="point")
p=p+rotate_x_text(60)
"p1=p+stat_compare_means(aes(group=Cluster),"
"                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c(""***"", ""**"", ""*"", "" "")),"
                        label = "p.signif")


"pdf(file=""boxplot.pdf"", width=7, height=5)"
print(p1)
dev.off()

####pca####

"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


library(limma)
library(ggplot2)

clusterFile="cluster.txt"      
setwd("F:/IPF1/Consesus Clustering")     


"rt=read.table(clusterFile, header=T, sep=""\t"", check.names=F, row.names=1)"
"data=rt[,1:(ncol(rt)-1),drop=F]"
"Cluster=as.vector(rt[,ncol(rt)])"


data.pca=prcomp(data)
pcaPredict=predict(data.pca)
"PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Cluster=Cluster)"
"PCA.mean=aggregate(PCA[,1:2], list(Cluster=PCA$Cluster), mean)"

"bioCol=c(""#0066FF"",""#FF0000"",""#FF9900"",""#6E568C"",""#7CC767"",""#223D6C"",""#D20A13"",""#FFD121"",""#088247"",""#11AA4D"")"
crgCluCol=bioCol[1:length(levels(factor(Cluster)))]



"veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {"
  theta <- (0:npoints) * 2 * pi/npoints
"  Circle <- cbind(cos(theta), sin(theta))"
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$Cluster))){
"  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$Cluster==g,],"
"                                                   veganCovEllipse(cov.wt(cbind(PC1,PC2),"
"                                                                          wt=rep(1/length(PC1),length(PC1)))$cov,"
"                                                                   center=c(mean(PC1),mean(PC2))))), Cluster=g))"
}


"pdf(file=""PCA.pdf"", width=6.5, height=5)"
"ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Cluster)) +"
"  scale_colour_manual(name=""Cluster"", values =crgCluCol)+"
  theme_bw()+
"  theme(plot.margin=unit(rep(1.5,4),'lines'))+"
"  geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=Cluster), size=1, linetype=2)+"
"  annotate(""text"",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$Cluster, cex=7)+"
"  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())"
dev.off()




####WGCNA####
"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
"#BiocManager::install(c(""GO.db"", ""preprocessCore"", ""impute"",""limma""))"

"#install.packages(c(""gplots"", ""matrixStats"", ""Hmisc"", ""foreach"", ""doParallel"", ""fastcluster"", ""dynamicTreeCut"", ""survival"")) "
#install.packages("WGCNA")

library(limma)
library(gplots)
library(WGCNA)
#BiocManager::install("WGCNA")

expFile="normalize.txt"       
clusterFile="cluster.txt"   
setwd("F:/IPF1/WGCNA")    

"rt=read.table(expFile, header=T, sep=""\t"", check.names=F)"
rt=as.matrix(rt)
"rownames(rt)=rt[,1]"
"exp=rt[,2:ncol(rt)]"
"dimnames=list(rownames(exp),colnames(exp))"
"data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)"
data=avereps(data)
"selectGenes=names(tail(sort(apply(data,1,sd)), n=round(nrow(data)*0.25)))"
data=as.data.frame(data)
"data=data[selectGenes,]"
data=as.matrix(data)

"cluster=read.table(clusterFile, header=T, sep=""\t"", check.names=F, row.names=1)"
"nameC1=row.names(cluster[cluster$Cluster==""C1"",,drop=F])"
"nameC2=row.names(cluster[cluster$Cluster==""C2"",,drop=F])"
"dataC1=data[,nameC1,drop=F]"
"dataC2=data[,nameC2,drop=F]"
conCount=ncol(dataC1)
treatCount=ncol(dataC2)
"data=cbind(dataC1, dataC2)"
datExpr0=t(data)

"gsg = goodSamplesGenes(datExpr0, verbose = 3)"
if (!gsg$allOK){
"  # Optionally, print the gene and sample names that were removed:"
  if (sum(!gsg$goodGenes)>0)
"    printFlush(paste(""Removing genes:"", paste(names(datExpr0)[!gsg$goodGenes], collapse = "", "")))"
  if (sum(!gsg$goodSamples)>0)
"    printFlush(paste(""Removing samples:"", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = "", "")))"
  # Remove the offending genes and samples from the data:
"  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]"
}


"sampleTree = hclust(dist(datExpr0), method = ""average"")"
"pdf(file = ""cluster01.sample_cluster.pdf"", width = 12, height = 9)"
par(cex = 0.6)
"par(mar = c(0,4,2,0))"
"plot(sampleTree, main = ""Sample clustering to detect outliers"", sub="""", xlab="""", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)"

"abline(h = 20000, col = ""red"")"
dev.off()


"clust=cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)"
table(clust)
keepSamples=(clust==1)
"datExpr0=datExpr0[keepSamples, ]"



"traitData=data.frame(C1=c(rep(1,conCount),rep(0,treatCount)),"
"                     C2=c(rep(0,conCount),rep(1,treatCount)))"
row.names(traitData)=colnames(data)
fpkmSamples=rownames(datExpr0)
traitSamples=rownames(traitData)
"sameSample=intersect(fpkmSamples,traitSamples)"
"datExpr0=datExpr0[sameSample,]"
"datTraits=traitData[sameSample,]"


"sampleTree2 = hclust(dist(datExpr0), method = ""average"")"
"traitColors = numbers2colors(datTraits, signed = FALSE)"
"pdf(file=""cluster02.sample_heatmap.pdf"",width=12,height=12)"
"plotDendroAndColors(sampleTree2, traitColors,"
"                    groupLabels = names(datTraits),"
                    main = "Sample dendrogram and trait heatmap")
dev.off()

###power
enableWGCNAThreads()   
powers = c(1:20)     
"sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)"
"pdf(file=""cluster03.scale_independence.pdf"",width=9,height=5)"
"par(mfrow = c(1,2))"
cex1 = 0.9

"plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],"
"     xlab=""Soft Threshold (power)"",ylab=""Scale Free Topology Model Fit,signed R^2"",type=""n"","
     main = paste("Scale independence"));
"text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],"
"     labels=powers,cex=cex1,col=""red"");"
"abline(h=0.90,col=""red"") "

"plot(sft$fitIndices[,1], sft$fitIndices[,5],"
"     xlab=""Soft Threshold (power)"",ylab=""Mean Connectivity"", type=""n"","
     main = paste("Mean connectivity"))
"text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col=""red"")"
dev.off()


sft 
softPower =sft$powerEstimate     
"adjacency = adjacency(datExpr0, power = softPower)"
softPower



TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


"geneTree = hclust(as.dist(dissTOM), method = ""average"");"
"pdf(file=""cluster04.gene_clustering.pdf"",width=12,height=9)"
"plot(geneTree, xlab="""", sub="""", main = ""Gene clustering on TOM-based dissimilarity"","
"     labels = FALSE, hang = 0.04)"
dev.off()


minModuleSize = 100      
"dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,"
"                            deepSplit = 2, pamRespectsDendro = FALSE,"
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
"pdf(file=""cluster05.Dynamic_Tree.pdf"",width=8,height=6)"
"plotDendroAndColors(geneTree, dynamicColors, ""Dynamic Tree Cut"","
"                    dendroLabels = FALSE, hang = 0.03,"
"                    addGuide = TRUE, guideHang = 0.05,"
                    main = "Gene dendrogram and module colors")
dev.off()



"MEList = moduleEigengenes(datExpr0, colors = dynamicColors)"
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
"METree = hclust(as.dist(MEDiss), method = ""average"")"
"pdf(file=""cluster06.Clustering_module.pdf"",width=7,height=6)"
"plot(METree, main = ""Clustering of module eigengenes"","
"     xlab = """", sub = """")"
dev.off()

moduleColors=dynamicColors
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
"select = sample(nGenes, size=1000)"
"selectTOM = dissTOM[select, select];"
"selectTree = hclust(as.dist(selectTOM), method=""average"")"
selectColors = moduleColors[select]
"#sizeGrWindow(9,9)"
plotDiss=selectTOM^softPower
diag(plotDiss)=NA
"myheatcol = colorpanel(250, ""red"", ""orange"", ""lemonchiffon"")    "
"pdf(file=""cluster07.TOMplot.pdf"", width=7, height=7)"
"TOMplot(plotDiss, selectTree, selectColors, main = ""Network heatmap plot, selected genes"", col=myheatcol)"
dev.off()



"moduleTraitCor = cor(MEs, datTraits, use = ""p"")"
"moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)"
"pdf(file=""cluster08.Module_trait.pdf"", width=6.5, height=5.5)"
"textMatrix = paste(signif(moduleTraitCor, 2), ""\n("","
"                   signif(moduleTraitPvalue, 1), "")"", sep = """")"
dim(textMatrix) = dim(moduleTraitCor)
"par(mar = c(3.5, 8, 3, 3))"
"labeledHeatmap(Matrix = moduleTraitCor,"
"               xLabels = names(datTraits),"
"               yLabels = names(MEs),"
"               ySymbols = names(MEs),"
"               colorLabels = FALSE,"
"               colors = blueWhiteRed(50),"
"               textMatrix = textMatrix,"
"               setStdMargins = FALSE,"
"               cex.text = 0.7,"
"               zlim = c(-1,1),"
               main = paste("Module-trait relationships"))
dev.off()



"modNames = substring(names(MEs), 3)"
"geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = ""p""))"
"MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))"
"names(geneModuleMembership) = paste(""MM"", modNames, sep="""")"
"names(MMPvalue) = paste(""p.MM"", modNames, sep="""")"
traitNames=names(datTraits)
"geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = ""p""))"
"GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))"
"names(geneTraitSignificance) = paste(""GS."", traitNames, sep="""")"
"names(GSPvalue) = paste(""p.GS."", traitNames, sep="""")"

trait="C2"
"traitColumn=match(trait,traitNames)  "
for (module in modNames){
"  column = match(module, modNames)"
  moduleGenes = moduleColors==module
"  if (nrow(geneModuleMembership[moduleGenes,]) > 1){"
"    outPdf=paste(""cluster09."", trait, ""_"", module,"".pdf"",sep="""")"
"    pdf(file=outPdf,width=7,height=7)"
"    par(mfrow = c(1,1))"
"    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),"
"                       abs(geneTraitSignificance[moduleGenes, traitColumn]),"
"                       xlab = paste(""Module Membership in"", module, ""module""),"
"                       ylab = paste(""Gene significance for "",trait),"
"                       main = paste(""Module membership vs. gene significance\n""),"
"                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)"
"    abline(v=0.8,h=0.5,col=""red"")"
    dev.off()
  }
}



probes = colnames(datExpr0)
"geneInfo0 = data.frame(probes= probes,"
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
"  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],"
"                         GSPvalue[, Tra])"
"  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],"
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
"  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],"
"                         MMPvalue[, mod])"
"  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],"
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
"geneInfo = geneInfo0[geneOrder, ]"
"write.table(geneInfo, file = ""cluster.GS_MM.xls"",sep=""\t"",row.names=F)"



for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
"  write.table(modGenes, file =paste0(""cluster.module_"",modules,"".txt""),sep=""\t"",row.names=F,col.names=F,quote=F)"
}


geneSigFilter=0.5         
moduleSigFilter=0.8       
"datMM=cbind(geneModuleMembership, geneTraitSignificance)"
"datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]"
for(mmi in colnames(datMM)[1:(ncol(datMM)-2)]){
"  dataMM2=datMM[abs(datMM[,mmi])>moduleSigFilter,]"
"  write.table(row.names(dataMM2), file =paste0(""cluster.hubGenes_"",mmi,"".txt""),sep=""\t"",row.names=F,col.names=F,quote=F)"
}


####lasso and training set####
#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")
"#rt=read.table(""./3/xiaochen/23/lasso/data.exp.txt"",header=T,sep=""\t"",row.names=1)          "

"group=sapply(strsplit(colnames(normalize),""\\_""), ""["", 2)"
"normalize=normalize[,group==""Treat""]"

geneFile="interGenes.txt"
setwd("F:/IPF1/Training Set")     
"geneRT=read.table(geneFile, header=F, sep=""\t"", check.names=F)"
"data=normalize[as.vector(geneRT[,1]),]"
"row.names(data)=gsub(""-"", ""_"", row.names(data))"
data <- t(data)
"write.csv(data,file=""./data2.csv"")"

rt$futime=rt$futime/365 

set.seed(3) 
"x=as.matrix(rt[,c(3:ncol(rt))])"
"y=data.matrix(Surv(rt$futime,rt$fustat))"
"fit=glmnet(x, y, family = ""cox"", maxit = 1000)"
"plot(fit, xvar = ""lambda"", label = TRUE)"

"cvfit = cv.glmnet(x, y, family=""cox"", maxit = 1000)"
plot(cvfit)

dev.off()


"coef=coef(fit, s = cvfit$lambda.min)"
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
"geneCoef=cbind(Gene=lassoGene,Coef=actCoef)"
geneCoef  

"FinalGeneExp = rt[,lassoGene]"
"myFun = function(x){crossprod(as.numeric(x),actCoef)}"
"riskScore = apply(FinalGeneExp,1,myFun)"
"outCol = c(""futime"", ""fustat"", lassoGene)"
"risk = as.vector(ifelse(riskScore > median(riskScore), ""high"", ""low""))"
"dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)"


#install.packages("ggpubr")
library(ggpubr)  
"p <- ggboxplot(dat, x = ""fustat"", y = ""riskScore"","
"               color = ""fustat"", palette = ""jco"","
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p  


#install.packages("ROCR")
library(ROCR)   
library(glmnet)
library(caret)
"pred <- prediction(dat$riskScore, dat$fustat)"
"perf <- performance(pred,""tpr"",""fpr"")"
"AUC <- performance(pred,""auc"")   "
"plot(perf,colorize=FALSE, col=""red"", print.auc =TRUE) "
"lines(c(0,1),c(0,1),col = ""gray"", lty = 4 )"
dev.off()


rt <- dat
color=as.vector(rt$fustat)
color[color==1]="indianred1"
color[color==0]="lightseagreen"
"plot(rt$futime, pch=19,"
"     xlab=""Patients (increasing risk socre)"", ylab=""Survival time (years)"","
     col=color)
"legend(""topleft"", c(""Dead"", ""Alive""),bty=""n"",pch=19,col=c(""indianred1"",""lightseagreen""),cex=1.2)"
"riskClass=rt[,""risk""]"
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
"abline(v=lowLength,lty=2)"
dev.off()
#y=riskscore
"rt <- rt[order(rt[,8]),]"
"riskClass=rt[,""risk""]"
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
"line=rt[,""riskScore""]"
line[line>10]=10
"plot(line, type=""p"", pch=20,"
"     xlab=""Patients (increasing risk socre)"", ylab=""Risk score"","
"     col=c(rep(""lightseagreen"",lowLength),rep(""indianred1"",highLength)) )"
"abline(h=median(rt$riskScore),v=lowLength,lty=2)"
"legend(""topleft"", c(""High risk"", ""low Risk""),bty=""n"",pch=19,col=c(""indianred1"",""lightseagreen""),cex=1.2)"
dev.off()
####timeROC####
library(ROCR)   
library(glmnet)
library(caret)
library(ggfortify)
library(DESeq2)
library(glmnet)
library(timeROC)
library(limma)
library(GEOquery)
library(survminer)
"fit <- survfit(Surv(futime,as.numeric(fustat))~risk, data=rt)"
"ggsurvplot(fit,data=rt,conf.int=T,pval=T,"
"           pval.size=5,legend.title=""Risk"","
"           legend.labs=c(""high risk"", ""low risk""),"
"           xlab=""Time(years)"",break.time.by=1,"
"           palette=c(""red"",""blue""),"
"           risk.table.title="""",risk.table.height=.25,"
           risk.table = TRUE)
"roc <- timeROC(T=rt$futime,delta=rt$fustat,marker=rt$riskScore,"
"               cause=1,weighting='aalen',times=c(1,2,3),ROC=T)"
"plot(roc,time=1,col='green',title=F,lwd=2)"
"plot(roc,time=2,col='blue',add=T,title=F,lwd=2)"
"plot(roc,time=3,col='red',add=T,title=F,lwd=2)"
"legend('bottomright',"
"       c(paste0('AUC at 1 years: ',sprintf(""%.03f"",roc$AUC[1])),"
"         paste0('AUC at 2 years: ',sprintf(""%.03f"",roc$AUC[2])),"
"         paste0('AUC at 3 years: ',sprintf(""%.03f"",roc$AUC[3]))),"
"       col=c(""green"",""blue"",""red""),lwd=2,bty='n')"
"write.table(rt, file = ""train_os_riskscore.txt"",sep=""\t"",row.names=T)"


####cox####

Coxoutput <- NULL
#surv.expr <- rt
#surv.expr$OS.time <- surv.expr$futime
#surv.expr$OS <- surv.expr$fustat
#surv.expr$OS.time=surv.expr$OS.time/365   
"#write.csv(surv.expr,file=""./data.csv"")"

for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
"  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) "
  coxSummary = summary(cox)
  
"  Coxoutput <- rbind.data.frame(Coxoutput,"
"                                data.frame(gene = g,"
"                                           HR = as.numeric(coxSummary$coefficients[,""exp(coef)""])[1],"
"                                           z = as.numeric(coxSummary$coefficients[,""z""])[1],"
"                                           pvalue = as.numeric(coxSummary$coefficients[,""Pr(>|z|)""])[1],"
"                                           lower = as.numeric(coxSummary$conf.int[,3][1]),"
"                                           upper = as.numeric(coxSummary$conf.int[,4][1]),"
"                                           stringsAsFactors = F),"
                                stringsAsFactors = F)
}


"write.table(Coxoutput, file = ""./cox results.txt"",sep = ""\t"",row.names = F,col.names = T,quote = F)"

pcutoff <- 0.05
"topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] "
"topgene <- topgene[1:5,]"
library(forestplot)


"tabletext <- cbind(c(""Gene"",topgene$gene),"
"                   c(""HR"",format(round(as.numeric(topgene$HR),3),nsmall = 3)),"
"                   c(""lower 95%CI"",format(round(as.numeric(topgene$lower),3),nsmall = 3)),"
"                   c(""upper 95%CI"",format(round(as.numeric(topgene$upper),3),nsmall = 3)),"
"                   c(""pvalue"",format(round(as.numeric(topgene$p),3),nsmall = 3)))"

"forestplot(labeltext=tabletext,"
"           mean=c(NA,as.numeric(topgene$HR)),"
"           lower=c(NA,as.numeric(topgene$lower)), "
"           upper=c(NA,as.numeric(topgene$upper)),"
"           graph.pos=5,"
"           graphwidth = unit(.25,""npc""),"
"           fn.ci_norm=""fpDrawDiamondCI"","
"           col=fpColors(box=""#00A896"", lines=""#02C39A"", zero = ""black""),"
           
"           boxsize=0.4,"
"           lwd.ci=1,"
"           ci.vertices.height = 0.1,ci.vertices=T,"
"           zero=1,"
"           lwd.zero=1.5,"
"           xticks = c(0.5,1,1.5),"
"           lwd.xaxis=2,"
"           xlab=""Hazard ratios"","
"           txt_gp=fpTxtGp(label=gpar(cex=1.2),"
"                          ticks=gpar(cex=0.85),"
"                          xlab=gpar(cex=1),"
"                          title=gpar(cex=1.5)),"
"           hrzl_lines=list(""1"" = gpar(lwd=2, col=""black""), "
"                           ""2"" = gpar(lwd=1.5, col=""black""), "
"                           ""7"" = gpar(lwd=2, col=""black"")), "
"           lineheight = unit(.75,""cm""),"
"           colgap = unit(0.3,""cm""),"
"           mar=unit(rep(1.5, times = 4), ""cm""),"
           new_page = F
)
dev.off()

library(dplyr)
library(survival)
library(survminer)
plot_data$OS=plot_data$OS/365

"fit.cox <- coxph(Surv(OS,OS.status) ~ Age + Gender + riskScore,data = plot_data)"
ggforest(fit.cox)
dev.off()

####nomogram####
library(rms)
clin_sig$OS.time <-clin_sig$OS.time*365
dd <- datadist(clin_sig)
options(datadist="dd")
setwd("F:/IPF1/lasso")
"f2 <- cph(Surv(OS.time,OS)~CXCL14+TM4SF1+CYTL1+SOD3+MMP10,"
"          data=clin_sig, x=T, y=T, surv = T, time.inc = 365)"
surv <- Survival(f2)
"nom <- nomogram(f2, fun=list(function(x) surv(365, x),"
"                             function(x) surv(730, x),"
"                             function(x) surv(1095, x)),"
"                funlabel=c(""OS > 1Y Probability"","
"                           ""OS > 2Y Probability"","
                           "OS > 3Y Probability"))
"plot(nom, xfrac=.5,cex.axis = .5,cex.var = .5)"
dev.off()

"f2 <- psm(Surv(OS.time,OS) ~ CXCL14+TM4SF1+CYTL1+SOD3+MMP10, data = clin_sig, x=T, y=T, dist='lognormal') "

"cal1 <- calibrate(f2, "
"                  cmethod='KM', "
"                  method=""boot"", "
"                  u=365, "
"                  m=37,"
                  B=1000) 
"plot(cal1,lwd=2,lty=1,"
"     conf.int=T,"
"     errbar.col=""blue"","
"     col=""red"", "
"     xlim=c(0,1),ylim=c(0,1),"
"     xlab=""Nomogram-Predicted Probability of 1-Year DFS"","
"     ylab=""Actual 1-Year DFS (proportion)"","
     subtitles = F)

"cal1 <- calibrate(f2, "
"                  cmethod='KM', "
"                  method=""boot"", "
"                  u=730, "
"                  m=37,"
                  B=1000) 
"plot(cal1,lwd=2,lty=1,"
"     conf.int=T,"
"     errbar.col=""blue"","
"     col=""red"", "
"     xlim=c(0,1),ylim=c(0,1),"
"     xlab=""Nomogram-Predicted Probability of 2-Year DFS"","
"     ylab=""Actual 2-Year DFS (proportion)"","
     subtitles = F)

"cal1 <- calibrate(f2, "
"                  cmethod='KM', "
"                  method=""boot"", "
"                  u=1095, "
"                  m=37,"
                  B=1000) 
"plot(cal1,lwd=2,lty=1,"
"     conf.int=T,"
"     errbar.col=""blue"","
"     col=""red"", "
"     xlim=c(0,1),ylim=c(0,1),"
"     xlab=""Nomogram-Predicted Probability of 3-Year DFS"","
"     ylab=""Actual 3-Year DFS (proportion)"","
     subtitles = F)

####test set####
#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")
"#rt=read.table(""./3/xiaochen/23/lasso/data.exp.txt"",header=T,sep=""\t"",row.names=1)          "

geneFile="gene.txt"
setwd("F:/IPF1/Validation Set/")    
"geneRT=read.table(geneFile, header=F, sep=""\t"", check.names=F)"
"data=normalize[as.vector(geneRT[,1]),]"
"row.names(data)=gsub(""-"", ""_"", row.names(data))"
data <- t(data)
"write.csv(data,file=""./data2.csv"")"

rt$futime=rt$futime/365 

"FinalGeneExp = rt[,lassoGene]"
"myFun = function(x){crossprod(as.numeric(x),actCoef)}"
"riskScore = apply(FinalGeneExp,1,myFun)"
"outCol = c(""futime"", ""fustat"", lassoGene)"
"risk = as.vector(ifelse(riskScore > median(riskScore), ""high"", ""low""))"
"dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)"


#install.packages("ggpubr")
library(ggpubr)  
"p <- ggboxplot(dat, x = ""fustat"", y = ""riskScore"","
"               color = ""fustat"", palette = ""jco"","
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   


library(ROCR)   
library(glmnet)
library(caret)
"pred <- prediction(dat$riskScore, dat$fustat)"
"perf <- performance(pred,""tpr"",""fpr"")"
"AUC <- performance(pred,""auc"")   "
"plot(perf,colorize=FALSE, col=""red"", print.auc =TRUE) "
"lines(c(0,1),c(0,1),col = ""gray"", lty = 4 )"
dev.off()


rt <- dat
color=as.vector(rt$fustat)
color[color==1]="indianred1"
color[color==0]="lightseagreen"
"plot(rt$futime, pch=19,"
"     xlab=""Patients (increasing risk socre)"", ylab=""Survival time (years)"","
     col=color)
"legend(""topleft"", c(""Dead"", ""Alive""),bty=""n"",pch=19,col=c(""indianred1"",""lightseagreen""),cex=1.2)"
"riskClass=rt[,""risk""]"
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
"abline(v=lowLength,lty=2)"
dev.off()
#y=riskscore
"rt <- rt[order(rt[,8]),]"
"riskClass=rt[,""risk""]"
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
"line=rt[,""riskScore""]"
line[line>10]=10
"plot(line, type=""p"", pch=20,"
"     xlab=""Patients (increasing risk socre)"", ylab=""Risk score"","
"     col=c(rep(""lightseagreen"",lowLength),rep(""indianred1"",highLength)) )"
"abline(h=median(rt$riskScore),v=lowLength,lty=2)"
"legend(""topleft"", c(""High risk"", ""low Risk""),bty=""n"",pch=19,col=c(""indianred1"",""lightseagreen""),cex=1.2)"
dev.off()
####timeROC####
library(ROCR)  
library(glmnet)
library(caret)
library(ggfortify)
library(DESeq2)
library(glmnet)
library(timeROC)
library(limma)
library(GEOquery)
library(survminer)
"fit <- survfit(Surv(futime,as.numeric(fustat))~risk, data=rt)"
"ggsurvplot(fit,data=rt,conf.int=T,pval=T,"
"           pval.size=5,legend.title=""Risk"","
"           legend.labs=c(""high risk"", ""low risk""),"
"           xlab=""Time(years)"",break.time.by=1,"
"           palette=c(""red"",""blue""),"
"           risk.table.title="""",risk.table.height=.25,"
           risk.table = TRUE)
"roc <- timeROC(T=rt$futime,delta=rt$fustat,marker=rt$riskScore,"
"               cause=1,weighting='aalen',times=c(1,2,3),ROC=T)"
"plot(roc,time=1,col='green',title=F,lwd=2)"
"plot(roc,time=2,col='blue',add=T,title=F,lwd=2)"
"plot(roc,time=3,col='red',add=T,title=F,lwd=2)"
"legend('bottomright',"
"       c(paste0('AUC at 1 years: ',sprintf(""%.03f"",roc$AUC[1])),"
"         paste0('AUC at 2 years: ',sprintf(""%.03f"",roc$AUC[2])),"
"         paste0('AUC at 3 years: ',sprintf(""%.03f"",roc$AUC[3]))),"
"       col=c(""green"",""blue"",""red""),lwd=2,bty='n')"
"write.table(rt, file = ""test_os_riskscore.txt"",sep=""\t"",row.names=T)"


####cox####

Coxoutput <- NULL
#surv.expr <- rt
#surv.expr$OS.time <- surv.expr$futime
#surv.expr$OS <- surv.expr$fustat
surv.expr$OS.time=surv.expr$OS.time/365   
"#write.csv(surv.expr,file=""./data.csv"")"

for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
"  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) "
  coxSummary = summary(cox)
  
"  Coxoutput <- rbind.data.frame(Coxoutput,"
"                                data.frame(gene = g,"
"                                           HR = as.numeric(coxSummary$coefficients[,""exp(coef)""])[1],"
"                                           z = as.numeric(coxSummary$coefficients[,""z""])[1],"
"                                           pvalue = as.numeric(coxSummary$coefficients[,""Pr(>|z|)""])[1],"
"                                           lower = as.numeric(coxSummary$conf.int[,3][1]),"
"                                           upper = as.numeric(coxSummary$conf.int[,4][1]),"
"                                           stringsAsFactors = F),"
                                stringsAsFactors = F)
}


"write.table(Coxoutput, file = ""./cox results.txt"",sep = ""\t"",row.names = F,col.names = T,quote = F)"

pcutoff <- 0.05
"topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] "
"topgene <- topgene[1:3,]"
library(forestplot)

"tabletext <- cbind(c(""Gene"",topgene$gene),"
"                   c(""HR"",format(round(as.numeric(topgene$HR),3),nsmall = 3)),"
"                   c(""lower 95%CI"",format(round(as.numeric(topgene$lower),3),nsmall = 3)),"
"                   c(""upper 95%CI"",format(round(as.numeric(topgene$upper),3),nsmall = 3)),"
"                   c(""pvalue"",format(round(as.numeric(topgene$p),3),nsmall = 3)))"

"forestplot(labeltext=tabletext,"
"           mean=c(NA,as.numeric(topgene$HR)),"
"           lower=c(NA,as.numeric(topgene$lower)), "
"           upper=c(NA,as.numeric(topgene$upper)),"
"           graph.pos=5,"
"           graphwidth = unit(.25,""npc""),"
"           fn.ci_norm=""fpDrawDiamondCI"","
"           col=fpColors(box=""#00A896"", lines=""#02C39A"", zero = ""black""),# box颜色"
           
"           boxsize=0.4,"
"           lwd.ci=1,"
"           ci.vertices.height = 0.1,ci.vertices=T,"
"           zero=1,"
"           lwd.zero=1.5,"
"           xticks = c(0.5,1,1.5),"
"           lwd.xaxis=2,"
"           xlab=""Hazard ratios"","
"           txt_gp=fpTxtGp(label=gpar(cex=1.2),"
"                          ticks=gpar(cex=0.85),"
"                          xlab=gpar(cex=1),"
"                          title=gpar(cex=1.5)),"
"           hrzl_lines=list(""1"" = gpar(lwd=2, col=""black""),"
"                           ""2"" = gpar(lwd=1.5, col=""black""), "
"                           ""5"" = gpar(lwd=2, col=""black"")), "
"           lineheight = unit(.75,""cm""),"
"           colgap = unit(0.3,""cm""),"
"           mar=unit(rep(1.5, times = 4), ""cm""),"
           new_page = F
)
dev.off()


library(dplyr)
library(survival)
library(survminer)
plot_data$OS=plot_data$OS/365
"fit.cox <- coxph(Surv(OS,OS.status) ~ Age + Gender + riskScore,data = plot_data)"
ggforest(fit.cox)
dev.off()

####model gene####
"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)             
expFile="normalize.txt"     
geneFile="gene.txt"         
setwd("F:/IPF1/Training Set")    

"rt=read.table(expFile, header=T, sep=""\t"", check.names=F)"
rt=as.matrix(rt)
"rownames(rt)=rt[,1]"
"exp=rt[,2:ncol(rt)]"
"dimnames=list(rownames(exp),colnames(exp))"
"data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)"
data=avereps(data)


"gene=read.table(geneFile, header=F, sep=""\t"", check.names=F)"
"sameGene=intersect(as.vector(gene[,1]), rownames(data))"
"geneExp=data[sameGene,]"


"out=rbind(ID=colnames(geneExp), geneExp)"
"write.table(out, file=""crgGeneExp.txt"", sep=""\t"", quote=F, col.names=F)"

"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("reshape2")
#install.packages("ggpubr")



library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

expFile="crgGeneExp.txt"     
setwd("F:/IPF1/Training Set")     


"rt=read.table(expFile, header=T, sep=""\t"", check.names=F)"
rt=as.matrix(rt)
"rownames(rt)=rt[,1]"
"exp=rt[,2:ncol(rt)]"
"dimnames=list(rownames(exp),colnames(exp))"
"data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)"
data=avereps(data)
exp=data


"Type=gsub(""(.*)\\_(.*)"", ""\\2"", colnames(data))"


sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
"  test=wilcox.test(data[i,] ~ Type)"
  pvalue=test$p.value
"  Sig=ifelse(pvalue<0.001,""***"",ifelse(pvalue<0.01,""**"",ifelse(pvalue<0.05,""*"","""")))"
  if(pvalue<0.05){
"    sigVec=c(sigVec, paste0(i, Sig))"
"    sigGeneVec=c(sigGeneVec, i)}"
}

"data=data[sigGeneVec,]"
"outTab=rbind(ID=colnames(data), data)"
"write.table(outTab, file=""diffGeneExp.txt"", sep=""\t"", quote=F, col.names=F)"
row.names(data)=sigVec


names(Type)=colnames(data)
Type=as.data.frame(Type)
"#pdf(file=""heatmap.pdf"", width=7, height=4.5)"
"pheatmap(data,"
"         annotation=Type,"
"         color = colorRampPalette(c(rep(""blue"",2), ""white"", rep(""red"",2)))(100),"
"         cluster_cols =F,"
"         cluster_rows =T,"
"         scale=""row"","
"         show_colnames=F,"
"         show_rownames=T,"
"         fontsize=7,"
"         fontsize_row=7,"
         fontsize_col=7)
dev.off()


exp=as.data.frame(t(exp))
"exp=cbind(exp, Type=Type)"
"data=melt(exp, id.vars=c(""Type""))"
"colnames(data)=c(""Type"", ""Gene"", ""Expression"")"


"p=ggboxplot(data, x=""Gene"", y=""Expression"", color = ""Type"", "
"            xlab="""","
"            ylab=""Gene expression"","
"            legend.title=""Type"","
"            palette = c(""blue"", ""red""),"
"            add=""point"","
            width=0.8)
p=p+rotate_x_text(60)
"p1=p+stat_compare_means(aes(group=Type),"
"                        method=""wilcox.test"","
"                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c(""***"", ""**"", ""*"", "" "")),"
                        label = "p.signif")


"pdf(file=""boxplot.pdf"", width=7, height=5)"
print(p1)
dev.off()

####cibersort####
#install.packages('e1071')

"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="normalize.txt"      
setwd("F:/IPF1/Training Set")      
source("geoCRG11.CIBERSORT.R")      

"outTab=CIBERSORT(""ref.txt"", inputFile, perm=1000, QN=T)"

"outTab=outTab[outTab[,""P-value""]<0.05,]"
"outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])"
"outTab=rbind(id=colnames(outTab),outTab)"
"write.table(outTab, file=""CIBERSORT-Results.txt"", sep=""\t"", quote=F, col.names=F)"


"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


library(limma)
library(reshape2)
library(ggpubr)

clusterFile="cluster.txt"          
immFile="CIBERSORT-Results.txt"     
setwd("F:/IPF1/cluster/Training Set")      

"immune=read.table(immFile, header=T, sep=""\t"", check.names=F, row.names=1)"

"group=gsub(""(.*)\\_(.*)"", ""\\2"", row.names(immune))"
"data=immune[group==""Treat"",,drop=F]"

"Cluster=read.table(clusterFile, header=T, sep=""\t"", check.names=F, row.names=1)"
"sameSample=intersect(row.names(data), row.names(Cluster))"
"rt=cbind(data[sameSample,,drop=F], Cluster[sameSample,""Cluster"",drop=F])"
"rt=rt[order(rt$Cluster, decreasing=F),]"
"conNum=nrow(rt[rt$Cluster==""C1"",])"
"treatNum=nrow(rt[rt$Cluster==""C2"",])"

"data=t(rt[,-ncol(rt)])"
"pdf(file=""barplot.pdf"", width=14.5, height=8)"
"col=rainbow(nrow(data), s=0.7, v=0.7)"
"par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)"
"a1=barplot(data, col=col, xaxt=""n"", yaxt=""n"", ylab=""Relative Percent"", cex.lab=1.8)"
"a2=axis(2,tick=F,labels=F)"
"axis(2,a2,paste0(a2*100,""%""))"
"par(srt=0,xpd=T)"
"rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col=""green"")"
"text(a1[conNum]/2,-0.035,""C1"",cex=2)"
"rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5 , ytop = -0.06,col=""red"")"
"text((a1[length(a1)]+a1[conNum])/2,-0.035,""C2"",cex=2)"
"ytick2 = cumsum(data[,ncol(data)])"
"ytick1 = c(0,ytick2[-length(ytick2)])"
"legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty=""n"",cex=1.3)"
dev.off()



data=rt
"data=melt(data, id.vars=c(""Cluster""))"
"colnames(data)=c(""Cluster"", ""Immune"", ""Expression"")"

group=levels(factor(data$Cluster))
"bioCol=c(""#0066FF"",""#FF0000"",""#6E568C"",""#7CC767"",""#223D6C"",""#D20A13"",""#FFD121"",""#088247"",""#11AA4D"")"
bioCol=bioCol[1:length(group)]
"boxplot=ggboxplot(data, x=""Immune"", y=""Expression"", color=""Cluster"","
"                  xlab="""","
"                  ylab=""Fraction"","
"                  legend.title=""Cluster"","
"                  add=""point"","
"                  width=0.8,"
                  palette=bioCol)+
  rotate_x_text(50)+
"  stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c(""***"", ""**"", ""*"", """")), label=""p.signif"")"
"pdf(file=""immune.diff.pdf"", width=8, height=6)"
print(boxplot)
dev.off()

####GSVA####
"#if (!requireNamespace(""BiocManager"", quietly = TRUE))"
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")



library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
#BiocManager::install("GSVA")

expFile="normalize.txt"             
clusterFile="cluster.txt"            
gmtFile="c5.all.v2023.1.Hs.symbols.gmt"     
setwd("F:/IPF1/cluster/Training Set")    


"rt=read.table(expFile, header=T, sep=""\t"", check.names=F)"
rt=as.matrix(rt)
"rownames(rt)=rt[,1]"
"exp=rt[,2:ncol(rt)]"
"dimnames=list(rownames(exp),colnames(exp))"
"data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)"
data=avereps(data)


"group=gsub(""(.*)\\_(.*)"", ""\\2"", colnames(data))"
"data=data[,group==""Treat"",drop=F]"


"geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())"


"ssgseaScore=gsva(data, geneSets, method='gsva')"

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)


"cluster=read.table(clusterFile, header=T, sep=""\t"", check.names=F, row.names=1)"
"nameC1=row.names(cluster[cluster$Cluster==""C1"",,drop=F])"
"nameC2=row.names(cluster[cluster$Cluster==""C2"",,drop=F])"
"dataC1=ssgseaScore[,nameC1,drop=F]"
"dataC2=ssgseaScore[,nameC2,drop=F]"
conNum=ncol(dataC1)
treatNum=ncol(dataC2)
"data=cbind(dataC1, dataC2)"
"Type=c(rep(""C1"",conNum), rep(""C2"",treatNum))"


outTab=data.frame()
for(i in row.names(data)){
"  test=t.test(data[i,] ~ Type)"
  pvalue=test$p.value
  t=test$statistic
  if(pvalue<0.05){
"    Sig=ifelse(pvalue>0.05, ""Not"", ifelse(t>0,""Up"",""Down""))"
"    outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))"
  }
}


termNum=10    
"outTab=outTab[order(outTab$t),]"
"outTab=outTab[c(1:termNum,(nrow(outTab)-termNum):nrow(outTab)),]"
"pdf(file=""barplot_kegg.pdf"", width=9, height=6)"
outTab$t=as.numeric(outTab$t)
"outTab$Sig=factor(outTab$Sig, levels=c(""Down"", ""Up""))"
"gg1=ggbarplot(outTab, x=""Pathway"", y=""t"", fill = ""Sig"", color = ""white"","
"              palette=c(""blue3"", ""red3""), sort.val = ""asc"", sort.by.groups = T,"
"              rotate=TRUE, legend=""right"", title="""","
"              xlab=""Term"", ylab=""t value of GSVA score, C2 vs C1"",  legend.title=""Group"", x.text.angle=60)"
print(gg1)
dev.off()

#install.packages("RCircos")


library(RCircos)     
#BiocManager::install("RCircos")
setwd("F:/IPF1/model gene/Training Set")     

"cytoBandIdeogram=read.table(""refer.txt"", header=T, sep=""\t"", check.names=F)"
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
"RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)"

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.8
rcircos.params$point.size=7
RCircos.Reset.Plot.Parameters(rcircos.params)

"pdf(file=""RCircos.pdf"", width=7, height=7)"
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

"RCircos.Gene.Label.Data=read.table(""Rcircos.geneLabel.txt"", header=T, sep=""\t"", check.names=F)"
name.col <- 4
side <- "in"
track.num <- 1
"RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)"
track.num <- 2
"RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)"
dev.off()

####model gene cor####
#install.packages("corrplot")
#install.packages("circlize")


library(corrplot)
library(circlize)

inputFile="diffGeneExp.txt"    
setwd("F:/IPF1/model gene/Training Set")     

"data=read.table(inputFile, header=T, sep=""\t"", check.names=F, row.names=1)"


"group=gsub(""(.*)\\_(.*)"", ""\\2"", colnames(data))"
"data=data[,group==""Treat"",drop=F]"
rt=t(data)


cor1=cor(rt)


"col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))"
cor1[cor1==1]=0
"c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))"
"col1 = matrix(c1,nc=ncol(rt))"


"pdf(file=""circos.pdf"", width=7, height=7)"
"par(mar=c(2,2,2,4))"
"circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)"
"chordDiagram(cor1, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)"
par(xpd=T)

"colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))"
dev.off()
circos.clear()


"pdf(file=""corrplot.pdf"", width=7, height=7)"
"corrplot(cor1,"
"         method = ""pie"","
"         order = ""hclust"","
"         type = ""upper"","
"         col=colorRampPalette(c(""green"", ""white"", ""red""))(50)"
)
dev.off()
