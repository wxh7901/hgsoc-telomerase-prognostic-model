#================加载R包===================================
library(R.utils)
library(stringr)
library(Seurat)
library(ggplot2)
library(scDblFinder)
library(ggraph)
library(clustree)
library(Rcpp)
library(harmony)
library(MetBrewer)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(AUCell) 
library(monocle)
library(Biobase)
library(proxy)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(CellChat)
library(patchwork)
library(rlang)
library(SingleR)
library(celldex)
library(pheatmap)
#=========解压缩====================
setwd("/work/")
dir.create('1.data')
untar("./1.data/GSE184880_RAW.tar",exdir = "./1.data/GSE184880_RAW")
filenames <- list.files(path = "./1.data/GSE184880_RAW/")
setwd('./1.data/GSE184880_RAW/')
sapply(filenames,gunzip)
setwd('../../')
###===========构建seurat对象====================

fs=list.files('./1.data/GSE184880_RAW/','^GSM')
fs

samples=str_split(fs,'_',simplify = T)[,1]
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste0("./1.data/GSE184880_RAW/", str_split(y[1],'_',simplify = T)[,1])
  dir.create(folder,recursive = T)
  file.rename(paste0("./1.data/GSE184880_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0("./1.data/GSE184880_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("./1.data/GSE184880_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
# unlink('./1.data/GSE184880_RAW/GSM6567167/',recursive = T)
# unlink('./1.data/GSE184880_RAW/GSM6567169/',recursive = T)
# unlink('./1.data/GSE184880_RAW/GSM6567170/',recursive = T)
samples=list.files("./1.data/GSE184880_RAW/")
samples
dir <- file.path('./1.data/GSE184880_RAW',samples)
names(dir) <- samples

scRNAlist <- list()
for(i in 1:length(dir)){
  print(i)
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=3,min.feature = 250,
                                       project = names(dir)[i]) #min.feature默认为0，早期的单数据可能没有进行过滤，需设置 = 200
}
saveRDS(scRNAlist,"./1.data/GSE184880_Ovaries_scRNAlist.rds")
##=================合并seurat对象并保存====================
scRNA <- merge(scRNAlist[[1]], y=unlist(scRNAlist[2:length(scRNAlist)]))
saveRDS(scRNA,"./1.data/GSE184880_Ovaries.rds")
dim(scRNA)   #查看基因数和细胞总数
table(scRNA@meta.data$orig.ident)  #查看每个样本的细胞数

##======================= meta.data中添加分组信息===============
scRNA$class <- ifelse(scRNA$orig.ident %in% c("GSM5599220","GSM5599221","GSM5599222","GSM5599223","GSM5599224"),
                      "control","disease")
saveRDS(scRNA,"./1.data/GSE184880_Ovaries.rds")
###===========================过滤=======================================
scRNA <- readRDS("./1.data/GSE184880_Ovaries.rds")
#计算线粒体read的比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^(MT|mt)(-|\\.)")
#计算每个UMI检测到的基因数目
scRNA$log10GenesPerUMI <- log10(scRNA$nFeature_RNA) / log10(scRNA$nCount_RNA)
#计算血红细胞read的比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA))
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
HB.genes <- intersect(rownames(scRNA@assays$RNA),HB.genes)
if (length(HB.genes) > 0) {
  scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
}
##过滤前小提琴图
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI", 
                            "percent.HB"), ncol = 4,pt.size = 0.1,combine = T) & 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45,hjust = 0.5,size = 5), 
        axis.text.y = element_text(size = 5), 
        text = element_text(size = 6), 
        axis.line = element_line(size = 0.5))
ggsave("2.filter/GSE184880_vlnplot_raw.pdf",width = 17,height = 8,units = "cm")
#查看nFeature_RNA\nCount_RNA数值分布
metadata <- scRNA@meta.data
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 8000)
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 400) +
  geom_vline(xintercept = 40000)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 10) +
  geom_vline(xintercept = 100)
metadata %>% 
  ggplot(aes(color=orig.ident, x=log10GenesPerUMI, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 0.7)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.HB, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 10)
scRNA1 <- subset(scRNA, subset = percent.mt < 35 &
                   log10GenesPerUMI > 0.75 &
                   percent.HB < 10)

##过滤后小提琴图
VlnPlot(scRNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI",
                             "percent.HB"), ncol = 4,pt.size = 0.1,combine = T) & 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45,hjust = 0.5,size = 5), 
        axis.text.y = element_text(size = 5), 
        text = element_text(size = 6), 
        axis.line = element_line(size = 0.5))
ggsave("2.filter/GSE184880_vlnplot_filtered.pdf",width = 17,height = 8,units = "cm")
#查看细胞数目
ncol(scRNA1) #57891个细胞
## 去双胞(这一步每次跑的结果都不一样)

sce = as.SingleCellExperiment(scRNA1)
sce <- scDblFinder(sce, samples="ident") #BPPARAM=MulticoreParam(10)
scRNA1$dbl.class = sce$scDblFinder.class
table(scRNA1$dbl.class) 
VlnPlot(scRNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), same.y.lims = F, split.by = "dbl.class", adjust = .8, split.plot = T, ncol = 3, pt.size = 0)
sc.sgl<- subset(scRNA1, dbl.class=="singlet")
scRNA1=sc.sgl
rm(sc.sgl)
rm(sce)
#查看细胞数目
ncol(scRNA1) #53671
##过滤双胞后小提琴图
VlnPlot(scRNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI",
                             "percent.HB"), ncol = 4,pt.size = 0.1,combine = T) & 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45,hjust = 0.5,size = 5), 
        axis.text.y = element_text(size = 5), 
        text = element_text(size = 6), 
        axis.line = element_line(size = 0.5))
ggsave("2.filter/GSE184880_vlnplot_filtered_singlet.pdf",width = 17,height = 8,units = "cm")

## 保存过滤后seurat对象
saveRDS(scRNA1,"./3.RDS/GSE184880_Ovaries_filtered.rds")
###===========================seurat标准流程=======================================
#seurat标准流程

scRNA <- readRDS("./3.RDS/GSE184880_Ovaries_filtered.rds")

##counts数据归一化，保存在scRNA@assays$RNA@data
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000) 

##高变异基因选择
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)
LabelPoints(plot = VariableFeaturePlot(scRNA), 
            points = top10, repel = T,cex=2)+
  theme(legend.key.height = unit(0.25,"cm"),legend.key.width = unit(0.1,"cm"),
        legend.text = element_text(size = 6),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "vertical",
        text = element_text(size = 7),
        axis.text = element_text(size = 6),axis.line = element_line(size = 0.6),
        axis.ticks = element_line(size = 0.6)
  )
ggsave("4.cluster/VariableFeatures_dot.pdf",width = 6,height = 9,units = "cm")

##Z-score标准化,给予基因同等的权重（以0为中心，方差为1），数据保存在scRNA@assays$RNA@scale.data
scRNA <- ScaleData(scRNA, features = rownames(scRNA))

##线性降维
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),npcs = 50)

##挑选主成分
ElbowPlot(scRNA,ndims = 50)+
  theme(text = element_text(size = 7),
        axis.text = element_text(size = 6),axis.line = element_line(size = 0.6),
        axis.ticks = element_line(size = 0.6))
ggsave("4.cluster/ElbowPlot.pdf",width = 17,height = 8,units = "cm")
pca.num = 1:17

##去批次

scRNA = scRNA %>% RunHarmony("class", plot_convergence = TRUE, project.dim = F)#按照样本去批次，也可以选择根据class/group

##非线性降维
#统一流形逼近与投影Uniform Manifold Approximation and Projection
scRNA <- RunUMAP(scRNA, dims = pca.num, reduction = "harmony")
#t分布随机邻域嵌入t-Distributed Stochastic Neighbor Embedding
scRNA <- RunTSNE(scRNA, dims = pca.num,check_duplicates = F, reduction = "harmony")

##细胞聚类
scRNA <- FindNeighbors(scRNA, dims = pca.num, reduction = "harmony")#通过KNN（最近邻算法）计算每个细胞之间的相互距离，依据每个细胞最近邻之间的邻域重叠构建SNN（共享最近邻图）
scRNA@meta.data[,grep("^RNA_snn_res.",rownames(scRNA@meta.data))] <- NULL #重聚类的时候要把之前的聚类信息删除
scRNA <- FindClusters(scRNA, resolution = c(0.2,0.3,0.4,0.5,0.6,0.8,1.2)) #通过基于共享最近邻居（SNN）模块化优化的聚类算法来识别聚类
#如果报错可以尝试更新R包ggraph
#clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
ggsave("clustree.pdf",width = 12,height = 10)
scRNA <- FindClusters(scRNA, resolution = 0.8)#默认0.8
saveRDS(scRNA,"3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution0.8.rds")


##===============二维聚类图：可视化class================
reduc <- "tsne" #tsne 或umap

DimPlot(scRNA, label = T,label.size = 2, 
        cols = met.brewer("Egypt",n=length(unique(scRNA$class))),
        reduction = reduc,group.by = "class",
        repel = T) + 
  #NoLegend()+ 
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 9,hjust = 0),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_blank(),
        legend.key.height = unit(0.8,"cm"), legend.key.width = unit(0.1,"cm"),
        legend.key = element_blank(), legend.margin = unit(0.1,"inches"),
        legend.background = element_blank(), 
        legend.text = element_text(size = 9),legend.title = element_blank())
ggsave("4.cluster/class_tsne2.pdf",width = 18.75,height = 15,units = "cm")
ggsave("4.cluster/class_tsne2.png",width = 18.75,height = 15,units = "cm")


##二维聚类图：可视化seurat_clusters
DimPlot(scRNA, label = T,label.size = 3, 
        cols = met.brewer("Egypt",n=length(unique(scRNA$seurat_clusters))),
        reduction = reduc,group.by = "seurat_clusters",
        repel = T)+
  #NoLegend()+
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 9,hjust = 0),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_blank(),
        legend.key.height = unit(0.8,"cm"), legend.key.width = unit(0.1,"cm"),
        legend.key = element_blank(), legend.margin = unit(0.1,"inches"),
        legend.background = element_blank(), 
        legend.text = element_text(size = 9),legend.title = element_blank())
ggsave("4.cluster/seurat_clusters_tsne2.pdf",width = 18.75,height = 15,units = "cm")
ggsave("4.cluster/seurat_clusters_tsne2.png",width = 18.75,height = 15,units = "cm")

##===========二维聚类图：分组可视化seurat_clusters=============
DimPlot(scRNA, label = T,label.size = 2,ncol = 5, 
        cols = met.brewer("Egypt",n=length(unique(scRNA$seurat_clusters))),
        reduction = reduc,group.by = "seurat_clusters",split.by = "class",
        repel = T)+
  #NoLegend()+
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 9,hjust = 0),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_blank(),
        legend.key.height = unit(0.8,"cm"), legend.key.width = unit(0.1,"cm"),
        legend.key = element_blank(), legend.margin = unit(0.1,"inches"),
        legend.background = element_blank(), 
        legend.text = element_text(size = 9),legend.title = element_blank())
ggsave("4.cluster/seurat_clusters_split.by.class_tsne2.pdf",width = 18.75,height = 15,units = "cm")
ggsave("4.cluster/seurat_clusters_split.by.class_tsne2.png",width = 18.75,height = 15,units = "cm")


# cluster差异分析并绘制热图 --------------------------------------------------------
## cluster差异分析 
scRNA<-readRDS("3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution0.8.rds")
Idents(scRNA) <- "seurat_clusters"
sce.markers1 <- FindAllMarkers(scRNA, test.use = "wilcox", #默认检验方法, 耗时30min-60min
                               logfc.threshold = 1, #两细胞群体差异倍数1以上
                               min.pct = 0.1, #检测的基因在两细胞群体具有10%以上细胞中表达
                               only.pos = F) #是否只保留上调基因
write.csv(sce.markers1,"4.cluster/DEG_per_clusters.csv")

##绘制热图
scRNA1 <- subset(scRNA,downsample=100)
deg <- sce.markers1 %>% group_by(cluster) %>% top_n(2,wt = avg_log2FC)
DoHeatmap(scRNA1, features = deg$gene,group.by = "seurat_clusters",angle = 45,raster = F,size = 0) +
  #scale_fill_gradientn(colors = c("blue", "white", "red"))+
  theme(axis.text = element_text(size = 6,face = "plain",colour = "black"),
        legend.text = element_text(size = 6,face = "plain",colour = "black"),
        legend.title = element_text(size = 6,face = "plain",colour = "black"),
        legend.key.height = unit(0.3,"cm"),legend.key.width = unit(0.1,"cm"))
ggsave("4.cluster/cluster_DEG_heatmap.pdf",width = 12,height = 12,units = "cm")


# singleR自动注释- HumanPrimaryCellAtlasData: the Human Primary Cell A --------


###===========================singleR自动注释=======================================


# 选择参考数据集
ref <- HumanPrimaryCellAtlasData()
pred.scRNA <- SingleR(test = scRNA@assays$RNA@data, ref = ref,labels = ref$label.fine,
                      clusters = scRNA@active.ident, fine.tune = TRUE )
pred.scRNA$pruned.labels




# 绘制热图看cluster注释情况
heatmap <- plotScoreHeatmap(pred.scRNA, clusters=pred.scRNA@rownames,
                            fontsize.row = 15,show_colnames = T)
ggsave("5.singleR/singleR_heatmap.pdf",plot = heatmap,width = 60,height = 10)

# 注释结果导入
Idents(scRNA) <- "seurat_clusters"
new.cluster.ids <- c('CD4+_effector_memory','Smooth_muscle_cells:bronchial','Osteoblasts','NK cell',
                     'NK cell','Macrophage:monocyte-derived','BMP2','CD4+_effector_memory',
                     'Chondrocytes:MSC-derived','Monocyte','Epithelial_cells:bladder','Epithelial_cells:bladder',
                     'Epithelial_cells:bladder','BMP2','Endothelial_cells:lymphatic','Epithelial_cells:bladder',
                     'gamma-delta','Endothelial_cells:lymphatic','Plasma_cell','Epithelial_cells:bladder',
                     'NA','Epithelial_cells:bladder','Smooth_muscle_cells:vascular:IL-17','BMP2',
                     'Macrophage:monocyte-derived','B_cell:Memory','MSC','CD4+_effector_memory',
                     'MSC','Osteoblasts','B_cell:Plasma_cell')

names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA@meta.data$celltype <- Idents(scRNA)

reduc <- "tsne"
# 二维聚类图：可视化celltype
DimPlot(scRNA, label = T,label.size = 2,
        cols = met.brewer("Egypt",n=length(unique(scRNA$celltype))),
        reduction = reduc,group.by = "celltype", #split.by = "class",
        repel = T) +
  #NoLegend()+
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 9,hjust = 0),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(),
        legend.key.height = unit(0.8,"cm"), legend.key.width = unit(0.1,"cm"),
        legend.key = element_blank(), legend.margin = unit(0.1,"inches"),
        legend.background = element_blank(),
        legend.text = element_text(size = 9),legend.title = element_blank())
ggsave(paste("5.singleR/celltype1_",reduc,".pdf"), width = 8.5, height = 6)

saveRDS(scRNA,"3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution0.8_anno.rds")

# ##cluster注释人工注释
#
# ###===========================人工注释v1=======================================
scRNA<-readRDS("3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution0.8.rds")
#找marker基因
marker<-read.csv('5.anno/marker_gene.csv',header = T)
marker <- marker%>% distinct(marker,.keep_all = T)
rownames(marker)<-marker$marker

genes <- intersect(rownames(scRNA@assays$RNA@data),unique(marker$marker))

##可视化marker基因在不同cluster中的表达，确定cluster属于哪种细胞类型
p <- DotPlot(scRNA, features = genes, group.by = "seurat_clusters") + RotatedAxis() +
  scale_x_discrete("") + scale_y_discrete("")

data<-p$data

data$celltype<-marker[as.character(data$features.plot),'celltype']

p1 <- ggplot(data, aes(id,features.plot, fill =avg.exp.scaled, size=pct.exp))+
  geom_point(shape = 21, colour="black")+
  theme(axis.text.x = element_text(angle = 45,vjust= 1, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        strip.background.x = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background.y = element_rect(color="black", size=0., linetype="solid"),
        axis.line = element_line(color="black",size=0.5),
        strip.text.y.left = element_text(angle=0, face="bold"),
        panel.spacing.x = unit(0.2, "cm"),
        panel.spacing.y = unit(0.2, "cm")
  )+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  facet_grid(celltype~., scales = "free", space="free", switch = "y")+
  labs(x = "Celltype", y = "Marker_genes")+
  scale_y_discrete(position = "right")
ggsave("5.anno/cellanno20.pdf")
print(p1)

##人工注释结果导入
#按照cluster顺序写入相应的注释结果
Idents(scRNA) <- "seurat_clusters"
new.cluster.ids <- c('T Cell','Cancer stem cell','Mesenchymal cell','T Cell','T Cell','Macrophage',
                     'Granulosa cell','Immune cell','Mesenchymal cell','Macrophage','Cancer cell','Cancer cell',
                     'Epithelial cell','Smooth muscle cell','Endothelial cell','Cancer cell','T Cell',
                     'Endothelial cell','B cell','Cancer cell','Endothelial cell',
                     'Cancer cell','Granulosa cell','Mesenchymal cell','Macrophage','Immune cell','Granulosa cell','T Cell','Mesenchymal cell','Fibroblast','B cell')#0-5,6-11,12-16,17-20
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA@meta.data$celltype <- Idents(scRNA)

# 二维聚类图：可视化celltype
DimPlot(scRNA, label = T,label.size = 2,
        cols = met.brewer("Egypt",n=length(unique(scRNA$celltype))),
        reduction = "tsne",group.by = "celltype", #split.by = "class",
        repel = T) +
  #NoLegend()+
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 9,hjust = 0),
        axis.text = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(),
        legend.key.height = unit(0.8,"cm"), legend.key.width = unit(0.1,"cm"),
        legend.key = element_blank(), legend.margin = unit(0.1,"inches"),
        legend.background = element_blank(),
        legend.text = element_text(size = 9),legend.title = element_blank())
ggsave("5.anno/celltype1_tsne4.pdf", width = 18.75, height = 15)

saveRDS(scRNA,"3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution_0.8_anno.rds")

#### 细胞类型marker小提琴图
## 选取 marker
marker <- read.csv("5.anno/marker_gene.csv",sep = ",",header = T)

marker<- c('CD69',
           'CXCR6',
           'ITGAE',
           'DCN',
           'CHD7',
           'KIT',
           'CD74',
           'CLDN4',
           'EPCAM',
           'FABP1',
           'KRT20',
           'LCN2',
           'PECAM1',
           'CALD1',
           'LGR5',
           'PROM1',
           'CD37',
           'MS4A1'
)

## 控制画图的细胞类型顺序
scRNA$celltype <- factor(scRNA$celltype, levels =  rev(c('Mesenchymal cell', 'Epithelial cell', 'Enterocyte',
                                                         'Macrophage', 'T cell' ,'B cell', 'Endothelial cell', 'Stromal cell',
                                                         'Fibroblast' ,'Mast cell')))

## 提取marker表达数据，导入meta.data

#marker$marker=as.character(marker$marker)

#marker<-marker[marker$marker%in%rownames(scRNA@assays$RNA@counts),]#可能表达矩阵不存在其中的marker基因

vln.df = as.data.frame(scRNA[["RNA"]]@data[marker,])
vln.df$gene = rownames(vln.df)
vln.df = melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)] = c("CB","exp")

anno = scRNA@meta.data
vln.df = merge(vln.df,anno,by.x="CB", by.y="row.names")

vln.df$gene = factor(vln.df$gene,levels = marker) #为了控制画图的基因顺序

## 绘制小提琴图--以细胞类型着色
vln.df %>% ggplot(aes(celltype,exp)) + geom_violin(aes(fill=celltype), size=0.3, scale = "width")+
  facet_grid(vln.df$gene~., scales = "free_y")+
  scale_fill_manual(values = met.brewer("Juarez",n=length(unique(scRNA$celltype))))+
  #scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("") + scale_y_continuous("")+
  theme_bw()+
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.5, colour = "grey40"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA, colour = "grey40"),
        legend.position = "none")
ggsave("Annotation_CRC_celltype_vln1.pdf", width = 11, height = 9, units = "cm")
#
#
## 绘制小提琴图--以基因着色
vln.df %>% ggplot(aes(celltype,exp)) + geom_violin(aes(fill=gene), size=0.3, scale = "width")+
  facet_grid(vln.df$gene~., scales = "free_y")+
  scale_fill_manual(values = met.brewer("Juarez", n=length(marker)))+
  #scale_fill_brewer(palette = "Set3", direction = 1)+
  scale_x_discrete("") + scale_y_continuous("")+
  theme_bw()+
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.5, colour = "grey40"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA, colour = "grey40"),
        legend.position = "none")
ggsave("Annotation_CRC_celltype_vln2.pdf", width = 11, height = 9, units = "cm")

###细胞类型marker气泡图
## 读取marker基因

marker <- read.csv("5.anno/marker_gene.csv",header = T)

## 提取画图数据
scRNA@active.ident <- scRNA$celltype
genes <- intersect(rownames(scRNA@assays$RNA@data),unique(marker$marker))
p <- DotPlot(scRNA, features = genes, assay="RNA", group.by="celltype") +
  theme(axis.text.x = element_text(angle=45, vjust=1))
data <- p$data
data <- filter(data, data$pct.exp != 0)
data[!complete.cases(data),]
data <- na.omit(data)
data <- merge(data,marker,by.x = "features.plot",by.y = "marker")

# 控制画图的细胞类型和基因顺序
data$celltype<-factor(data$celltype,levels = levels(data$id))
data$features.plot<-as.character(data$features.plot)
data<-data[order(data$celltype),]
data$features.plot<-as.factor(data$features.plot)

# ## 读取marker基因
marker <- read.table("marker.txt",sep = "\t",header = T)

## 提取画图数据
scRNA@active.ident <- scRNA$celltype
genes <- intersect(rownames(scRNA@assays$RNA@data),unique(marker$marker))
p <- DotPlot(scRNA, features = genes, assay="RNA", group.by="celltype") +
  theme(axis.text.x = element_text(angle=45, vjust=1))
data <- p$data
data <- dplyr::filter(data, data$pct.exp != 0)
data[!complete.cases(data),]
data <- na.omit(data)
data <- merge(data,marker,by.x = "features.plot",by.y = "marker")

# 控制画图的细胞类型和基因顺序
data$celltype <- factor(data$celltype,levels = levels(data$id))
data$features.plot <- as.character(data$features.plot)
data <- data[order(data$celltype),]
data$features.plot <- as.factor(data$features.plot)


# marker点图--形式1 -----------------------------------------------------------


## 绘制气泡图
pdf("5.anno/Annotation_XX_celltype_dotplot4.pdf",width = 10,height = 7)
p1 <- ggplot(data, aes(id,features.plot, fill =avg.exp.scaled, size=pct.exp))+
  geom_point(shape = 21, colour="black")+
  theme(axis.text.x = element_text(angle = 45,vjust= 1, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        strip.background.x = element_rect(color="black", size=1.5, linetype="solid"),
        strip.background.y = element_rect(color="black", size=1.5, linetype="solid"),
        axis.line = element_line(color="black",size=0.5),
        strip.text.y.left = element_text(angle=0, face="bold"),
        panel.spacing.x = unit(0.2, "cm"),
        panel.spacing.y = unit(0.2, "cm")
  )+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  facet_grid(celltype~., scales = "free", space="free", switch = "y")+
  labs(x = "Celltype", y = "Marker_genes")+
  scale_y_discrete(position = "right")
print(p1)
while (!is.null(dev.list()))  dev.off()

#
#
# # marker点图--形式2 -----------------------------------------------------------
# ## 绘制气泡图
p <- ggplot(data, aes(x=features.plot, y = as.numeric(id), size = pct.exp, color = avg.exp.scaled))+
  geom_point() +
  scale_size("% detected", range = c(0,3)) + #调整绘图点的相对大小
  scale_color_gradientn(colours = viridis::viridis(20),
                        guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"),
                        name = "Average\nexpression") +
  cowplot::theme_cowplot() +ylab("") + xlab("Markers") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(data$id)), labels = levels(data$id), sec.axis = dup_axis())+ #复制 y轴 代替边框效果
  facet_grid(~celltype, scales="free_x", space = "free")+
  theme_classic() +
  theme(axis.text.x = element_text(size=6, angle=45, hjust=1, vjust=1, color="black", face="plain"),#x轴标签样式
        axis.text.y = element_text(size=7, color="black",face="plain"),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),#坐标轴刻度
        axis.text.y.right = element_blank(),#坐标轴标签隐藏
        axis.ticks.x = element_line(size = 0.5,colour = 'grey30'),
        axis.line = element_line(colour = 'grey30',size = 0.2), #坐标轴轴线样式
        panel.spacing=unit(0, "mm"), #分面图图块间距
        legend.text = element_text(size=6, color="black",face="plain"),
        legend.title = element_text(size=6, color="black",face="plain"),
        legend.key.width = unit(0.1,"cm"), legend.key.height = unit(0.4,"cm"),
        legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"), legend.box.spacing = unit(3,"pt"),
        legend.background = element_blank(),legend.key = element_blank(),
        strip.text.x = element_text(size=6, face="plain", color = "#FFFFFF",
                                    vjust = 0.5,margin = margin(b = 3,t=3)),#分面标签样式
        strip.background = element_rect(colour="grey30", fill="grey60", size = 0.1)
  )
p
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
cols <- met.brewer("Juarez", n=length(unique(scRNA$celltype))) %>% as.vector()
cols <- brewer.pal(length(unique(scRNA$celltype)),"Set3") %>% as.vector()
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
plot(g)
ggsave("5.anno/Annotation_OV_celltype_dotplot4.pdf",g,width = 19,height = 8,units = "cm")
#
#
# #####marker表达特征图
# ## 需绘图的基因
marker<- c('PRF1',
           'CD3D',
           'CD69',
           'EPHX2',
           'SAMHD1',
           'CHD7',
           'MZB1',
           'LCN2',
           'KRT20',
           'DCN',
           'CALD1',
           'EPCAM',
           'KRT8',
           'FABP1',
           'TSC2',
           'CD74',
           'CPA3',
           'IGFBP7'
)

## 绘制marker表达特征图
for (gene in marker) {
  p <- FeaturePlot(scRNA,features = gene,reduction = "tsne")
  data <- p$data
  colnames(data)[4] <- "geneX"
  module_plot <- ggplot(data,aes(x = tSNE_1,y =tSNE_2,colour = geneX))+
    geom_point(size = 0.1)+
    ggtitle(gene)+
    scale_color_gradientn(values = seq(0,1,0.2),colours = rev(brewer.pal(11,"Spectral")))+
    theme(plot.title = element_text(size = 14,face = "plain",hjust = 0.5),
          plot.subtitle = element_text(size = 12,face = "plain",hjust = 0.5),
          plot.background = element_blank(),plot.margin = margin(t=4,r=4,b=4,l=4,unit="pt"),
          panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.background = element_blank(),legend.key = element_blank(),
          legend.title = element_blank(),
          legend.position = "right",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),
          legend.box.spacing = unit(4,"pt"))
  pdf(paste(gene,"_tSNE.pdf",sep = ""),width = 6,height = 6)
  plot(module_plot)
  while (!is.null(dev.list()))  dev.off()
}


##细胞类型占比统计

###===========================细胞数目与百分比统计======================================
data <- dplyr::mutate(as.data.frame(table(data.frame(Samples=scRNA$class,Celltype=scRNA$celltype))))
colnames(data) <- c("Samples","Celltype","cell_num")
groupnum <- dplyr::mutate(as.data.frame(table(Groupnum=scRNA$class)))
colnames(groupnum) <- c("Samples","Groupnum")
data <- merge(data,groupnum,by.x ="Samples",by.y="Samples")
data$Cell_ratio <- data$cell_num/data$Groupnum
data$Cell_percent <- scales::percent(data$Cell_ratio, 0.01)

write.csv(data,"6.Cell_statistics/cell_percent.csv",quote = F,row.names = F)

cols <- met.brewer("Egypt",n=length(unique(scRNA$celltype)))

## ========================饼图============================


## 对照组饼图
data1 <- subset(data,subset = Samples %in% c("control"))
pdf("6.Cell_statistics/Ctrl_cell_percent_pieplot2.pdf",width = 6,height = 6)
pie(data1$cell_num,labels = data1$Cell_percent, border="white",col =cols,radius=0.8)
while (!is.null(dev.list()))  dev.off()

## 实验组饼图
data2 <- subset(data,subset = Samples %in% c("disease"))
pdf("6.Cell_statistics/Case_cell_percent_pieplot2.pdf",width = 6,height = 6)
pie(data2$cell_num,labels = data2$Cell_percent, border="white",col =cols,radius=0.8)
while (!is.null(dev.list()))  dev.off()


#===========柱状图================================
## 无标签柱状图
p1 <- ggplot(data,aes(x=Samples,y=cell_num,fill=Celltype))+
  geom_bar(position = "fill",stat="identity",alpha = 0.8)+
  scale_fill_manual(values=met.brewer("Juarez",n=length(unique(scRNA$celltype))))+
  guides(fill=guide_legend(title=NULL))+
  theme(plot.title = element_blank(),plot.subtitle = element_blank(),
        plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
        panel.border = element_rect(size = 0.5,fill = NA,colour = "grey40"),
        panel.background = element_blank(),panel.grid = element_blank(),
        axis.title = element_blank(),axis.ticks = element_line(size = 0.5,colour = "grey40"),
        axis.line = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "right",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),
        legend.box.spacing = unit(3,"pt"))

## 有标签柱状图
p2 <- ggplot(data0,aes(x=Samples,y=Cell_ratio,fill=Celltype))+
  geom_bar(stat ="identity",alpha=0.8)+
  geom_text(aes(label=Cell_percent),position = position_stack(vjust =0.5),size=3)+
  scale_fill_manual(values=met.brewer("Juarez",n=length(unique(scRNA$celltype))))+
  labs(y="Percentage")+
  theme(plot.title = element_blank(),plot.subtitle = element_blank(),
        plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
        panel.border = element_rect(size = 0.5,fill = NA,colour = "grey40"),
        panel.background = element_blank(),panel.grid = element_blank(),
        axis.title = element_blank(),axis.ticks = element_line(size = 0.5,colour = "grey40"),
        axis.line = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "right",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),
        legend.box.spacing = unit(3,"pt"))

pdf("6.Cell_statistics/cell_percentage2.pdf",width = 6,height = 4)
print(p1)
while (!is.null(dev.list()))  dev.off()

##============AUCell计算细胞的表型活性分数========================================

scRNA <- readRDS("3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution0.8_anno.rds")

## 准备基因集
genelist <- read.table("端粒相关基因列表_TelNet人类.txt", check.names = FALSE,sep = "\t",header = T)[,1] # 2089
genes <- unique(genelist)
genes <- capitalize(genes)
genes <- intersect(genes,rownames(scRNA))

###============AUCell: 计算单细胞转录组的每个细胞中特定基因集的活性程度======================


## 1.AUCell_buildRankings: 对每一个细胞的每一个基因进行排名，得到一个排名的矩阵
cells_rankings <- AUCell_buildRankings(scRNA@assays$RNA@data)
cells_rankings

## 2.AUCell_calcAUC: 为了确定基因集是否在每个细胞的基因排名顶部富集，
cells_AUC <- AUCell_calcAUC(genes,cells_rankings,aucMaxRank = nrow(cells_rankings)*0.1)

## 3.AUCell_exploreThresholds：区分较高和较低的AUC值，也就是自动划定双峰分布的阈值
cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist = T,nCores = 1,assignCells = T)
saveRDS(cells_AUC,"7.AUCell/AUCell_calcAUC.rds")

#查看阈值,selected为选中的阈值
cells_assignment$geneSet$aucThr

selectedThresholds <- cells_assignment$geneSet$aucThr$thresholds['Global_k1','threshold']
# 基因集名称
geneSetName <- rownames(cells_AUC)[grep("geneSet", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=selectedThresholds)
abline(v=selectedThresholds)

## 满足新阈值的细胞id
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,] > selectedThresholds))
length(newSelectedCells)#13071
head(newSelectedCells)
table(scRNA@meta.data[newSelectedCells,]$celltype)
write.csv(data.frame(cells = newSelectedCells),paste("7.AUCell/newSelectedCells",selectedThresholds,".csv",sep = ""))

## 调整阈值并查看新阈值下的细胞分布与基因集活性
if (exists("selectedThresholds") == "FALSE") {
  selectedThresholds <- getThresholdSelected(cells_assignment)
}
names(selectedThresholds) <- "geneSet"
# #查看在所选阈值下处于on状态的细胞数目和细胞类型
# length(cells_assignment$geneSet$assignment) 
# table(scRNA@meta.data[cells_assignment$geneSet$assignment,]$celltype)
# 
# 
# # 基因集名称
# geneSetName <- rownames(cells_AUC)[grep("geneSet", rownames(cells_AUC))]
# 
# ## 满足阈值的细胞id
# SelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,] > getThresholdSelected(cells_assignment)))
# length(SelectedCells)
# head(SelectedCells)
# write.csv(data.frame(cells = SelectedCells),
#           paste("SelectedCells",getThresholdSelected(cells_assignment),".csv",sep = ""))


#获取每个细胞的AUC值,导入元数据
aucs <- as.numeric(getAUC(cells_AUC))
scRNA$AUC <- aucs



# 细胞坐标
cellsTsne <- scRNA@reductions$tsne@cell.embeddings
#plot(cellsTsne, pch=16, cex=.3)
cellsUMAP <- scRNA@reductions$umap@cell.embeddings


if (reduc == "tsne") {
  ###========================tsne===============================
  # 获取细胞的二维坐标信息
  df<- data.frame(scRNA@meta.data, scRNA@reductions$tsne@cell.embeddings)
  
  #细胞类型中心位置，为画图的标签做准备
  class_avg <- df %>%
    group_by(celltype) %>%
    summarise(
      tSNE_1 = median(tSNE_1),
      tSNE_2 = median(tSNE_2)
    )
  
  #绘制geneSet活性二维图
  p <- ggplot(df, aes(tSNE_1, tSNE_2))  +
    geom_point(aes(colour = AUC),size=0.3) + viridis::scale_color_viridis(option="A") +
    theme(plot.background = element_blank(),panel.background = element_blank(),
          legend.position = "none",legend.background = element_blank(), 
          legend.text = element_text(size = 7),legend.title = element_text(size = 7),
          axis.title = element_text(size = 9,hjust = 0.5),
          axis.text = element_text(size = 7,hjust = 0),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5)) + 
    ggrepel::geom_label_repel(aes(label = celltype),
                              data = class_avg,
                              size = 3,
                              family = "sans",
                              segment.color = NA)
  pdf("7.AUCell/genelist_AUC_score_tsne.pdf",width = 5,height = 4.5)
  print(p)
  while (!is.null(dev.list()))  dev.off()
  
  ## 绘制直方图、过阈值的细胞分布图、geneSet活性图
  pdf(paste("7.AUCell/geneSet_AUC_score",round(getThresholdSelected(cells_assignment),2),"barplot-tsne.pdf",sep = "_"),width = 13,height = 4.5)
  par(mfrow=c(1,3))
  AUCell_plotTSNE(tSNE=cellsTsne, exprMat=scRNA@assays$RNA@data,
                  cellsAUC=cells_AUC, 
                  thresholds=getThresholdSelected(cells_assignment))
  while (!is.null(dev.list()))  dev.off()
  
}else{
  ###========================UMAP===============================
  # 获取细胞的二维坐标信息
  df<- data.frame(scRNA@meta.data, scRNA@reductions$umap@cell.embeddings)
  
  #细胞类型中心位置，为画图的标签做准备
  class_avg <- df %>%
    group_by(celltype) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  
  #绘制geneSet活性二维图
  p <- ggplot(df, aes(UMAP_1, UMAP_2))  +
    geom_point(aes(colour = AUC),size=0.3) + viridis::scale_color_viridis(option="A") +
    theme(plot.background = element_blank(),panel.background = element_blank(),
          legend.position = "none",legend.background = element_blank(), 
          legend.text = element_text(size = 7),legend.title = element_text(size = 7),
          axis.title = element_text(size = 9,hjust = 0.5),
          axis.text = element_text(size = 7,hjust = 0),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5)) + 
    ggrepel::geom_label_repel(aes(label = celltype),
                              data = class_avg,
                              size = 3,
                              family = "sans",
                              segment.color = NA)
  pdf("7.AUCell/genelist_AUC_score_umap.pdf",width = 5,height = 4.5)
  print(p)
  while (!is.null(dev.list()))  dev.off()
}

selectedThresholds <- cells_assignment$geneSet$aucThr$thresholds['Global_k1','threshold']
df$cellColor<-ifelse(df$AUC>=selectedThresholds,'on','off')
df$cellColor<-factor(df$cellColor,levels = c('on','off'))

p1<-ggplot(df, aes(tSNE_1, tSNE_2))  +
  geom_point(aes(colour = cellColor),size=0.3) + 
  theme(plot.background = element_blank(),panel.background = element_blank(),
        legend.background = element_blank(), 
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),
        axis.title = element_text(size = 9,hjust = 0.5),
        axis.text = element_text(size = 7,hjust = 0),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_line(size = 0.5)) + 
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 3,
                            family = "sans",
                            segment.color = NA)

pdf(paste0("7.AUCell/genelist_on_",reduc,".pdf"),width = 5,height = 4.5)
p1
dev.off()

ggplot(df,aes(x = celltype,y=AUC,color=celltype))+geom_boxplot()+theme_bw()+
  #geom_jitter(aes(fill=celltype),width =0.2,shape = 21,size=2.5)+
  labs(y='AUCell score')+
  theme(axis.text.x = element_text(angle = 90,hjust = 0.5))
ggsave('7.AUCell/genelist_score_boxplot.pdf',width = 6,height = 4.5)

# 手动调阈值 -------------------------------------------------------------------

###================手动调整阈值tsne(可选)================================
## 绘制特定基因集的条型图，设定新阈值
selectedThresholds <- cells_assignment$geneSet$aucThr$thresholds['Global_k1','threshold']
# 基因集名称
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=selectedThresholds)
abline(v=selectedThresholds)

## 满足新阈值的细胞id
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,] > selectedThresholds))
length(newSelectedCells)
head(newSelectedCells)
write.csv(data.frame(cells = newSelectedCells),paste("7.AUCell/newSelectedCells",selectedThresholds,".csv",sep = ""))



## 调整阈值并查看新阈值下的细胞分布与基因集活性
if (exists("selectedThresholds") == "FALSE") {
  selectedThresholds <- getThresholdSelected(cells_assignment)
}
names(selectedThresholds) <- "geneSet"
pdf(paste("7.AUCell/geneSet_AUC_score",selectedThresholds,"barplot-tsne.pdf",sep = "_"),width = 13,height = 4.5)
par(mfrow=c(1,3))
AUCell_plotTSNE(tSNE=cellsTsne, exprMat=scRNA@assays$RNA@data,
                cellsAUC=cells_AUC, 
                thresholds=selectedThresholds)
while (!is.null(dev.list()))  dev.off()



#===========================================拟时序分析==================================================
scRNA <- readRDS("3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution0.8_anno.rds")
SelectedCells <- read.csv("7.AUCell/newSelectedCells0.224142264724898.csv",row.names = 1)[,1]
scRNA1 <- subset(scRNA,cells = SelectedCells)
scRNA1$celltype <- as.character(scRNA1$celltype)
table(scRNA1$celltype)
sp <- "org.Hs.eg.db" #人类："org.Hs.eg.db"；小鼠：org.Mm.eg.db

###============================拟时序分析monocle2=========================================

####一、导入数据，创建cds对象####
#表达矩阵
Mono_matrix <- as(as.matrix(GetAssayData(scRNA1, slot = "counts")), 'sparseMatrix')
dim(Mono_matrix)
#基因注释（基因的各种名字）
feature_ann <- data.frame(gene_id=rownames(Mono_matrix), gene_short_name = rownames(Mono_matrix))
rownames(feature_ann) <- rownames(Mono_matrix)
Mono_fd <- new("AnnotatedDataFrame", data = feature_ann)
#细胞信息（细胞类型、病人分组等）
sample_ann <- scRNA1@meta.data
Mono_pd <- new("AnnotatedDataFrame", data = sample_ann)
#创建对象
Mono.cds <- newCellDataSet(Mono_matrix,phenoData = Mono_pd, featureData = Mono_fd, 
                           expressionFamily = negbinomial.size())
head(pData(Mono.cds))
head(fData(Mono.cds))

####二、归一化####
Mono.cds <- estimateSizeFactors(Mono.cds)#size facotr帮助标准化细胞之间的mRNA的差异
Mono.cds <- estimateDispersions(Mono.cds)#离散度值帮助进行后续的差异分析

##计算每个基因在多少细胞中表达
Mono.cds=detectGenes(Mono.cds,min_expr = 0.1)#这一操作会在fData(Mono.cds)中添加一列num_cells_expressed
nrow(fData(Mono.cds)) 
nrow(pData(Mono.cds)) 
saveRDS(Mono.cds,"8.monocle/Mono.cds.rds")
expressed_genes <- row.names(subset(fData(Mono.cds),num_cells_expressed >= 10))
length(expressed_genes)

####三、轨迹定义基因选择（两种方法）####
dir.create("8.monocle/monocle_celltype_DEG")
dir.create("8.monocle/monocle_dispersion_gene")

## 1.选择celltype差异表达基因
Mono.cds <- readRDS("8.monocle/Mono.cds.rds")
celltype_DEG_genes <- differentialGeneTest(Mono.cds[expressed_genes,],fullModelFormulaStr = '~celltype',cores = 8)
write.csv(celltype_DEG_genes,"8.monocle/monocle_celltype_DEG/celltype_DEG_genes.csv")
ordering_genes <- row.names(celltype_DEG_genes)[order(celltype_DEG_genes$qval)][1:1000] #排序后取前1000个基因
mycds <- setOrderingFilter(Mono.cds,ordering_genes = ordering_genes)
plot_ordering_genes(mycds)#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)
ggsave("8.monocle/monocle_celltype_DEG/ordering_genes_plot.pdf",width = 7,height = 5)
saveRDS(mycds,"8.monocle/monocle_celltype_DEG/mycds.rds")

## 2.选择离散度高的基因
Mono.cds <- readRDS("8.monocle/Mono.cds.rds")
disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
ordering_genes <- unsup_clustering_genes$gene_id
mycds <- setOrderingFilter(Mono.cds, ordering_genes = ordering_genes)
plot_ordering_genes(mycds)
ggsave("8.monocle/monocle_dispersion_gene/ordering_genes_plot.pdf",width = 7,height = 5)
saveRDS(mycds,"8.monocle/monocle_dispersion_gene/mycds.rds")

####四、用DDrtree方法进行降维并在拟时间内排列细胞####
f <- c("8.monocle/monocle_celltype_DEG/",
       "8.monocle/monocle_dispersion_gene/")
for (r in f) {
  mycds <- readRDS(paste(r,"mycds.rds",sep = ""))
  # residualModelFormulaStr减少因素class的影响
  mycds1 <- reduceDimension(mycds,reduction_method = "DDRTree",max_components = 2,
                            norm_method = 'log',
                            residualModelFormulaStr = "~class")#如果没有class，可以改成orig.ident
  Mono.cds1 <- orderCells(mycds1)
  saveRDS(Mono.cds1, paste(r,"Mono.cds_rmf.class.rds",sep = ""))
  
  # 不设置residualModelFormulaStr
  mycds2 <- reduceDimension(mycds,reduction_method = "DDRTree",max_components = 2,
                            norm_method = 'log')
  Mono.cds2 <- orderCells(mycds2)
  saveRDS(Mono.cds2,paste(r,"Mono.cds_no.rmf.rds",sep = ""))
}

# 调整根部
# Mono.cds1 <- orderCells(Mono.cds1,root_state = 3)


#==================拟时序可视化=========================
## 遍历4个cds对象，绘制轨迹图和柱状图
f <- c("8.monocle/monocle_celltype_DEG/",
       "8.monocle/monocle_dispersion_gene/")
m <- c("Mono.cds_rmf.class.rds",
       "Mono.cds_no.rmf.rds")

for (r in f) {
  for (n in m) {
    dirn <- paste(r,strsplit(n,".rds"),sep = "")
    dir.create(dirn)
    Mono.cds1 <- readRDS(paste(r,n,sep = ""))
    ## 1.不同着色绘图####
    plot_cell_trajectory(Mono.cds1, color_by = "Pseudotime",
                         cell_size = 0.8)
    ggsave(paste(dirn,"monocle_Pseudotime.pdf",sep = "/"),width = 5,height = 5)
    
    plot_cell_trajectory(Mono.cds1, color_by = "class",show_backbone=FALSE,
                         cell_size = 0.8)+ ggsci::scale_color_npg()
    ggsave(paste(dirn,"monocle_class.pdf",sep = "/"),width = 5,height = 5)
    
    plot_cell_trajectory(Mono.cds1, color_by = "State",show_backbone = FALSE,
                         cell_size = 0.8) + ggsci::scale_color_npg()
    ggsave(paste(dirn,"monocle_State.pdf",sep = "/"),width = 5,height = 5)
    library(ggsci)
    plot_cell_trajectory(Mono.cds1, color_by = "celltype",show_backbone = FALSE,
                         cell_size = 0.8) + 
      scale_color_manual(values = colorRampPalette(pal_npg("nrc")(10))(length(unique(Mono.cds1$celltype))))
    ggsave(paste(dirn,"monocle_celltype.pdf",sep = "/"),width = 8,height = 8)
    
    
    ## 2.不同状态下细胞的数量和比例柱状图####
    ph <- Mono.cds1@phenoData@data
    ph$celltype <- as.character(ph$celltype)
    data1 <- ph[,c("State","class")]
    data2=melt(data.frame(table(data1)))
    ggplot(data2,aes(class,value,fill=State))+
      geom_bar(stat="identity",position="fill")+ #position不设置时：同列显示频数；position="fill":同列显示比例；position=position_dodge(0.9)：分列显示频数或频率
      scale_fill_manual(values=met.brewer("Juarez",n=length(unique(Mono.cds1$State))))+
      xlab("State") + ylab("")+ labs(fill = "State")+ 
      theme_classic()
    ggsave(paste(dirn,'class_State_freq_bar.pdf',sep = "/"),width = 5,height = 4) 
    
    data1 <- ph[,c("celltype","State")]
    data2=melt(data.frame(table(data1)))
    ggplot(data2,aes(State,value,fill=celltype))+
      geom_bar(stat="identity",position="fill")+ #position不设置时：同列显示频数；position="fill":同列显示比例；position=position_dodge(0.9)：分列显示频数或频率
      scale_fill_manual(values=met.brewer("Juarez",n=length(unique(Mono.cds1$celltype))))+
      xlab("State") + ylab("")+ labs(fill = "celltype")+ 
      theme_classic()
    ggsave(paste(dirn,'State_celltype_freq_bar.pdf',sep = "/"),width = 5,height = 4)
  } 
}

# 拟时间相关差异基因 ---------------------------------------------------------------
# 从以下四个对象中选择一个拟时序轨迹和细胞占比状况较好的cds对象，进行分析：
# ./monocle_cluster_DEG/Mono.cds_rmf.class.rds (或 ./monocle_celltype_DEG/Mono.cds_rmf.class.rds)
# ./monocle_cluster_DEG/Mono.cds_no.rmf.rds (或 ./monocle_celltype_DEG/Mono.cds_no.rmf.rds)
# ./monocle_dispersion_gene/Mono.cds_rmf.class.rds 
# ./monocle_dispersion_gene/Mono.cds_no.rmf.rds
## 根据前面可视化内容，选择一个cds对象，进行后续拟时序相关分析

Mono.cds1  <- readRDS("8.monocle/monocle_celltype_DEG/Mono.cds_no.rmf.rds")
celltype_DEG_genes <- read.csv('8.monocle/monocle_celltype_DEG/celltype_DEG_genes.csv',row.names = 1)
ordering_genes <- row.names(celltype_DEG_genes)[order(celltype_DEG_genes$qval)][1:1000] #排序后取前1000个基因
## 3.探索拟时序上基因的显著性####
pseudotime_de <- differentialGeneTest(Mono.cds1[expressed_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)") 
write.csv(pseudotime_de,"8.monocle/monocle_celltype_DEG/Mono.cds_no.rmf/pseudotime_DEG_genes.csv")
pseudotime_de1 <- subset(pseudotime_de, qval < 1e-4) %>% pull(gene_short_name) %>% as.character()#取q值小于0.0001的基因
#pseudotime_de1 <- top_n(pseudotime_de, n = 20, desc(qval)) %>% pull(gene_short_name) %>% as.character()

###a.热图####
p = plot_pseudotime_heatmap(Mono.cds1[pseudotime_de1,],return_heatmap=T,
                            hmcols = colorRampPalette(rev(brewer.pal(11, "Spectral")))(70),
                            num_clusters = 2,cluster_rows = T,
                            show_rownames = F)
pdf("8.monocle/monocle_celltype_DEG/Mono.cds_no.rmf/pseudotime_DEG_heatmap.pdf",width = 5,height = 5)
print(p)
while (!is.null(dev.list()))  dev.off()
#提取聚类信息
clusters <- cutree(p$tree_row, k = 2) %>% data.frame()
clusters[,1] <- as.character(clusters[,1])
colnames(clusters) <- "Gene_Clusters"
clusters$gene <- rownames(clusters)
table(clusters$Gene_Clusters)
write.csv(clusters,"8.monocle/monocle_celltype_DEG/Mono.cds_no.rmf/pseudotime_de_cluster_gene.csv",row.names = T)

###b.GO分析####
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
allcluster_go=data.frame()
for (i in unique(clusters$Gene_Clusters)) {
  small_gene_group=filter(clusters,clusters$Gene_Clusters==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=sp)
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = sp,
                 keyType       = 'ENTREZID',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
write.csv(allcluster_go,"8.monocle/monocle_celltype_DEG/Mono.cds_no.rmf/GO_pseudotime_de_cluster.csv")

# BEAM分析 ------------------------------------------------------------------
## 根据前面可视化内容，选择一个cds对象，进行后续拟时序相关分析
Mono.cds1  <- readRDS("8.monocle/monocle_celltype_DEG/Mono.cds_no.rmf.rds")

## 4.BEAM分析：与结点相连的分支之间基因变化情况####
BEAM_res <- BEAM(Mono.cds1[ordering_genes,], branch_point = 1, cores = 2)#选取了前面筛选的基因做BEAM分析
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "8.monocle/BEAM/BEAM_res.rds")

###a.热图####
pdf("8.monocle/BEAM/BEAM_branchPoint1_heatmap.pdf",width = 5,height = 5)
tmp1=plot_genes_branched_heatmap(
  Mono.cds1[row.names(subset(BEAM_res,qval < 1e-4)),],
  branch_point = 1,#绘制的是哪个分支
  num_clusters = 4, #分成几个cluster，根据需要调整
  cores = 1,
  branch_labels = c("Cell fate 1", "Cell fate 2"),#分支
  hmcols = colorRampPalette(rev(brewer.pal(11, "Spectral")))(70),
  branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
  use_gene_short_name = T,
  show_rownames = F,#基因很多
  return_heatmap = T)
while (!is.null(dev.list()))  dev.off()
#提取聚类信息
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)
write.csv(gene_group,"8.monocle/BEAM/BEAM_cluster_gene.csv",row.names = T)

###b.计算热图中的基因富集通路(GO)####
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=sp)
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = sp,
                 keyType       = 'ENTREZID',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
write.csv(allcluster_go,"8.monocle/BEAM/GO_branchPoint1_cluster.csv")

allcluster_go <- read.csv('8.monocle/BEAM/GO_branchPoint1_cluster.csv')
i = 1
i = 2
i = 3
i = 4
print(i)
go <- allcluster_go[allcluster_go$cluster == i,c("ID","Description","p.adjust","cluster")]
go <- go[order(go$p.adjust),]
View(go)




## 5.基因按照热图顺序排序####
hp.genes <- tmp1$ph_res$tree_row$labels[tmp1$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(BEAM_sig, "8.monocle/BEAM/BEAM_sig.csv", row.names = F)







###===========================细胞通讯分析=================================
table(scRNA$celltype,scRNA$class)


##=============单样本分析（可以是同一组的多个生物学重复一起分析）==========
scRNA <- readRDS("3.RDS/GSE184880_Ovaries_filtered_combine_pca17_resolution_0.8_anno.rds") 
sp <- "human" #human/mouse
Idents(scRNA) <- "class"
cl <- c("control","disease")#左边对照组，右边实验组

###===========================细胞通讯分析==================================

### 一、创建cellchat对象####
sc_disease <- subset(scRNA,idents=cl[2])
sc_disease$celltype<-as.character(sc_disease$celltype)
cc_disease <- createCellChat(sc_disease@assays$RNA@data,meta = sc_disease@meta.data,group.by = "celltype")
cellchat <- cc_disease

cellchat <- setIdent(cellchat, ident.use = "celltype")
groupSize <- as.numeric(table(cellchat@idents))

## 导入配体受体数据库
# CellChatDB 是一个手动整理的文献支持的配体受体在人和小鼠中的交互数据库
# Secreted Signaling:自分泌/旁分泌信号相互作用;
# ECM:细胞外基质;
# Cell-Cell Contact细胞-细胞接触相互作用
if (sp == "human") {
  CellChatDB <- CellChatDB.human
}else{
  CellChatDB <- CellChatDB.mouse
}
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

### 二、预处理用于细胞通信分析的表达数据####
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
if (sp == "human") {
  cellchat <- projectData(cellchat, PPI.human)
}else{
  cellchat <- projectData(cellchat, PPI.mouse)
}

### 三、推断细胞通讯网络####
# 计算细胞群间的配受体水平的通信概率/强度   
cellchat <- computeCommunProb(cellchat)  #<<<耗时较长1-2h>>>
# 过滤细胞数很少的细胞群通讯
cellchat <- filterCommunication(cellchat, min.cells = 10) #Fibrolast cell
df.net <- subsetCommunication(cellchat,slot.name = "net")
write.csv(df.net,"9.cellchat/net_lr.csv")
# 汇总所有相关配体/受体，计算信号通路水平上的通信概率
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.netp,"9.cellchat/net_pathway.csv")
# 整合通讯网络
cellchat <- aggregateNet(cellchat)

### 四、识别全局通信模式####
# 计算网络中心性得分，以识别所有推断的通信网络中的主要发送者、接收者、媒介和影响者
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat,"9.cellchat/cellchat.rds")

######==============================可视化============================
#看不同细胞类群间的互作数量和强度 
cols <- met.brewer("Juarez",n=length(unique(cellchat@meta$celltype)))
###==============================circleplot=====================================
pdf("9.cellchat/Circle_celltype_count.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= T,color.use = cols, title.name = "Number of interactions")
dev.off()
pdf("9.cellchat/Circle_celltype_weight.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F,color.use = cols, title.name = "Interaction weights/strength")
dev.off()

dir.create('./9.cellchat/source')
for (i in unique(row.names(cellchat@net$prob))) {
  pdf(paste0("9.cellchat/source/Circle_source_",i,".pdf"),width = 4,height = 5)
  netVisual_circle(cellchat@net$weight,sources.use = i, vertex.weight = groupSize, color.use = cols,weight.scale = T, label.edge= F, title.name = paste0("Source：",i," Interaction weights"))
  dev.off()
}



#看各个细胞类群的Signaling pathway情况
#============1. 通过heatmap看outgoing和incoming的communication probability========
cellchat<-netAnalysis_computeCentrality(cellchat)

pdf("9.cellchat/Heatmap_celltype_outgoing2.pdf")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",color.use = cols,font.size = 5.6)
dev.off()

pdf("9.cellchat/Heatmap_celltype_incoming2.pdf")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.use = cols,font.size = 5.6)
dev.off()



# 2. 看特定细胞群体之间的互作 ####
levels(cellchat@idents)
su <- c(8)
tu <- c(1,2,3,4,5,6,7,9,10)
# a. 通过气泡图看互作————选择好sources.use和targets.use，展示对应细胞群体间所有互作
pdf("9.cellchat/Bubble_celltype.pdf",width = 5,height = 12)
netVisual_bubble(cellchat, sources.use = su, targets.use = tu, remove.isolate = FALSE,
                 font.size = 10,font.size.title = 10)
while (!is.null(dev.list()))  dev.off()

su <- c(8)
tu <- c(1,2,3,4,5,6,7,9,10)
# b. 通过弦图看互作————targets.use可以设置为多个receiver
pdf("9.cellchat/Chord_celltype.pdf",width = 10,height = 10)
netVisual_chord_gene(cellchat, sources.use = su, targets.use = tu, lab.cex = 0.5,
                     color.use = cols,legend.pos.y = 10,legend.pos.x = 10)
while (!is.null(dev.list()))  dev.off()


# 3.选择特定pathway在特定细胞群体之间的互作####
df.net <- subsetCommunication(cellchat,slot.name = "net")
# 填写选中的配体
choice_pathway <- df.net$pathway_name[which(df.net$ligand %in% c("NAMPT","SPP1"))] %>% unique()

pdf("9.cellchat/Chord_celltype_pathway.pdf")
netVisual_chord_gene(cellchat, sources.use = su, targets.use = tu,lab.cex = 0.8, 
                     small.gap = 2.5,
                     big.gap = 10,
                     signaling = choice_pathway,
                     color.use = cols, legend.pos.y = 10,legend.pos.x = 10) 
while (!is.null(dev.list()))  dev.off()


# 4.看特定pathway中相关基因的表达
for (i in choice_pathway) {
  p <- plotGeneExpression(cellchat, signaling = i,color.use = cols)
  pdf(paste0("9.cellchat/express/Vlnplot_celltype_",i,".pdf"),width = 6,height = 5)
  print(p)
  while (!is.null(dev.list()))  dev.off()
}




# 取交集基因 -------------------------------------------------------------------



# 2. 基因集活化on和off的细胞群差异基因 与 疾病与对照的差异基因  取交集 --------------------------------

m <- scRNA@meta.data
SelectedCells <- read.csv("7.AUCell/newSelectedCells0.224142264724898.csv",row.names = 1)[,1]
m$class2 <- ifelse(rownames(m) %in% SelectedCells,"on","off")
table(m$class2)
scRNA$class2 <- m$class2

## 基因集活性on和off的细胞群差异分析 
Idents(scRNA) <- "class2"
sce.markers4 <- FindMarkers(scRNA, ident.1 = "on",ident.2 = "off",
                            test.use = "wilcox", #默认检验方法
                            logfc.threshold = 0.25, #两细胞群体差异倍数0.25以上
                            min.pct = 0.1, #检测的基因在两细胞群体具有10%以上细胞中表达
                            only.pos = F) #是否只保留上调基因
sce.markers4$cluster <- ifelse(sce.markers4$avg_log2FC>0,"on","off")
sce.markers4$gene <- rownames(sce.markers4)
write.csv(sce.markers4,"10.intersect/DEG_of_AUC_on.off.csv")
write.table(data.frame(gene_name = sce.markers4$gene),"10.intersect/on.off_gene.txt",row.names = F)
##绘制热图
scRNA1 <- subset(scRNA,downsample=100)
deg1 <- sce.markers4 %>% group_by(cluster) %>% top_n(10,wt = abs(avg_log2FC))
DoHeatmap(scRNA1, features = deg1$gene,group.by = "class2",angle = 45,raster = F,size = 0,combine = T) +
  scale_fill_gradientn(colors = c("#4DBBD5", "black", "#F39B7F"))+
  theme(axis.text = element_text(size = 6,face = "plain",colour = "black"),
        legend.text = element_text(size = 6,face = "plain",colour = "black"),
        legend.title = element_text(size = 6,face = "plain",colour = "black"),
        legend.key.height = unit(0.3,"cm"),legend.key.width = unit(0.1,"cm"))
ggsave("10.intersect/DEG_of_AUC_on.off_heatmap7.pdf",width = 8,height = 8,units = "cm")

## 疾病与对照差异分析 
Idents(scRNA) <- "class"
sce.markers3 <- FindMarkers(scRNA, ident.1 = "disease",ident.2 = "control",
                            test.use = "wilcox", #默认检验方法
                            logfc.threshold = 0.25, #两细胞群体差异倍数0.25以上
                            min.pct = 0.1, #检测的基因在两细胞群体具有10%以上细胞中表达
                            only.pos = F) #是否只保留上调基因
sce.markers3$cluster <- ifelse(sce.markers3$avg_log2FC>0,"disease","control")
sce.markers3$gene <- rownames(sce.markers3)
write.csv(sce.markers3,"10.intersect/DEG_of_case.ctrl.csv")
##绘制热图
scRNA2 <- subset(scRNA,downsample=100)
deg2 <- sce.markers3 %>% group_by(cluster) %>% top_n(10,wt = abs(avg_log2FC))
DoHeatmap(scRNA2, features = deg2$gene,group.by = "class",angle = 45,raster = F,size = 0,combine = T) +
  scale_fill_gradientn(colors = c("#4DBBD5", "black", "#F39B7F"))+
  theme(axis.text = element_text(size = 6,face = "plain",colour = "black"),
        legend.text = element_text(size = 6,face = "plain",colour = "black"),
        legend.title = element_text(size = 6,face = "plain",colour = "black"),
        legend.key.height = unit(0.3,"cm"),legend.key.width = unit(0.1,"cm"))
ggsave("10.intersect/DEG_of_case.ctrl_heatmap.pdf",width = 8,height = 8,units = "cm")
dev.off()
## 取交集

# 韦恩图 ---------------------------------------------------------------------
#绘制成 pdf 格式
#将 filename 设置成 NULL
#将画出来的图先保存到 venn.plot中
library(futile.logger)
library(VennDiagram)
library(grDevices)
genes<-read.csv('10.intersect/DEG_of_AUC_on.off.csv',row.names = 1)
genes1<-read.csv('10.intersect/DEG_of_case.ctrl.csv',row.names = 1)
venn.plot = venn.diagram(
  x = list(
    on.off = genes$gene,
    case.ctrl = genes1$gene
  ),  
  cat.col = c("#F39B7F",
              "#4DBBD5"
  ),
  fill = c("#F39B7F",
           "#4DBBD5"
  ),
  filename = NULL
)
#将 venn.plot 通过 grid.draw 画到pdf文件中
pdf("10.intersect/venn.pdf")
grid.draw(venn.plot)
dev.off()

genes <- intersect(rownames(sce.markers4),rownames(sce.markers3))
write.table(data.frame(gene_name = genes),"10.intersect/intersect_gene.txt",row.names = F)
write.table(data.frame(gene_name = genes),"10.intersect/intersect_gene2.txt",row.names = F,quote = FALSE)






# 交集基因GO KEGG富集分析 ---------------------------------------------------------


library(ggplot2)
library(clusterProfiler)
library(GOplot) 
library(dplyr)
library(ggrepel)

#genes为待做富集分析的基因集

sp <- "org.Hs.eg.db"
genes <- read.table('10.intersect/intersect_gene.txt',header = T)[,1]
#ID转换
genes_id <- genes
genes_id <- bitr(genes_id,fromType = "SYMBOL",toType = "ENTREZID",
                 OrgDb = sp,drop = T)
dir.create('10.intersect/GOKEGG')
####GO富集####
ego_n <- enrichGO(gene = genes_id$ENTREZID,
                  OrgDb = sp,
                  keyType = 'ENTREZID',
                  ont = "ALL", 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = T)
ego_list <- data.frame(ego_n)
write.csv(ego_list, "10.intersect/GOKEGG/enrichGO.csv")

####KEGG富集####
library(tidygraph)
kegg_name<-read.table('../转录组部分/prepare/KEGG_pathway_name.txt',sep='\t',header = T)
rownames(kegg_name)<-kegg_name$path_id
genelist <- pull(genes_id,ENTREZID)
ekegg <- enrichKEGG(gene = genelist, 
                    organism = 'hsa',
                    keyType = "kegg",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    use_internal_data = T)
ekegg_list <- data.frame(ekegg)
ekegg_list$Description<-kegg_name[ekegg_list$ID,'path_name']
entrezid=strsplit(ekegg_list$geneID,"/")
symbol=sapply(entrezid,function(x){
  y=bitr(x, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  #一对多，取第一个
  y=y[!duplicated(y$ENTREZID),-1]
  y=paste(y,collapse = "/")
})
ekegg_list$geneID = symbol
write.csv(ekegg_list,"10.intersect/GOKEGG/enrichKEGG.csv")


colors = c("#4DBBD5", "#E64B35", "#00A087")
#GO富集棒棒糖图，展示BP、CC、MF的显著性前三条目
library(stringr)
library(ggpubr)
g<-read.csv("10.intersect/GOKEGG/enrichGO.csv")
g$pvalue <- -log10(g$pvalue)
g1=g %>% group_by(ONTOLOGY) %>% top_n(3,wt = pvalue)
g1plot<-ggdotchart(g1, x="Description", y="pvalue", color = "ONTOLOGY",
                   palette = colors,
                   sorting = "descending",   #上升排序，区别于desc
                   add = "segments",    #增加线段
                   xlab = 'GO',ylab = "-log10(pvalue)",
                   rotate = T,       #横向显示
                   dot.size = 7,        #圆圈大小
                   ggtheme = theme_pubr()+
                     theme(axis.text.y = element_text(size = 12))
)+scale_x_discrete(labels=function(x) str_wrap(x, width=30))

ggsave("10.intersect/GOKEGG/GO.pdf",height = 10,width = 15,units = "cm")

#KEGG富集棒棒糖图，展示显著性前9通路
e<-read.csv("10.intersect/GOKEGG/enrichKEGG.csv")
e$pvalue <- -log10(e$pvalue)
n<-ifelse(nrow(e)>9,9,nrow(e))
e1plot<-ggdotchart(e[1:n,], x="Description", y="pvalue", color = "Description",
                   palette = colorRampPalette(colors = c("#4DBBD5", "#E64B35", "#00A087"))(n),     #配色方案
                   sorting = "descending",   #上升排序，区别于desc
                   add = "segments",    #增加线段
                   add.params = list(color = colorRampPalette(colors = c("#4DBBD5", "#E64B35", "#00A087"))(n), size = 1),
                   xlab = 'KEGG',ylab = "-log10(pvalue)",
                   rotate = T,       #横向显示
                   dot.size = 7,        #圆圈大小
                   ggtheme = theme_pubr()+
                     theme(axis.text.y = element_text(size = 12))
)+NoLegend()+scale_x_discrete(labels=function(x) str_wrap(x, width=35))

ggsave("10.intersect/GOKEGG/KEGG.pdf",height = 10,width = 18,units = "cm")

#基因logFC数据整理
genes1<-read.csv('10.intersect/DEG_of_case.ctrl.csv',row.names = 1)
genes1<-genes1[genes,c("gene","avg_log2FC")]
colnames(genes1) <- c('ID','logFC')
#富集结果整理
ego_list1 <- ego_list[,c('ONTOLOGY','ID','Description','Count','p.adjust','geneID')]
ekegg_list$ONTOLOGY <- "KEGG"
ekegg_list1 <- ekegg_list[,c('ONTOLOGY','ID','Description','Count','p.adjust','geneID')]
go_kegg <- rbind(ego_list1,ekegg_list1)
colnames(go_kegg) <- c('category','ID','term','count','adj_pval','genes')
go_kegg$genes <- gsub("/",",",go_kegg$genes)

####计算z-score####
circ <- circle_dat(go_kegg, genes1)
write.csv(circ,"10.intersect/GOKEGG/go_kegg_zscore.csv")

#选取子集
x <- circ[,c(2,8)]
x <- x[!duplicated(x),]
circ1 <- merge(go_kegg,x,by = "ID")
circ2 <- circ1 %>% arrange(adj_pval) %>%  group_by(category) %>% do(head(.,n = 3))

#柱状图
#GOBar(subset(circ2, category == 'BP'))
GOBar(circ2,display = 'single',order.by.zscore = T)#display = 'multiple'
ggsave("10.intersect/GOKEGG/GOBar.pdf",width = 5.5,height = 7)

#气泡图
source("../转录组部分/prepare/myGOBubble.R") #myGOBubble见后面模块
environment(myGOBubble) <- environment(GOBubble)
#reduced_circ <- reduce_overlap(circ, overlap = 0.75)
myGOBubble(circ, title = '', colour = brewer.pal(4, "Set3"), display = 'single', labels = 3,table.legend = F)  
ggsave("10.intersect/GOKEGG/GOBubble.pdf",width = 6,height = 5)

#圈图
source("../转录组部分/prepare/myGOCircle.R")#myGOCircle见后面模块
environment(myGOCircle) <- environment(GOCircle)
pdf("10.intersect/GOKEGG/GOCircle.pdf",width = 6.7,height = 5)
myGOCircle(circ, table.legend = T,nsub = circ2$ID,label.fontface = "plain",label.size = 1.25)
dev.off()

#弦图
#chord <- chord_dat(data = circ, genes = genes1, process = circ2$ID)
circ3<-circ[circ$ID%in%circ2$ID,]
top2 = circ3 %>% group_by(ID) %>% top_n(2,abs(logFC))
if(length(unique(top2$genes))<10){top2 = circ3 %>% group_by(ID) %>% top_n(3,abs(logFC))
print('genes=3')}
genes3<-genes1[genes1$ID%in%top2$genes,]
chord <- chord_dat(data = circ3, genes = genes3, process = circ2$ID)

source("../转录组部分/prepare/myGOChord.R")#myGOChord见后面模块
environment(myGOChord) <- environment(GOChord)
myGOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,
          limit = c(0, 0),ribbon.col=colorRampPalette(brewer.pal(12, "Set3"))(length(circ2$ID)))
ggsave("10.intersect/GOKEGG/GOChord.pdf",width = 9.2,height = 10.5)





