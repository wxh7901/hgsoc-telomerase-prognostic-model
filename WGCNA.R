####WGCNA####
dir.create('WGCNA')
combat_expr <- read.csv(tumor_exp_path,row.names = 1)
colnames(combat_expr)<-gsub('[.]','-',colnames(combat_expr))
keygene<-read.table(keygene_path,header = T)
phenotype<-colnames(keygene)
exp1=as.matrix(combat_expr)
pheno_symbol<-list(keygene[,1])
names(pheno_symbol)<-phenotype
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)

type <- "unsigned"
corType <- "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

#导入表型数据#
#ssGSEA分组打分
library(GSVA)
min.sz = 2 # 单个基因集基因数下限
max.sz = 10000 # 基因数上限
parallel.sz= 10 # 并行处理cpu数量
mx.diff= T # ES 计算为最大正和负随机游走偏差之间的幅度差。
tau=0.25 # 默认1，tau=1 when method="gsva" and tau=0.25 when method="ssgsea" 
method='ssgsea' # 默认
kcdf = "Gaussian" # 默认"Gaussian"，当输入数值为整数时，设置为"Poisson"

ssgsea = gsva(expr=exp1, gset.idx.list=pheno_symbol, 
              method=method, kcdf=kcdf, min.sz=min.sz, max.sz=max.sz, 
              parallel.sz=parallel.sz, mx.diff=mx.diff,tau=0.25 )

write.csv(ssgsea,'WGCNA/ssGSEA_pheno.csv')
ssgsea_group=t(ssgsea)

dataExpr <- exp1  
dim(dataExpr)        
datTraits <- data.frame(gsm=rownames(ssgsea_group),group=ssgsea_group) #group是表型打分用作分组信息
colnames(datTraits)=c('gsm','group')

#数据筛选#
#选择挑选方差前5000的基因进行后续分析
m.mad <- apply(dataExpr,1,var)    
dataExprVar <- dataExpr[order(m.mad,decreasing = T)[1:5000],]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))   #方差前5000的基因

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK          #TRUE
#如果gsg$allOK的结果为TRUE，证明没有缺失值，可以直接下一步。如果为FALSE，则需要用以下函数进行删除缺失值。

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
head(dataExpr)[,1:8]

## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                colors = c("blue","red"),signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(sampleTree, sample_colors,
                    groupLabels = "Group",
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

#软阈值筛选#
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)
pdf("WGCNA/软阈值筛选.pdf",height = 5.5,width = 8)
par(mfrow = c(1,2))
cex1 = 0.9

# 横轴是Soft threshold(power),纵轴是无标度网络的评估参数,数值越高,网络越符合无标度特征(non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")# 筛选标准 R-square=0.85

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")

dev.off()
power = sft$powerEstimate   #设置软阈值beta值
power

write.csv(sft$fitIndices,"WGCNA/sft_fitIndices.csv",row.names = F,quote = F)

##无满足条件的power时选用经验power
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

#网络构建#
net = blockwiseModules(dataExpr, power = power, 
                       maxBlockSize = nGenes,#计算机能处理的最大模块的基因数量
                       TOMType = type, minModuleSize = 30,     #模块数量太多时，提高模块内最低基因数值
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25, #合并模块的阈值，越大模块越少
                       numericLabels = TRUE, #返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
                       pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("exprMat.tom"),
                       verbose = 3)
table(net$colors)#0(grey)表示未分入任何模块的基因      

f=0
while (ncol(net$MEs)>16){
  f=f+1
  net = blockwiseModules(dataExpr, power = power, 
                         maxBlockSize = nGenes,#计算机能处理的最大模块的基因数量
                         TOMType = type, minModuleSize = 30+10*f,     #模块数量太多时，提高模块内最低基因数值
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25, #合并模块的阈值，越大模块越少
                         numericLabels = TRUE, #返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
                         pamRespectsDendro = FALSE,
                         saveTOMs=TRUE, corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = paste0("exprMat.tom"),
                         verbose = 3)
}


write.csv(data.frame(table(net$colors)),'WGCNA/Module_Nofgenes.csv')
saveRDS(net,"WGCNA/net.rds")

#模块可视化#
##层级聚类树
pdf("WGCNA/基因聚类模块.pdf",height = 6,width = 10)
moduleLabels = net$colors     
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()
##绘制模块之间相关性图
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf("WGCNA/模块间的相关性热图.pdf",height = 6,width = 6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(1,3,2,4),
                      marHeatmap = c(3,3,1,2), # marDendro/marHeatmap设置下左上右的边距
                      plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()

#可视化基因网络 (TOM plot)#
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
#抽样
nSelect <- 400
set.seed(10)
select <- sample(nGenes, size = nSelect)
head(select)
selectTOM <- dissTOM[select, select]
selectTOM[1:6,1:6]

# 根据TOM矩阵重新对基因聚类；
selectTree <- hclust(as.dist(selectTOM), method = "average")
#提取选择基因的模块颜色；
selectColors <- moduleColors[select]
# 对TOM矩阵进行指数转化，使热图能展示更丰富的信息；
plotDiss <- selectTOM^7
diag(plotDiss) <- NA
#绘制热图，默认数值越大（越接近1）颜色越深（底色）；
pdf("WGCNA/基因网络TOM.pdf",height = 6,width = 6)
TOMplot((1-plotDiss), selectTree, selectColors,
        main = "Network heatmap plot, selected genes")
dev.off()

#关联表型数据#
traitData=model.matrix(~0+ datTraits$group)
colnames(traitData)=phenotype

### 模块与表型数据关联
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## 模块与表型相关性热图
textMatrix = paste(signif(modTraitCor, 4), "\n(", signif(modTraitP, 1), ")", 
                   sep = "")# signif表示保留几位小数
dim(textMatrix) = dim(modTraitCor)
pdf("WGCNA/模块和性状的关系.pdf",width = 8,height = 8)
par(mar = c(3,6,4,3))
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = colnames(traitData),xLabelsAngle = 0,
               yLabels = colnames(MEs_col),
               cex.lab = 0.8,
               ySymbols = colnames(MEs_col), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.8, textAdj = c(0.5, 0.5),
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


## 计算模块与基因的相关性矩阵
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

## 计算性状与基因的相关性矩阵
if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

## 两个相关性矩阵联合起来,指定感兴趣模块进行分析
modTraitCor
modTraitP
write.csv(modTraitCor,'WGCNA/modTraitCor.csv')
write.csv(modTraitP,'WGCNA/modTraitP.csv')

#筛选关联性最强的模块名字
#先找出p值大于0.05的模块,再找出关联值最大的模块名字
module_names=data.frame()
for (c in 1:nrow(modTraitCor)){
  if (modTraitP[c]<0.05){
    module_names=rbind(module_names,data.frame(row.names=rownames(modTraitCor)[c],modCor=modTraitCor[c]))
  }
}
module=rownames(module_names)[which(abs(module_names$modCor)==max(abs(module_names$modCor)))]
module=gsub('ME','',module)         #module_trait relationship 中感兴趣的模块名字
pheno = phenotype            #module_trait relationship 中感兴趣的性状名字
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

pdf("WGCNA/模块中基因相关性散点图.pdf",height = 6,width = 6)
#sizeGrWindow(7, 7)
par(mfrow = c(1,1),mar=c(6,6,6,6))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#提取指定模块内的基因名#
gene <- rownames(geneTraitCor)[moduleGenes]    
write.table(gene,"WGCNA/Module_genes.txt",sep='\t',quote = F,row.names = F,col.names = F)
write.csv(gene,"WGCNA/Module_genes.csv")

#方差前5000筛选的与表型最相关模块内的基因，与差异基因取交集
DEG<-read.csv(paste0('DEG/',test,'-',control,'_DEG.csv'),header = T,row.names = 1)
genes_id<-DEG$gene_symbol[DEG$group!='no']
gene1=intersect(gene,genes_id)
write.table(gene1,"WGCNA/关键表型相关基因.txt",sep='\t',quote = F,row.names = F,col.names = F)
write.csv(gene1,'WGCNA/关键表型相关基因.csv')




####GO/KEGG####
dir.create('DEG/GO_KEGG')

library(clusterProfiler)
library(Seurat)
#ID转换，通常采用ENTREZID进行功能富集分析
genes_id <- bitr(gene1,fromType = "SYMBOL",toType = "ENTREZID",
                 OrgDb = 'org.Hs.eg.db',drop = T)
####GO富集
ego_n <- enrichGO(gene = genes_id$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'ENTREZID',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = T)
ego_list <- data.frame(ego_n)
write.csv(ego_list, "DEG/GO_KEGG/DEG_keygene_enrichGO.csv")

####KEGG富集
#KEGG富集没有readable = T参数，需要自行将结果中的ENTREZID换回SYMBOL。
library(tidygraph)
genelist <- pull(genes_id,ENTREZID)
ekegg <- enrichKEGG(gene = genelist,
                    organism = 'hsa',
                    keyType = "kegg",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)
ekegg_list <- data.frame(ekegg)
if(nrow(ekegg_list)!=0){entrezid=strsplit(ekegg_list$geneID,"/")
symbol=sapply(entrezid,function(x){
  y=bitr(x, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  y=y[!duplicated(y$ENTREZID),-1]
  y=paste(y,collapse = "/")
})
ekegg_list$geneID = symbol
write.csv(ekegg_list,"DEG/GO_KEGG/DEG_keygene_enrichKEGG.csv")
}
#GO富集棒棒糖图，展示BP、CC、MF的显著性前三条目
library(stringr)
library(ggpubr)
g<-read.csv("DEG/GO_KEGG/DEG_keygene_enrichGO.csv")
g$pvalue <- -log10(g$pvalue)
g1=g %>% group_by(ONTOLOGY) %>% top_n(3,wt = pvalue)
g1plot<-ggdotchart(g1, x="Description", y="pvalue", color = "ONTOLOGY",
                   sorting = "descending",   #上升排序，区别于desc
                   add = "segments",    #增加线段
                   xlab = 'GO',ylab = "-log10(pvalue)",
                   rotate = T,       #横向显示
                   dot.size = 7,        #圆圈大小
                   ggtheme = theme_pubr()+
                     theme(axis.text.y = element_text(size = 12))
)+scale_x_discrete(labels=function(x) str_wrap(x, width=35))

while (!is.null(dev.list()))  dev.off()
pdf("DEG/GO_KEGG/DEG_keygene_GO.pdf",height = 8,width = 6)
print(g1plot)
dev.off()
#KEGG富集棒棒糖图，展示显著性前9通路
if(nrow(ekegg_list)!=0){e<-read.csv("DEG/GO_KEGG/DEG_keygene_enrichKEGG.csv")
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

while (!is.null(dev.list()))  dev.off()
pdf("DEG/GO_KEGG/DEG_keygene_KEGG.pdf",height = 8,width = 6)
print(e1plot)
dev.off()}

library(GOplot)
library(ggrepel)
#基因logFC数据整理
genes<-DEG[DEG$group!='no',]
genes1 <- genes[,c('gene_symbol','logFC')]
colnames(genes1) <- c('ID','logFC')
#富集结果整理
ego_list1 <- ego_list[,c('ONTOLOGY','ID','Description','Count','p.adjust','geneID')]
if(nrow(ekegg_list)!=0){ekegg_list$ONTOLOGY <- "KEGG"
ekegg_list1 <- ekegg_list[,c('ONTOLOGY','ID','Description','Count','p.adjust','geneID')]
go_kegg <- rbind(ego_list1,ekegg_list1)}else{go_kegg<-ego_list1}
colnames(go_kegg) <- c('category','ID','term','count','adj_pval','genes')
go_kegg$genes <- gsub("/",",",go_kegg$genes)

####计算z-score
circ <- circle_dat(go_kegg, genes1)
write.csv(circ,"DEG/GO_KEGG/DEG_keygene_go_kegg_zscore.csv")

x <- circ[,c(2,8)]
x <- x[!duplicated(x),]
circ1 <- merge(go_kegg,x,by = "ID")

#选取子集
circ2 <- circ1 %>% arrange(adj_pval) %>%  group_by(category) %>% do(head(.,n = 3))

#柱状图
#GOBar(subset(circ2, category == 'BP'))
while (!is.null(dev.list()))  dev.off()
pdf("DEG/GO_KEGG/DEG_keygene_Bar.pdf",width = 5.5,height = 7)
print(GOBar(circ2,display = 'single',order.by.zscore = T))#display = 'multiple'
dev.off()
#气泡图
source("prepare/myGOBubble.R")
#myGOBubble <- edit(GOBubble)
environment(myGOBubble) <- environment(GOBubble)
#reduced_circ <- reduce_overlap(circ, overlap = 0.75)
while (!is.null(dev.list()))  dev.off()
pdf("DEG/GO_KEGG/DEG_keygene_Bubble.pdf",width = 6,height = 5)
print(myGOBubble(circ, title = '', colour = brewer.pal(4, "Set3"), display = 'single', labels = 10,table.legend = F))  
dev.off()
#圈图
source("prepare/myGOCircle.R")
#myGOCircle <- edit(GOCircle)
environment(myGOCircle) <- environment(GOCircle)
while (!is.null(dev.list()))  dev.off()
pdf("DEG/GO_KEGG/DEG_keygene_Circle.pdf",width = 6,height = 3.4)
print(myGOCircle(circ, table.legend = T,nsub = circ2$ID,label.fontface = "plain",label.size = 2))
dev.off()
#弦图
circ3<-circ[circ$ID%in%circ2$ID,]
top2 = circ3 %>% group_by(ID) %>% top_n(2,abs(logFC))
if(length(unique(top2$genes))<10){top2 = circ3 %>% group_by(ID) %>% top_n(3,abs(logFC))
print('genes=3')}
genes3<-genes1[genes1$ID%in%top2$genes,]
chord <- chord_dat(data = circ3, genes = genes3, process = circ2$ID)
source("prepare/myGOChord.R")
#myGOChord <- edit(GOChord)
environment(myGOChord) <- environment(GOChord)
while (!is.null(dev.list()))  dev.off()
pdf("DEG/GO_KEGG/DEG_keygene_Chord.pdf",width = 6,height = 7)
print(myGOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,
                limit = c(0, 0),ribbon.col=colorRampPalette(brewer.pal(12, "Set3"))(length(circ2$ID))))
dev.off()
