
dir.create('risk_high_vs_low')
tumor_exp<-read.csv(tumor_exp_path,row.names = 1)
colnames(tumor_exp)<-gsub('[.]','-',colnames(tumor_exp))
total_risk<-read.table('prognosis/totalRisk.txt',sep = '\t',header = T)
total_risk<-total_risk[colnames(tumor_exp),]

####DEG####
dir.create('risk_high_vs_low/DEG')
exp<-tumor_exp
group <- total_risk$risk
group<-factor(group,levels = c('High','Low'))

## 实验设计矩阵
design <- model.matrix(~ 0 + group)
rownames(design) <- colnames(exp)
colnames(design) <- levels(group)
library(limma)
## 线性建模
fit <- lmFit(exp,design)
cont.matrix <- makeContrasts(contrasts = paste0("High",'-',"Low"), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
## 经验贝叶斯调整
fit2 <- eBayes(fit2)
## 筛选差异基因
dif <- topTable(fit2, coef = 1, n = Inf)
dif <- na.omit(dif)
dif$gene_symbol <- row.names(dif)
library(dplyr)
if(P=='p'){
  dif.up <- dif %>%
    filter(logFC > logFC_value & P.Value < P_value)
  dif.down <- dif %>%
    filter(logFC < (-logFC_value) & P.Value < P_value)
  dif <- dif %>% mutate(group = case_when(
    gene_symbol %in% dif.up$gene_symbol ~ 'up',
    gene_symbol %in% dif.down$gene_symbol ~ 'down',
    TRUE ~ "no"
  ))}else{
    dif.up <- dif %>%
      filter(logFC > logFC_value & adj.P.Val < P_value)
    dif.down <- dif %>%
      filter(logFC < (-logFC_value) & adj.P.Val < P_value)
    dif <- dif %>% mutate(group = case_when(
      gene_symbol %in% dif.up$gene_symbol ~ 'up',
      gene_symbol %in% dif.down$gene_symbol ~ 'down',
      TRUE ~ "no"
    ))
  }

DEG <- dif
print(table(DEG$group))
write.csv(DEG,'risk_high_vs_low/DEG/high_low_DEG.csv')


#热图
top10_padj = DEG[DEG$group!='no',] %>% group_by(group) %>% top_n(5,-log10(adj.P.Val))#选择P值最小的差异基因上下调各五个
top10_p = DEG[DEG$group!='no',] %>% group_by(group) %>% top_n(5,-log10(P.Value))
if(P=="p")top10=top10_p else top10=top10_padj
library(pheatmap)
annotate<-as.character(group)
annotate<-as.data.frame(annotate)
annotate$sample<-colnames(exp)
annotate<-annotate[order(annotate$annotate),]
mat<-exp[top10$gene_symbol,annotate$sample]
annotate<-data.frame(sample=annotate$annotate)
rownames(annotate)<-colnames(mat)
annotate_color<-list(sample=c(colors[2],colors[1]))
names(annotate_color$sample)<-levels(group)

mat=t(scale(t(mat)))
mat[mat > 3] = 3 #限定上限，使表达量大于3的等于3
mat[mat < -3] = -3 #限定下限，使表达量小于-3的等于-3

p<-pheatmap(mat,annotation_col  = annotate,scale = 'none',
            show_rownames = T,annotation_colors = annotate_color,
            show_colnames = F,cluster_cols = F,cluster_rows = T,color = colorRampPalette(colors = c(colors[1],"white",colors[2]))(50))

while (!is.null(dev.list()))  dev.off()
pdf('risk_high_vs_low/DEG/high_low_DEG_heatmap.pdf',width = 6,height = 6)
print(p)
dev.off()
#火山图
library(ggplot2)
if(P=="p")label=top10_p else label=top10_padj#选择P值最小的差异基因上下调各五个
DEG$Label=ifelse(DEG$gene_symbol%in%label$gene_symbol,DEG$gene_symbol,NA)#添加一列存放选取出来的差异基因
Significant=factor(DEG$group,levels = c('up','down',"no"))
dp<-ifelse(P=='p',4,5)
library(ggplot2)
library(ggrepel)
p1<-ggplot(DEG, aes(logFC, -log10(DEG[,dp])))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c(colors[2],colors[1],"grey"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-logFC_value,logFC_value), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  str(DEG, max.level = c(-1, 1))+theme_bw()+geom_text_repel(aes(label=Label))#添加注释文本

while (!is.null(dev.list()))  dev.off()
pdf('risk_high_vs_low/DEG/high_low_DEG_volplot.pdf',width = 6,height = 6)
print(p1)
dev.off()

####GO/KEGG####
dir.create('risk_high_vs_low/GO_KEGG')
DEG<-read.csv('risk_high_vs_low/DEG/high_low_DEG.csv',header = T,row.names = 1)
genes_id<-DEG$gene_symbol[DEG$group!='no']
library(clusterProfiler)
library(Seurat)
#ID转换，通常采用ENTREZID进行功能富集分析
genes_id <- bitr(genes_id,fromType = "SYMBOL",toType = "ENTREZID",
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
write.csv(ego_list, "risk_high_vs_low/GO_KEGG/high_low_enrichGO.csv")

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
write.csv(ekegg_list,"risk_high_vs_low/GO_KEGG/high_low_enrichKEGG.csv")}

#GO富集棒棒糖图，展示BP、CC、MF的显著性前三条目
library(stringr)
library(ggpubr)
g<-read.csv("risk_high_vs_low/GO_KEGG/high_low_enrichGO.csv")
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
pdf("risk_high_vs_low/GO_KEGG/high_low_GO.pdf",height = 8,width = 6)
print(g1plot)
dev.off()
#KEGG富集棒棒糖图，展示显著性前9通路
if(nrow(ekegg_list)!=0){e<-read.csv("risk_high_vs_low/GO_KEGG/high_low_enrichKEGG.csv")
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
pdf("risk_high_vs_low/GO_KEGG/high_low_KEGG.pdf",height = 8,width = 6)
print(e1plot)
dev.off()
}
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
write.csv(circ,"risk_high_vs_low/GO_KEGG/high_low_go_kegg_zscore.csv")

x <- circ[,c(2,8)]
x <- x[!duplicated(x),]
circ1 <- merge(go_kegg,x,by = "ID")

#选取子集
circ2 <- circ1 %>% arrange(adj_pval) %>%  group_by(category) %>% do(head(.,n = 3))

#柱状图
#GOBar(subset(circ2, category == 'BP'))
while (!is.null(dev.list()))  dev.off()
pdf("risk_high_vs_low/GO_KEGG/high_low_Bar.pdf",width = 5.5,height = 7)
print(GOBar(circ2,display = 'single',order.by.zscore = T))#display = 'multiple'
dev.off()
#气泡图
source("prepare/myGOBubble.R")
#myGOBubble <- edit(GOBubble)
environment(myGOBubble) <- environment(GOBubble)
#reduced_circ <- reduce_overlap(circ, overlap = 0.75)
while (!is.null(dev.list()))  dev.off()
pdf("risk_high_vs_low/GO_KEGG/high_low_Bubble.pdf",width = 6,height = 5)
print(myGOBubble(circ, title = '', colour = brewer.pal(4, "Set3"), display = 'single', labels = 10,table.legend = F))  
dev.off()
#圈图
source("prepare/myGOCircle.R")
#myGOCircle <- edit(GOCircle)
environment(myGOCircle) <- environment(GOCircle)
while (!is.null(dev.list()))  dev.off()
pdf("risk_high_vs_low/GO_KEGG/high_low_Circle.pdf",width = 6,height = 3.4)
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
pdf("risk_high_vs_low/GO_KEGG/high_low_Chord.pdf",width = 6,height = 7)
print(myGOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,
          limit = c(0, 0),ribbon.col=colorRampPalette(brewer.pal(12, "Set3"))(length(circ2$ID))))
dev.off()
####GSEA####
dir.create('risk_high_vs_low/GSEA')
DEG<-read.csv('risk_high_vs_low/DEG/high_low_DEG.csv',header = T,row.names = 1)
#GSEA 需要样本中全部的基因及其logFC（差异和非差异都需要）
library(GSEABase)
library(enrichplot)
library(cowplot)
library(msigdbr)
ge = DEG$logFC
names(ge) = DEG$gene_symbol
ge = sort(ge,decreasing = T)
head(ge)
msgdC2 = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")

pathway_gene=cbind(data.frame(msgdC2$gs_name),data.frame(msgdC2$gene_symbol))  
colnames(pathway_gene)=c('term','gene')

gsea<-GSEA(ge, TERM2GENE =pathway_gene,
           nPerm=1000,
           verbose=FALSE,by="fgsea",
           pAdjustMethod="BH",
           pvalueCutoff=0.05)
res<-gsea@result
write.csv(res,"risk_high_vs_low/GSEA/high_vs_low_GSEA.csv")
a<-order(gsea$NES,decreasing = T)[c(1:3,(nrow(res)-2):nrow(res))]
titles =str_sub(gsea$Description[a],6)
titles2<-c()
for(i in titles){
  titles_d=str_split(i,'_')
  for(j in 1:length(titles_d[[1]])){
  titles1<-ifelse(j==1,titles_d[[1]][j],ifelse(j==4,paste(titles1,titles_d[[1]][j],sep = '\n'),paste(titles1,titles_d[[1]][j],sep = ' ')))
  }
  titles2<-c(titles2,titles1)
}

titles =titles2
source('prepare/myGSEA.R')
environment(myGSEA) <- environment(gseaplot2)
gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]] <- myGSEA(gsea, geneSetID = a[i], title =titles[i],pvalue_table = F,base_size = 6,subplots = 1:3,color = colors[1])
}

gsea_gseplot<-plot_grid(gseaplot[[1]],gseaplot[[2]],gseaplot[[3]],
              gseaplot[[4]],gseaplot[[5]],gseaplot[[6]],
              nrow = 2,labels=c('A','B','C','D','E','F'),
              rel_widths=c(1,1),label_size =15,label_fontface = 'plain')+
  theme(plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm"))

while (!is.null(dev.list()))  dev.off()
pdf('risk_high_vs_low/GSEA/high_vs_low_GSEA.pdf',width = 17/2.54,height = 11.4/2.54)
print(gsea_gseplot)
dev.off()


####GSVA####
dir.create('risk_high_vs_low/GSVA')
exp2<-as.matrix(tumor_exp)
#需要表达矩阵
library(msigdbr)
library(GSVA)
##参数设置
min.sz = 2 # 单个基因集基因数下限
max.sz = 10000 # 基因数上限
parallel.sz= 10 # 并行处理cpu数量
mx.diff= T # ES 计算为最大正和负随机游走偏差之间的幅度差。
tau=1 # 默认1，tau=1 when method="gsva" and tau=0.25 when method="ssgsea" 
method='gsva' # 默认
kcdf = "Gaussian" # 默认"Gaussian"，当输入数值为整数时，设置为"Poisson"
msgdC2 = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
gene_sets = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)

gsva_scores = gsva(expr=exp2, gset.idx.list=gene_sets, 
                   method=method, kcdf=kcdf, min.sz=min.sz, max.sz=max.sz, 
                   parallel.sz=parallel.sz, mx.diff=mx.diff )
## 导出数据
gsva_scores = as.data.frame(gsva_scores)
gsva_scores$geneset = rownames(gsva_scores)
gsva_scores = gsva_scores %>% dplyr::select(geneset,everything())
write.csv(gsva_scores,file = "risk_high_vs_low/GSVA/high_vs_low_GSVA.csv", row.names = F)
#差异通路热图
library(limma)

exp<-read.csv('risk_high_vs_low/GSVA/high_vs_low_GSVA.csv',row.names = 1)
colnames(exp)<-gsub('[.]','-',colnames(exp))
group <- total_risk$risk
group<-factor(group,levels = c('High','Low'))

## 实验设计矩阵
design <- model.matrix(~ 0 + group)
rownames(design) <- colnames(exp)
colnames(design) <- levels(group)
library(limma)
## 线性建模
fit <- lmFit(exp,design)
cont.matrix <- makeContrasts(contrasts = paste0("High",'-',"Low"), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
## 经验贝叶斯调整
fit2 <- eBayes(fit2)
## 筛选差异通路
dif <- topTable(fit2, coef = 1, n = Inf)
dif <- na.omit(dif)
dif$pathway <- row.names(dif)
library(dplyr)
if(P=='p'){
  dif.up <- dif %>%
    filter(logFC > 0 & P.Value < P_value)
  dif.down <- dif %>%
    filter(logFC < (-0) & P.Value < P_value)
  dif <- dif %>% mutate(group = case_when(
    pathway %in% dif.up$pathway ~ 'up',
    pathway %in% dif.down$pathway ~ 'down',
    TRUE ~ "no"
  ))}else{
    dif.up <- dif %>%
      filter(logFC > 0 & adj.P.Val < P_value)
    dif.down <- dif %>%
      filter(logFC < (-0) & adj.P.Val < P_value)
    dif <- dif %>% mutate(group = case_when(
      pathway %in% dif.up$pathway ~ 'up',
      pathway %in% dif.down$pathway ~ 'down',
      TRUE ~ "no"
    ))
  }

DEG <- dif
print(table(DEG$group))
write.csv(DEG,'risk_high_vs_low/GSVA/high_low_DEG_pathway.csv')


#热图
top10_padj = DEG[DEG$group!='no',] %>% group_by(group) %>% top_n(5,-log10(adj.P.Val))#选择P值最小的差异通路上下调各五个
top10_p = DEG[DEG$group!='no',] %>% group_by(group) %>% top_n(5,-log10(P.Value))
if(P=="p")top10=top10_p else top10=top10_padj
library(pheatmap)
annotate<-as.character(group)
annotate<-as.data.frame(annotate)
annotate$sample<-colnames(exp)
annotate<-annotate[order(annotate$annotate),]
mat<-exp[top10$pathway,annotate$sample]
annotate<-data.frame(sample=annotate$annotate)
rownames(annotate)<-colnames(mat)
annotate_color<-list(sample=c(colors[2],colors[1]))
names(annotate_color$sample)<-levels(group)

mat=t(scale(t(mat)))
mat[mat > 3] = 3 #限定上限，使表达量大于3的等于3
mat[mat < -3] = -3 #限定下限，使表达量小于-3的等于-3


p<-pheatmap(mat,annotation_col  = annotate,scale = 'none',
            show_rownames = T,annotation_colors = annotate_color,
            show_colnames = F,cluster_cols = F,cluster_rows = T,color = colorRampPalette(colors = c(colors[1],"white",colors[2]))(50))

while (!is.null(dev.list()))  dev.off()
pdf('risk_high_vs_low/GSVA/high_low_DEG_pathway_heatmap.pdf',width = 17,height = 9.6)
print(p)
dev.off()

####ssGSEA####
dir.create('risk_high_vs_low/ssGSEA')
exp<-as.matrix(tumor_exp)#行为基因列为样本的矩阵
load('prepare/celltype28.rda')
##参数设置
min.sz = 2 # 单个基因集基因数下限
max.sz = 10000 # 基因数上限
parallel.sz= 10 # 并行处理cpu数量
mx.diff= T # ES 计算为最大正和负随机游走偏差之间的幅度差。
tau=0.25 # 默认1，tau=1 when method="gsva" and tau=0.25 when method="ssgsea" 
method='ssgsea' # 默认
kcdf = "Gaussian" # 默认"Gaussian"，当输入数值为整数时，设置为"Poisson"
#ssGSEA生成普通转录组中单细胞的细胞类型的细胞丰度
ssgsea = gsva(expr=exp, gset.idx.list=celltype28, 
              method=method, kcdf=kcdf, min.sz=min.sz, max.sz=max.sz, 
              parallel.sz=parallel.sz, mx.diff=mx.diff )

#绘制热图查看细胞丰度的大致情况
pheatmap::pheatmap(ssgsea,scale = 'none',show_colnames = F)

#min-max标准化
ssgsea.1 <- ssgsea
ssgsea.2<-ssgsea.1
for (i in colnames(ssgsea.1)) { 
  ssgsea.2[,i] <- (ssgsea.1[,i] -min(ssgsea.1[,i]))/(max(ssgsea.1[,i] )-min(ssgsea.1[,i] ))
}

write.csv(ssgsea.2,'risk_high_vs_low/ssGSEA/ssGSEA.csv')

#读取
TME.results <- read.csv("risk_high_vs_low/ssGSEA/ssGSEA.csv",row.names = 1)
colnames(TME.results)<-gsub('[.]','-',colnames(TME.results))
re <- as.data.frame(t(TME.results))

#堆叠图
re1<-apply(re,1,function(x){y=x/sum(x);return(y)})
re1<-as.data.frame(t(re1))
library(RColorBrewer)
library(tibble)
dat <- re1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

dat$risk<-ifelse(dat$Sample%in%rownames(total_risk)[total_risk$risk=='High'],'High','Low')

identitplot<-ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "", y = "Estiamted Proportion") + 
  #facet_grid(. ~ risk,scales = 'free_x',space = 'free_x')+#按照高低风险组分面，risk是分组列的列名
  scale_fill_manual(values = colorRampPalette(colors)(28))+
  theme(plot.title = element_blank(),plot.subtitle = element_blank(),
        plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
        panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.title = element_text(size = 7,colour = "black"),
        axis.text.y = element_text(size = 6,colour = "black"),
        legend.background = element_blank(),legend.key = element_blank(),
        legend.key.width = unit(0.15,"cm"),legend.key.height = unit(0.15,"cm"),
        legend.title = element_blank(),legend.text = element_text(size = 6),
        legend.position = "bottom",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))#+coord_flip() #旋转箱线图
while (!is.null(dev.list()))  dev.off()
pdf("risk_high_vs_low/ssGSEA/stacked_barplot.pdf",width = 17/2.54,height = 7/2.54)
print(identitplot)
dev.off()
#热图
library(pheatmap)
k <- apply(re,2,function(x) {sum(x == 0) < nrow(TME.results)/2})#在一半以上样本里丰度为0的免疫细胞，就不展示在热图里了
table(k)

re2 <- as.data.frame(t(re[,k]))

Group = total_risk$risk
annotation = data.frame(group = Group,
                        row.names = colnames(re2))
annotation_color<-list(c(colors[2],colors[1]))
names(annotation_color[[1]])<-unique(annotation$group)[order(unique(annotation$group))]
names(annotation_color)<-'group'

mat=t(scale(t(re2)))
mat[mat > 3] = 3 #限定上限，使表达量大于3的等于3
mat[mat < -3] = -3 #限定下限，使表达量小于-3的等于-3

pheatmap(mat[,order(annotation$group)],
         scale = "none",
         #cellwidth = 5, cellheight = 25,## 设置热图块大小
         show_colnames = F, 
         cluster_cols = F,cluster_rows = T,
         annotation_col = annotation,
         annotation_colors = annotation_color,
         fontsize = 7,
         colorRampPalette(colors = c(colors[1],"white",colors[2]))(50))


#箱线图
library(ggpubr)

#读取免疫浸润分析结果，并转化形式
#re为全部免疫细胞在样本中的丰度，也可用提取在超一半样本中丰度不为零的免疫细胞，即t(re2)
dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

dat$group <- ifelse(dat$Sample%in%rownames(total_risk)[total_risk$risk=='High'],'High','Low')

#绘图
boxplot1<-ggplot(dat,aes(Cell_type,Proportion,fill = group)) + 
  geom_boxplot(outlier.shape = 21,outlier.color = NA,outlier.size = 0.8,
               #notch=TRUE,notchwidth=0.8, #在箱子上生成槽口,notchwidth越小则越往里凹
               alpha = 0.6) + 
  # stat_summary(fun.y="mean",geom="point",shape=23,size=1,fill="white",
  #              position = position_dodge(0.8))+ #添加均值标记
  labs(x = "", y = "Estimated Proportion") +
  scale_fill_manual(values = c(colors[2],colors[1]))+
  stat_compare_means(aes(group = group,label = ..p.signif..),bracket.size = 0.6, size = 3,
                     label.y = max(dat$Proportion)-0.02, #p值位置
                     hide.ns = T,#隐藏ns
                     method = ifelse(length(unique(dat$group))==2,"wilcox.test","kruskal.test"))+
  theme(plot.title = element_blank(),plot.subtitle = element_blank(),
        plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
        panel.border = element_rect(size = 0.5,fill = NA,colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size = 7,colour = "black"),
        axis.title = element_text(size = 7,colour = "black"),
        axis.text.y = element_text(size = 6,colour = "black"),
        legend.background = element_blank(),legend.key = element_blank(),
        legend.title = element_text(size = 7),legend.text = element_text(size = 7),
        legend.position = "top",legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),
        legend.box.spacing = unit(3,"pt"))#+coord_flip() #旋转箱线图
while (!is.null(dev.list()))  dev.off()
pdf("risk_high_vs_low/ssGSEA/immune_cells_boxplot.pdf",width = 17/2.54,height = 10/2.54)  
print(boxplot1)
dev.off()
#
re_high<-re[rownames(total_risk)[total_risk$risk=='High'],]
re_low<-re[rownames(total_risk)[total_risk$risk=='Low'],]
for(i in colnames(re)){
  res<-wilcox.test(re_high[,i],re_low[,i])
  if(i == colnames(re[1]))
    wilcoxtest<-data.frame(celltype=i,p.value=res$p.value) else{
      tmp<-data.frame(celltype=i,p.value=res$p.value)
      wilcoxtest<-rbind(wilcoxtest,tmp)
    }
}
write.csv(wilcoxtest,'risk_high_vs_low/ssGSEA/wilcoxtest.csv')
cells<-wilcoxtest$celltype[wilcoxtest$p.value<0.05]
cells<-na.omit(cells)

#带显著性的相关性热图
library(tidyverse)
library(reshape2)
library(Hmisc)

exp<-t(re2)
exp=cbind(exp,total_risk$riskScore)
colnames(exp)[ncol(exp)]<-"risk"
table(rownames(exp)==rownames(total_risk))

exp_high<-exp[total_risk$risk=='High',]
exp_low<-exp[total_risk$risk=='Low',]
###high###
rcorr_list<-rcorr(exp_high,type = "pearson")#"pearson","spearman"
cor_ma<-rcorr_list$r
cor_ma1<-cor_ma[-nrow(cor_ma),-ncol(cor_ma)]
cor_risk<-cor_ma[,'risk']

p_ma<-rcorr_list$P
p_ma1<-p_ma[-nrow(p_ma),-ncol(p_ma)]
p_risk<-p_ma[,'risk']
  
data_cor<-melt(get_upper_tri(cor_ma1))
data_p<-melt(get_upper_tri(p_ma1))

data_all<-merge(data_cor,data_p,by=c('Var1','Var2'))
colnames(data_all)<-c("x","y","cor","p")

data_all$p_text<-ifelse(data_all$p<0.0001,"****",
                        ifelse(data_all$p<0.001,"***",
                               ifelse(data_all$p<0.01,"**",
                                      ifelse(data_all$p<0.05,"*",""))))

data_dot<-data.frame(x=names(cor_risk),y=names(cor_risk),cor=cor_risk,p=p_risk)
data_dot<-rbind(data_dot,data_dot)
data_dot<-data_dot[data_dot$x!='risk',]
data_dot[(length(unique(data_all$x))+1):(length(unique(data_all$x))*2),'x']<-levels(data_all$x)[length(unique(data_all$x))]
data_dot[(length(unique(data_all$x))+1):(length(unique(data_all$x))*2),'y']<-levels(data_all$x)[round(length(unique(data_all$x))/2)]

data_dot$group<-rep(1:length(unique(data_all$x)),2)

a<-c()
for(i in 1:nrow(data_all)){
  if(!data_all$x[i]==data_all$y[i]){
    if(sum(is.na(data_all[i,]))>0)a=c(a,i)
  }
}
data_all<-data_all[-a,]
p<-ggplot(data_all, aes(x, y)) + 
  geom_tile(aes(fill = cor), colour = "white", size = 0.1)+
  geom_text(aes(label=p_text),col ="black",size = 2) +
  geom_point(data = data_dot,aes(x=x,y=y),size=3,color="#B59283")+
  geom_line(data = data_dot,aes(x=x,y=y,group=group,color=cor),linetype="dashed",lwd=1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + 
  scale_color_gradient2(low = "#5C5DAF",mid = "#FCB886",high = "#EA2E2D") + 
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_blank(), # 调整x轴文字
        axis.text.y = element_text(size = 6, face = "plain"),#调整y轴文字
        legend.text = element_text(size = 7, face = "plain"),
        legend.title = element_text(size = 7, face = "plain"),
        panel.grid.major =  element_blank(),
        legend.position = 'left') + 
  labs(fill =paste0("* p < 0.05","\n\n","** p <0.01","\n\n","*** p <0.001","\n\n","**** p <0.0001","\n\n","Correlation"))# 修改 legend 内容
while (!is.null(dev.list()))  dev.off()
pdf('risk_high_vs_low/ssGSEA/immune_cells_cor_heatmap_high.pdf',width =8 ,height = 6)
print(p)
dev.off()

###low###
rcorr_list<-rcorr(exp_low,type = "pearson")#"pearson","spearman"
cor_ma<-rcorr_list$r
cor_ma1<-cor_ma[-nrow(cor_ma),-ncol(cor_ma)]
cor_risk<-cor_ma[,'risk']

p_ma<-rcorr_list$P
p_ma1<-p_ma[-nrow(p_ma),-ncol(p_ma)]
p_risk<-p_ma[,'risk']

data_cor<-melt(get_upper_tri(cor_ma1))
data_p<-melt(get_upper_tri(p_ma1))
data_all<-merge(data_cor,data_p,by=c('Var1','Var2'))
colnames(data_all)<-c("y","x","cor","p")

data_all$p_text<-ifelse(data_all$p<0.0001,"****",
                        ifelse(data_all$p<0.001,"***",
                               ifelse(data_all$p<0.01,"**",
                                      ifelse(data_all$p<0.05,"*",""))))

data_dot<-data.frame(x=names(cor_risk),y=names(cor_risk),cor=cor_risk,p=p_risk)
data_dot<-rbind(data_dot,data_dot)
data_dot<-data_dot[data_dot$x!='risk',]
data_dot[(length(unique(data_all$x))+1):(length(unique(data_all$x))*2),'x']<-levels(data_all$x)[1]
data_dot[(length(unique(data_all$x))+1):(length(unique(data_all$x))*2),'y']<-levels(data_all$x)[round(length(unique(data_all$x))/2)]

data_dot$group<-rep(1:length(unique(data_all$x)),2)

a<-c()
for(i in 1:nrow(data_all)){
  if(!data_all$x[i]==data_all$y[i]){
    if(sum(is.na(data_all[i,]))>0)a=c(a,i)
  }
}
data_all<-data_all[-a,]
p<-ggplot(data_all, aes(x, y)) + 
  geom_tile(aes(fill = cor), colour = "white", size = 0.1)+
  geom_text(aes(label=p_text),col ="black",size = 2) +
  geom_point(data = data_dot,aes(x=x,y=y),size=3,color="#B59283")+
  geom_line(data = data_dot,aes(x=x,y=y,group=group,color=cor),linetype="dashed",lwd=1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + 
  scale_color_gradient2(low = "#5C5DAF",mid = "#FCB886",high = "#EA2E2D") + 
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_blank(), # 调整x轴文字
        axis.text.y = element_text(size = 6, face = "plain"),#调整y轴文字
        legend.text = element_text(size = 7, face = "plain"),
        legend.title = element_text(size = 7, face = "plain"),
        panel.grid.major =  element_blank(),
        legend.position = 'right') + 
  scale_y_discrete(position = 'right')+
  labs(fill =paste0("* p < 0.05","\n\n","** p <0.01","\n\n","*** p <0.001","\n\n","**** p <0.0001","\n\n","Correlation"))# 修改 legend 内容
while (!is.null(dev.list()))  dev.off()
pdf('risk_high_vs_low/ssGSEA/immune_cells_cor_heatmap_low.pdf',width =8 ,height = 6)
print(p)
dev.off()

#相关性散点图
library(dplyr)
library(ggplot2)
library(ggExtra)
genes <- read.csv('prognosis/prognosis_coefficients.csv')[,1]
exp <- t(tumor_exp[genes,])
wilcoxtest<-read.csv('risk_high_vs_low/ssGSEA/wilcoxtest.csv',row.names = 1)
cells<-wilcoxtest$celltype[wilcoxtest$p.value<0.05]
cells<-na.omit(cells)

dir.create("risk_high_vs_low/ssGSEA/cor_dot")
plot_sort<-data.frame()
for (c in 1:length(cells)) {
  for (g in 1:length(genes)) {
    data <- data.frame(gene = as.numeric(exp[,genes[g]]),#这里如果用的是已经log转换过的表达矩阵就不用再log转换了
                       cell = re[,cells[c]])
    data <- data[which(data$cell != 0),]
    cor <- cor.test(data$gene, data$cell, method = "pearson")
    p1<-ggplot(data,aes(x = cell, y = gene))+ 
      xlim(min(data$cell),max(data$cell))+ylim(min(data$gene),max(data$gene))+
      labs(x=cells[c], y=paste("Expression level of",genes[g])) +
      geom_point(size=0.5)+ #shape=21空心圆圈
      geom_smooth(method=lm,color=ifelse(cor$estimate < 0,colors[1],colors[2]))+
      theme(plot.title = element_blank(),plot.subtitle = element_blank(),
            plot.background = element_blank(),plot.margin = margin(t=2,r=2,b=2,l=2,unit="pt"),
            panel.border = element_rect(size = 0.45,fill = NA,colour = "black"),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.text.x = element_blank(),axis.ticks.x = element_blank(),
            axis.title = element_text(size = 6,colour = "black"),
            axis.text.y = element_text(size = 5,colour = "black"))+
      annotate("text", 0.85*(max(data$cell)-min(data$cell))+min(data$cell), y = 0.9*(max(data$gene)-min(data$gene))+min(data$gene),
               label=paste0('R=',round(cor$estimate,4),'\n',ifelse(cor$p.value < 0.001, "p < 0.001",paste0('p=',round(cor$p.value,4)))),
               colour="black",size=2.5)
    p2 <- ggMarginal(p1, type="histogram",colour=ifelse(cor$estimate < 0,colors[1],colors[2]),fill = ifelse(cor$estimate < 0,colors[1],colors[2]))
    if (cor$p.value < 0.05) {
      while (!is.null(dev.list()))  dev.off()
      pdf(paste0("risk_high_vs_low/ssGSEA/cor_dot/",genes[g],"_",gsub("[ ]",".",cells[c]),"_","cor_dot.pdf"),width = 5/2.54,height = 5/2.54)
      print(p2)
      dev.off()
      plot_sort<-rbind(plot_sort,data.frame(plot=paste0("risk_high_vs_low/ssGSEA/cor_dot/",genes[g],"_",gsub("[ ]",".",cells[c]),"_","cor_dot.pdf"),
                                            R=cor$estimate))
    }
  }
}

dir.create('risk_high_vs_low/ssGSEA/cor_dot_use')
plot_sort<-plot_sort[order(abs(plot_sort$R),decreasing = T),]
sapply(plot_sort$plot[1:9],function(x){file.copy(from = x,to ='risk_high_vs_low/ssGSEA/cor_dot_use')})

####TMB####
if(TCGA){
  library(maftools)
  library(ggplot2)
  
  dir.create('risk_high_vs_low/TMB')
  
  filepath = dir(path =paste0(maf_path,'/',project,'/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/'),
                 pattern = "maf.gz$",
                 full.names = T,
                 recursive = T)
  
  x = lapply(filepath, data.table::fread, skip = "Hugo_Symbol")
  x = data.table::rbindlist(l = x, use.names = TRUE, fill = TRUE)
  luad = maftools::read.maf(maf = x)
  
  luad # 直接输入MAF对象可以查看MAF文件的基本信息
  getSampleSummary(luad) # 显示样品的统计
  getGeneSummary(luad) # 显示基因的统计
  getClinicalData(luad) # 显示样品关联的临床数据
  getFields(luad) # 显示MAF文件中的所有字段
  write.mafSummary(maf=luad, basename=paste0(maf_path,"/luad")) # 将详细统计结果输出到文件
  while (!is.null(dev.list()))  dev.off()
  pdf('risk_high_vs_low/TMB/mafSummary.pdf',width = 10,height = 10)
  plotmafSummary(maf=luad, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)
  dev.off()
  
  maf<-read.maf(paste0(maf_path,'/luad_maftools.maf'))
  dat<-maf@data
  sample1=unique(as.character(dat$Tumor_Sample_Barcode))
  sample2=substr(sample1,1,12)
  sample1=sample1[!duplicated(sample2)]
  sample2=sample2[!duplicated(sample2)]
  hname=rownames(total_risk)[total_risk$risk=='High']
  lname=rownames(total_risk)[total_risk$risk=='Low']
  
  samp=intersect(sample2,hname)
  name1=sample1[which(sample2 %in% samp)]
  samp=intersect(sample2,lname)
  name2=sample1[which(sample2 %in% samp)]
  
  high_maf=subsetMaf(maf,tsb = name1)
  low_maf=subsetMaf(maf,tsb = name2)
  maf1=subsetMaf(maf,tsb = c(name1,name2))
  
  save(maf1,file = "risk_high_vs_low/TMB/maf1.rda")
  save(high_maf,low_maf,file="risk_high_vs_low/TMB/maf_hl.rda")
  
  a<-maf1@data
  b<-unique(a[,c('Hugo_Symbol','Tumor_Sample_Barcode')])
  times<-as.data.frame(table(b$Hugo_Symbol))
  times<-times[order(times$Freq,decreasing = T),]
  top20_gene<-as.character(times$Var1[1:20])
  while (!is.null(dev.list()))  dev.off()
  pdf('risk_high_vs_low/TMB/maf.pdf',width = 7,height = 7)
  oncoplot(maf=maf1, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
  dev.off()
  while (!is.null(dev.list()))  dev.off()
  pdf('risk_high_vs_low/TMB/maf_High_risk.pdf',width = 7,height = 7)
  oncoplot(maf=high_maf, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
  dev.off()
  while (!is.null(dev.list()))  dev.off()
  pdf('risk_high_vs_low/TMB/maf_Low_risk.pdf',width = 7,height = 7)
  oncoplot(maf=low_maf, borderCol=NULL,gene_mar = 8,genes = top20_gene,keepGeneOrder = T)
  dev.off()
  
  #肿瘤突变负荷计算
  #maftools版本2.10.05
  tmb1<-tmb(maf = maf1)
  
  tmb1$risk<-ifelse(tmb1$Tumor_Sample_Barcode%in%name1,'High','Low')
  write.csv(tmb1,'risk_high_vs_low/TMB/tmb.csv')
  
  tmbplot<-ggplot(data = tmb1, aes(x = risk, y = total_perMB_log, fill = risk))+ 
    #scale_fill_manual(values = c("class1"="#02786A", "class2"="#EB4B17")) +
    geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                size = 0.8, color="black") +
    geom_boxplot(notch = TRUE, outlier.size = -1, 
                 color="black", lwd=0.8, alpha = 0.7) +
    geom_point(shape = 21, size=2,
               position = position_jitterdodge(), # 让点散开
               color="black", alpha = 1) +
    theme_classic() + 
    ylab("log10 (Tumor mutation burden)") +
    xlab("") +
    theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "none",
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12)) +
    stat_compare_means(aes(group = risk,label = ..p.signif..),bracket.size = 0.6, size = 3,
                       label.x = 1.5,
                       label.y = max(tmb1$total_perMB_log)*0.95+min(tmb1$total_perMB_log)*0.05, #p值位置
                       hide.ns = T,#隐藏ns
                       method = ifelse(length(unique(tmb1$risk))==2,"wilcox.test","kruskal.test")) 
  while (!is.null(dev.list()))  dev.off()
  pdf("risk_high_vs_low/TMB/TMB-class_box.pdf",width = 6,height = 6)
  print(tmbplot)
  dev.off()
}

####药物敏感性####

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
testPtype <- fread(paste0(drug_path,'/DrugPredictions.csv'), data.table = F) %>% column_to_rownames("V1")

data <- merge(total_risk,testPtype,by = "row.names") %>% column_to_rownames("Row.names")

##绘图
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
sample_size = data %>% group_by(risk) %>% dplyr::summarise(num=n())

class1 <- data[which(data$risk == "High"),]
class2 <- data[which(data$risk == "Low"),]

dir.create("risk_high_vs_low/drug")

DrugPredictions<-data.frame()

for (i in colnames(testPtype)) {
  test <- wilcox.test(class1[,i],class2[,i])
  if (test$p.value < 0.05) {
    P <- data %>%
      left_join(sample_size) %>%
      mutate(myaxis = paste0(risk, "\n", "n=", num)) %>%
      ggplot(aes(x=myaxis, y=data[,i], fill=risk)) +
      geom_violin(width=1) +
      geom_boxplot(width=0.1, color="grey", alpha=0.2) +
      scale_fill_manual(values = c(colors[2],colors[1])) +
      theme(
        legend.position="none",
        plot.title = element_text(size=8),
        panel.background = element_blank(),
        axis.line = element_line(size = 0.5),axis.text = element_text(size=7),
        axis.title = element_text(size=8)
      ) +
      labs(x="",y=i)+
      stat_compare_means(aes(group = risk),method = "wilcox.test",size = 3)
    while (!is.null(dev.list()))  dev.off()
    pdf(paste0("risk_high_vs_low/drug/",i,"_vln.pdf"),width = 3,height = 3)
    print(P)
    dev.off()
    
    DrugPredictions<-rbind(DrugPredictions,data.frame(path=paste0("risk_high_vs_low/drug/",i,"_vln.pdf"),
                                                      drug=i,
                                                      P=test$p.value,
                                                      mp=(median(class1[,i])-median(class2[,i]))/(max(data[,i])-min(data[,i])),
                                                      median=median(data[,i])))
  }
}

DrugPredictions<-DrugPredictions[order(DrugPredictions$median,decreasing = F),]
DrugPredictions1<-DrugPredictions[abs(DrugPredictions$mp)>0.025&DrugPredictions$median<7,]
if(nrow(DrugPredictions1)<1)DrugPredictions1<-DrugPredictions[abs(DrugPredictions$mp)>0.01&DrugPredictions$median<7,]
DrugPredictions1<-DrugPredictions1[1:ifelse(nrow(DrugPredictions1)>9,9,nrow(DrugPredictions1)),]

dir.create('risk_high_vs_low/drug/vlnplot_use')
sapply(DrugPredictions1$path,function(x){file.copy(from = x,to ='risk_high_vs_low/drug/vlnplot_use')})


