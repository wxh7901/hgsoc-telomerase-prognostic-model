
####DEG####
dir.create('DEG')
exp<-read.csv(exp_path,row.names = 1)
colnames(exp)<-gsub('[.]','-',colnames(exp))
group <- read.csv(group_path,row.names = 1)
group<-factor(group$group,levels = c(test,control))

## 实验设计矩阵
design <- model.matrix(~ 0 + group)
rownames(design) <- colnames(exp)
colnames(design) <- levels(group)
library(limma)
## 线性建模
fit <- lmFit(exp,design)
cont.matrix <- makeContrasts(contrasts = paste0(test,'-',control), levels = design)
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
write.csv(DEG,paste0('DEG/',test,'-',control,'_DEG.csv'))


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
pdf(paste0('DEG/',test,'-',control,'_DEG_heatmap.pdf'),width = 6,height = 6)
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
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(plot.title = element_blank(),plot.subtitle = element_blank(),
        plot.background = element_blank(),plot.margin = margin(t=1,r=1,b=1,l=1,unit="pt"),
        panel.border = element_rect(size = 0.6,fill = NA,colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6,colour = "black"),axis.title=element_text(size=7),
        legend.background = element_blank(),legend.key = element_blank(),
        legend.title = element_text(size = 7),legend.text = element_text(size = 7),
        legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),legend.box.spacing = unit(3,"pt"))+
  str(DEG, max.level = c(-1, 1))+geom_text_repel(aes(label=Label),size=2.5)#添加注释文本

while (!is.null(dev.list()))  dev.off()
pdf(paste0('DEG/',test,'-',control,'_DEG_volplot.pdf'),width = 8/2.54,height = 8/2.54)
print(p1)
dev.off()

