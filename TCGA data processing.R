####TCGA####
projects<-‘TCGA-OV’
#加载 R 包 TCGAbiolinks
library(TCGAbiolinks)

#选择研究的癌症项目，下载 Counts 数据
query_exp = GDCquery(project = project,
                     legacy = FALSE, 
                     data.category = "Transcriptome Profiling",#选择转录组数据
                     data.type = "Gene Expression Quantification",#选择基因表达量，这个选项和 GDC 官网上看到的一样
                     workflow.type = "STAR - Counts")#数据类型选择 STAR - Counts格式

GDCdownload(query = query_exp)

pre_exp = GDCprepare(query = query_exp)
library(SummarizedExperiment)

# mRNA的counts矩阵
expr_counts_mrna <- assay(pre_exp,"unstranded")

# mRNA的TPM矩阵
expr_tpm_mrna <- assay(pre_exp,"tpm_unstrand")


#添加gene_symbol也就非常简单了，只要提取gene_name这一列，然后和原来的表达矩阵合并即可！
#先提取 gene_name
symbol_mrna <- rowData(pre_exp)$gene_name

#和基因表达矩阵合并就行了：
symbol_mrna_symbol = cbind(data.frame(symbol_mrna),
                           as.data.frame(expr_counts_mrna))

#使用aggregate根据symbol列中的相同基因进行合并
counts = aggregate(.~symbol_mrna, mean, data = symbol_mrna_symbol)

symbol_mrna_symbol = cbind(data.frame(symbol_mrna),
                           as.data.frame(expr_tpm_mrna))

#使用aggregate根据symbol列中的相同基因进行合并
tpms = aggregate(.~symbol_mrna, mean, data = symbol_mrna_symbol)

#保存该基因表达矩阵文件
dir.create('data')
write.csv(counts, file = paste0('data/',project,'_count.csv'),row.names = F)
write.csv(tpms, file = paste0('data/',project,'_tpm.csv'),row.names = F)

#log2转换
tpms_log<-tpms[,-1]
rownames(tpms_log)<-tpms$symbol_mrna
tpms_log<-log2(tpms_log+1)
write.csv(tpms_log,paste0('data/',project,'_tpm_log.csv'))
#分组信息
group<-str_sub(colnames(tpms_log),14,14)
group<-ifelse(group=='0','tumor','control')
group<-data.frame(sample=colnames(tpms_log),group=group)
print(table(group$group))
write.csv(group,paste0('data/',project,'_group.csv'))

####下载临床数据####
library(stringr)
clinical <- GDCquery_clinic(project = project, type = "clinical")
write.csv(clinical,file = paste0('GDCdata/',project,'/clinical.csv'))

####提取几列临床信息####
cli <- read.csv(paste0('GDCdata/',project,'/clinical.csv'), row.names = 1)

TCGA_cli<-read.csv('prepare/TCGA_cli.csv',row.names = 1)

TCGA_cli=TCGA_cli[TCGA_cli$TCGA_cli%in%colnames(cli),]
cli <- cli[,TCGA_cli$TCGA_cli]
colnames(cli) <- TCGA_cli$cli

#去除没有生存状态的样本
cli <- cli[!is.na(cli$vital_status),]
status<-cli$vital_status=='Alive'
cli$fustat<-ifelse(status,0,1)
cli$futime<-ifelse(status,cli$follow_up,cli$days_to_death)
cli$vital_status <- NULL
cli$follow_up <- NULL
cli$days_to_death <- NULL
cli$fustat <- as.numeric(cli$fustat)
cli$futime <- as.numeric(cli$futime)

cli<-na.omit(cli)
cli<-cli[cli$futime>0,]
cli$futime<-cli$futime/365

sample<-group
sample$sample1<-str_sub(sample$sample,1,12)
sample<-sample[!duplicated(sample$sample1),]
rownames(sample)<-sample$sample1
sample1<-intersect(sample$sample1,cli$id)
sample<-sample[sample1,]

cli<-cli[!duplicated(cli$id),]
rownames(cli)<-cli$id
cli<-cli[sample1,]
write.csv(cli,paste0('data/',project,'_clinical.csv'),row.names = F)

tpms_log_tumor<-tpms_log[,sample$sample]
colnames(tpms_log_tumor)<-sample$sample1
write.csv(tpms_log_tumor,paste0('data/',project,'_tpm_log_tumor.csv'))

####SNV突变数据####
library(TCGAbiolinks)
query_exp = GDCquery(project = project,
                     data.category = "Simple Nucleotide Variation",
                     data.type="Masked Somatic Mutation")

GDCdownload(query = query_exp,directory = "GDCdata/SNV")
####CNV突变数据####

print('数据处理完成')
}else{print('请输入正确的TCGA')}
