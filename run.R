####运行环境####
setwd('/work/')
source('prepare/R包加载.R')
#颜色主题,去转录组分析模块代码里复制新的颜色主题
colors <- c('#0073C2','#EFC000','#868686','#CD534C','#7AA6DC','#003C67','#8F7700','#3B3B3B','#A73030','#4A6990')
####1.数据处理（TCGA）####
####参数####
#TCGA项目的名称
project="TCGA-OV"
####运行####
source('prepare/TCGA数据处理.R')
getwd()
####2.差异分析####
####参数####

#log2转换之后的表达矩阵和分组信息路径
exp_path='data/TCGA-OV_tpm_log.csv'
group_path='data/TCGA-OV_group.csv'
keygene_path='../单细胞部分/10.intersect/intersect_gene.txt'

#分组信息文件中的分组名,test是实验组，control是对照组
test="tumor"
control="control"

#差异分析阈值
P="p.adj"
#可选"p","p.adj"
P_value=0.05
#可选0.05,0.01
logFC_value=0.5
#可选0,0.5,1,2
####运行####
source('prepare/差异分析.R')

####3.构建预后模型####
####参数1####
#仅肿瘤样本的表达矩阵和临床信息以及感兴趣的基因路径
tumor_exp_path='data/TCGA-OV_tpm_log_tumor.csv'
tumor_cli_path='data/TCGA-OV_clinical.csv'

#随机数种子取值范围
seed_range=1:1000

#ROC、列线图、校准曲线用到的时间
years=c(1,2,3)
#可选c(1,3,5),c(1,2,3)
####运行####
source('prepare/cox_lasso.R')#最主流的预后模型构建方法



####参数2####
#构建预后模型采用的方法
prognosis_method="cox_lasso"


#根据前面的结果挑选合适的随机数种子
prognosis_seed=436

####运行####
source('prepare/prognosis.R')

####4.构建列线图模型####
source('prepare/nomogram.R')
#如果还报错将years从c(1,3,5)改成c(1,2,3)

####5.高低风险组分析汇总####
####参数####
TCGA=TRUE 
####运行####
source('prepare/高低风险组分析汇总.R')
