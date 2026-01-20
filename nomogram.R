dir.create('nomogram')
library(rms)
library(survival)
library(survminer)
library(glmnet)
library(GEOquery)
library(nomogramEx)
library(ggplotify)
library(ggplot2)

total_risk<-read.table('prognosis/totalRisk.txt',header = T)
clinical<-read.csv(tumor_cli_path,row.names = 1)
TCGA_cli<-read.csv('prepare/TCGA_cli.csv',row.names = 1)

cli<-intersect(colnames(clinical),TCGA_cli$cli)
data<-clinical
for(i in cli){
  if(i=='age'){data[,i]=ifelse(data[,i]>65,'>65','<=65')}
  if(i=='gender'){data[,i]=ifelse(data[,i]%in%c('F','female','Female'),'Female','Male')}
  if(i=='stage'){data[,i]=ifelse(grepl('X',data[,i]),NA,
                                 ifelse(grepl('IV|III',data[,i]),'StageIII&IV','StageI&II'))}
  if(i=='Tstage'){data[,i]=ifelse(grepl('X',data[,i]),NA,
                                  ifelse(grepl('3|4',data[,i]),'T3&4','T1&2'))}
  if(i=='Nstage'){data[,i]=ifelse(grepl('X',data[,i]),NA,
                                  ifelse(grepl('0',data[,i]),'N0','N1&2&3'))}
  if(i=='Mstage'){data[,i]=ifelse(grepl('X',data[,i]),NA,
                                  ifelse(grepl('0',data[,i]),'M0','M1'))}
}

data<-na.omit(data)

data$risk<-total_risk[rownames(data),'risk']
data$risk<-factor(data$risk,levels = c('Low','High'))
cli<-c(cli,'risk')

y<- Surv(time = data$futime,event = data$fustat==1)#1为感兴趣事件
#批量单因素回归模型建立：Uni_cox_model
Uni_cox_model<-
  function(x){
    FML <- as.formula(paste0 ("y~",x))
    cox<- coxph(FML,data=data)
    cox1<-summary(cox)
    HR <- round(cox1$coefficients[,2],2)#提取HR值，保留2位小数
    PValue <- ifelse(cox1$coefficients[,5] < 0.05,"<0.05",round(cox1$coefficients[,5],3))#提取p值，保留3位小数
    CI5 <-round(cox1$conf.int[,3],2)#提取CI，保留2位小数
    CI95 <-round(cox1$conf.int[,4],2)
    #将提取到的信息放入表格中（Uni_cox_model）
    Uni_cox_model<- data.frame(
      names <-rownames(cox1$conf.int),#第1列为亚变量名
      'HR' = HR,#第2列为HR值
      'CI5' = CI5,#第3列为95%ci下区间
      'CI95' = CI95,#第4列为95%ci上区间
      'P' = PValue)#第5列为P值
    return(Uni_cox_model)#返回，开始，进行循环
  }  
#查看原始数据变量的名字
names(data)

variable.names<- cli

Uni_cox <- lapply(variable.names, Uni_cox_model)
library(plyr)
Uni_cox <- ldply(Uni_cox,data.frame)
#将95%CI连接起来
Uni_cox$HR.CI95 <- paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")");Uni_cox
#第一列列名为'Characteristics'
colnames(Uni_cox)[1] <- 'Characteristics'
result <- Uni_cox[,c(1:4,6,5)]

####表格整理
#删除部分变量名，只保留亚变量
result$Characteristics<-str_remove(result$Characteristics,"age|gender|stage|Tstage|Nstage|Mstage|risk")
result

#给参考变量插入空行
ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
#插入空行，形成一个新表
result1<-data.frame("Characteristics", NA, NA, NA, "HR(95%CI)","p")
colnames(result1)<-colnames(result)

for(i in variable.names){
  a<-as.data.frame(table(data[,i]))
  result1<-rbind(result1,
                 ins(i),
                 ins(as.character(a$Var1[!a$Var1%in%result$Characteristics])),
                 result[which(variable.names==i),])
}

result1<-rbind(result1,
               c(NA, NA, NA, NA, NA,NA))

rownames(result1)<-1:nrow(result1)

myVars <- variable.names
catVars <-  variable.names
library(tableone)
table1<- print(CreateTableOne(vars=myVars,
                              data = data,
                              factorVars = catVars),
               showAllLevels=TRUE)

N<-data.frame(c(NA,NA),
              c(NA,NA))
colnames(N)<-colnames(table1)
for(i in 1:length(variable.names)){
  N=rbind(N,
          table1[(i*2):(i*2+1),],
          c(NA,NA))
}

N<-N[,-1]
N<-data.frame(N)
result2<-cbind(result1,N)
result2<-result2[,c(1,7,2:6)]

result2[1,]<-c("Characteristics","Number (%)",NA,NA,NA,"HR (95%CI)","P.value")

####图形进行美化
library(forestplot)
hrzl_lines<-list(gpar(lty=1,lwd=2),#表头上方添加实线
                 gpar(lty=2),#表头下方添加虚线
                 gpar(lwd=2,lty=1,columns=c(1:4)))#最后一行下方添加实线

names(hrzl_lines)<-c('1','2',nrow(result2)+1)

is_summary<-is.na(result2[,2]) #按顺序指定每行是否加粗，T为加粗，F不加粗
is_summary[1]=TRUE

fig2<- forestplot(result2[,c(1,2,6,7)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                  mean=as.numeric(result2[,3]),   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                  lower=as.numeric(result2[,4]),  #告诉函数表格第3列为5%CI，
                  upper=as.numeric(result2[,5]),  #表格第5列为95%CI，它俩要化作线段，穿过方块
                  zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                  boxsize=0.5,       #设置小黑块的大小
                  graph.pos= "right" ,
                  hrzl_lines=hrzl_lines, #最后一行下方添加实线
                  graphwidth = unit(.25,"npc"),
                  #xticks=c(-1,1,2,3,4,5) , #森林图刻度
                  is.summary=is_summary, #按顺序指定每行是否加粗，T为加粗，F不加粗
                  txt_gp=fpTxtGp(
                    label=gpar(cex=0.5),
                    ticks=gpar(cex=0.5), 
                    xlab=gpar(cex=0.5), 
                    title=gpar(cex=0.5)),
                  lwd.zero=1,
                  lwd.ci=1.5,
                  lwd.xaxis=2, 
                  lty.ci=1.5,
                  ci.vertices =T,
                  ci.vertices.height=0.2, 
                  clip=c(0.1,8),
                  ineheight=unit(9, 'mm'), 
                  line.margin=unit(9, 'mm'),
                  colgap=unit(1.7, 'mm'),
                  fn.ci_norm="fpDrawDiamondCI", 
                  title="",
                  col=fpColors(box =colors, 
                               lines =colors, 
                               zero = "black"))       #森林图应插在图形第2列
while (!is.null(dev.list()))  dev.off()
pdf("nomogram/单因素Cox回归森林图.pdf",height = 8/2.54,width = 8.5/2.54)
print(fig2)
dev.off()

####多因素####

Uni_cox_1<-Uni_cox
Uni_cox_1$Char<-variable.names
Uni_cox_1<-Uni_cox_1[Uni_cox_1$P=='<0.05',]

variable.names_1<-Uni_cox_1$Char

formula<-'Surv(time=futime,event=fustat)~'
for(i in variable.names_1)formula<-ifelse(i==variable.names_1[1],paste(formula,i,sep = ''),paste(formula,i,sep = '+'))

mul_cox <- coxph(as.formula(formula), #对单因素COX中p<0.05的因素进行多因素cox分析
                 data=data)
mul_cox1 <- summary(mul_cox)
multi1 <- as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
#multi2：提取：HR(95%CI)和P
multi2<-ShowRegTable(mul_cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#将两次提取结果合并成表；取名result
result <-cbind(multi1,multi2)
#行名转为表格第一列，并给予命名"Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics")
result
####表格整理
#删除部分变量名，只保留亚变量
result$Characteristics<-str_remove(result$Characteristics,"age|gender|stage|Tstage|Nstage|Mstage|risk")
result

#给参考变量插入空行
ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
#插入空行，形成一个新表
result1<-data.frame("Characteristics", NA, NA, NA, "HR(95%CI)","p")
colnames(result1)<-colnames(result)

for(i in variable.names_1){
  a<-as.data.frame(table(data[,i]))
  result1<-rbind(result1,
                 ins(i),
                 ins(as.character(a$Var1[!a$Var1%in%result$Characteristics])),
                 result[which(variable.names_1==i),])
}

result1<-rbind(result1,
               c(NA, NA, NA, NA, NA,NA))

rownames(result1)<-1:nrow(result1)

myVars <- variable.names_1
catVars <-  variable.names_1
library(tableone)
table1<- print(CreateTableOne(vars=myVars,
                              data = data,
                              factorVars = catVars),
               showAllLevels=TRUE)

N<-data.frame(c(NA,NA),
              c(NA,NA))
colnames(N)<-colnames(table1)
for(i in 1:length(variable.names_1)){
  N=rbind(N,
          table1[(i*2):(i*2+1),],
          c(NA,NA))
}

N<-N[,-1]
N<-data.frame(N)
result2<-cbind(result1,N)
result2<-result2[,c(1,7,2:6)]

result2[1,]<-c("Characteristics","Number (%)",NA,NA,NA,"HR (95%CI)","P.value")

####图形进行美化
library(forestplot)
hrzl_lines<-list(gpar(lty=1,lwd=2),#表头上方添加实线
                 gpar(lty=2),#表头下方添加虚线
                 gpar(lwd=2,lty=1,columns=c(1:4)))#最后一行下方添加实线

names(hrzl_lines)<-c('1','2',nrow(result2)+1)

is_summary<-is.na(result2[,2]) #按顺序指定每行是否加粗，T为加粗，F不加粗
is_summary[1]=TRUE

fig2<- forestplot(result2[,c(1,2,6,7)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                  mean=as.numeric(result2[,3]),   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                  lower=as.numeric(result2[,4]),  #告诉函数表格第3列为5%CI，
                  upper=as.numeric(result2[,5]),  #表格第5列为95%CI，它俩要化作线段，穿过方块
                  zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                  boxsize=0.5,       #设置小黑块的大小
                  graph.pos= "right" ,
                  hrzl_lines=hrzl_lines, #最后一行下方添加实线
                  graphwidth = unit(.25,"npc"),
                  #xticks=c(-1,1,2,3,4,5) , #森林图刻度
                  is.summary=is_summary, #按顺序指定每行是否加粗，T为加粗，F不加粗
                  txt_gp=fpTxtGp(
                    label=gpar(cex=0.5),
                    ticks=gpar(cex=0.5), 
                    xlab=gpar(cex=0.5), 
                    title=gpar(cex=0.5)),
                  lwd.zero=1,
                  lwd.ci=1.5,
                  lwd.xaxis=2, 
                  lty.ci=1.5,
                  ci.vertices =T,
                  ci.vertices.height=0.2, 
                  clip=c(0.1,8),
                  ineheight=unit(9, 'mm'), 
                  line.margin=unit(9, 'mm'),
                  colgap=unit(1.7, 'mm'),
                  fn.ci_norm="fpDrawDiamondCI", 
                  title="",
                  col=fpColors(box =colors, 
                               lines =colors, 
                               zero = "black"))       #森林图应插在图形第2列
while (!is.null(dev.list()))  dev.off()
pdf("nomogram/多因素Cox回归森林图.pdf",height = 8/2.54,width = 8.5/2.54)
print(fig2)
dev.off()

####画nomogram图####

multi3<-data.frame(Char=variable.names_1,p=multi2[,2])
multi3$p_1<-ifelse(multi3$p=='<0.001','<0.001',ifelse(as.numeric(multi3$p)<0.05,'<0.05',NA))
multi3<-na.omit(multi3)

vl<-multi3$Char
formula<-'Surv(time=futime,event=fustat)~'
for(i in vl)formula<-ifelse(i==vl[1],paste(formula,i,sep = ''),paste(formula,i,sep = '+'))

####二分类变量####
dir.create('nomogram/categorical')
dd=datadist(data)#data必须是data.frame
options(datadist="dd")
fit<-cph(as.formula(formula),data=data,x=TRUE,y=TRUE,surv=TRUE)
mode(data)

survival<-Survival(fit)
survival1<-function(x)survival(years[1],x)
survival2<-function(x)survival(years[2],x)
survival3<-function(x)survival(years[3],x)

nomo<-nomogram(fit,fun=list(survival1,survival2,survival3),
               funlabel =paste0(years," Year Survival"))

#C-index检验预测效果
B<-rcorrcens(Surv(futime,fustat)~ predict(fit),data=data)
C_Index=(1-B[1])
print(C_Index)#0.7767172


save(nomo,C_Index,file = "nomogram/categorical/nomplot.rda")
while (!is.null(dev.list()))  dev.off()
pdf("nomogram/categorical/nomplot.pdf",width=6,height=6)
par(family="sans",ps="6") #sans为Arial字体，serif为新罗马字体 ,ps为字体磅值
plot(nomo, lplabel="Linear Predictor",
     xfrac=.25, #左边标题与右边图形间隔
     label.every = 1,col.grid = gray(c(0.8, 0.95)), #对应上方points的线条
     cex.var = 1.2,cex.axis = 1,cex.lab = 1.2, #cex文本属性
     lty=1, #lty指定线型
     lwd=5 #lwd改变线条粗细
)
title(main="")
dev.off()

####画校准曲线calibration curve
fit1<-cph(as.formula(formula),data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[1])
fit2<-cph(as.formula(formula),data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[2])
fit3<-cph(as.formula(formula),data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[3])

cal1<-calibrate(fit1, cmethod="KM",method="boot", u=years[1], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

cal2<-calibrate(fit2, cmethod="KM",method="boot", u=years[2], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

cal3<-calibrate(fit3, cmethod="KM",method="boot", u=years[3], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

while (!is.null(dev.list()))  dev.off()
pdf("nomogram/categorical/Calibration.pdf",width=4,height=4.5)
par(family="sans",ps="6",mar=c(6,4,1,1))
plot(cal1,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[1],col=colors[1],lwd=1,
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
#lines(cal1[,c("mean.predicted","KM")],lwd=1,col=colors[1],type= "b",cex=1,pch=15)#设置校准线

par(new=TRUE,family="sans",ps="0")
plot(cal2,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[2],col=colors[2],lwd=1,axes=FALSE,
     xlab="",ylab="", add=T)
#lines(cal2[,c("mean.predicted","KM")],lwd=1,col=colors[2],type= "b",cex=1.5,pch=16)#设置校准线粗细

par(new=TRUE,family="sans",ps="6")
plot(cal3,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[3],col=colors[3],lwd=1,axes=FALSE,
     xlab="",ylab="", add=T)  #errbar.col定义误差线的颜色，col定义校准曲线的颜色
#lines(cal3[,c("mean.predicted","KM")],lwd=1,col=colors[3],type= "b",cex=1.5,pch=17)#设置校准线

#abline(0,1,lty=3,lwd=1,col=c(rgb(0,118,192,maxColorValue=255)))

legend("topleft", inset=0.05,paste0(years,'-year'), col=colors[1:3], lty =1,lwd =1, bty = "n")

dev.off()

####用ROC曲线验证
library(rmda)

nomoRisk=predict(fit, data=data, type="lp")
data$nomoRisk=nomoRisk

rt <- data
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
while (!is.null(dev.list()))  dev.off()
pdf(file="nomogram/categorical/ROC.pdf",width=5,height=5)
plot(ROC_rt,time=1,col=colors[3],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=colors[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=colors[1],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c(colors[3],colors[2],colors[1]),lwd=2,bty = 'n')
dev.off()

####连续变量####
dir.create('nomogram/continuous')
data1<-clinical[rownames(data),]
data1$risk<-total_risk[rownames(data),'riskScore']

dd=datadist(data1)#data必须是data.frame
options(datadist="dd")
fit<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE)
mode(data1)

survival<-Survival(fit)
survival1<-function(x)survival(years[1],x)
survival2<-function(x)survival(years[2],x)
survival3<-function(x)survival(years[3],x)

nomo<-nomogram(fit,fun=list(survival1,survival2,survival3),
               funlabel = paste0(years," Year Survival"))

#C-index检验预测效果
B<-rcorrcens(Surv(futime,fustat)~ predict(fit),data=data1)
C_Index=(1-B[1])
print(C_Index)#0.7767172


save(nomo,C_Index,file = "nomogram/continuous/nomplot.rda")
while (!is.null(dev.list()))  dev.off()
pdf("nomogram/continuous/nomplot.pdf",width=6,height=6)
par(family="sans",ps="6") #sans为Arial字体，serif为新罗马字体 ,ps为字体磅值
plot(nomo, lplabel="Linear Predictor",
     xfrac=.25, #左边标题与右边图形间隔
     label.every = 1,col.grid = gray(c(0.8, 0.95)), #对应上方points的线条
     cex.var = 1.2,cex.axis = 1,cex.lab = 1.2, #cex文本属性
     lty=1, #lty指定线型
     lwd=5 #lwd改变线条粗细
)
title(main="")
dev.off()

####画校准曲线calibration curve
fit1<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[1])
fit2<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[2])
fit3<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[3])

cal1<-calibrate(fit1, cmethod="KM",method="boot", u=years[1], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

cal2<-calibrate(fit2, cmethod="KM",method="boot", u=years[2], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

cal3<-calibrate(fit3, cmethod="KM",method="boot", u=years[3], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

while (!is.null(dev.list()))  dev.off()
pdf("nomogram/continuous/Calibration.pdf",width=4,height=4.5)
par(family="sans",ps="6",mar=c(6,4,1,1))
plot(cal1,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[1],col=colors[1],lwd=1,
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
#lines(cal1[,c("mean.predicted","KM")],lwd=1,col=colors[1],type= "b",cex=1,pch=15)#设置校准线

par(new=TRUE,family="sans",ps="0")
plot(cal2,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[2],col=colors[2],lwd=1,axes=FALSE,
     xlab="",ylab="", add=T)
#lines(cal2[,c("mean.predicted","KM")],lwd=1,col=colors[2],type= "b",cex=1.5,pch=16)#设置校准线粗细

par(new=TRUE,family="sans",ps="6")
plot(cal3,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[3],col=colors[3],lwd=1,axes=FALSE,
     xlab="",ylab="", add=T)  #errbar.col定义误差线的颜色，col定义校准曲线的颜色
#lines(cal3[,c("mean.predicted","KM")],lwd=1,col=colors[3],type= "b",cex=1.5,pch=17)#设置校准线

#abline(0,1,lty=3,lwd=1,col=c(rgb(0,118,192,maxColorValue=255)))

legend("topleft", inset=0.05,paste0(years,'-year'), col=colors[1:3], lty =1,lwd =1, bty = "n")

dev.off()

####ROC曲线
library(rmda)

nomoRisk=predict(fit, data=data1, type="lp")
data1$nomoRisk=nomoRisk

rt <- data1
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
while (!is.null(dev.list()))  dev.off()
pdf(file="nomogram/continuous/ROC.pdf",width=5,height=5)
plot(ROC_rt,time=1,col=colors[3],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=colors[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=colors[1],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c(colors[3],colors[2],colors[1]),lwd=2,bty = 'n')
dev.off()


####连续变量仅risk####
dir.create('nomogram/continuous_only_risk')
data1<-data
data1$risk<-total_risk[rownames(data),'riskScore']

dd=datadist(data1)#data必须是data.frame
options(datadist="dd")
fit<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE)
mode(data1)

survival<-Survival(fit)
survival1<-function(x)survival(years[1],x)
survival2<-function(x)survival(years[2],x)
survival3<-function(x)survival(years[3],x)

nomo<-nomogram(fit,fun=list(survival1,survival2,survival3),
               funlabel = paste0(years," Year Survival"))

#C-index检验预测效果
B<-rcorrcens(Surv(futime,fustat)~ predict(fit),data=data1)
C_Index=(1-B[1])
print(C_Index)#0.7767172


save(nomo,C_Index,file = "nomogram/continuous_only_risk/nomplot.rda")
while (!is.null(dev.list()))  dev.off()
pdf("nomogram/continuous_only_risk/nomplot.pdf",width=6,height=6)
par(family="sans",ps="6") #sans为Arial字体，serif为新罗马字体 ,ps为字体磅值
plot(nomo, lplabel="Linear Predictor",
     xfrac=.25, #左边标题与右边图形间隔
     label.every = 1,col.grid = gray(c(0.8, 0.95)), #对应上方points的线条
     cex.var = 1.2,cex.axis = 1,cex.lab = 1.2, #cex文本属性
     lty=1, #lty指定线型
     lwd=5 #lwd改变线条粗细
)
title(main="")
dev.off()


####画校准曲线calibration curve
fit1<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[1])
fit2<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[2])
fit3<-cph(as.formula(formula),data=data1,x=TRUE,y=TRUE,surv=TRUE,time.inc=years[3])

cal1<-calibrate(fit1, cmethod="KM",method="boot", u=years[1], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

cal2<-calibrate(fit2, cmethod="KM",method="boot", u=years[2], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

cal3<-calibrate(fit3, cmethod="KM",method="boot", u=years[3], m=50, B=1000)#u为时间，B为重复的次数，m一般为50、100，和总样本量有关

while (!is.null(dev.list()))  dev.off()
pdf("nomogram/continuous_only_risk/Calibration.pdf",width=4,height=4.5)
par(family="sans",ps="6",mar=c(6,4,1,1))
plot(cal1,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[1],col=colors[1],lwd=1,
     xlab="Nomogram-Predicted OS",ylab="Observed OS(%)")
#lines(cal1[,c("mean.predicted","KM")],lwd=1,col=colors[1],type= "b",cex=1,pch=15)#设置校准线

par(new=TRUE,family="sans",ps="0")
plot(cal2,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[2],col=colors[2],lwd=1,axes=FALSE,
     xlab="",ylab="", add=T)
#lines(cal2[,c("mean.predicted","KM")],lwd=1,col=colors[2],type= "b",cex=1.5,pch=16)#设置校准线粗细

par(new=TRUE,family="sans",ps="6")
plot(cal3,xlim = c(0,1),ylim= c(0,1),errbar.col=colors[3],col=colors[3],lwd=1,axes=FALSE,
     xlab="",ylab="", add=T)  #errbar.col定义误差线的颜色，col定义校准曲线的颜色
#lines(cal3[,c("mean.predicted","KM")],lwd=1,col=colors[3],type= "b",cex=1.5,pch=17)#设置校准线

#abline(0,1,lty=3,lwd=1,col=c(rgb(0,118,192,maxColorValue=255)))

legend("topleft", inset=0.05,paste0(years,'-year'), col=colors[1:3], lty =1,lwd =1, bty = "n")

dev.off()

####ROC曲线
library(rmda)

nomoRisk=predict(fit, data=data1, type="lp")
data1$nomoRisk=nomoRisk

rt <- data1
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
while (!is.null(dev.list()))  dev.off()
pdf(file="nomogram/continuous_only_risk/ROC.pdf",width=5,height=5)
plot(ROC_rt,time=1,col=colors[3],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=colors[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=colors[1],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c(colors[3],colors[2],colors[1]),lwd=2,bty = 'n')
dev.off()

