
total_tumor<-read.csv(tumor_exp_path,row.names = 1)
total_tumor_cli<-read.csv(tumor_cli_path,row.names = 1)
colnames(total_tumor)<-gsub('[.]','-',colnames(total_tumor))

gene<-read.csv('WGCNA/关键表型相关基因.csv',row.names = 1)[,1]

total<-cbind(t(total_tumor[gene,]),total_tumor_cli[,c('fustat','futime')])

library(survival)
library(survminer)
sur.res <- matrix(nrow = ncol(total)-2, ncol = 6)
for (j in 1:nrow(sur.res)) { #
  Bcox<-coxph(Surv(futime, fustat)~as.numeric(as.character(total[,j]))>median(as.numeric(as.character(total[,j]))),data=total)
  summcph<-summary(Bcox)
  sur.res[j,1]<-summcph$conf.int[1]
  sur.res[j,2]<-summcph$conf.int[3]
  sur.res[j,3]<-summcph$conf.int[4]
  sur.res[j,4]<-as.matrix(summcph$logtest)[3]
  sur.res[j,5]<-as.matrix(summcph$sctest)[3]
  sur.res[j,6]<-summcph$coefficients[5]
}
rownames(sur.res)=colnames(total)[1:(nrow(sur.res))]#
colnames(sur.res)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")
sur.res<-as.data.frame(sur.res)
sur.res1<-sur.res[sur.res$p_value<0.05,]
sur.res1<-na.omit(sur.res1)

cox_gene<-rownames(sur.res1)
dir.create('prognosis')
write.csv(sur.res,'prognosis/cox.csv')
write.csv(cox_gene,'prognosis/cox_gene.csv')

library(timeROC)
seed<-data.frame(seed=0,p=0,p1=0,train1=0,train2=0,train3=0,varify1=0,varify2=0,varify3=0)
pb <- txtProgressBar(style=3)
for(i in seed_range){
  set.seed(i)
  ind <- sample(2, nrow(total), replace = T, prob = c(0.7, 0.3))
  train <- total[ind==1,]
  varify <- total[ind==2,]
  if(length(cox_gene)>2){
    ####lasso####
    #对与生存显著相关的基因进行LASSO回归构建预后模型
    library(glmnet)
    x<-as.matrix(train[,cox_gene])
    y<-train[,c("fustat","futime")]
    y <- as.matrix(survival::Surv(y$futime, y$fustat))
    
    fit <-glmnet(x,y,family = "cox",alpha = 1)
    #主要在做交叉验证,lasso
    fitcv <- cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)
    coefficient <- coef(fitcv, s="lambda.min")
    Active.Index <- which(as.numeric(coefficient) != 0)
    if(length(Active.Index)>1){
      lasso_coefficients <- as.numeric(coefficient)[Active.Index]
      lasso_gene <- rownames(coefficient)[Active.Index]
      names(lasso_coefficients)<-lasso_gene
      
      lasso_train<-train[,c(lasso_gene,'fustat','futime')]
      lasso_varify<-varify[,c(lasso_gene,'fustat','futime')]
      
      x_risk<-function(x){
        a<-as.data.frame(lapply(lasso_gene,function(i){x[,i]*lasso_coefficients[i]}))
        a$risk<-apply(a,1,sum)
        return(a$risk)
      }
      
      lasso_train$risk<-x_risk(lasso_train)
      lasso_varify$risk<-x_risk(lasso_varify)
      #采用surv_cutpoint函数选择最佳截断值区分高低风险组
      cutoff<-surv_cutpoint(lasso_train,time = "futime",event = "fustat",variables = "risk",minprop = 0.3)
      cut=cutoff$cutpoint$cutpoint
      lasso_train$class<-ifelse(lasso_train$risk>=cut,'High','Low')
      #生存分析
      diff<-survdiff(Surv(futime,fustat)~lasso_train$class,data = lasso_train)
      p.val <- 1 - pchisq(diff$chisq, length(diff$n) - 1) 
      
      #测试集数据进行验证
      lasso_varify$class<-ifelse(lasso_varify$risk>=cut,'High','Low')
      diff<-survdiff(Surv(futime,fustat)~lasso_varify$class,data = lasso_varify)
      p.val1 <- 1 - pchisq(diff$chisq, length(diff$n) - 1) 
      
      if(p.val<0.05&p.val1<0.05){
        ROC_rt1=timeROC(T=lasso_train$futime, delta=lasso_train$fustat,
                        marker=lasso_train$risk, cause=1,
                        weighting='aalen',
                        times=years, ROC=TRUE)
        ROC_rt2=timeROC(T=lasso_varify$futime, delta=lasso_varify$fustat,
                        marker=lasso_varify$risk, cause=1,
                        weighting='aalen',
                        times=years, ROC=TRUE)
        if(ROC_rt1$AUC[1]>0.6&ROC_rt1$AUC[2]>0.6&ROC_rt1$AUC[3]>0.6){
          if(ROC_rt2$AUC[1]>0.6&ROC_rt2$AUC[2]>0.6&ROC_rt2$AUC[3]>0.6)seed=rbind(seed,
                                                                              data.frame(seed=i,
                                                                                         p=p.val,
                                                                                         p1=p.val1,
                                                                                         train1=ROC_rt1$AUC[1],
                                                                                         train2=ROC_rt1$AUC[2],
                                                                                         train3=ROC_rt1$AUC[3],
                                                                                         varify1=ROC_rt2$AUC[1],
                                                                                         varify2=ROC_rt2$AUC[2],
                                                                                         varify3=ROC_rt2$AUC[3]))
        }
      }
    }}
  setTxtProgressBar(pb, (i - seed_range[1]+1)/(seed_range[length(seed_range)]-seed_range[1]+1))}
close(pb)
seed=seed[-1,]
write.csv(seed,'lasso_cox_seed.csv')
View(seed)
