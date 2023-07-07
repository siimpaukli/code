setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical')

###############################


GSE10846_clinical_plot<-merge(GSE10846_lar_score_final,clinical_GSE10846,by='sample')
GSE10846_clinical_plot[,6:7]<-NULL
colnames(GSE10846_clinical_plot)[2:3]<-c("status","time")

GSE10846_clinical_plot$age<-ifelse(GSE10846_clinical_plot$age <=60 ,"<=60",">60")


GSE181063_clinical_1<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical.txt',header = T,check.names = F,sep = '\t')
GSE181063_lar_score_final

GSE181063_clinical_plot<-merge(GSE181063_lar_score_final,GSE181063_clinical_1,by='sample')
GSE181063_clinical_plot$age<-ifelse(GSE181063_clinical_plot$age <=60 ,"<=60",">60")


library(tidyr)
GSE10846_clinical_plot_high<-GSE10846_clinical_plot[which(GSE10846_clinical_plot$riskgroup == "High"),]
GSE10846_clinical_plot_high[,1:4]<-NULL

GSE10846_clinical_plot_low<-GSE10846_clinical_plot[which(GSE10846_clinical_plot$riskgroup == "Low"),]
GSE10846_clinical_plot_low[,1:4]<-NULL

GSE10846_clinical_high<-gather(GSE10846_clinical_plot_high,Clinial,type,gender:stage)
GSE10846_clinical_low<-gather(GSE10846_clinical_plot_low,Clinial,type,gender:stage)



GSE10846_clinical_plot[,1:4]<-NULL
GSE10846_clinical_plot_final<-gather(GSE10846_clinical_plot,Clinial,type,gender:stage)

clinical_type<-as.matrix(as.data.frame(table(GSE10846_clinical_plot_final$Clinial))[,1])

j<-1
for (i in 1:4)
{
  clinical<-clinical_type[1,i]
  temp<-GSE10846_clinical_plot_final[which(GSE10846_clinical_plot_final$Clinial %in%clinical),]
  if(clinical == "age")
  {
    my_comparisons<-list(c("<=60",">60"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = c("aquamarine3","lightcoral"),scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("<=60",">60"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  if(clinical == "gender")
  {
    my_comparisons<-list(c("Male","Female"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = c("aquamarine3","lightcoral"),scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("Female","Male"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  if(clinical == "ecog")
  {
    my_comparisons<-list(c("0","1"),c("0","2"),c("0","3"),c("0","4"),c("1","2"),c("1","3"),c("1","4"),c("2","3"),c("2","4"),c("3","4"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = "jco",scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("0","1","2","3","4"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  if(clinical == "stage")
  {
    my_comparisons<-list(c("1","2"),c("1","3"),c("1","4"),c("2","3"),c("2","4"),c("3","4"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = "jco",scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("1","2","3","4"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  assign(paste0("p", j),p)
  j<-j+1
}


p131 <- cowplot::plot_grid(p1,p3,p2,p4,ncol = 4,nrow = 1,labels = 'AUTO',rel_widths = c(2,2,5,4),label_size = 20,scale = c(0.9,0.9,0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/LARrisk_clinical_GSE10846.pdf' ,
                 plot = p131,width = 15,height = 5)



###############################
library(tidyr)


GSE181063_clinical_plot[,1:4]<-NULL
GSE181063_clinical_plot_final<-gather(GSE181063_clinical_plot,Clinial,type,age:ipiScore)

clinical_type<-as.matrix(as.data.frame(table(GSE181063_clinical_plot_final$Clinial))[,1])
GSE181063_clinical_plot_final<-na.omit(GSE181063_clinical_plot_final)

j<-1
for (i in 1:4)
{
  clinical<-clinical_type[i,1]
  temp<-GSE181063_clinical_plot_final[which(GSE181063_clinical_plot_final$Clinial %in%clinical),]
  if(clinical == "age")
  {
    my_comparisons<-list(c("<=60",">60"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = c("aquamarine3","lightcoral"),scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("<=60",">60"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  if(clinical == "gender")
  {
    my_comparisons<-list(c("M","F"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = c("aquamarine3","lightcoral"),scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("F","M"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  if(clinical == "ipiScore")
  {
    my_comparisons<-list(c("ipi_score=0","ipi_score=1"),c("ipi_score=0","ipi_score=2"),c("ipi_score=0","ipi_score=3"),c("ipi_score=0","ipi_score=4"),c("ipi_score=1","ipi_score=2"),c("ipi_score=1","ipi_score=3"),c("ipi_score=1","ipi_score=4"),c("ipi_score=2","ipi_score=3"),c("ipi_score=2","ipi_score=4"),c("ipi_score=3","ipi_score=4"),c("ipi_score=0","ipi_score=5"),c("ipi_score=1","ipi_score=5"),c("ipi_score=2","ipi_score=5"),c("ipi_score=3","ipi_score=5"),c("ipi_score=4","ipi_score=5"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = "jco",scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("ipi_score=0","ipi_score=1","ipi_score=2","ipi_score=3","ipi_score=4","ipi_score=5"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  if(clinical == "stage")
  {
    my_comparisons<-list(c("Stage I","Stage II"),c("Stage I","Stage III"),c("Stage I","Stage IV"),c("Stage II","Stage III"),c("Stage II","Stage IV"),c("Stage III","Stage IV"))
    p<-ggboxplot(temp,x="type" , y="riskscore" ,fill= "type", palette = "jco",scales = "free_x",xlab=F,ylab="LARscore",width = 0.5,order=c("Stage I","Stage II","Stage III","Stage IV"))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(clinical)
  }
  assign(paste0("p", j),p)
  j<-j+1
}


p131 <- cowplot::plot_grid(p1,p2,p3,p4,ncol = 4,nrow = 1,labels = 'AUTO',rel_widths = c(2,2,6,4),label_size = 20,scale = c(0.9,0.9,0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/LARrisk_clinical_GSE181063.pdf' ,
                 plot = p131,width = 15,height = 5)







###############################

GSE10846_clinical_multicox<-merge(GSE10846_lar_score_final,clinical_GSE10846,by='sample')
GSE10846_clinical_multicox[,6:7]<-NULL
colnames(GSE10846_clinical_multicox)[2:3]<-c("status","time")

GSE10846_clinical_multicox$age<-ifelse(GSE10846_clinical_multicox$age <=60 ,"<=60",">60")
GSE10846_clinical_multicox$riskscore<-NULL
GSE10846_clinical_multicox$age<-factor(GSE10846_clinical_multicox$age,levels = c("<=60",">60"))
GSE10846_clinical_multicox$age=relevel(GSE10846_clinical_multicox$age, ref = "<=60")



GSE10846_clinical_multicox$stage<-ifelse(GSE10846_clinical_multicox$stage == "1" | GSE10846_clinical_multicox$stage == "2","I+II","III+IV")
GSE10846_clinical_multicox$stage<-factor(GSE10846_clinical_multicox$stage,levels = c("I+II","III+IV"))
GSE10846_clinical_multicox$stage=relevel(GSE10846_clinical_multicox$stage, ref = "I+II")


GSE10846_clinical_multicox$ecog<-ifelse(GSE10846_clinical_multicox$ecog == "0" | GSE10846_clinical_multicox$ecog == "1","0+1","2+3+4")
GSE10846_clinical_multicox$ecog<-factor(GSE10846_clinical_multicox$ecog,levels = c("0+1","2+3+4"))
GSE10846_clinical_multicox$ecog=relevel(GSE10846_clinical_multicox$ecog, ref = "0+1")


GSE10846_clinical_multicox$gender<-factor(GSE10846_clinical_multicox$gender,levels = c("Female","Male"))
GSE10846_clinical_multicox$gender=relevel(GSE10846_clinical_multicox$gender, ref = "Female")


GSE10846_clinical_multicox$riskgroup<-factor(GSE10846_clinical_multicox$riskgroup,levels = c("Low","High"))
GSE10846_clinical_multicox$riskgroup=relevel(GSE10846_clinical_multicox$riskgroup, ref = "Low")



library("survival")
library("survminer") 
library("survivalROC") 
library("pROC")
library("timeROC")

covariates <- c(colnames(GSE10846_clinical_multicox)[4:8])
#covariates<-as.factor(covariates)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = GSE10846_clinical_multicox)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         c_index <-x$concordance[1]
                         p.value<-signif(x$wald["pvalue"], digits=6)
                         wald.test<-signif(x$wald["test"], digits=6)
                         beta<-signif(x$coef[1], digits=6);
                         HR <-signif(x$coef[2], digits=6);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 6)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],6)
                         CI_for_HR<-paste0(" (", 
                                           HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper,CI_for_HR, wald.test, p.value,c_index)
                         names(res)<-c("beta", "HR","CI_lower","CI_upper", "(95%_CI_for_HR)", "wald.test", 
                                       "p.value","C-index")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res_tcga<- t(as.data.frame(univ_results, check.names = FALSE))
res_tcga<-as.data.frame(res_tcga)
res_tcga$p.value<-as.numeric(res_tcga$p.value)
res_tcga<-res_tcga[order(res_tcga$p.value),]

write.table(res_tcga,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_clinical_single_cox_results.txt',row.names=T,col.names=T,quote=F,sep="\t")

res<-res_tcga
rownames(res)
rownames(res)[1]<-"ecog (2+3+4 vs 0+1)"
rownames(res)[4]<-"stage (III+IV vs I+II)"
rownames(res)[5]<-"gender (Male vs Female)"
rownames(res)[3]<-"age (>60 vs <=60)"
rownames(res)[2]<-"LARscore (High vs Low)"

library(forestplot)
forest_table<-data.frame(Risk_factors=rownames(res),HR=res$HR,`(95%CI)`=res$`(95%_CI_for_HR)`, pvalue=res$p.value,check.names = F)
forest_table1<-data.frame(Risk_factors=NA,HR=NA,`(95%CI)`=NA, pvalue=NA,check.names = F)
forest_table2<-data.frame(Risk_factors="Risk_factors",HR="HR",`(95%CI)`="(95%CI)", pvalue="p-value",check.names = F)
tabletext<-rbind(forest_table2,forest_table1,forest_table)
tabletext$sig<-NA
for(i in 3:7)
{
  if(as.numeric(tabletext$pvalue[i])< 0.05 & as.numeric(tabletext$pvalue[i]) > 0.01 )
  {
    tabletext$sig[i]<-"*"
  }
  else if(as.numeric(tabletext$pvalue[i])< 0.01 & as.numeric(tabletext$pvalue[i])> 0.001)
  {
    tabletext$sig[i]<-"**"
  }
  else if (as.numeric(tabletext$pvalue[i]) < 0.001)
  {
    tabletext$sig[i]<-"***"
  }
  else
  {
    tabletext$sig[i]<-" "
  }
}
tabletext$sig[1]<-"Significant"
forest_stastic<-data.frame(mean=as.numeric(res$HR),lower=as.numeric(res$CI_lower),upper=as.numeric(res$CI_upper))
forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
cochrane_from_rmeta<-rbind(forest_stastic1,forest_stastic)


pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_clinical_single_forest.pdf',width = 8,height = 4)
forestplot(tabletext,  #显示的文本
           mean=c(NA,cochrane_from_rmeta$mean), #误差条的均值(此处为差值的中值)
           lower=c(NA,cochrane_from_rmeta$lower), #误差条的下界(此处为差值的25%分位数)
           upper=c(NA,cochrane_from_rmeta$upper), #误差条的上界(此处为差值的75%分位数)
           graph.pos=2,
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.1, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 1,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.7), cex = 0.7), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="Hazard ratio", #x轴的标题
           title="Univariate Cox Regression (GSE10846)"
)
dev.off()




fmla2 <- as.formula(paste0("Surv(time, status)~",paste0(colnames(GSE10846_clinical_multicox)[4:8],collapse = '+')))
multi.res2 <- coxph(fmla2, data = GSE10846_clinical_multicox)


MultiName <- as.character(colnames(GSE10846_clinical_multicox)[4:8])
Multisum <- summary(multi.res2)
MHR <- round(Multisum$coefficients[, 2], 2)
MPValue <- round(Multisum$coefficients[, 5], 3)
MCIL <- paste0(round(Multisum$conf.int[, 3], 2))
MCIU <- paste0(round(Multisum$conf.int[, 4], 2))
MCI <- paste0(MCIL, "-", MCIU)
Multicox <- data.frame("Factor" = MultiName,
                       "HR" = MHR,
                       "CI95" = MCI,
                       "PValue" = MPValue,
                       "CI_lower"= MCIL,
                       "CI_upper"=MCIU)

Multicox<-Multicox[order(Multicox$PValue),]

write.table(Multicox,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_clinical_multi_cox_results.txt',row.names=T,col.names=T,quote=F,sep="\t")

library(forestplot)

Multicox$Factor[2]<-"ecog (2+3+4 vs 0+1)"
Multicox$Factor[4]<-"stage (III+IV vs I+II)"
Multicox$Factor[5]<-"gender (Male vs Female)"
Multicox$Factor[3]<-"age (>60 vs <=60)"
Multicox$Factor[1]<-"LARscore (High vs Low)"




forest_table<-data.frame(Risk_factors=Multicox$Factor,HR=Multicox$HR,`(95%CI)`=Multicox$CI95, pvalue=Multicox$PValue,check.names = F)
forest_table1<-data.frame(Risk_factors=NA,HR=NA,`(95%CI)`=NA, pvalue=NA,check.names = F)
forest_table2<-data.frame(Risk_factors="Risk_factors",HR="HR",`(95%CI)`="(95%CI)", pvalue="p-value",check.names = F)
tabletext<-rbind(forest_table2,forest_table1,forest_table)
tabletext$sig<-NA
for(i in 3:7)
{
  if(as.numeric(tabletext$pvalue[i])< 0.05 & as.numeric(tabletext$pvalue[i]) > 0.01 )
  {
    tabletext$sig[i]<-"*"
  }
  else if(as.numeric(tabletext$pvalue[i])< 0.01 & as.numeric(tabletext$pvalue[i])> 0.001)
  {
    tabletext$sig[i]<-"**"
  }
  else if (as.numeric(tabletext$pvalue[i]) < 0.001)
  {
    tabletext$sig[i]<-"***"
  }
  else
  {
    tabletext$sig[i]<-" "
  }
}
tabletext$sig[1]<-"Significant"
forest_stastic<-data.frame(mean=as.numeric(Multicox$HR),lower=as.numeric(Multicox$CI_lower),upper=as.numeric(Multicox$CI_upper))
forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
cochrane_from_rmeta<-rbind(forest_stastic1,forest_stastic)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_clinical_multi_forest.pdf',width = 8,height = 4)
forestplot(tabletext,  #显示的文本
           mean=c(NA,cochrane_from_rmeta$mean), #误差条的均值(此处为差值的中值)
           lower=c(NA,cochrane_from_rmeta$lower), #误差条的下界(此处为差值的25%分位数)
           upper=c(NA,cochrane_from_rmeta$upper), #误差条的上界(此处为差值的75%分位数)
           graph.pos=2,
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.1, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 1,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.7), cex = 0.7), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="Hazard ratio",#x轴的标题
           title='Multivariate Cox regression (GSE10846)'
)
dev.off()


###############################
GSE181063_clinical_1<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical.txt',header = T,check.names = F,sep = '\t')
GSE181063_clinical_multicox<-merge(GSE181063_lar_score_final,GSE181063_clinical_1,by='sample')

GSE181063_clinical_multicox$riskscore<-NULL



GSE181063_clinical_multicox$age<-ifelse(GSE181063_clinical_multicox$age <=60 ,"<=60",">60")
GSE181063_clinical_multicox$age<-factor(GSE181063_clinical_multicox$age,levels = c("<=60",">60"))
GSE181063_clinical_multicox$age=relevel(GSE181063_clinical_multicox$age, ref = "<=60")



GSE181063_clinical_multicox$ipiScore<-ifelse(GSE181063_clinical_multicox$ipiScore == "ipi_score=0" | GSE181063_clinical_multicox$ipiScore == "ipi_score=1" | GSE181063_clinical_multicox$ipiScore == "ipi_score=2","ipi_score_0+1+2","ipi_score_3+4+5")
GSE181063_clinical_multicox$ipiScore<-factor(GSE181063_clinical_multicox$ipiScore,levels = c("ipi_score_0+1+2","ipi_score_3+4+5"))
GSE181063_clinical_multicox$ipiScore=relevel(GSE181063_clinical_multicox$ipiScore, ref = "ipi_score_0+1+2")


GSE181063_clinical_multicox$stage<-ifelse(GSE181063_clinical_multicox$stage == "Stage I" | GSE181063_clinical_multicox$stage == "Stage II","I+II","III+IV")
GSE181063_clinical_multicox$stage<-factor(GSE181063_clinical_multicox$stage,levels = c("I+II","III+IV"))
GSE181063_clinical_multicox$stage=relevel(GSE181063_clinical_multicox$stage, ref = "I+II")


GSE181063_clinical_multicox$gender<-factor(GSE181063_clinical_multicox$gender,levels = c("F","M"))
GSE181063_clinical_multicox$gender=relevel(GSE181063_clinical_multicox$gender, ref = "F")


GSE181063_clinical_multicox$riskgroup<-factor(GSE181063_clinical_multicox$riskgroup,levels = c("Low","High"))
GSE181063_clinical_multicox$riskgroup=relevel(GSE181063_clinical_multicox$riskgroup, ref = "Low")



library("survival")
library("survminer") 
library("survivalROC") 
library("pROC")
library("timeROC")

covariates <- c(colnames(GSE181063_clinical_multicox)[4:8])
#covariates<-as.factor(covariates)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = GSE181063_clinical_multicox)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         c_index <-x$concordance[1]
                         p.value<-signif(x$wald["pvalue"], digits=6)
                         wald.test<-signif(x$wald["test"], digits=6)
                         beta<-signif(x$coef[1], digits=6);
                         HR <-signif(x$coef[2], digits=6);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 6)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],6)
                         CI_for_HR<-paste0(" (", 
                                           HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper,CI_for_HR, wald.test, p.value,c_index)
                         names(res)<-c("beta", "HR","CI_lower","CI_upper", "(95%_CI_for_HR)", "wald.test", 
                                       "p.value","C-index")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res_tcga<- t(as.data.frame(univ_results, check.names = FALSE))
res_tcga<-as.data.frame(res_tcga)
res_tcga$p.value<-as.numeric(res_tcga$p.value)
res_tcga<-res_tcga[order(res_tcga$p.value),]

write.table(res_tcga,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical_single_cox_results.txt',row.names=T,col.names=T,quote=F,sep="\t")

res<-res_tcga
rownames(res)
rownames(res)[1]<-"ipiScore (3+4+5 vs 0+1+2)"
rownames(res)[4]<-"stage (III+IV vs I+II)"
rownames(res)[5]<-"gender (Male vs Female)"
rownames(res)[2]<-"age >60 vs <=60)"
rownames(res)[3]<-"LARscore (High vs Low)"

library(forestplot)
forest_table<-data.frame(Risk_factors=rownames(res),HR=res$HR,`(95%CI)`=res$`(95%_CI_for_HR)`, pvalue=res$p.value,check.names = F)
forest_table1<-data.frame(Risk_factors=NA,HR=NA,`(95%CI)`=NA, pvalue=NA,check.names = F)
forest_table2<-data.frame(Risk_factors="Risk_factors",HR="HR",`(95%CI)`="(95%CI)", pvalue="p-value",check.names = F)
tabletext<-rbind(forest_table2,forest_table1,forest_table)
tabletext$sig<-NA
for(i in 3:7)
{
  if(as.numeric(tabletext$pvalue[i])< 0.05 & as.numeric(tabletext$pvalue[i]) > 0.01 )
  {
    tabletext$sig[i]<-"*"
  }
  else if(as.numeric(tabletext$pvalue[i])< 0.01 & as.numeric(tabletext$pvalue[i])> 0.001)
  {
    tabletext$sig[i]<-"**"
  }
  else if (as.numeric(tabletext$pvalue[i]) < 0.001)
  {
    tabletext$sig[i]<-"***"
  }
  else
  {
    tabletext$sig[i]<-" "
  }
}
tabletext$sig[1]<-"Significant"
forest_stastic<-data.frame(mean=as.numeric(res$HR),lower=as.numeric(res$CI_lower),upper=as.numeric(res$CI_upper))
forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
cochrane_from_rmeta<-rbind(forest_stastic1,forest_stastic)


pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical_single_forest.pdf',width = 8,height = 4)
forestplot(tabletext,  #显示的文本
           mean=c(NA,cochrane_from_rmeta$mean), #误差条的均值(此处为差值的中值)
           lower=c(NA,cochrane_from_rmeta$lower), #误差条的下界(此处为差值的25%分位数)
           upper=c(NA,cochrane_from_rmeta$upper), #误差条的上界(此处为差值的75%分位数)
           graph.pos=2,
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.1, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 1,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.7), cex = 0.7), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="Hazard ratio", #x轴的标题
           title="Univariate Cox Regression (GSE181063)"
)
dev.off()




fmla2 <- as.formula(paste0("Surv(time, status)~",paste0(colnames(GSE181063_clinical_multicox)[4:8],collapse = '+')))
multi.res2 <- coxph(fmla2, data = GSE181063_clinical_multicox)


MultiName <- as.character(colnames(GSE181063_clinical_multicox)[4:8])
Multisum <- summary(multi.res2)
MHR <- round(Multisum$coefficients[, 2], 2)
MPValue <- round(Multisum$coefficients[, 5], 3)
MCIL <- paste0(round(Multisum$conf.int[, 3], 2))
MCIU <- paste0(round(Multisum$conf.int[, 4], 2))
MCI <- paste0(MCIL, "-", MCIU)
Multicox <- data.frame("Factor" = MultiName,
                       "HR" = MHR,
                       "CI95" = MCI,
                       "PValue" = MPValue,
                       "CI_lower"= MCIL,
                       "CI_upper"=MCIU)

Multicox<-Multicox[order(Multicox$PValue),]

write.table(Multicox,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical_multi_cox_results.txt',row.names=T,col.names=T,quote=F,sep="\t")

library(forestplot)

Multicox$Factor[4]<-"ipiScore (3+4+5 vs 0+1+2)"
Multicox$Factor[5]<-"stage (III+IV vs I+II)"
Multicox$Factor[3]<-"gender (Male vs Female)"
Multicox$Factor[2]<-"age >60 vs <=60)"
Multicox$Factor[1]<-"LARscore (High vs Low)"




forest_table<-data.frame(Risk_factors=Multicox$Factor,HR=Multicox$HR,`(95%CI)`=Multicox$CI95, pvalue=Multicox$PValue,check.names = F)
forest_table1<-data.frame(Risk_factors=NA,HR=NA,`(95%CI)`=NA, pvalue=NA,check.names = F)
forest_table2<-data.frame(Risk_factors="Risk_factors",HR="HR",`(95%CI)`="(95%CI)", pvalue="p-value",check.names = F)
tabletext<-rbind(forest_table2,forest_table1,forest_table)
tabletext$sig<-NA
for(i in 3:7)
{
  if(as.numeric(tabletext$pvalue[i])< 0.05 & as.numeric(tabletext$pvalue[i]) > 0.01 )
  {
    tabletext$sig[i]<-"*"
  }
  else if(as.numeric(tabletext$pvalue[i])< 0.01 & as.numeric(tabletext$pvalue[i])> 0.001)
  {
    tabletext$sig[i]<-"**"
  }
  else if (as.numeric(tabletext$pvalue[i]) < 0.001)
  {
    tabletext$sig[i]<-"***"
  }
  else
  {
    tabletext$sig[i]<-" "
  }
}
tabletext$sig[1]<-"Significant"
forest_stastic<-data.frame(mean=as.numeric(Multicox$HR),lower=as.numeric(Multicox$CI_lower),upper=as.numeric(Multicox$CI_upper))
forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
cochrane_from_rmeta<-rbind(forest_stastic1,forest_stastic)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical_multi_forest.pdf',width = 8,height = 4)
forestplot(tabletext,  #显示的文本
           mean=c(NA,cochrane_from_rmeta$mean), #误差条的均值(此处为差值的中值)
           lower=c(NA,cochrane_from_rmeta$lower), #误差条的下界(此处为差值的25%分位数)
           upper=c(NA,cochrane_from_rmeta$upper), #误差条的上界(此处为差值的75%分位数)
           graph.pos=2,
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.1, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 1,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.7), cex = 0.7), #文本大小设置
           lineheight = "auto", #线的高度 
           xlab="Hazard ratio",#x轴的标题
           title='Multivariate Cox regression (GSE181063)'
)
dev.off()


###############################

library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2) 
library(rms)

GSE10846_clinical_multicox<-merge(GSE10846_lar_score_final,clinical_GSE10846,by='sample')
GSE10846_clinical_multicox[,6:7]<-NULL
colnames(GSE10846_clinical_multicox)[2:3]<-c("status","time")

GSE10846_clinical_multicox$age<-ifelse(GSE10846_clinical_multicox$age <=60 ,"<=60",">60")
GSE10846_clinical_multicox$age<-factor(GSE10846_clinical_multicox$age,levels = c("<=60",">60"))
GSE10846_clinical_multicox$age=relevel(GSE10846_clinical_multicox$age, ref = "<=60")



GSE10846_clinical_multicox$stage<-ifelse(GSE10846_clinical_multicox$stage == "1" | GSE10846_clinical_multicox$stage == "2","I+II","III+IV")
GSE10846_clinical_multicox$stage<-factor(GSE10846_clinical_multicox$stage,levels = c("I+II","III+IV"))
GSE10846_clinical_multicox$stage=relevel(GSE10846_clinical_multicox$stage, ref = "I+II")


GSE10846_clinical_multicox$ecog<-ifelse(GSE10846_clinical_multicox$ecog == "0" | GSE10846_clinical_multicox$ecog == "1","0+1","2+3+4")
GSE10846_clinical_multicox$ecog<-factor(GSE10846_clinical_multicox$ecog,levels = c("0+1","2+3+4"))
GSE10846_clinical_multicox$ecog=relevel(GSE10846_clinical_multicox$ecog, ref = "0+1")


GSE10846_clinical_multicox$gender<-factor(GSE10846_clinical_multicox$gender,levels = c("Female","Male"))
GSE10846_clinical_multicox$gender=relevel(GSE10846_clinical_multicox$gender, ref = "Female")

GSE10846_clinical_multicox$riskgroup<-NULL


gse_data<-GSE10846_clinical_multicox[,2:8]
dd=datadist(gse_data)
options(datadist="dd") 
f2 <- psm(Surv(time, status) ~ riskscore  + gender + age + ecog + stage, data =  gse_data)
med <- Quantile(f2) # 计算中位生存时间
surv <- Survival(f2) # 构建生存概率函数
nom <- nomogram(f2, fun=list(function(x) surv(1, x),
                             function(x) surv(3, x),
                             function(x) surv(5, x),
                             function(x) surv(10, x)),
                funlabel=c("1-year Survival Probability",
                           "3-year Survival Probability",
                           "5-year Survival Probability",
                           "10-year Survival Probability"))
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_clinical_nomogram_forest.pdf',width = 20,height = 10)
plot(nom, xfrac=.6)
dev.off()


fcox1 <- cph(Surv(time, status)~ riskscore ,surv=T,x=T, y=T,time.inc = 1,data=gse_data)
cal1 <- calibrate(fcox1, cmethod="KM", method="boot", u=1, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_cali_curve_1y.pdf',width = 6,height = 6)
plot(cal1,col="blue",xlim=c(0.5,1),ylim=c(0.5,1),xlab="Predicted Probability 1-Year OS",ylab="Actucal 1-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()


fcox3 <- cph(Surv(time, status)~ riskscore ,surv=T,x=T, y=T,time.inc = 3,data=gse_data)
cal3 <- calibrate(fcox3, cmethod="KM", method="boot", u=3, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_cali_curve_3y.pdf',width = 6,height = 6)
plot(cal3,col="blue",xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability 3-Year OS",ylab="Actucal 3-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()

fcox5 <- cph(Surv(time, status)~ riskscore  ,surv=T,x=T, y=T,time.inc = 5,data=gse_data)
cal5 <- calibrate(fcox5, cmethod="KM", method="boot", u=5, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_cali_curve_5y.pdf',width = 6,height = 6)
plot(cal5,col="blue",xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability 5-Year OS",ylab="Actucal 5-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()


fcox10 <- cph(Surv(time, status)~ riskscore  ,surv=T,x=T, y=T,time.inc = 10,data=gse_data)
cal10 <- calibrate(fcox10, cmethod="KM", method="boot", u=10, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_cali_curve_10y.pdf',width = 6,height = 6)
plot(cal10,col="blue",xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability 10-Year OS",ylab="Actucal 10-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()




###############################
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2) 
library(rms)

GSE181063_clinical_1<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical.txt',header = T,check.names = F,sep = '\t')
GSE181063_clinical_multicox<-merge(GSE181063_lar_score_final,GSE181063_clinical_1,by='sample')




GSE181063_clinical_multicox$age<-ifelse(GSE181063_clinical_multicox$age <=60 ,"<=60",">60")
GSE181063_clinical_multicox$age<-factor(GSE181063_clinical_multicox$age,levels = c("<=60",">60"))
GSE181063_clinical_multicox$age=relevel(GSE181063_clinical_multicox$age, ref = "<=60")



GSE181063_clinical_multicox$ipiScore<-ifelse(GSE181063_clinical_multicox$ipiScore == "ipi_score=0" | GSE181063_clinical_multicox$ipiScore == "ipi_score=1" | GSE181063_clinical_multicox$ipiScore == "ipi_score=2","ipi_score_0+1+2","ipi_score_3+4+5")
GSE181063_clinical_multicox$ipiScore<-factor(GSE181063_clinical_multicox$ipiScore,levels = c("ipi_score_0+1+2","ipi_score_3+4+5"))
GSE181063_clinical_multicox$ipiScore=relevel(GSE181063_clinical_multicox$ipiScore, ref = "ipi_score_0+1+2")


GSE181063_clinical_multicox$stage<-ifelse(GSE181063_clinical_multicox$stage == "Stage I" | GSE181063_clinical_multicox$stage == "Stage II","I+II","III+IV")
GSE181063_clinical_multicox$stage<-factor(GSE181063_clinical_multicox$stage,levels = c("I+II","III+IV"))
GSE181063_clinical_multicox$stage=relevel(GSE181063_clinical_multicox$stage, ref = "I+II")


GSE181063_clinical_multicox$gender<-factor(GSE181063_clinical_multicox$gender,levels = c("F","M"))
GSE181063_clinical_multicox$gender=relevel(GSE181063_clinical_multicox$gender, ref = "F")

GSE181063_clinical_multicox$riskgroup<-NULL

gse_data<-GSE181063_clinical_multicox[,2:8]
dd=datadist(gse_data)
options(datadist="dd") 
f2 <- psm(Surv(time, status) ~ riskscore  + gender + age + ipiScore + stage, data =  gse_data)
med <- Quantile(f2) # 计算中位生存时间
surv <- Survival(f2) # 构建生存概率函数
nom <- nomogram(f2, fun=list(function(x) surv(1, x),
                             function(x) surv(3, x),
                             function(x) surv(5, x),
                             function(x) surv(10, x)),
                funlabel=c("1-year Survival Probability",
                           "3-year Survival Probability",
                           "5-year Survival Probability",
                           "10-year Survival Probability"))
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_clinical_nomogram_forest.pdf',width = 20,height = 10)
plot(nom, xfrac=.6)
dev.off()


fcox1 <- cph(Surv(time, status)~ riskscore ,surv=T,x=T, y=T,time.inc = 1,data=gse_data)
cal1 <- calibrate(fcox1, cmethod="KM", method="boot", u=1, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_cali_curve_1y.pdf',width = 6,height = 6)
plot(cal1,col="blue",xlim=c(0.5,1),ylim=c(0.5,1),xlab="Predicted Probability 1-Year OS",ylab="Actucal 1-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()


fcox3 <- cph(Surv(time, status)~ riskscore ,surv=T,x=T, y=T,time.inc = 3,data=gse_data)
cal3 <- calibrate(fcox3, cmethod="KM", method="boot", u=3, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_cali_curve_3y.pdf',width = 6,height = 6)
plot(cal3,col="blue",xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability 3-Year OS",ylab="Actucal 3-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()

fcox5 <- cph(Surv(time, status)~ riskscore  ,surv=T,x=T, y=T,time.inc = 5,data=gse_data)
cal5 <- calibrate(fcox5, cmethod="KM", method="boot", u=5, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_cali_curve_5y.pdf',width = 6,height = 6)
plot(cal5,col="blue",xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability 5-Year OS",ylab="Actucal 5-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()


fcox10 <- cph(Surv(time, status)~ riskscore  ,surv=T,x=T, y=T,time.inc = 10,data=gse_data)
cal10 <- calibrate(fcox10, cmethod="KM", method="boot", u=10, m=100, B=500)
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE181063_cali_curve_10y.pdf',width = 6,height = 6)
plot(cal10,col="blue",xlim=c(0,1),ylim=c(0,1),xlab="Predicted Probability 10-Year OS",ylab="Actucal 10-Year OS (proportion)",cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.3)
dev.off()



#########################

GSE10846_clinical_multicox_2<-merge(GSE10846_lar_score_final,clinical_GSE10846,by='sample')
GSE10846_clinical_multicox_2[,6:7]<-NULL
colnames(GSE10846_clinical_multicox_2)[2:3]<-c("status","time")


GSE10846_clinical_multicox_2$stage<-factor(GSE10846_clinical_multicox_2$stage)
levels(GSE10846_clinical_multicox_2$stage)

GSE10846_clinical_multicox_2$ecog<-factor(GSE10846_clinical_multicox_2$ecog)
levels(GSE10846_clinical_multicox_2$ecog)

GSE10846_clinical_multicox_3<-GSE10846_clinical_multicox_2[,c(2,5,8,9)]
GSE10846_clinical_multicox_3[,c(1,2,3,4)]<-lapply(GSE10846_clinical_multicox_3[,c(1,2,3,4)], as.numeric)



rocobj_riskscore<- roc( GSE10846_clinical_multicox_3$status,GSE10846_clinical_multicox_3$riskscore,plot=TRUE, print.auc=TRUE)
rocobj_riskscore$auc
rocobj_stage<- roc( GSE10846_clinical_multicox_3$status,GSE10846_clinical_multicox_3$stage, plot=TRUE, print.auc=TRUE)
rocobj_stage$auc
rocobj_ecog<- roc(GSE10846_clinical_multicox_3$status,GSE10846_clinical_multicox_3$ecog, plot=TRUE, print.auc=TRUE)
rocobj_ecog$auc


pdf('//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_clinical_ROC.pdf',width = 5,height = 5)
ggroc(list(Risk_Score=rocobj_riskscore, Stage=rocobj_stage, ecog= rocobj_ecog))+annotate("text",x=c(0.5,0.5,0.5),y=c(0.05,0.15,0.25),label=c("Risk_Score AUC=0.674","Stage AUC=0.663","Grade AUC=0.667"))
dev.off()
#########################

GSE10846_clinical_multicox_2<-merge(GSE10846_lar_score_final,clinical_GSE10846,by='sample')
GSE10846_clinical_multicox_2[,6:7]<-NULL
colnames(GSE10846_clinical_multicox_2)[2:3]<-c("status","time")


GSE10846_clinical_multicox_2$stage<-factor(GSE10846_clinical_multicox_2$stage)
levels(GSE10846_clinical_multicox_2$stage)

GSE10846_clinical_multicox_2$ecog<-factor(GSE10846_clinical_multicox_2$ecog)
levels(GSE10846_clinical_multicox_2$ecog)

GSE10846_clinical_multicox_3<-GSE10846_clinical_multicox_2[,c(2,5,8,9)]
GSE10846_clinical_multicox_3[,c(1,2,3,4)]<-lapply(GSE10846_clinical_multicox_3[,c(1,2,3,4)], as.numeric)



rocobj_riskscore<- roc( GSE10846_clinical_multicox_3$status,GSE10846_clinical_multicox_3$riskscore,plot=TRUE, print.auc=TRUE)
rocobj_riskscore$auc
rocobj_stage<- roc( GSE10846_clinical_multicox_3$status,GSE10846_clinical_multicox_3$stage, plot=TRUE, print.auc=TRUE)
rocobj_stage$auc
rocobj_ecog<- roc(GSE10846_clinical_multicox_3$status,GSE10846_clinical_multicox_3$ecog, plot=TRUE, print.auc=TRUE)
rocobj_ecog$auc


pdf('//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/GSE10846_clinical_ROC.pdf',width = 5,height = 5)
ggroc(list(Risk_Score=rocobj_riskscore, Stage=rocobj_stage, ecog= rocobj_ecog))+annotate("text",x=c(0.5,0.5,0.5),y=c(0.05,0.15,0.25),label=c("Risk_Score AUC=0.674","Stage AUC=0.663","Grade AUC=0.667"))
dev.off()
###############################

save.image('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/05_LAR_clinical/05_LAR_clinical.RData')




