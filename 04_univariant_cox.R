setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/03_univariant_cox')
#####################################
univariate_Cox <- function(exp = NULL, clinical = NULL) {
  library(tidyverse)
  library(survival)
  library(survminer)
  exp <- as.data.frame(t(exp))
  gene_name <- colnames(exp)
  clinical <- clinical[na.omit(match(rownames(exp), clinical[, 1])), ]
  exp <- exp[na.omit(match(clinical[, 1], rownames(exp))), ]
  
  i <- 1
  a <- tibble()
  b <- tibble()
  result <- data.frame()
  result2 <- sapply(gene_name, function(gene_name) {
    result <- coxph(
      formula = Surv(clinical$time, clinical$status) ~ exp[, gene_name],
      data = exp
    )
    
    result <- summary(result)
    coef <- result$coef[1]
    Hazard_Ratio <- result$coef[2]
    p.value <- result$wald["pvalue"]
    lower_.95 <- result$conf.int[, "lower .95"]
    upper_.95 <- result$conf.int[, "upper .95"]
    logrank_pvalue <- result$sctest["pvalue"]
    wald_pvalue <- result$waldtest["pvalue"]
    Likelihood_pvalue <- result$logtest["pvalue"]
    
    a <- as_tibble(
      cbind(gene_name[i], coef, p.value, Hazard_Ratio, lower_.95, upper_.95, logrank_pvalue, wald_pvalue, Likelihood_pvalue)
    )
    b <- bind_rows(b, a)
    return(b)
  })

  result2 <- as.data.frame(t(result2))
  result2[, 2:ncol(result2)] <- apply(result2 %>% dplyr::select(2:ncol(result2)), 2, as.numeric)
  
  colnames(result2)[1] <- "gene"
  result2 <- result2 %>%
    filter(p.value < 0.05)
  result2$gene <- unlist(result2$gene)
  
  result2 <- result2 %>%
    mutate(HR = paste0(
      round(Hazard_Ratio, 2),
      "(",
      round(lower_.95, 2),
      "-",
      round(upper_.95, 2),
      ")"
    ))
  
  str(result2)
  result2<-result2[order(result2$p.value),]
  return(result2)
}


GSE181063_DLBC_epr_lactic<-GSE181063_DLBC_epr[which(rownames(GSE181063_DLBC_epr) %in% gene$Tags),]



GSE181063_univcox<-univariate_Cox(exp= GSE181063_DLBC_epr_lactic,clinical=clinical_GSE181063_DLBC)

GSE181063_univcox<-GSE181063_univcox[order(GSE181063_univcox$p.value),]
write.table(GSE181063_univcox,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/03_univariant_cox/GSE181063_univcox.txt',row.names=F,col.names=T,quote=F,sep="\t")




unicox_forest_plot<- function(unicox_re=NULL,plot_file_name=NULL)
{
  library(forestplot)
  forest_table<-data.frame(Risk_factors=rownames(unicox_re),HR=unicox_re$HR, pvalue=unicox_re$p.value,check.names = F)
  forest_table1<-data.frame(Risk_factors=NA,HR=NA, pvalue=NA,check.names = F)
  forest_table2<-data.frame(Risk_factors="Risk_factors",HR="HR", pvalue="p-value",check.names = F)
  tabletext<-rbind(forest_table2,forest_table1,forest_table)
  tabletext$sig<-NA
  for(i in 3:dim(tabletext)[1])
  {
    if(as.numeric(tabletext$pvalue[i])< 0.05 & as.numeric(tabletext$pvalue[i]) > 0.01 )
    {
      tabletext$sig[i]<-"*"
      tabletext$pvalue[i]<-round(as.numeric(tabletext$pvalue[i]),6)
    }
    else if(as.numeric(tabletext$pvalue[i])< 0.01 & as.numeric(tabletext$pvalue[i])> 0.001)
    {
      tabletext$sig[i]<-"**"
      tabletext$pvalue[i]<-round(as.numeric(tabletext$pvalue[i]),6)
    }
    else if (as.numeric(tabletext$pvalue[i]) < 0.001)
    {
      tabletext$sig[i]<-"***"
      tabletext$pvalue[i]<-round(as.numeric(tabletext$pvalue[i]),6)
      if( tabletext$pvalue[i] == 0)
      {
        tabletext$pvalue[i]<- "<1e-06"
      }
    }
    else
    {
      tabletext$sig[i]<-" "
      tabletext$pvalue[i]<-round(as.numeric(tabletext$pvalue[i]),6)
    }
  }
  tabletext$sig[1]<-"Significant"
  forest_stastic<-data.frame(mean=as.numeric(unicox_re$Hazard_Ratio),lower=as.numeric(unicox_re$lower_.95),upper=as.numeric(unicox_re$upper_.95))
  forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
  cochrane_from_rmeta<-rbind(forest_stastic1,forest_stastic)
  

  pdf(plot_file_name,width = 8,height = dim(tabletext)[1]*0.2)
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
  
  result=list(tabletext=tabletext,forest_stastic=forest_stastic)
  return(result)
}


GSE181063_univcox_plot<-unicox_forest_plot(unicox_re=GSE181063_univcox,plot_file_name='/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/03_univariant_cox/clinical_single_forest.pdf')


#####################################

GSE181063_univcox_hr_over1<-GSE181063_univcox[which(GSE181063_univcox$Hazard_Ratio >1),]
GSE181063_univcox_hr_less1<-GSE181063_univcox[which(GSE181063_univcox$Hazard_Ratio <1),]

########### HR over 1 #####
GSE181063_univcox_hr_over1_top_10<-GSE181063_univcox_hr_over1[1:6,]


GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1<-GSE181063_DLBC_epr_lactic[which(rownames(GSE181063_DLBC_epr_lactic) %in% GSE181063_univcox_hr_over1_top_10$gene),]



GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1<-as.data.frame(t(GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1))

GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1$sample<-rownames(GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1)

GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical<-merge(clinical_GSE181063_DLBC,GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1,by='sample')

colnames(GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical)





library("survival")
library("survminer") 
library("survivalROC")
library("timeROC")
library(ggpubr)

j<-1
for(i in 1:dim(GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1)[2])
{
  gene_name<-colnames(GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical)[i+8]
  
  temp<-data.frame(time=GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical$time,status=GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical$status,gene=GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical[,i+8])
  #colnames(temp)[3]<-gene_name
  
  
  temp$gene<-ifelse(temp$gene <= median(temp$gene),"Low","High")
  
  fit_validation <- survfit( Surv(time, status) ~ gene,data = temp )
  
  p<-ggsurvplot(fit_validation, data = temp,
                conf.int = TRUE,
                pval = TRUE,
                fun = "pct",
                risk.table = TRUE,
                size = 1,
                linetype = "strata",
                palette = c("red","#2E9FDF"),
                legend = "bottom",
                legend.title = paste0(gene_name," expression"),
                legend.labs = c("High","Low"),
                surv.median.line = "hv",
                xlab = "Time (Years)",
                tables.height = 0.1)
  
  q1 <- ggarrange(p$plot, p$table, ncol = 1, align = 'v',heights = c(0.75, 0.3),legend="top")
  
  
  temp1<-data.frame(time=GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical$time,status=GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical$status,gene=GSE181063_DLBC_epr_lactic_unicox_top5_hr_over1_clinical[,i+8])
  temp1$expression<-ifelse(temp1$gene <= median(temp1$gene),"Low","High")
  
  my_comparisons<-list(c("Low","High"))
  
  q2<- ggboxplot(temp1,x="expression" , y="gene" ,fill= "expression", palette = c("#2E9FDF","red"),scales = "free_x",xlab=F,ylab=paste0(gene_name," expression (log2(TPM+1))"),width = 0.5,order=c("Low","High"))+
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position='top')+theme(legend.text = element_text( size = 8))+guides(fill=guide_legend(title=NULL))

  
  s<-ggarrange(q2, q1,ncol = 2,nrow = 1,widths=c(2,4))
  
  assign(paste0("s", j),s)
  j<-j+1
  
}

p131 <- cowplot::plot_grid(s1,s2,s3,s4,s5,s6,ncol = 3,nrow = 2,labels = 'AUTO',rel_heights = c(1,1,1),label_size = 20,scale = c(0.9,0.9,0.9,0.9,0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/03_univariant_cox/gene_expression_survival_HR_over1.pdf' ,
                 plot = p131,width = 21,height = 14)


########### HR less 1 #####


GSE181063_univcox_hr_less1_top_10<-GSE181063_univcox_hr_less1[1:4,]


GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1<-GSE181063_DLBC_epr_lactic[which(rownames(GSE181063_DLBC_epr_lactic) %in% GSE181063_univcox_hr_less1_top_10$gene),]



GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1<-as.data.frame(t(GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1))

GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1$sample<-rownames(GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1)

GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical<-merge(clinical_GSE181063_DLBC,GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1,by='sample')

colnames(GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical)





library("survival")
library("survminer") 
library("survivalROC")
library("timeROC")
library(ggpubr)

j<-1
for(i in 1:dim(GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1)[2])
{
  gene_name<-colnames(GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical)[i+8]
  
  temp<-data.frame(time=GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical$time,status=GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical$status,gene=GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical[,i+8])
  #colnames(temp)[3]<-gene_name
  
  
  temp$gene<-ifelse(temp$gene <= median(temp$gene),"Low","High")
  
  fit_validation <- survfit( Surv(time, status) ~ gene,data = temp )
  
  p<-ggsurvplot(fit_validation, data = temp,
                conf.int = TRUE,
                pval = TRUE,
                fun = "pct",
                risk.table = TRUE,
                size = 1,
                linetype = "strata",
                palette = c("red","#2E9FDF"),
                legend = "bottom",
                legend.title = paste0(gene_name," expression"),
                legend.labs = c("High","Low"),
                surv.median.line = "hv",
                xlab = "Time (Years)",
                tables.height = 0.1)
  
  q1 <- ggarrange(p$plot, p$table, ncol = 1, align = 'v',heights = c(0.75, 0.3),legend="top")
  
  
  temp1<-data.frame(time=GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical$time,status=GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical$status,gene=GSE181063_DLBC_epr_lactic_unicox_top5_hr_less1_clinical[,i+8])
  temp1$expression<-ifelse(temp1$gene <= median(temp1$gene),"Low","High")
  
  my_comparisons<-list(c("Low","High"))
  
  q2<- ggboxplot(temp1,x="expression" , y="gene" ,fill= "expression", palette = c("#2E9FDF","red"),scales = "free_x",xlab=F,ylab=paste0(gene_name," expression (log2(TPM+1))"),width = 0.5,order=c("Low","High"))+
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position='top')+theme(legend.text = element_text( size = 8))+guides(fill=guide_legend(title=NULL))
  
  s<-ggarrange(q2, q1,ncol = 2,nrow = 1,widths=c(2,4))
  
  assign(paste0("s", j),s)
  j<-j+1
  
}

p131 <- cowplot::plot_grid(s1,s2,s3,s4,ncol = 3,nrow = 2,labels = 'AUTO',rel_heights = c(1,1,1),label_size = 20,scale = c(0.9,0.9,0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/03_univariant_cox/gene_expression_survival_hr_less1.pdf' ,
                 plot = p131,width = 21,height = 14)








#####################################

save.image('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/03_univariant_cox/03_univariant_cox.RData')











