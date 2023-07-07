setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular')
library(tidyr)
#############################  cibersort ##############
source('/Volumes/Work/genome/supportFunc_cibersort.R')


GSE181063_result <- CIBERSORT('/Volumes/Work/genome/LM22.txt','/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE181063_DLBC_epr.txt', perm = 1000, QN = TRUE)


GSE10846_result <- CIBERSORT('/Volumes/Work/genome/LM22.txt','/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE10846_epr.txt', perm = 1000, QN = TRUE)


GSE181063_ciber<-GSE181063_result[,1:22]
GSE10846_ciber<-GSE10846_result[,1:22]



GSE181063_lar_score_final

GSE181063_ciber_plot<-as.data.frame(GSE181063_ciber)
GSE181063_ciber_plot$sample<-rownames(GSE181063_ciber_plot)
GSE181063_ciber_plot_final<-gather(GSE181063_ciber_plot,cell_type,cell_proportion,`B cells naive`:Neutrophils)
GSE181063_ciber_plot_final$LARgroup<-ifelse(GSE181063_ciber_plot_final$sample %in% GSE181063_lar_score_final[which(GSE181063_lar_score_final$riskgroup == "High"),]$sample,"High","Low")

GSE10846_lar_score_final

GSE10846_ciber_plot<-as.data.frame(GSE10846_ciber)
GSE10846_ciber_plot$sample<-rownames(GSE10846_ciber_plot)
GSE10846_ciber_plot_final<-gather(GSE10846_ciber_plot,cell_type,cell_proportion,`B cells naive`:Neutrophils)
GSE10846_ciber_plot_final$LARgroup<-ifelse(GSE10846_ciber_plot_final$sample %in% GSE10846_lar_score_final[which(GSE10846_lar_score_final$riskgroup == "High"),]$sample,"High","Low")


p1<- ggboxplot(GSE181063_ciber_plot_final, x="cell_type", y="cell_proportion", fill = "LARgroup", 
               palette = "jama",short.panel.labs = T,xlab=F)+stat_compare_means(aes(group=LARgroup), label = "p.signif", method = "wilcox.test")+theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle("GSE181063")


p2<-  ggboxplot(GSE10846_ciber_plot_final, x="cell_type", y="cell_proportion", fill = "LARgroup", 
                palette = "jama",short.panel.labs = T,xlab=F)+stat_compare_means(aes(group=LARgroup), label = "p.signif", method = "wilcox.test")+theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle("GSE10846")


p131 <- cowplot::plot_grid(p1,p2,ncol = 1,nrow = 2,labels = 'AUTO',rel_heights = c(1,1),label_size = 20,scale = c(0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/LARrisk_ciber.pdf' ,
                 plot = p131,width = 15,height = 10)





#############################  xCell ##############
library(xCell)


GSE181063_xCell  = xCellAnalysis(GSE181063_DLBC_epr, rnaseq = F)
GSE10846_xCell  = xCellAnalysis(GSE10846_epr, rnaseq = F)



GSE181063_xCell_results<-as.data.frame(t(GSE181063_xCell[1:64,]))
GSE10846_xCell_results<-as.data.frame(t(GSE10846_xCell[1:64,]))



GSE181063_xcell_plot<-as.data.frame(GSE181063_xCell_results)
GSE181063_xcell_plot$sample<-rownames(GSE181063_xcell_plot)
GSE181063_xcell_plot_final<-gather(GSE181063_xcell_plot,cell_type,xCell,aDC:Tregs)
GSE181063_xcell_plot_final$LARgroup<-ifelse(GSE181063_xcell_plot_final$sample %in% GSE181063_lar_score_final[which(GSE181063_lar_score_final$riskgroup == "High"),]$sample,"High","Low")

GSE10846_lar_score_final

GSE10846_xcell_plot<-as.data.frame(GSE10846_xCell_results)
GSE10846_xcell_plot$sample<-rownames(GSE10846_xcell_plot)
GSE10846_xcell_plot_final<-gather(GSE10846_xcell_plot,cell_type,xCell,aDC:Tregs)
GSE10846_xcell_plot_final$LARgroup<-ifelse(GSE10846_xcell_plot_final$sample %in% GSE10846_lar_score_final[which(GSE10846_lar_score_final$riskgroup == "High"),]$sample,"High","Low")


p3<- ggboxplot(GSE181063_xcell_plot_final, x="cell_type", y="xCell", fill = "LARgroup", 
               palette = "jama",short.panel.labs = T,xlab=F)+stat_compare_means(aes(group=LARgroup), label = "p.signif", method = "wilcox.test")+theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle("GSE181063")


p4<-  ggboxplot(GSE10846_xcell_plot_final, x="cell_type", y="xCell", fill = "LARgroup", 
                palette = "jama",short.panel.labs = T,xlab=F)+stat_compare_means(aes(group=LARgroup), label = "p.signif", method = "wilcox.test")+theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle("GSE10846")


p131 <- cowplot::plot_grid(p3,p4,ncol = 1,nrow = 2,labels = 'AUTO',rel_heights = c(1,1),label_size = 20,scale = c(0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/LARrisk_xcell.pdf' ,
                 plot = p131,width = 25,height = 15)




#############################  ssGSEA  ##############

immune_gene<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/PMID28052254-ssGSEA 细胞对应的基因集.tsv',header = T,check.names = F,sep = '\t')

immune_gene_list<- split(as.matrix(immune_gene)[,1], immune_gene[,2])


library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)

GSE181063_ssGSEA<- gsva(as.matrix(GSE181063_DLBC_epr), immune_gene_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

GSE10846_ssGSEA<- gsva(as.matrix(GSE10846_epr), immune_gene_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)



GSE181063_ssGSEA_results<-as.data.frame(t(GSE181063_ssGSEA))
GSE10846_ssGSEA_results<-as.data.frame(t(GSE10846_ssGSEA))


GSE181063_ssGSEA_plot<-as.data.frame(GSE181063_ssGSEA_results)
GSE181063_ssGSEA_plot$sample<-rownames(GSE181063_ssGSEA_plot)
GSE181063_ssGSEA_plot_final<-gather(GSE181063_ssGSEA_plot,cell_type,ssGSEA,`B cells memory`:`T cells regulatory (Tregs)`)
GSE181063_ssGSEA_plot_final$LARgroup<-ifelse(GSE181063_ssGSEA_plot_final$sample %in% GSE181063_lar_score_final[which(GSE181063_lar_score_final$riskgroup == "High"),]$sample,"High","Low")

GSE10846_lar_score_final

GSE10846_ssGSEA_plot<-as.data.frame(GSE10846_ssGSEA_results)
GSE10846_ssGSEA_plot$sample<-rownames(GSE10846_ssGSEA_plot)
GSE10846_ssGSEA_plot_final<-gather(GSE10846_ssGSEA_plot,cell_type,ssGSEA,`B cells memory`:`T cells regulatory (Tregs)`)
GSE10846_ssGSEA_plot_final$LARgroup<-ifelse(GSE10846_ssGSEA_plot_final$sample %in% GSE10846_lar_score_final[which(GSE10846_lar_score_final$riskgroup == "High"),]$sample,"High","Low")


p5<- ggboxplot(GSE181063_ssGSEA_plot_final, x="cell_type", y="ssGSEA", fill = "LARgroup", 
               palette = "jama",short.panel.labs = T,xlab=F)+stat_compare_means(aes(group=LARgroup), label = "p.signif", method = "wilcox.test")+theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle("GSE181063")


p6<-  ggboxplot(GSE10846_ssGSEA_plot_final, x="cell_type", y="ssGSEA", fill = "LARgroup", 
                palette = "jama",short.panel.labs = T,xlab=F)+stat_compare_means(aes(group=LARgroup), label = "p.signif", method = "wilcox.test")+theme(axis.text.x = element_text(angle = 30, hjust = 1))+ggtitle("GSE10846")


p131 <- cowplot::plot_grid(p5,p6,ncol = 1,nrow = 2,labels = 'AUTO',rel_heights = c(1,1),label_size = 20,scale = c(0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/LARrisk_ssGSEA.pdf' ,
                 plot = p131,width = 15,height = 10)







#############################  ICP  ##############

ICP_gene<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/immue_chekpoint_PMID32814346.txt')

genename_all <- as.character(rownames(GSE181063_DLBC_epr)) 

gene_map_all <- select(org.Hs.eg.db, keys=genename_all, keytype="SYMBOL", columns=c("ENTREZID"))

GSE181063_DLBC_epr_use<-GSE181063_DLBC_epr
GSE181063_DLBC_epr_use$SYMBOL<-rownames(GSE181063_DLBC_epr_use)

GSE181063_DLBC_epr_use_entzid<-base::merge(gene_map_all,GSE181063_DLBC_epr_use,by='SYMBOL')

ICP_gene_1<-as.character(ICP_gene$V1)
gene_map_icp <- select(org.Hs.eg.db, keys=ICP_gene_1, keytype="SYMBOL", columns=c("ENTREZID"))

gene_map_icp<-na.omit(gene_map_icp)

GSE181063_ICP<-GSE181063_DLBC_epr_use_entzid[which(GSE181063_DLBC_epr_use_entzid$ENTREZID %in% gene_map_icp$ENTREZID),]

rownames(GSE181063_ICP)<-GSE181063_ICP$SYMBOL
GSE181063_ICP$ENTREZID<-NULL
GSE181063_ICP$SYMBOL<-NULL


GSE181063_ICP<-as.data.frame(t(GSE181063_ICP))
GSE181063_ICP$sample<-rownames(GSE181063_ICP)

GSE181063_ICP<-base::merge(GSE181063_ICP,GSE181063_lar_score_final,by='sample')
rownames(GSE181063_ICP)<-GSE181063_ICP$sample

colnames(GSE181063_ICP)

wilcox_p_ICP<-data.frame(matrix(data = 1,nrow = 58,ncol=3))
colnames(wilcox_p_ICP)<-c("gene",'pvalue','type')

for(i in 2:59)
{
  temp<-GSE181063_ICP[,c(i,62)]
  temp1<-temp[which(temp$riskgroup == 'Low'),]
  temp2<-temp[which(temp$riskgroup == 'High'),]
  re<-wilcox.test(temp1[,1],temp2[,1])
  wilcox_p_ICP[i-1,1]<-colnames(temp)[1]
  wilcox_p_ICP[i-1,2]<-re$p.value
  delta<-median(temp1[,1])-median(temp2[,1])
  if(delta <0)
  {
    wilcox_p_ICP[i-1,3]<-"High over Low"
  }
  if(delta > 0)
  {
    wilcox_p_ICP[i-1,3]<-"Low over High"
  }
}

write.table(wilcox_p_ICP,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/GSE181063_ICP_MPIrisk_group_wilcox_test.txt',row.names=F,col.names=TRUE,quote=FALSE,sep="\t")

wilcox_p_ICP_sig<-wilcox_p_ICP[which(wilcox_p_ICP$pvalue <0.01),]

GSE181063_ICP_sig<-GSE181063_ICP[,which(colnames(GSE181063_ICP) %in% wilcox_p_ICP_sig$gene)]

GSE181063_ICP_sig$sample<-rownames(GSE181063_ICP_sig)
GSE181063_ICP_sig<-base::merge(GSE181063_ICP_sig,GSE181063_lar_score_final,by='sample')


colnames(GSE181063_ICP_sig)[47]<-'LAR Group'



table(wilcox_p_ICP_sig$type)

low_ocer_high_gene<-wilcox_p_ICP_sig$gene[which(wilcox_p_ICP_sig$type == "Low over High")]#37
high_ocer_low_gene<-wilcox_p_ICP_sig$gene[which(wilcox_p_ICP_sig$type == "High over Low")]#6


j<-1
m<-1
for(i in 2:44)
{
  temp<-GSE181063_ICP_sig[,c(i,47)]
  type<-colnames(temp)[1]
  if(type %in% low_ocer_high_gene)
  {
    my_comparisons<-list(c("High","Low"))
    p<-ggboxplot(temp,x="LAR Group" , y=type ,fill= "LAR Group", palette = c("lightgreen","pink"),scales = "free_x",xlab=F,ylab="Expression",width = 0.5,order=c("Low","High"))+
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+ggtitle(type)
    
    assign(paste0("p", j),p)
    j<-j+1
  }
  if(type %in% high_ocer_low_gene)
  {
    my_comparisons<-list(c("High","Low"))
    q<-ggboxplot(temp,x="LAR Group" , y=type ,fill= "LAR Group", palette = c("lightgreen","pink"),scales = "free_x",xlab=F,ylab="Expression",width = 0.5,order=c("Low","High"))+
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+ggtitle(type)
    
    assign(paste0("q", m),q)
    m<-m+1
  }
}

p131 <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,ncol = 7,nrow =6,labels = NULL,rel_heights = c(1,1,1),label_size = 20,scale = c(0.9,0.9),align='hv')


p132 <- cowplot::plot_grid(q1,q2,q3,q4,q5,q6,ncol = 6,nrow =1,labels = NULL,rel_heights = c(1,1,1),label_size = 20,scale = c(1,1,1,1,1,1,1,1),align='hv')


p133 <- cowplot::plot_grid(p131,p132,ncol = 1,nrow =2,labels = 'AUTO',rel_heights = c(7,1),label_size = 20,scale = c(0.9,0.9),align='hv')


cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/GSE181063_LARrisk_ICP.pdf' ,
                 plot = p133,width = 20,height = 35)


#################################
genename_all <- as.character(rownames(GSE10846_epr)) 

gene_map_all <- select(org.Hs.eg.db, keys=genename_all, keytype="SYMBOL", columns=c("ENTREZID"))

GSE10846_DLBC_epr_use<-GSE10846_epr
GSE10846_DLBC_epr_use$SYMBOL<-rownames(GSE10846_DLBC_epr_use)

GSE10846_DLBC_epr_use_entzid<-base::merge(gene_map_all,GSE10846_DLBC_epr_use,by='SYMBOL')


GSE10846_ICP<-GSE10846_DLBC_epr_use_entzid[which(GSE10846_DLBC_epr_use_entzid$ENTREZID %in% gene_map_icp$ENTREZID),]

rownames(GSE10846_ICP)<-GSE10846_ICP$SYMBOL
GSE10846_ICP$ENTREZID<-NULL
GSE10846_ICP$SYMBOL<-NULL

#######################################

GSE10846_ICP<-as.data.frame(t(GSE10846_ICP))
GSE10846_ICP$sample<-rownames(GSE10846_ICP)

GSE10846_ICP<-base::merge(GSE10846_ICP,GSE10846_lar_score_final,by='sample')
rownames(GSE10846_ICP)<-GSE10846_ICP$sample

colnames(GSE10846_ICP)

wilcox_p_ICP<-data.frame(matrix(data = 1,nrow = 56,ncol=3))
colnames(wilcox_p_ICP)<-c("gene",'pvalue','type')

for(i in 2:57)
{
  temp<-GSE10846_ICP[,c(i,60)]
  temp1<-temp[which(temp$riskgroup == 'Low'),]
  temp2<-temp[which(temp$riskgroup == 'High'),]
  re<-wilcox.test(temp1[,1],temp2[,1])
  wilcox_p_ICP[i-1,1]<-colnames(temp)[1]
  wilcox_p_ICP[i-1,2]<-re$p.value
  delta<-median(temp1[,1])-median(temp2[,1])
  if(delta <0)
  {
    wilcox_p_ICP[i-1,3]<-"High over Low"
  }
  if(delta > 0)
  {
    wilcox_p_ICP[i-1,3]<-"Low over High"
  }
}

write.table(wilcox_p_ICP,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/GSE10846_ICP_MPIrisk_group_wilcox_test.txt',row.names=F,col.names=TRUE,quote=FALSE,sep="\t")

wilcox_p_ICP_sig<-wilcox_p_ICP[which(wilcox_p_ICP$pvalue <0.01),]

GSE10846_ICP_sig<-GSE10846_ICP[,which(colnames(GSE10846_ICP) %in% wilcox_p_ICP_sig$gene)]

GSE10846_ICP_sig$sample<-rownames(GSE10846_ICP_sig)
GSE10846_ICP_sig<-base::merge(GSE10846_ICP_sig,GSE10846_lar_score_final,by='sample')


colnames(GSE10846_ICP_sig)[31]<-'LAR Group'



table(wilcox_p_ICP_sig$type)

low_ocer_high_gene<-wilcox_p_ICP_sig$gene[which(wilcox_p_ICP_sig$type == "Low over High")]#20
high_ocer_low_gene<-wilcox_p_ICP_sig$gene[which(wilcox_p_ICP_sig$type == "High over Low")]#7


j<-1
m<-1
for(i in 2:28)
{
  temp<-GSE10846_ICP_sig[,c(i,31)]
  type<-colnames(temp)[1]
  if(type %in% low_ocer_high_gene)
  {
    my_comparisons<-list(c("High","Low"))
    p<-ggboxplot(temp,x="LAR Group" , y=type ,fill= "LAR Group", palette = c("lightgreen","pink"),scales = "free_x",xlab=F,ylab="Expression",width = 0.5,order=c("Low","High"))+
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+ggtitle(type)
    
    assign(paste0("p", j),p)
    j<-j+1
  }
  if(type %in% high_ocer_low_gene)
  {
    my_comparisons<-list(c("High","Low"))
    q<-ggboxplot(temp,x="LAR Group" , y=type ,fill= "LAR Group", palette = c("lightgreen","pink"),scales = "free_x",xlab=F,ylab="Expression",width = 0.5,order=c("Low","High"))+
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
      theme(plot.title = element_text(hjust = 0.5))+ggtitle(type)
    
    assign(paste0("q", m),q)
    m<-m+1
  }
}

p131 <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,ncol = 5,nrow =4,labels = NULL,rel_heights = c(1,1,1),label_size = 20,scale = c(0.9,0.9),align='hv')


p132 <- cowplot::plot_grid(q1,q2,q3,q4,q5,q6,ncol = 5,nrow =2,labels = NULL,rel_heights = c(1,1,1),label_size = 20,scale = c(0.9,0.9),align='hv')


p133 <- cowplot::plot_grid(p131,p132,ncol = 1,nrow =2,labels = 'AUTO',rel_heights = c(4,2),label_size = 20,scale = c(0.9,0.9),align='hv')


cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/GSE10846_LARrisk_ICP.pdf' ,
                 plot = p133,width = 15,height = 30)




################################# hallmark #######
t<-read.gmt("/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/h.all.v7.4.symbols.gmt") #读gmt文件

hm_gene_list<- split(as.matrix(hm_gmt)[,2], hm_gmt[,1])



library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)

GSE181063_ssGSEA_hm<- gsva(as.matrix(GSE181063_DLBC_epr), hm_gene_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

GSE10846_ssGSEA_hm<- gsva(as.matrix(GSE10846_epr), hm_gene_list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

GSE181063_ssGSEA_hm_results<-as.data.frame(t(GSE181063_ssGSEA_hm))
GSE10846_ssGSEA_hm_results<-as.data.frame(t(GSE10846_ssGSEA_hm))

GSE181063_ssGSEA_hm_results$sample<-rownames(GSE181063_ssGSEA_hm_results)
GSE10846_ssGSEA_hm_results$sample<-rownames(GSE10846_ssGSEA_hm_results)

GSE181063_hm_results<-merge(GSE181063_ssGSEA_hm_results,GSE181063_lar_score_final,by="sample")
GSE10846_hm_results<-merge(GSE10846_ssGSEA_hm_results,GSE10846_lar_score_final,by="sample")

GSE181063_hm_results<-GSE181063_hm_results[,match(colnames(GSE10846_hm_results),colnames(GSE181063_hm_results))]


###
corr_p<-data.frame(matrix(data=1,ncol = 50,nrow = 2))
colnames(corr_p)<-colnames(GSE181063_hm_results)[2:51]
rownames(corr_p)<-c("GSE181063","GSE10846")


corr_corr<-data.frame(matrix(data=0,ncol = 50,nrow = 2))
colnames(corr_corr)<-colnames(GSE181063_hm_results)[2:51]

rownames(corr_corr)<-c("GSE181063","GSE10846")

for(i in 2:51)
{
  M=cor.test(GSE181063_hm_results[,i],GSE181063_hm_results[,55],method = "spearman")
  if (M$p.value ==0)
  {
    corr_p[1,i-1]<-0.000000000000000001
  }
  else
  {
    corr_p[1,i-1]<-M$p.value
  }
  corr_corr[1,i-1]<-M$estimate
  M1=cor.test(GSE10846_hm_results[,i],GSE10846_hm_results[,55],method = "spearman")
  if (M1$p.value ==0)
  {
    corr_p[2,i-1]<-0.000000000000000001
  }
  else
  {
    corr_p[2,i-1]<-M1$p.value
  }
  corr_corr[2,i-1]<-M1$estimate
 
}



library(lattice)
library("circlize")
library("ComplexHeatmap")


plotMutiHeatmap2=function(up,down,up_break,up_colors,down_break,down_colors,title){
  UpColor <- colorRamp2(breaks = up_break, colors = up_colors)
  DnColor <- colorRamp2(breaks = down_break, colors = down_colors)
  
  
  
  DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                   unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                   unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
      if(up[i, j]>=1.3){
        txt="***"
        if(up[i, j]>=1.3&up[i, j]<2){
          txt='*'
        }else if(up[i, j]>=2&up[i, j]<3){
          txt='**'
        }else if(up[i, j]>=3&up[i, j]<4){
          txt='***'
        }
        grid.text(label=txt,x=(x + 0.5*width),
                  y=(y+ 0.5*height),just = c('right','top'))
      }
      if(down[i, j]>0){
      }
    }
  }
  
  p1 <- Heatmap(up, column_title = title
                , rect_gp = gpar(type = "none")
                , show_heatmap_legend = F
                , cluster_rows = F
                , cluster_columns = F, 
                #top_annotation = row_an, 
                #width=unit(20, "cm"),
                cell_fun = DiagFunc(up = up, down = down) 
  ) 
  col_fun = colorRamp2(down_break, down_colors) 
  lgd <- Legend(title = "log2FC", 
                col_fun = col_fun, 
                at = c(-0.4,0,0.4), 
                labels = c("-0.4","0","0.4"),  
                direction = "horizontal" 
  )
  col_fun2 = colorRamp2(up_break, up_colors) 
  lgd2 <- Legend(title = "-log10(p value)", 
                 col_fun = col_fun2, 
                 at = c(0,1,2,3,4,5), 
                 labels = c('0',"1","2","3","4",">5"),  
                 direction = "horizontal"
  )
  
  draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom"
       ,heatmap_legend_side = "bottom", merge_legend = TRUE)
}

up_break=c(0, 5)
down_break=c(-0.4, 0,0.4)
up_colors=c("#FFFFFF","#6f9a8d")
down_colors=c("blue",'white',"red")


up=-log10(corr_p)
down=corr_corr
title='Correlation between 50 hallmark ssGSEA Score and LARscore'
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/heatmap_corr.pdf',width = 20,height = 4)
plotMutiHeatmap2(up,down,up_break,up_colors,down_break,down_colors,title)
dev.off()



################################# estimate ####
library(estimate)
filterCommonGenes(input.f="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE181063_DLBC_epr.txt", 
                  output.f="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE181063_finaldata.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE181063_finaldata.gct",
              output.ds="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE181063_finaldata_estimate_score.gct", 
              platform="illumina")
GSE181063_scores=read.table("/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE181063_finaldata_estimate_score.gct",skip = 2,header = T)
rownames(GSE181063_scores)=GSE181063_scores[,1]
GSE181063_scores=t(GSE181063_scores[,3:ncol(GSE181063_scores)])
GSE181063_scores
write.table(GSE181063_scores,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/GSE181063_estimat_scorese.txt',row.names=T,col.names=T,quote=F,sep="\t")




filterCommonGenes(input.f="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE10846_epr.txt", 
                  output.f="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE10846_finaldata.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE10846_finaldata.gct",
              output.ds="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE10846_finaldata_estimate_score.gct", 
              platform="illumina")
GSE10846_scores=read.table("/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/estimate/GSE10846_finaldata_estimate_score.gct",skip = 2,header = T)
rownames(GSE10846_scores)=GSE10846_scores[,1]
GSE10846_scores=t(GSE10846_scores[,3:ncol(GSE10846_scores)])
GSE10846_scores
write.table(GSE10846_scores,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/GSE10846_estimat_scorese.txt',row.names=T,col.names=T,quote=F,sep="\t")

GSE181063_scores<-as.data.frame(GSE181063_scores)

GSE10846_scores<-as.data.frame(GSE10846_scores)


GSE181063_scores$sample<-rownames(GSE181063_scores)
GSE10846_scores$sample<-rownames(GSE10846_scores)




GSE181063_estimate<-merge(GSE181063_scores,GSE181063_lar_score_final,by="sample")
GSE10846_estimate<-merge(GSE10846_scores,GSE10846_lar_score_final,by="sample")

library(ggExtra)
colnames(GSE181063_estimate)[7]<-'LARrisk Group'
q1<-ggscatter(GSE181063_estimate,x = "riskscore", #x变量
          y = "ImmuneScore",#y变量
          add = "reg.line",##拟合曲线
          #conf.int = TRUE,##置信区间阴影带
          cor.coef = TRUE, ##系数
          color =  "LARrisk Group",
          palette = c("pink","lightgreen"),
          add.params = list(color = "black",fill = "lightgray"),
          cor.method = "spearman",
          xlab = "LAR score", ## x轴
          ylab = 'ImmuneScore',
          title="GSE181063")+ theme(legend.position="bottom")

q4<-ggMarginal(q1, type = "boxplot",groupColour = TRUE, groupFill = TRUE,margins="x",size=10)

my_comparisons<-list(c("Low","High"))
q2<- ggboxplot(GSE181063_estimate,x="LARrisk Group" , y="ImmuneScore" ,fill= "LARrisk Group", palette = c("lightgreen","pink"),scales = "free_x",xlab=F,ylab=F,width = 0.5,order = c("Low","High"))+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+scale_x_discrete(breaks = NULL)+scale_y_continuous(breaks = NULL)+ theme(legend.position="bottom")

p1 <- cowplot::insert_yaxis_grob(q4, q2, grid::unit(.2, "null"), position = "right")






colnames(GSE10846_estimate)[7]<-'LARrisk Group'
q2<-ggscatter(GSE10846_estimate,x = "riskscore", #x变量
              y = "ImmuneScore",#y变量
              add = "reg.line",##拟合曲线
              #conf.int = TRUE,##置信区间阴影带
              cor.coef = TRUE, ##系数
              color =  "LARrisk Group",
              palette = c("pink","lightgreen"),
              add.params = list(color = "black",fill = "lightgray"),
              cor.method = "spearman",
              xlab = "LAR score", ## x轴
              ylab = 'ImmuneScore',
              title="GSE10846")+ theme(legend.position="bottom")

q5<-ggMarginal(q2, type = "boxplot",groupColour = TRUE, groupFill = TRUE,margins="x",size=10)

my_comparisons<-list(c("Low","High"))
q3<- ggboxplot(GSE10846_estimate,x="LARrisk Group" , y="ImmuneScore" ,fill= "LARrisk Group", palette = c("lightgreen","pink"),scales = "free_x",xlab=F,ylab=F,width = 0.5,order = c("Low","High"))+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position='none')+scale_x_discrete(breaks = NULL)+scale_y_continuous(breaks = NULL)+ theme(legend.position="bottom")

p2 <- cowplot::insert_yaxis_grob(q5, q3, grid::unit(.2, "null"), position = "right")


p131 <- cowplot::plot_grid(p1,p2,ncol = 2,nrow = 1,labels = 'AUTO',rel_heights = c(1,1),label_size = 20,scale = c(0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/immu_larscore_correlation.pdf' ,
                 plot = p131,width = 10,height = 5)



#################################


save.image('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/06_LAR_molecular/06_LAR_molecular.RData')



#################################



















