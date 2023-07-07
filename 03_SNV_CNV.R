setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV')


library(maftools)


DLBC<-read.maf("/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/TCGA/DLBC.maf")

DLBC_data<-DLBC@data

table(DLBC_data$Tumor_Sample_Barcode)


DLBC_data$Tumor_Sample_Barcode<-substr(DLBC_data$Tumor_Sample_Barcode,1,12)


DLBC_data_lactic<-DLBC_data[which(DLBC_data$Hugo_Symbol %in% gene$Tags),]

lactic_mutant_sample<-DLBC_data_lactic$Tumor_Sample_Barcode[!duplicated(DLBC_data_lactic$Tumor_Sample_Barcode)]




#################################
TCGA_clinical<-read.table("/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/TCGA/TCGA_DLBC_clinical.txt",header = T,check.names = F,sep = '\t')

DLBC_LACTIC<-read.maf(maf=DLBC_data_lactic,isTCGA=T,clinicalData = TCGA_clinical)

table(DLBC_data_lactic$Hugo_Symbol) #30 gene

mut_gene_lactic<-DLBC_data_lactic$Hugo_Symbol[!duplicated(DLBC_data_lactic$Hugo_Symbol)]


oncoplot(maf = DLBC_LACTIC, top = 50, fontSize = 0.8 ,showTumorSampleBarcodes = F )


for( abbr in mut_gene_lactic)
{
  
  pdf(paste0('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/01_gene_survival/',abbr,'_survival.pdf'),width = 5,height = 5)
  
  mafSurvival(maf = DLBC_LACTIC, genes = abbr, time = 'time', Status = 'status', isTCGA = TRUE)
  
  dev.off()
}


#################################


TCGA_clinical$lactic_mut<-ifelse(TCGA_clinical$sample %in% lactic_mutant_sample,"Mut","Wild_type")


library("survival")
library("survminer") 
library("survivalROC")
library("timeROC")


fit_validation <- survfit( Surv(time, status) ~ lactic_mut,data = TCGA_clinical )
pdf(file='/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/lactic_mutation_wildtype_survival.pdf',width = 7,height = 7)
ggsurvplot(fit_validation, data = TCGA_clinical,
           conf.int = TRUE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           linetype = "strata",
           palette = c("red","#2E9FDF"),
           legend = "bottom",
           legend.title = "Lactic genes mutation",
           legend.labs = c("Mut","Wild_type"))
dev.off()


#################################

library(dplyr)
library(TCGAbiolinks)
setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV')
query <- GDCquery(project = "TCGA-DLBC", 
                  data.category = "Copy Number Variation", 
                  data.type = "Masked Copy Number Segment")

GDCdownload(query, method = "api", files.per.chunk = 100)
segment_dat <- GDCprepare(query = query,save=TRUE, save.filename = "DLBC_CNV_download.rad")


A=load("/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/DLBC_CNV_download.rad")

tumorCNV<- eval(parse(text=A))

tumorCNV=tumorCNV[,2:7]
tumorCNV=tumorCNV[,c('Sample','Chromosome','Start','End','Num_Probes','Segment_Mean')]
write.table(tumorCNV,"DLBC_CNV.txt",sep="\t",
            quote = F,col.names = F,row.names = F)

tumorCNV$Sample <- substring(tumorCNV$Sample,1,16)
tumorCNV <- grep("01A$",tumorCNV$Sample) %>% 
  tumorCNV[.,]
tumorCNV[,1] <- tumorCNV$Sample
#tumorCNV <- tumorCNV[,-7]

write.table(tumorCNV,"DLBC_MaskedCopyNumberSegment.txt",sep="\t",
            quote = F,col.names = F,row.names = F)



################################# CNV

library(maftools)
DLBC.gistic <- readGistic(gisticAllLesionsFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/all_lesions.conf_99.txt", gisticAmpGenesFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/amp_genes.conf_99.txt", gisticDelGenesFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/del_genes.conf_99.txt", gisticScoresFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/scores.gistic", isTCGA=TRUE)




gisticOncoPlot(gistic=DLBC.gistic, sortByAnnotation=TRUE, top=10)

oncoplot(maf = DLBC_LACTIC, top = 30, fontSize = 0.8 ,showTumorSampleBarcodes = F )


DLBC.plus.gistic <- read.maf(maf=DLBC_data_lactic, 
                             gisticAllLesionsFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/all_lesions.conf_99.txt", 
                             gisticAmpGenesFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/amp_genes.conf_99.txt", 
                             gisticDelGenesFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/del_genes.conf_99.txt", 
                             gisticScoresFile="/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/scores.gistic")

DLBC.plus.gistic = read.maf(
  maf = DLBC_data_lactic,
  gisticAllLesionsFile = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/all_lesions.conf_99.txt",
  gisticAmpGenesFile = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/amp_genes.conf_99.txt",
  gisticDelGenesFile = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/del_genes.conf_99.txt",
  gisticScoresFile = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/scores.gistic",
  isTCGA = TRUE,
  verbose = FALSE, 
  clinicalData = TCGA_clinical
)


oncoplot(maf=DLBC.plus.gistic, borderCol=NULL, top=50)

#################################

CNV_gistic<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/GisticResults/all_thresholded.by_genes.txt',header = T,check.names = F,sep = '\t',row.names = 1)

CNV_gistic$`Locus ID`<-NULL
CNV_gistic$Cytoband<-NULL


CNV_gistic_lactic<-CNV_gistic[which(rownames(CNV_gistic) %in% gene$Tags),]


dim(CNV_gistic_lactic) #196  47


CNV_plot_data<-as.data.frame(matrix(data = 0,nrow = dim(CNV_gistic_lactic)[1],ncol = 4))
colnames(CNV_plot_data)<-c("gene","loss","none","gain")
CNV_plot_data$gene<-rownames(CNV_gistic_lactic)

for(i in 1:dim(CNV_gistic_lactic)[1])
{
  temp<-CNV_gistic_lactic[i,]
  loss<-dim(temp[,which(temp <0)])[2]
  none<-dim(temp[,which(temp ==0)])[2]
  gain<-dim(temp[,which(temp >0)])[2]
  if (is.null(gain))
  {
    CNV_plot_data[i,4]<-0
  }
  else
  {
    CNV_plot_data[i,4]<-gain/dim(CNV_gistic_lactic)[2]
  }
  
  if (is.null(loss))
  {
    CNV_plot_data[i,2]<-0
  }
  else
  {
    CNV_plot_data[i,2]<-loss/dim(CNV_gistic_lactic)[2]
  }
  
  if (is.null(none))
  {
    CNV_plot_data[i,3]<-0
  }
  else
  {
    CNV_plot_data[i,3]<-none/dim(CNV_gistic_lactic)[2]
  }

}








#################################

library(reshape2)

CNV_plot_data_use<-melt(CNV_plot_data)



plot5<-ggplot(CNV_plot_data_use,mapping = aes(gene,value,fill=variable))+geom_bar(stat='identity',position='fill') +labs(x = 'Gene',y = 'frequnency') +theme(axis.title =element_text(size = 5),axis.text =element_text(size = 5, color = 'black'))+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(name = "variable", values=c("loss"="aquamarine3","none" = "lightskyblue","gain" = "lightcoral"))



cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/CNV_ratio.pdf' ,
                 plot = plot5,width = 20,height = 5)



CNV_plot_data_use_snv<-CNV_plot_data_use[which(CNV_plot_data_use$gene %in% mut_gene_lactic),]


plot6<-ggplot(CNV_plot_data_use_snv,mapping = aes(gene,value,fill=variable))+geom_bar(stat='identity',position='fill') +labs(x = 'Gene',y = 'frequnency') +theme(axis.title =element_text(size = 15),axis.text =element_text(size = 15, color = 'black'))+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(name = "variable", values=c("loss"="aquamarine3","none" = "lightskyblue","gain" = "lightcoral"))+coord_flip() 



cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/CNV_ratio_mut_gene.pdf' ,
                 plot = plot6,width = 5,height = 10)



#################################

CNV_plot_data<-CNV_plot_data[order(CNV_plot_data$gain),]


CNV_loss_top15<-CNV_plot_data[c(1:15),]
CNV_gain_top15<-CNV_plot_data[c(182:196),]
CNV_Top<-rbind(CNV_gain_top15,CNV_loss_top15)

CNV_gistic_lactic_top<-CNV_gistic_lactic[which(rownames(CNV_gistic_lactic) %in% CNV_Top$gene),]



CNV_gene_sample<-as.data.frame(matrix(data = 0,nrow = dim(CNV_gistic_lactic_top)[1],ncol = dim(CNV_gistic_lactic_top)[2]))
colnames(CNV_gene_sample)<-colnames(CNV_gistic_lactic_top)
rownames(CNV_gene_sample)<-rownames(CNV_gistic_lactic_top)

for(i in 1:dim(CNV_gistic_lactic_top)[1])
{
  for (j in 1:dim(CNV_gistic_lactic_top)[2])
  {
    if(CNV_gistic_lactic_top[i,j]==0)
    {
      CNV_gene_sample[i,j]<-"Normal"
    }
    if(CNV_gistic_lactic_top[i,j]>0)
    {
      CNV_gene_sample[i,j]<-"Gain"
    }
    if(CNV_gistic_lactic_top[i,j]<0)
    {
      CNV_gene_sample[i,j]<-"Loss"
    }
  }
}
colnames(CNV_gene_sample)<-substr(colnames(CNV_gene_sample),1,15)

CNV_gene_sample<-CNV_gene_sample[order(rownames(CNV_gene_sample)),]



TCGA_DLBC_exp_tpm<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/TCGA/Merge_TCGA-DLBC_TPM.txt',header = T,row.names = 1,check.names = F,sep = '\t')

TCGA_DLBC_exp_tpm<-TCGA_DLBC_exp_tpm[,which(colnames(TCGA_DLBC_exp_tpm) %in% colnames(CNV_gene_sample))]
TCGA_DLBC_exp_tpm<-TCGA_DLBC_exp_tpm[,match(colnames(CNV_gene_sample),colnames(TCGA_DLBC_exp_tpm))]

TCGA_DLBC_exp_tpm<-TCGA_DLBC_exp_tpm[which(rownames(TCGA_DLBC_exp_tpm) %in% rownames(CNV_gene_sample)),]


TCGA_DLBC_exp_tpm<-TCGA_DLBC_exp_tpm[match(rownames(CNV_gene_sample),rownames(TCGA_DLBC_exp_tpm)),]


TCGA_DLBC_exp_tpm<-log2(TCGA_DLBC_exp_tpm+1)


j<-1
for(i in 1:dim(TCGA_DLBC_exp_tpm)[1])
{
  gene_name<-rownames(TCGA_DLBC_exp_tpm)[i]
  temp<-data.frame(CNV=t(CNV_gene_sample[i,]),Expression=t(TCGA_DLBC_exp_tpm[i,]))
  colnames(temp)<-c("CNV","Expression")
  

  if(dim(as.data.frame(table(temp$CNV)))[1] != 3)
  {
    temp[48:50,1]<-"Gain"
    temp[48:50,2]<-0
  }

  my_comparisons<-list(c("Normal","Gain"),c("Normal","Loss"),c("Loss","Gain"))
  p<- ggboxplot(temp,x="CNV" , y="Expression" ,fill= "CNV", palette = c("aquamarine3","lightskyblue","lightcoral"),scales = "free_x",xlab=F,ylab="log2(TPM+1)",width = 0.5,order=c("Loss","Normal","Gain"))+
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif")+
    theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="bottom")+ggtitle(gene_name)
  
  assign(paste0("p", j),p)
  j<-j+1

}

p131 <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,ncol = 6,nrow = 5,labels = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD"),rel_heights = c(1,1,1,1,1),label_size = 20,scale = c(0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/CNV_top_gene_expression.pdf' ,
                 plot = p131,width = 30,height = 25)



#################################


save.image('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/02_CNV_SNV/02_CNV_SNV.RData')












