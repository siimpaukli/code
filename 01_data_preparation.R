setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes')


############################
library(Biobase)
library(GEOquery)

GSE56315 <- getGEO('GSE56315', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL570	 <- getGEO('GPL570', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")



GSE32018 <- getGEO('GSE32018', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL6480 <- getGEO('GPL6480', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")


GSE56315 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE56315_series_matrix.txt.gz') 
GPL570 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL570.soft')

GSE32018 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE32018_series_matrix.txt.gz') 
GPL6480 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL6480.soft')

GSE181063 <- getGEO('GSE181063', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL14951 <- getGEO('GPL14951', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")

GSE181063 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE181063_series_matrix.txt.gz') 
GPL14951 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL14951.soft')


############################ GSE56315 ###########

setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO')


colnames(Table(GPL570))
Table(GPL570)[1:10,1:28]
write.csv(Table(GPL570)[,c(1,11)],"GPL570.csv", row.names = F)
genename = read.csv("GPL570.csv")

exprSet <- as.data.frame(exprs(GSE56315))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]
GSE56315_epr = express[!duplicated(express[,89]),]

GSE56315_epr<-as.data.frame(GSE56315_epr)

GSE56315_epr$Gene.Symbol<-strsplit(GSE56315_epr$Gene.Symbol, "///")
a<-dim(GSE56315_epr)
GSE56315_epr$Symbol<-"a"

for (i in 1:a[1])
{
  if(!is.na(GSE56315_epr$Gene.Symbol[i][[1]][1]))
  {
    GSE56315_epr$Symbol[i]<-GSE56315_epr$Gene.Symbol[i][[1]][1]
  }
}

GSE56315_epr = GSE56315_epr[!duplicated(GSE56315_epr$Symbol),]

GSE56315_epr<-GSE56315_epr[-1,]
rownames(GSE56315_epr) = GSE56315_epr$Symbol
GSE56315_epr<-GSE56315_epr[,-90]
GSE56315_epr<-GSE56315_epr[,-89]


pdata = pData(GSE56315)
clinical_GSE56315<-pdata[,c(2,10)]

colnames(clinical_GSE56315)<-c("sample","tissue")
clinical_GSE56315$type<-ifelse(clinical_GSE56315$tissue == "tissue: human healthy tonsils","Normal","Tumor")

write.table(clinical_GSE56315,'//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE56315.txt',col.names=T,quote=F,sep="\t")
write.table(GSE56315_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE56315_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")




############################ GSE32018 ###########

colnames(Table(GPL6480))
Table(GPL6480)[1:10,1:13]
write.csv(Table(GPL6480)[,c(1,7)],"GPL6480.csv", row.names = F)
genename = read.csv("GPL6480.csv")

exprSet <- as.data.frame(exprs(GSE32018))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]

GSE32018_epr = express[!duplicated(express[,128]),]


rownames(GSE32018_epr) = GSE32018_epr$GENE_SYMBOL
GSE32018_epr<-GSE32018_epr[-16,]
GSE32018_epr<-GSE32018_epr[,-128]



pdata = pData(GSE32018)
table(pdata[,8])


clinical_GSE32018<-pdata[,c(2,8)]

colnames(clinical_GSE32018)<-c("sample","tissue")
clinical_GSE32018_DLBC<-clinical_GSE32018[which(clinical_GSE32018$tissue == "Diffuse Large B cell Lymphoma patient" | clinical_GSE32018$tissue == "Lymph-node"),]
clinical_GSE32018_DLBC$type<-ifelse(clinical_GSE32018_DLBC$tissue == "Lymph-node","Normal","Tumor")


write.table(clinical_GSE32018_DLBC,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE32018_DLBC.txt',row.names=F,col.names=T,quote=F,sep="\t")

write.table(GSE32018_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE32018_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")

GSE32018_DLBCL_epr<-GSE32018_epr[,which(colnames(GSE32018_epr) %in% clinical_GSE32018_DLBC$sample)]

write.table(GSE32018_DLBCL_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE32018_DLBCL_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")




############################ GSE181063 ###########
colnames(Table(GPL14951))
Table(GPL14951)[1:10,1:13]
write.csv(Table(GPL14951)[,c(1,12)],"GPL14951.csv", row.names = F)
genename = read.csv("GPL14951.csv")

exprSet <- as.data.frame(exprs(GSE181063))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]

GSE181063_epr = express[!duplicated(express[,1311]),]

rownames(GSE181063_epr) = GSE181063_epr$Symbol
GSE181063_epr<-GSE181063_epr[,-1311]


clinical_GSE181063<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE181063_clinical.txt',header = T,check.names = F,sep = '\t')

rownames(clinical_GSE181063)<-clinical_GSE181063$sample

clinical_GSE181063_DLBC<-clinical_GSE181063[which(clinical_GSE181063$type == "DLBCL"),]


GSE181063_DLBC_epr<-GSE181063_epr[,which(colnames(GSE181063_epr) %in% clinical_GSE181063_DLBC$sample)]
write.table(GSE181063_DLBC_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE181063_DLBC_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")




########################### combat ###################


library(sva)

batch1 <- c(rep('GSE32018',29),rep('GSE56315',88))


GSE181063_DLBC_epr_use<-GSE181063_DLBC_epr
GSE181063_DLBC_epr_use$gene<-rownames(GSE181063_DLBC_epr_use)




GSE56315_epr_use<-GSE56315_epr
GSE56315_epr_use$gene<-rownames(GSE56315_epr_use)


GSE32018_DLBCL_epr_use<-GSE32018_DLBCL_epr
GSE32018_DLBCL_epr_use$gene<-rownames(GSE32018_DLBCL_epr_use)


GSE_data<-merge(GSE32018_DLBCL_epr_use,GSE56315_epr_use,by='gene')
rownames(GSE_data)<-GSE_data$gene
GSE_data$gene<-NULL

x_merge1 <- as.data.frame(GSE_data)
class(x_merge1)
x_merge2 <- as.matrix(x_merge1)

design <- model.matrix(~0 + batch1)


combat_edata <- ComBat(dat = x_merge2, batch = batch1, ref.batch="GSE56315")


GSE_data_clinical<-rbind(clinical_GSE32018_DLBC,clinical_GSE56315)

GSE_data_clinical<-GSE_data_clinical[order(GSE_data_clinical$type),]

GSE_data_normalized<-combat_edata[,match(GSE_data_clinical$sample,colnames(combat_edata))]


write.table(GSE_data_normalized,'//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/GSE_data_normalized.txt',row.names=T,col.names=T,quote=F,sep="\t")

write.table(GSE_data_clinical,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/GSE_data_clinical.txt',row.names=F,col.names=T,quote=F,sep="\t")



############################### GSE10846 #######

GSE10846 <- getGEO('GSE10846', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL570 <- getGEO('GPL570', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")

GSE10846 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE10846_series_matrix.txt.gz') 
GPL570 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL570.soft')


colnames(Table(GPL570))

Table(GPL570)[1:10,1:28]
write.csv(Table(GPL570)[,c(1,11)],"GPL570.csv", row.names = F)
genename = read.csv("GPL570.csv")

exprSet <- as.data.frame(exprs(GSE10846))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]
GSE10846_epr = express[!duplicated(express[,421]),]

GSE10846_epr<-as.data.frame(GSE10846_epr)

GSE10846_epr$Gene.Symbol<-strsplit(GSE10846_epr$Gene.Symbol, "///")
a<-dim(GSE10846_epr)
GSE10846_epr$Symbol<-"a"

for (i in 1:a[1])
{
  if(!is.na(GSE10846_epr$Gene.Symbol[i][[1]][1]))
  {
    GSE10846_epr$Symbol[i]<-GSE10846_epr$Gene.Symbol[i][[1]][1]
  }
}

GSE10846_epr = GSE10846_epr[!duplicated(GSE10846_epr$Symbol),]

GSE10846_epr<-GSE10846_epr[-1,]
rownames(GSE10846_epr) = GSE10846_epr$Symbol
GSE10846_epr<-GSE10846_epr[,-422]
GSE10846_epr<-GSE10846_epr[,-421]


pdata = pData(GSE10846)
clinical_GSE10846<-pdata[,c(2,17,18,10,11,20,21)]

colnames(clinical_GSE10846)<-c("sample","status","time","gender","age","ecog","stage")

clinical_GSE10846$status<-ifelse(clinical_GSE10846$status== 'Clinical info: Follow up status: DEAD',1,0)
clinical_GSE10846$time<-as.numeric(gsub("Clinical info: Follow up years\\: ","",clinical_GSE10846$time))

clinical_GSE10846$gender<-ifelse(clinical_GSE10846$gender== 'Gender: male',"Male","Female")

clinical_GSE10846$age<-as.numeric(gsub("Age\\: ","",clinical_GSE10846$age))
clinical_GSE10846$ecog<-as.numeric(gsub("Clinical info: ECOG performance status\\: ","",clinical_GSE10846$ecog))

clinical_GSE10846$stage<-as.numeric(gsub("Clinical info\\: Stage\\: ","",clinical_GSE10846$stage))



write.table(clinical_GSE10846,'//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE10846.txt',col.names=T,quote=F,sep="\t")
write.table(GSE10846_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE10846_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")

 


#################### GSE69053 ########
GSE69053 <- getGEO('GSE69053', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL14951 <- getGEO('GPL14951', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")

GSE69053 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE69053-GPL14951_series_matrix.txt.gz') 
GPL14951 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL14951.soft')


colnames(Table(GPL14951))

write.csv(Table(GPL14951)[,c(1,12)],"GPL14951.csv", row.names = F)
genename = read.csv("GPL14951.csv")

exprSet <- as.data.frame(exprs(GSE69053))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]

GSE69053_epr = express[!duplicated(express[,251]),]

rownames(GSE69053_epr) = GSE69053_epr$Symbol
GSE69053_epr<-GSE69053_epr[,-251]




clinical_GSE69053<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE69053.txt',header = T,check.names = F,sep = '\t')

rownames(clinical_GSE69053)<-clinical_GSE69053$sample


GSE69053_DLBC_epr<-GSE69053_epr[,which(colnames(GSE69053_epr) %in% clinical_GSE69053$sample)]


write.table(clinical_GSE69053,'//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE69053.txt',col.names=T,quote=F,sep="\t")
write.table(GSE69053_DLBC_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE69053_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")

#################### GSE23501 ########

GSE23501 <- getGEO('GSE23501', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL570 <- getGEO('GPL570', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")

GSE23501 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE23501_series_matrix.txt.gz') 
GPL570 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL570.soft')


colnames(Table(GPL570))

Table(GPL570)[1:10,1:28]
write.csv(Table(GPL570)[,c(1,11)],"GPL570.csv", row.names = F)
genename = read.csv("GPL570.csv")

exprSet <- as.data.frame(exprs(GSE23501))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]
GSE23501_epr = express[!duplicated(express[,70]),]

GSE23501_epr<-as.data.frame(GSE23501_epr)

GSE23501_epr$Gene.Symbol<-strsplit(GSE23501_epr$Gene.Symbol, "///")
a<-dim(GSE23501_epr)
GSE23501_epr$Symbol<-"a"

for (i in 1:a[1])
{
  if(!is.na(GSE23501_epr$Gene.Symbol[i][[1]][1]))
  {
    GSE23501_epr$Symbol[i]<-GSE23501_epr$Gene.Symbol[i][[1]][1]
  }
}

GSE23501_epr = GSE23501_epr[!duplicated(GSE23501_epr$Symbol),]

GSE23501_epr<-GSE23501_epr[-1,]
rownames(GSE23501_epr) = GSE23501_epr$Symbol
GSE23501_epr<-GSE23501_epr[,-71]
GSE23501_epr<-GSE23501_epr[,-70]



clinical_GSE23501<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE23501.txt',header = T,check.names = F,sep = '\t')


write.table(GSE23501_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE23501_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")


############################### GSE87371 #####

GSE87371 <- getGEO('GSE87371', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL570 <- getGEO('GPL570', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")

GSE87371 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE87371_series_matrix.txt.gz') 
GPL570 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL570.soft')


colnames(Table(GPL570))

Table(GPL570)[1:10,1:28]
write.csv(Table(GPL570)[,c(1,11)],"GPL570.csv", row.names = F)
genename = read.csv("GPL570.csv")

exprSet <- as.data.frame(exprs(GSE87371))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]

GSE87371_epr = express[!duplicated(express[,224]),]

GSE87371_epr<-as.data.frame(GSE87371_epr)

GSE87371_epr$Gene.Symbol<-strsplit(GSE87371_epr$Gene.Symbol, "///")
a<-dim(GSE87371_epr)
GSE87371_epr$Symbol<-"a"

for (i in 1:a[1])
{
  if(!is.na(GSE87371_epr$Gene.Symbol[i][[1]][1]))
  {
    GSE87371_epr$Symbol[i]<-GSE87371_epr$Gene.Symbol[i][[1]][1]
  }
}

GSE87371_epr = GSE87371_epr[!duplicated(GSE87371_epr$Symbol),]

GSE87371_epr<-GSE87371_epr[-1,]
rownames(GSE87371_epr) = GSE87371_epr$Symbol
GSE87371_epr<-GSE87371_epr[,-225]
GSE87371_epr<-GSE87371_epr[,-224]


clinical_GSE87371<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE87371.txt',header = T,check.names = F,sep = '\t')

clinical_GSE87371$time<-clinical_GSE87371$time/12


write.table(clinical_GSE87371,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE87371.txt',col.names=T,quote=F,sep="\t")
write.table(GSE87371_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE87371_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")



#################### GSE32918 #########

GSE32918 <- getGEO('GSE32918', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO") 
GPL8432 <- getGEO('GPL8432', destdir = "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO")

GSE32918 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE32918_series_matrix.txt.gz') 
GPL8432 <- getGEO(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GPL8432.soft')


colnames(Table(GPL8432))

write.csv(Table(GPL14951)[,c(1,12)],"GPL8432.csv", row.names = F)
genename = read.csv("GPL8432.csv")

exprSet <- as.data.frame(exprs(GSE32918))

exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))
express = express[order(rowMeans,decreasing = T),]

GSE32918_epr = express[!duplicated(express[,250]),]
GSE32918_epr<-na.omit(GSE32918_epr)
rownames(GSE32918_epr) = GSE32918_epr$Symbol

GSE32918_epr<-GSE32918_epr[,-250]


clinical_GSE32918<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE32918.txt',header = T,check.names = F,sep = '\t')

rownames(clinical_GSE32918)<-clinical_GSE32918$sample


write.table(clinical_GSE32918,'//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/clinical_GSE32918.txt',col.names=T,quote=F,sep="\t")
write.table(GSE32918_epr,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/GSE32918_epr.txt',row.names=T,col.names=T,quote=F,sep="\t")

##################################################

save(GSE181063_DLBC_epr,clinical_GSE181063_DLBC,GSE181063_epr,clinical_GSE181063,GSE32018_DLBCL_epr,clinical_GSE32018_DLBC,GSE32018_epr,clinical_GSE32018,GSE56315_epr,clinical_GSE56315,GSE_data_normalized,GSE_data_clinical,GSE69053_DLBC_epr,clinical_GSE69053,clinical_GSE23501,GSE23501_epr,clinical_GSE10846,GSE10846_epr,clinical_GSE87371,GSE87371_epr,GSE32918_epr,clinical_GSE32918,file = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/expression.RData')


##################################################
save.image('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GEO/01_data_preparation.RData')










