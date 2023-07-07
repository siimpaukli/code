library(inSilicoMerging)
library(estimate)
library(inSilicoDb)

setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/01_N_T')

load('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/expression.RData')

GSE_data_normalized<-as.data.frame(GSE_data_normalized)

gene<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/genes.txt',header = T,sep = '\t')



################################## DEGs N vs T #########

library(limma)



GSE_data_normalized


GSE_data_clinical

table(GSE_data_clinical$type)


design1 <- model.matrix(~0+ factor(c(rep("Normal",40),rep("Tumor",77))))
colnames(design1) <- c("Normal", "Tumor")
row.names(design1)<-colnames(GSE_data_normalized)


Contrasts.matrix<-makeContrasts('Tumor-Normal',levels = design1)

es1<-as.matrix(GSE_data_normalized)

fit <- lmFit(es1, design1)
fit2<-contrasts.fit(fit,Contrasts.matrix)
fit2 <- eBayes(fit2)
allGeneSets <- topTable(fit2, coef=1, number=Inf)


allGeneSets$type <- ifelse(allGeneSets$adj.P.Val< 0.01 & abs(allGeneSets$logFC) > 0.5849625,
                           ifelse(allGeneSets$logFC > 0.5849625,'Up_in_Tumor','Up_in_Normal'),'Others')

write.table(allGeneSets,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/01_N_T/T_N_lastic_DEGS.txt',quote=F,sep="\t")

table(allGeneSets$type)


GSE_degs_sig<-allGeneSets[which(allGeneSets$type != "Others"),]


GSE_degs_sig_lactic<-GSE_degs_sig[which(rownames(GSE_degs_sig) %in% gene$Tags),]

write.table(GSE_degs_sig_lactic,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/01_N_T/T_N_lastic_DEGS_lac.txt',quote=F,sep="\t")

################################## 

deg_plot<-allGeneSets

deg_plot$log.adj.P.Val= -log10(deg_plot$adj.P.Val)

deg_plot$Description<-NA

deg_plot$Description<-ifelse(rownames(deg_plot) %in% rownames(GSE_degs_sig_lactic),rownames(GSE_degs_sig_lactic),NA)

library(ggpubr)
library(ggrepel)
library(ggplot2)
library(Cairo)

title <- "Lactic-DEGs between Tumor and Normal in DLBC"

plot4<-ggscatter(deg_plot, 
                 x = "logFC", 
                 y = "log.adj.P.Val", 
                 color = "type",
                 size = 1,
                 label = "Description", 
                 repel = T,
                 label.rectangle = T,
                 palette = c("grey", "blue","red"),
                 font.label = c(3, "plain"),
)+xlim(c(-7,7))+geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.6)+
        geom_vline(xintercept = 0.5849625,lty=4,col="black",lwd=0.6)+
        geom_vline(xintercept = -0.5849625,lty=4,col="black",lwd=0.6)+
        theme(legend.position="bottom",
              panel.grid=element_blank(),
              legend.title = element_blank(),
              legend.text= element_text(face="bold", color="black",family = "Times", size=8),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(face="bold", color="black", size=12),
              axis.text.y = element_text(face="bold",  color="black", size=12),
              axis.title.x = element_text(face="bold", color="black", size=12),
              axis.title.y = element_text(face="bold",color="black", size=12))+
        labs(x="logFC",y="-log10(adj.P.Val)",title= title)


################################## 


DEG <- GSE_degs_sig_lactic
Exp <- GSE_data_normalized[which(rownames(GSE_data_normalized) %in% rownames(DEG)),]

# group_list <- as.data.frame(c(rep("case", 5), rep("control", 5)))
# colnames(group_list) = "group"
# rownames(group_list) <- colnames(Exp[,2:ncol(Exp)])


group_list <- GSE_data_clinical


Exp_scale<-t(scale(t(Exp)))

library("ComplexHeatmap")
library("circlize")
library('dplyr')
heat_colors2 <- colorRamp2(c(-3,0,3), c( "springgreen4","white","red"))

col_an <-  HeatmapAnnotation(type =group_list$type, 
                             show_annotation_name = F, 
                             col = list(type = c("Tumor" = "orange", "Normal" = "lightblue")), 
                             show_legend = T,  
                             annotation_legend_param = list(title = "Type",nrow = 2), 
                             which = "col" )

plot3<-Heatmap(Exp_scale,
        name="Expression\n(Z-score)",
        #left_annotation = row_an,
        top_annotation = col_an,
        border='grey',
        rect_gp = gpar(col = NA),
        cluster_rows = T,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        col=heat_colors2,
        color_space = "RGB",
        cluster_columns = F,
        show_column_dend=F,
        show_row_dend=T,
        #row_order=data_subgroup_gene$Writer,
        #column_order=data_subgroup_subgroup$sample,
        show_column_names = F,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 4),
        gap = unit(2, "mm"),
        column_title = "",
        column_title_gp = gpar(fontsize = 6),
        #width=unit(12, "cm"),
        #heatmap_width=unit(20, "cm"),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 6, fontface = "bold"))) 

pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/01_N_T/heatmap_N_T.pdf',width = 10,height = 5)
draw(plot3)
dev.off()







################################## 

pheno<-GSE_data_clinical

GSE_data_normalized_lac<-GSE_data_normalized[which(rownames(GSE_data_normalized) %in% gene$Tags),]
data.pca <-prcomp(t(GSE_data_normalized_lac), center= F,scale.= F)

library(ggplot2)

PCi<-data.frame(data.pca$x,Type=pheno$type)

plot1<-ggplot(PCi,aes(x=PC1,y=PC2,col=Type))+
        geom_point(size=3,alpha=0.5)+ #Size and alpha just for fun
        scale_color_manual(values = c("lightblue","orange"))+ #your colors here
        theme_classic()+ggtitle("PCA anaylisis (Lactic genes)")

################################## 

library(patchwork)


p131 <- cowplot::plot_grid(plot4,plot1,ncol = 2,nrow = 1,labels = 'AUTO',rel_heights = c(1,1),label_size = 20,scale = c(0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/01_N_T/deg_pca.pdf' ,
                 plot = p131,width = 10,height = 5)


################################## 
setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment')




Gene <- data.frame(symbol = unique(gene$Tags))  # 1038 基因

# 本地分析
# 转id
library("org.Hs.eg.db")
keytypes(org.Hs.eg.db)
ensids <- Gene$symbol  # esembl ENSG00000223972
ensids <- as.character(ensids)
cols <- c("SYMBOL","ENTREZID")
resu <- select(org.Hs.eg.db, keys=ensids,columns=cols,keytype="SYMBOL")
gene_id <- unique(resu$ENTREZID)  # KIF21B, 23046

library(clusterProfiler)
KEGG <- enrichKEGG(gene = gene_id, 
                   #organism = "mmu", 
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 0.05, pAdjustMethod = "holm", minGSSize = 10,
                   qvalueCutoff = 0.2, use_internal_data = F)
KEGG_Pathway <- as.data.frame(KEGG @ result)

KEGG_Pathway_sig<-KEGG_Pathway[which(KEGG_Pathway$p.adjust <0.05),]

write.table(KEGG_Pathway, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment/KEGG_Pathway.txt",row.names = F,quote = F,sep="\t")
write.table(KEGG_Pathway_sig, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment/KEGG_Pathway_sig.txt",row.names = F,quote = F,sep="\t")

# GO分析
library(clusterProfiler)
library(org.Hs.eg.db)
ego_CC <- enrichGO(gene = gene_id,
                   #OrgDb=org.Mm.eg.db,
                   OrgDb=org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
CC = as.data.frame(ego_CC @ result)
CC_select<-CC[which(CC$p.adjust <0.05),]
write.table(CC, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichmentGO_CC.txt",row.names = F,quote = F,sep="\t")
write.table(CC_select, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment/CC_select.txt",row.names = F,quote = F,sep="\t")



ego_BP <- enrichGO(gene = gene_id,
                   #OrgDb=org.Mm.eg.db,
                   OrgDb=org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
BP = as.data.frame(ego_BP @ result)
BP_select<-BP[which(BP$p.adjust <0.05),]
write.table(BP, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment/GO_BP.txt" ,row.names = F,quote = F,sep="\t")
write.table(CC_select, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment/BP_select.txt",row.names = F,quote = F,sep="\t")


ego_MF <- enrichGO(gene = gene_id,
                   #OrgDb=org.Mm.eg.db,
                   OrgDb=org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
MF = as.data.frame(ego_MF @ result)
MF_select<-MF[which(MF$p.adjust <0.05),]
write.table(MF, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment/GO_MF.txt" ,row.names = F,quote = F,sep="\t")
write.table(CC_select, "/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/02_enrichment/MF_select.txt",row.names = F,quote = F,sep="\t")


library(stringr)

dotplot(ego_CC)+scale_y_discrete(labels=function(x) str_wrap(x, width=30))
dotplot(ego_BP)+scale_y_discrete(labels=function(x) str_wrap(x, width=30))
dotplot(ego_MF)+scale_y_discrete(labels=function(x) str_wrap(x, width=30))



# 和旋图-CC
library(clusterProfiler)
library(GOplot)
library(stringr)
library(org.Hs.eg.db)
ego <- data.frame(KEGG_Pathway)
ego$geneID <- str_replace_all(ego$geneID, "/", ", ")
ego$Category <- "KEGG"
ego <- ego[order(ego$Count,decreasing = T), ]
ego<-ego[order(ego$qvalue),]
ego <- ego[1:5, c(10,1,2,8,5)]
colnames(ego) <- c("Category", "ID", "Term", "Genes", "adj_pval")
head(ego)



genes <-allGeneSets[which(rownames(allGeneSets) %in% gene$Tags),]


genes$SYMBOL<-rownames(genes)
genes<-genes[,c(8,1,4,5)]

genes<-base::merge(resu,genes,by='SYMBOL')
colnames(genes)[2]<-'ID'
genes$SYMBOL<-NULL
colnames(genes)<-c('ID','logFC','P.Value','adj.P.Val')

circ <- circle_dat(ego, genes)  # 如果这里报错，看看下面的Load the included datasetego$Term,重点是colnames要一样
# 深入体会，这里其实只要差异基因的第一列DEG的colnames是ID就行，后面多少信息都可以带着，与上面的那个ID（基因名字）一致就可以
# 把esembl id换成基因名字,并且构建成"Category", "ID", "Term", "Genes", "adj_pval"这种格式
#circ<-circ[which(!is.na(circ$zscore)),]

circ2 <- merge(resu, circ, by.x = "ENTREZID", by.y = "genes", all = F)
#circ2 <- circ2[, c(4,5,6,2,9)]
colnames(circ2)[2] <-  "Genes"

temp_gene<-circ2$Genes[!duplicated(circ2$Genes)]
temp_Term<-circ2$term[!duplicated(circ2$term)]

chord <- chord_dat(data = circ2, genes =temp_gene, process = temp_Term)
chord<-as.data.frame(chord)
chord$Genes<-rownames(chord)
circ2_1<-circ2[!duplicated(circ2$Genes),]
circ2_1<-circ2_1[,c(2,7)]

chord<-base::merge(chord,circ2_1,by='Genes')
chord<-na.omit(chord)
rownames(chord)<-chord$Genes
chord$Genes<-NULL
chord<-as.matrix(chord)

GOChord(chord, space = 0.05, gene.order = "logFC", gene.space = 0.25, 
        gene.size = 2, 
        border.size = 0,
        lfc.col=c('red','black','cyan'),
        ribbon.col="")



#Load the included dataset
data(EC)
library(GOplot)
#Building the circ object
temp<-EC$david
temp1<-EC$genelist
circ<-circle_dat(temp, temp1)  # 如果上面报错，看这里的数据样式即可

################################## 

save.image('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/01_lactic_acis_genes/01_N_T.RData')












