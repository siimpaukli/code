setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR')


compute_lasso_coef <- function(clinical = NULL, exp = NULL, unicox_res = NULL, seed = 123) {
  library(tidyverse)
  library(survival)
  library(ggpubr)
  var_name <- deparse(substitute(exp))
  
  unicox_exp <- exp[unicox_res$gene, ]
  clinical <- clinical[na.omit(match(colnames(unicox_exp), clinical[, 1])), ]
  unicox_exp <- unicox_exp[, na.omit(match(clinical[, 1], colnames(unicox_exp)))]
  # dim(unicox_exp)
  # dim(clinical)
  
  lasso_res <<- list()
  library(glmnet)
  y <- data.matrix(Surv(clinical$time, clinical$status))
  x <- unicox_exp %>%
    t() %>%
    as.data.frame() %>%
    as.matrix()
  
  set.seed(seed)
  fit <- glmnet(x, y, family = "cox", alpha = 1)
  lasso_res[[2]] <<- fit
  plot(fit)
  
  set.seed(seed)
  cvfit <- cv.glmnet(x, y, family = "cox")
  lasso_res[[3]] <<- cvfit
  plot(cvfit)
  
  tmp <- coef(object = cvfit, s = "lambda.min")
  Signature_Coef_min_os <- tmp@Dimnames[[1]][tmp %>%
                                               summary() %>%
                                               as.data.frame() %>%
                                               .[, 1]] %>%
    as.data.frame()
  
  Signature_Coef_min_os[, 2] <- tmp %>%
    summary() %>%
    as.data.frame() %>%
    .[, 3] %>%
    as.numeric()
  
  colnames(Signature_Coef_min_os) <- c("symbol", "coef")
  
  str(Signature_Coef_min_os)
  print(paste0("Lasso-cox gets ", dim(Signature_Coef_min_os)[1], " genes"))
  print(paste0("Univ_cox has ", nrow(unicox_exp), " genes"))
  
  lasso_res[[1]] <- Signature_Coef_min_os
  lasso_res[[1]][,1] <- as.character(lasso_res[[1]][,1])
  
  return(lasso_res)
}

multivariate_Cox <- function(exp = NULL, clinical = NULL) {
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
      formula =  as.formula(paste0("Surv(clinical$time, clinical$status)~",paste0("`",gene_name[1:length(gene_name)],"`",collapse = '+'))),
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
  
  return(result2)
}


lasso <- function(clinical = NULL, exp = NULL, unicox_res = NULL, seed = 123, saveplot = TRUE) {
  var_name <- deparse(substitute(exp))
  
  lasso_res <- compute_lasso_coef(
    clinical = clinical, exp = exp,
    unicox_res = unicox_res, seed = seed
  )
  
  temp_coef <- lasso_res[[1]] %>%
    mutate(
      group = ifelse(.$coef > 0, 1, 2),
      min_coef = min(coef), max_coef = max(coef)
    )
  
  coef_bar_p <- ggplot(
    data = temp_coef,
    mapping = aes(x = coef, y = factor(reorder(symbol, coef)), fill = factor(group))
  ) +
    geom_bar(stat = "identity", color = "black", width = .8, size = .6) +
    theme_pubr() +
    guides(fill = "none") +
    ylab("Gene") +
    xlab("Coeffecients") +
    # theme(
    #   axis.ticks.y = element_blank(),
    #   axis.line.y = element_blank()
    # ) +
    # scale_x_continuous(
    #   expand = c(0, 0),
    #   breaks = seq(round(temp_coef$min_coef[1] - 0.1, 1), round(temp_coef$max_coef[1] + 0.1, 1), 0.4)
    # ) +
    scale_fill_manual(values = c("1" = "#F4A01A", "2" = "#437BAA"))+
    theme(axis.text.y=element_text(vjust=1,size=10))
  
  lasso_res[[4]] <- coef_bar_p
  
  p <- cowplot::plot_grid(
    ~ plot(lasso_res[[2]]),
    ~ plot(lasso_res[[3]]),
    coef_bar_p + theme(plot.margin = unit(c(0.5, 0.55, 0.25, 0.25), "cm")), # 调整画布大小,
    labels = "AUTO",
    label_size = 20, nrow = 1
  )
  
  print(p)
  
  if (saveplot) {
    if (!dir.exists("./lasso")) {
      dir.create("./lasso")
    }
    
    if (!dir.exists("./lasso/lasso_res")) {
      dir.create("./lasso/lasso_res", mode = "7777", showWarnings = F)
    } else {
      (
        print("Dir is ready.")
      )
    }
    
    cowplot::ggsave2(
      filename = sprintf("./lasso/lasso_res/lasso_%s.pdf", var_name),
      plot = p, width = 12, height = 4
    )
    cowplot::ggsave2(
      filename = sprintf("./lasso/lasso_res/lasso_%s.tiff", var_name),
      plot = p, width = 12, height = 4, dpi = 300
    )
    cowplot::ggsave2(
      filename = sprintf("./lasso/lasso_res/lasso_%s_72dpi.tiff", var_name),
      plot = p, width = 12, height = 4, dpi = 72
    )
  } else {
    print("Save plot by Parameter 'saveplot = T'.")
  }
  return(lasso_res)
}

compute_score <- function(exp = NULL, model_coef = NULL, clinical = NULL) {
  var_name <- deparse(substitute(exp))
  clinical_name <- deparse(substitute(clinical))
  
  library(tidyverse)
  
  if(!exists("var_name")){
    var_name <- deparse(substitute(exp))
  }
  
  
  if (length(clinical)!=0) {
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ] %>%
      filter(time > 0 & time != "NA")
    exp <- exp[,na.omit(match(clinical[, 1], colnames(exp))) ]
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ]
  }
  
  exp <- exp[model_coef[, 1], ] %>%
    t() %>%
    as.data.frame()
  
  miss_gene <- model_coef[, 1][which(is.na(exp[1, ]))]
  miss_gene_symbol <- paste0(miss_gene, collapse = " ")
  
  exp[is.na(exp)] <- 0
  
  print(dim(exp))
  print(dim(model_coef[, 2, drop = F]))
  
  riskscore_res <- as.matrix(exp) %*% as.matrix(model_coef[, 2, drop = F]) %>%
    as.data.frame() %>%
    rename(riskscore = 1)
  
  if (length(miss_gene) == 1) {
    print(sprintf("%s is missd in %s", miss_gene_symbol, var_name))
  } else if (length(miss_gene) > 1) {
    print(sprintf("%s are missd in %s", miss_gene_symbol, var_name))
  } else {
    print(sprintf("All model_gene is mapped in %s.", var_name))
  }
  
  print(head(riskscore_res,5))
  return(riskscore_res)
} # 封装score计算函数
score_res <- function(riskscore_res = NULL, exp = NULL, clinical = NULL) {
  library(tidyverse)
  
  if (length(exp) != 0 & length(clinical) != 0) {
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ] %>%
      filter(time > 0 & time != "NA")
    exp <- exp[, na.omit(match(clinical[, 1], colnames(exp)))]
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ]
  } # 仅用于确认，在临床信息表中的样本，是同时具有临床信息和表达数据的样本
  
  score_res <- riskscore_res %>%
    rownames_to_column(var = "sample") %>%
    mutate(riskgroup = ifelse(riskscore > median(riskscore), "High", "Low")) %>% 
    inner_join(clinical, .) %>%
    mutate(time = as.numeric(.$time)) %>%
    na.omit() %>%
    dplyr::select(sample, status, time, riskgroup, riskscore)
  
  print(head(score_res, 5))
  print(str(score_res))
  return(score_res)
}
KM_ROC<- function(score_res = NULL, surtime_unit = 1) {
  library(tidyverse)
  library(survival)
  library(survminer)
  score_res <- score_res %>%
    mutate(HR_group = ifelse(.$riskscore > median(.$riskscore), 1, 0))
  
  coxtmp <- summary(coxph(Surv(time, status) ~ riskscore, data = score_res))
  
  HR <- coxtmp$coefficients[2]
  logrank_pvalue <- coxtmp$sctest["pvalue"]
  lower_.95 <- coxtmp$conf.int[, "lower .95"]
  upper_.95 <- coxtmp$conf.int[, "upper .95"]
  C <- coxtmp$concordance[1]
  
  test<-survdiff(Surv(time, status) ~ HR_group, data = score_res)
  # HR = (test$obs[2]/test$exp[2])/(test$obs[1]/test$exp[1])
  # HR = (test$obs[1]/test$exp[1])/(test$obs[2]/test$exp[2])
  # upper_.95 = exp(log(HR) + qnorm(0.975)*sqrt(1/test$exp[2]+1/test$exp[1]))
  # lower_.95 = exp(log(HR) - qnorm(0.975)*sqrt(1/test$exp[2]+1/test$exp[1]))
  # logrank_pvalue <- 1 - pchisq(test$chisq, length(test$n) -1)
  
  
  
  p_chara <- paste0(
    ifelse(logrank_pvalue < 0.001, "P < 0.001", paste0("P = ", round(logrank_pvalue, 3))),
    "  ",
    "HR = ", round(HR, 2),
    "\n95% CI = ", round(lower_.95, 2), " - ", round(upper_.95, 2),
    "  C-index = ", round(C, 2)
  )
  
  fit <- surv_fit(Surv(time, status) ~ riskgroup, data = score_res)
  # KM <- ggsurvplot(fit,
  #                  data = score_res,
  #                  surv.median.line = "hv",
  #                  legend.title = "Risk Score",
  #                  legend.labs = c("High", "Low"),
  #                  palette = "Set1",
  #                  ggtheme = theme_survminer(),
  #                  pval = p_chara,
  #                  pval.size = 4.5,
  #                  xlab = "Time (Days)",
  #                  tables.height = 0.28,
  #                  risk.table = T,
  # )
  
  
  KM <-ggsurvplot(fit, data = score_res,
                  conf.int = TRUE,
                  pval = p_chara,
                  pval.size = 3,
                  fun = "pct",
                  risk.table = TRUE,
                  #size = 1,
                  linetype = "strata",
                  palette = "Set1",
                  #legend = "bottom",
                  legend.title = "Risk Score",
                  legend.labs = c("High","Low"),
                  surv.median.line = "hv",
                  xlab = "Time (Years)",
                  tables.height = 0.1)
  
  
  library(timeROC)
  set.seed(1110)
  ROC_data <- timeROC(
    T = score_res$time,
    delta = score_res$status,
    marker = score_res$riskscore,
    cause = 1,
    weighting = "marginal",
    times = c(1 * surtime_unit, 3 * surtime_unit, 5 * surtime_unit),
    ROC= TRUE,
    iid = TRUE
  )
  
  rocFORplot <- tibble(
    year1x <- ROC_data$FP[, 1],
    year1y <- ROC_data$TP[, 1],
    year3x <- ROC_data$FP[, 2],
    year3y <- ROC_data$TP[, 2],
    year5x <- ROC_data$FP[, 3],
    year5y <- ROC_data$TP[, 3],
  )
  colnames(rocFORplot) <- c("year1x", "year1y", "year3x", "year3y", "year5x", "year5y")
  
  AUC_1year <- paste0("AUC at 1 years = ", sprintf('%.3f',ROC_data$AUC[[1]]))
  AUC_3year <- paste0("AUC at 3 years = ", sprintf('%.3f',ROC_data$AUC[[2]]))
  AUC_5year <- paste0("AUC at 5 years = ", sprintf('%.3f',ROC_data$AUC[[3]]))
  
  ROC_1 <- ggplot(data = rocFORplot) +
    geom_line(aes(x = year1x, y = year1y), size = 1.2, color = "#E83D3F") +
    geom_line(aes(x = year3x, y = year3y), size = 1.2, color = "#377EB8") +
    geom_line(aes(x = year5x, y = year5y), size = 1.2, color = "#4DAF4A") +
    geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
    egg::theme_article() +
    annotate(geom = "line", x = c(.43, 0.49), y = .17, colour = "#E83D3F", size = 1.2) +
    annotate("text", x = 0.5, y = .17, size = 4.35, label = AUC_1year, color = "black",hjust = "left") +
    annotate(geom = "line", x = c(.43, 0.49), y = .11, colour = "#377EB8", size = 1.2) +
    annotate("text", x = 0.5, y = .11, size = 4.35, label = AUC_3year, color = "black",hjust = "left") +
    annotate(geom = "line", x = c(.43, 0.49), y = .05, colour = "#4DAF4A", size = 1.2) +
    annotate("text", x = 0.5, y = .05, size = 4.35, label = AUC_5year, color = "black",hjust = "left") +
    labs(x = "1-Specificity", y = "Sensitivity") +
    theme(
      axis.text.x = element_text(face = "plain", size =12 , color = "black"),
      axis.text.y = element_text(face = "plain", size = 12, color = "black"),
      axis.title.x = element_text(face = "plain", size = 14, color = "black"),
      axis.title.y = element_text(face = "plain", size = 14, color = "black")
    )
  res <- list()
  res[[1]] <- KM
  res[[2]] <- ROC_1
  
  library(patchwork)
  p <- ggarrange(res[[1]]$plot, res[[1]]$table, ncol = 1, align = 'v',heights = c(0.75, 0.3))
  print(cowplot::plot_grid(p, res[[2]], rel_widths = c(0.95, 1.05),scale = c(1,.97)))
  
  return(res)
} # 封装KM、ROC曲线函数
KM_ROC_curve <- function(model_coef = NULL, exp = NULL, clinical = NULL, surtime_unit = 1,saveplot = TRUE) {
  riskscore_res <- compute_score(model_coef = model_coef, exp = exp, clinical = clinical)
  score_res <- score_res(riskscore_res = riskscore_res, clinical = clinical, exp = exp)
  KM_ROC_curve <- KM_ROC(score_res = score_res, surtime_unit = surtime_unit)
  KM_ROC_curve[[3]] <- score_res
  
  var_name <- deparse(substitute(exp))
  clinical_name <- deparse(substitute(clinical))
  
  if (saveplot) {
    
    if(!dir.exists('./lasso')){
      dir.create('./lasso')
    }
    
    if (!dir.exists(sprintf("./lasso/KM_ROC_curve_%s", var_name))) {
      dir.create(sprintf("./lasso/KM_ROC_curve_%s", var_name))
    } else {
      print(sprintf("Dir './lasso/KM_ROC_curve_%s' is existed.", var_name))
    }
    
    p <- ggarrange(KM_ROC_curve[[1]]$plot, KM_ROC_curve[[1]]$table, ncol = 1,align = 'v', heights = c(0.8, 0.28))
    
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_KM.pdf", var_name,var_name) ,
                     plot = p,width = 5.2,height = 5.5)
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_KM.tiff", var_name,var_name) ,
                     plot = p,width = 5.2,height = 5.5,dpi = 300)
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_KM_dpi72.tiff", var_name,var_name) ,
                     plot = p,width = 5.2,height = 5.5,dpi = 72)
    
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_ROC.pdf", var_name,var_name) ,
                     plot = KM_ROC_curve[[2]],width = 5.5,height = 5.2)
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_ROC.tiff", var_name,var_name) ,
                     plot = KM_ROC_curve[[2]],width = 5.5,height = 5.2,dpi = 300)
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_ROC_dpi72.tiff", var_name,var_name) ,
                     plot = KM_ROC_curve[[2]],width = 5.5,height = 5.2,dpi = 72)
    
    p2 <- cowplot::plot_grid(p, KM_ROC_curve[[2]], 
                             nrow = 1,labels = 'AUTO',rel_widths = c(0.95,1.05),label_size = 20,scale = c(1,.97))
    
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_KM-ROC.pdf", var_name,var_name) ,
                     plot = p2,width = 10,height = 5.2)
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_KM-ROC.tiff", var_name,var_name) ,
                     plot = p2,width = 10,height = 5.2,dpi = 300)
    cowplot::ggsave2(filename = sprintf("./lasso/KM_ROC_curve_%s/%s_KM-ROC_dpi72.tiff", var_name,var_name) ,
                     plot = p2,width = 10,height = 5.2,dpi = 72)
  }
  
  return(KM_ROC_curve)
}



compute_risk_score_clinical <- function(exp = NULL, model_coef = NULL, clinical = NULL) {
  var_name <- deparse(substitute(exp))
  clinical_name <- deparse(substitute(clinical))
  
  library(tidyverse)
  
  if(!exists("var_name")){
    var_name <- deparse(substitute(exp))
  }
  
  
  if (length(clinical)!=0) {
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ] %>%
      filter(time > 0 & time != "NA")
    exp <- exp[,na.omit(match(clinical[, 1], colnames(exp))) ]
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ]
  }
  
  exp <- exp[model_coef[, 1], ] %>%
    t() %>%
    as.data.frame()
  
  
  miss_gene <- rownames(model_coef)[which(is.na(exp[1, ]))]
  miss_gene_symbol <- paste0(miss_gene, collapse = " ")
  
  exp[is.na(exp)] <- 0
  
  print(dim(exp))
  print(dim(model_coef[, 2, drop = F]))
  
  riskscore_res <- as.matrix(exp) %*% as.matrix(model_coef[, 2, drop = F]) %>%
    as.data.frame() #%>%
  colnames(riskscore_res)<-'riskscore'
  #rename(riskscore = 1)
  
  if (length(miss_gene) == 1) {
    print(sprintf("%s is missd in %s", miss_gene_symbol, var_name))
  } else if (length(miss_gene) > 1) {
    print(sprintf("%s are missd in %s", miss_gene_symbol, var_name))
  } else {
    print(sprintf("All model_gene is mapped in %s.", var_name))
  }
  
  print(head(riskscore_res,5))
  return(riskscore_res)
} # 封装score计算函数
risk_score_res_clinical <- function(riskscore_res = NULL, exp = NULL, clinical = NULL) {
  library(tidyverse)
  
  if (length(exp) != 0 & length(clinical) != 0) {
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ] %>%
      filter(time > 0 & time != "NA")
    exp <- exp[, na.omit(match(clinical[, 1], colnames(exp)))]
    clinical <- clinical[na.omit(match(colnames(exp), clinical[, 1])), ]
  } # 仅用于确认，在临床信息表中的样本，是同时具有临床信息和表达数据的样本
  
  score_res <- riskscore_res %>%
    rownames_to_column(var = "sample") %>%
    mutate(riskgroup = ifelse(riskscore > median(riskscore), "High", "Low")) %>% 
    inner_join(clinical, .) %>%
    mutate(time = as.numeric(.$time)) %>%
    na.omit() %>%
    dplyr::select(sample, status, time, riskgroup, riskscore)
  
  print(head(score_res, 5))
  print(str(score_res))
  return(score_res)
}



###############################

GSE181063_univcox_sig_0.01<-GSE181063_univcox[which(GSE181063_univcox$p.value < 0.01),]



clinical_GSE181063_DLBC<-clinical_GSE181063_DLBC[which(clinical_GSE181063_DLBC$time >0),]


GSE181063_multivcox<-multivariate_Cox(exp= GSE181063_DLBC_epr_lactic,clinical=clinical_GSE181063_DLBC)
write.table(GSE181063_multivcox,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/GSE181063_multivcox.txt',row.names=F,col.names=T,quote=F,sep="\t")
GSE181063_multivcox_sig_0.01<-GSE181063_multivcox[which(GSE181063_multivcox$p.value < 0.01),]


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
GSE181063_multivcox<-GSE181063_multivcox[order(GSE181063_multivcox$p.value),]

GSE181063_multivcox_plot<-unicox_forest_plot(unicox_re=GSE181063_multivcox,plot_file_name='/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/clinical_multi_forest.pdf')



GSE181063_univcox_multicox<-GSE181063_univcox_sig_0.01[which(rownames(GSE181063_univcox_sig_0.01) %in% rownames(GSE181063_multivcox_sig_0.01)),]



GSE181063_lasso<-lasso(clinical = clinical_GSE181063_DLBC, exp = GSE181063_DLBC_epr_lactic, unicox_res = GSE181063_univcox_multicox, seed = 123, saveplot = TRUE)

#lasso gene 27

GSE181063_lasso[[1]]
write.table(GSE181063_lasso[[1]],'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/lasso/GSE181063_lasso_results.txt',quote=F,sep="\t",row.names = F)

############################### GSE181063 #######


a <- KM_ROC_curve(model_coef =  GSE181063_lasso[[1]],exp = GSE181063_DLBC_epr_lactic,clinical = clinical_GSE181063_DLBC,surtime_unit = 1)




############################### validation ########
############################### TCGA ####### no sig ######



TCGA_clinical_lasso<-TCGA_clinical
TCGA_clinical_lasso$time<-TCGA_clinical_lasso$time/365


TCGA_DLBC_tpm<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/TCGA/Merge_TCGA-DLBC_TPM.txt',header = T,row.names = 1,check.names = F,sep = '\t')

colnames(TCGA_DLBC_tpm)<-substr(colnames(TCGA_DLBC_tpm),1,12)
TCGA_DLBC_log2tpm<-log2(TCGA_DLBC_tpm+1)

b <- KM_ROC_curve(model_coef =  GSE181063_lasso[[1]],exp = TCGA_DLBC_log2tpm,clinical = TCGA_clinical_lasso,surtime_unit = 1)
############################### GSE10846 ####### sig ######

c <- KM_ROC_curve(model_coef =  GSE181063_lasso[[1]],exp = GSE10846_epr,clinical = clinical_GSE10846,surtime_unit = 1)


#################### GSE69053 ####### sig ######

d <- KM_ROC_curve(model_coef =  GSE181063_lasso[[1]],exp = GSE69053_DLBC_epr,clinical = clinical_GSE69053,surtime_unit = 1)
#################### GSE23501 #######   no sig ######

e <- KM_ROC_curve(model_coef =  GSE181063_lasso[[1]],exp = GSE23501_epr,clinical = clinical_GSE23501,surtime_unit = 1)

############################### GSE87371 no sig ############# 
f <- KM_ROC_curve(model_coef =  GSE181063_lasso[[1]],exp = GSE87371_epr,clinical = clinical_GSE87371,surtime_unit = 1)


############################### GSE32918 ####### sig ######
g <- KM_ROC_curve(model_coef =  GSE181063_lasso[[1]],exp = GSE32918_epr,clinical = clinical_GSE32918,surtime_unit = 1)




############################### LARscore #########

GSE181063_lar_score<-compute_risk_score_clinical(model_coef =  GSE181063_lasso[[1]],exp = GSE181063_DLBC_epr_lactic,clinical = clinical_GSE181063_DLBC)
GSE181063_lar_score_final<-risk_score_res_clinical(riskscore_res = GSE181063_lar_score, exp = GSE181063_DLBC_epr_lactic, clinical = clinical_GSE181063_DLBC)


GSE10846_lar_score<-compute_risk_score_clinical(model_coef =  GSE181063_lasso[[1]],exp = GSE10846_epr,clinical = clinical_GSE10846)
GSE10846_lar_score_final<-risk_score_res_clinical(riskscore_res = GSE10846_lar_score, exp = GSE10846_epr, clinical = clinical_GSE10846)

GSE69053_lar_score<-compute_risk_score_clinical(model_coef =  GSE181063_lasso[[1]],exp = GSE69053_DLBC_epr,clinical = clinical_GSE69053)
GSE69053_lar_score_final<-risk_score_res_clinical(riskscore_res = GSE69053_lar_score, exp = GSE69053_DLBC_epr, clinical = clinical_GSE69053)

GSE32918_lar_score<-compute_risk_score_clinical(model_coef =  GSE181063_lasso[[1]],exp = GSE32918_epr,clinical = clinical_GSE32918)
GSE32918_lar_score_final<-risk_score_res_clinical(riskscore_res = GSE32918_lar_score, exp = GSE32918_epr, clinical = clinical_GSE32918)



GSE181063_DLBC_epr_lactic_1<-as.data.frame(t(GSE181063_DLBC_epr_lactic))
GSE181063_DLBC_epr_lactic_1<-GSE181063_DLBC_epr_lactic_1[,which(colnames(GSE181063_DLBC_epr_lactic_1) %in% GSE181063_lasso[[1]]$symbol)]
GSE181063_DLBC_epr_lactic_1$sample<-rownames(GSE181063_DLBC_epr_lactic_1)
GSE181063_lar_score_final_1<-merge(GSE181063_lar_score_final,GSE181063_DLBC_epr_lactic_1,by='sample')


GSE10846_epr_1<-as.data.frame(t(GSE10846_epr))
GSE10846_epr_1<-GSE10846_epr_1[,which(colnames(GSE10846_epr_1) %in% GSE181063_lasso[[1]]$symbol)]
GSE10846_epr_1$sample<-rownames(GSE10846_epr_1)
GSE10846_lar_score_final_1<-merge(GSE10846_lar_score_final,GSE10846_epr_1,by='sample')

GSE69053_DLBC_epr_1<-as.data.frame(t(GSE69053_DLBC_epr))
GSE69053_DLBC_epr_1<-GSE69053_DLBC_epr_1[,which(colnames(GSE69053_DLBC_epr_1) %in% GSE181063_lasso[[1]]$symbol)]
GSE69053_DLBC_epr_1$sample<-rownames(GSE69053_DLBC_epr_1)
GSE69053_lar_score_finall_1<-merge(GSE69053_lar_score_final,GSE69053_DLBC_epr_1,by='sample')

GSE32918_epr_1<-as.data.frame(t(GSE32918_epr))
GSE32918_epr_1<-GSE32918_epr_1[,which(colnames(GSE32918_epr_1) %in% GSE181063_lasso[[1]]$symbol)]
GSE32918_epr_1$sample<-rownames(GSE32918_epr_1)
GSE32918_lar_score_final_1<-merge(GSE32918_lar_score_final,GSE32918_epr_1,by='sample')


Score<-list(GSE181063=GSE181063_lar_score_final_1,GSE32918=GSE32918_lar_score_final_1,GSE10846=GSE10846_lar_score_final_1,GSE69053=GSE69053_lar_score_finall_1)



###############################

path_plot<-'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/lasso/'


library(lattice)
library("circlize")
library("ComplexHeatmap")


for(i in 1:4)
{
  risk<-Score[[i]]
  var_name <-names(Score)[i]
  risk<-risk[order(risk$riskscore),]
  risk$patient<-c(1:dim(risk)[1])
  p1<-xyplot(riskscore~patient,data=risk,groups=riskgroup,auto.key=list(corner=c(1,1)))
  pdf(paste0(path_plot,var_name,'_riskplot.pdf'),width = 6,height = 3)
  print(p1)
  dev.off()
  
  risk$status<-ifelse(risk$status == 0 ,"Alive","Dead")
  p2<-xyplot(riskscore~time,data=risk,groups=status,auto.key=list(corner=c(1,1)))
  pdf(paste0(path_plot,var_name,'_liveplot.pdf'),width = 6,height = 3)
  print(p2)
  dev.off()
  
  heat_colors2 <- colorRamp2(c(-2, 0, 2), c("midnightblue", "white", "red"))
  col_an <-  HeatmapAnnotation(type = risk$riskgroup, 
                               show_annotation_name = F, 
                               col = list(type = c("Low" = "Blue","High"="red")), 
                               show_legend = T,  
                               annotation_legend_param = list(title = "Type",nrow = 1), 
                               which = "col" )
  
  
  
  data_ccp_heatmap_2<-as.data.frame(t(scale(risk[,6:25])))
  data_ccp_heatmap_2<-data_ccp_heatmap_2[which(rownames(data_ccp_heatmap_2) %in% GSE181063_lasso[[1]]$symbol),]
  data_ccp_heatmap_2<-as.matrix(data_ccp_heatmap_2)
  p3<-Heatmap(data_ccp_heatmap_2,
              name="Subgroup",
              #left_annotation = row_an,
              top_annotation = col_an,
              border='grey',
              rect_gp = gpar(col = NA),
              cluster_rows = T,
              col=heat_colors2,
              color_space = "RGB",
              cluster_columns = F,
              show_column_dend=F,
              show_row_dend=T,
              #row_order=data_subgroup_gene$Writer,
              #column_order=data_subgroup_subgroup$sample,
              show_column_names = F,
              show_row_names = T,
              row_names_gp = gpar(fontsize = 12),
              gap = unit(2, "mm"),
              column_title = "",
              column_title_gp = gpar(fontsize = 6),
              #width=unit(12, "cm"),
              #heatmap_width=unit(10, "cm"),
              show_heatmap_legend = TRUE,
              heatmap_legend_param=list(labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize = 6, fontface = "bold")))
  pdf(paste0(path_plot,var_name,'_heatmap.pdf'),width = 6,height = 6)
  draw(p3)
  dev.off()
}

############################### TIDE calculation #####

library(clusterProfiler)
library(org.Hs.eg.db)

gene.tcga <- bitr(rownames(GSE181063_DLBC_epr), fromType = "SYMBOL", toType ="ENTREZID", OrgDb = org.Hs.eg.db)

GSE181063_tide<-GSE181063_DLBC_epr

GSE181063_tide$SYMBOL<-rownames(GSE181063_tide)

GSE181063_tide_final<-base::merge(gene.tcga,GSE181063_tide,by='SYMBOL')
GSE181063_tide_final$SYMBOL<-NULL
rownames(GSE181063_tide_final)<-GSE181063_tide_final$ENTREZID
GSE181063_tide_final$ENTREZID<-NULL
GSE181063_tide_final<-log2(GSE181063_tide_final+1)
write.table(GSE181063_tide_final,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE181063_tide.txt',row.names=T,col.names=T,quote=F,sep="\t")



gene.tcga <- bitr(rownames(GSE10846_epr), fromType = "SYMBOL", toType ="ENTREZID", OrgDb = org.Hs.eg.db)

GSE10846_tide<-GSE10846_epr

GSE10846_tide$SYMBOL<-rownames(GSE10846_tide)

GSE10846_tide_final<-base::merge(gene.tcga,GSE10846_tide,by='SYMBOL')
GSE10846_tide_final$SYMBOL<-NULL
rownames(GSE10846_tide_final)<-GSE10846_tide_final$ENTREZID
GSE10846_tide_final$ENTREZID<-NULL
GSE10846_tide_final<-log2(GSE10846_tide_final+1)
write.table(GSE10846_tide_final,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE10846_tide.txt',row.names=T,col.names=T,quote=F,sep="\t")



gene.tcga <- bitr(rownames(GSE32918_epr), fromType = "SYMBOL", toType ="ENTREZID", OrgDb = org.Hs.eg.db)

GSE32918_tide<-GSE32918_epr

GSE32918_tide$SYMBOL<-rownames(GSE32918_tide)

GSE32918_tide_final<-base::merge(gene.tcga,GSE32918_tide,by='SYMBOL')
GSE32918_tide_final$SYMBOL<-NULL
rownames(GSE32918_tide_final)<-GSE32918_tide_final$ENTREZID
GSE32918_tide_final$ENTREZID<-NULL
GSE32918_tide_final<-log2(GSE32918_tide_final+1)
write.table(GSE32918_tide_final,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE32918_tide.txt',row.names=T,col.names=T,quote=F,sep="\t")


gene.tcga <- bitr(rownames(GSE69053_DLBC_epr), fromType = "SYMBOL", toType ="ENTREZID", OrgDb = org.Hs.eg.db)

GSE69053_tide<-GSE69053_DLBC_epr

GSE69053_tide$SYMBOL<-rownames(GSE69053_tide)

GSE69053_tide_final<-base::merge(gene.tcga,GSE69053_tide,by='SYMBOL')
GSE69053_tide_final$SYMBOL<-NULL
rownames(GSE69053_tide_final)<-GSE69053_tide_final$ENTREZID
GSE69053_tide_final$ENTREZID<-NULL
GSE32918_tide_final<-log2(GSE32918_tide_final+1)
write.table(GSE69053_tide_final,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE69053_tide.txt',row.names=T,col.names=T,quote=F,sep="\t")


system("tidepy /Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE10846_tide.txt -o /Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE10846_tide_tide.txt -c Other")


GSE69053_tide<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE69053_tide.txt.tide.txt',header = T,row.names = 1,check.names = F,sep = '\t')
GSE10846_tide<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE10846_tide.txt.tide.txt',header = T,row.names = 1,check.names = F,sep = '\t')
GSE32918_tide<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE32918_tide.txt.tide.txt',header = T,row.names = 1,check.names = F,sep = '\t')
GSE181063_tide<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/GSE181063_tide.txt.tide.txt',header = T,row.names = 1,check.names = F,sep = '\t')


###############################

TIDE_results<-list(GSE181063=GSE181063_tide,GSE32918=GSE32918_tide,GSE10846=GSE10846_tide,GSE69053=GSE69053_tide)



path_plot<-'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/'


library(lattice)
library("circlize")
library("ComplexHeatmap")

pmat<-data.frame(matrix(data = 1,nrow = 2,ncol = 4))
colnames(pmat)<-names(TIDE_results)
rownames(pmat)<-c("LARscore Low","LARscore High")

coormat<-data.frame(matrix(data = 1,nrow = 2,ncol = 4))
colnames(coormat)<-names(TIDE_results)
rownames(coormat)<-c("LARscore Low","LARscore High")


j<-1
for(i in 1:4)
{
  TIDE<-TIDE_results[[i]]
  var_name <-names(TIDE_results)[i]
  TIDE$sample<-rownames(TIDE)
  TIDE<-TIDE[,c(14,3)]
  
  risk<-Score[[i]]
  risk<-risk[,c(1,4,5)]
  temp<-merge(risk,TIDE,by='sample')
  colnames(temp)<-c("sample","LAR group","LARscore","TIDE")
  low_group<-temp[which(temp$`LAR group` =="Low"),]
  high_group<-temp[which(temp$`LAR group` =="High"),]
  
  corr_re_low<-cor.test(low_group$LARscore,low_group$TIDE,method = 'spearman')
  corr_re_high<-cor.test(high_group$LARscore,high_group$TIDE,method = 'spearman')
  
  pmat[1,i]<-corr_re_low$p.value
  pmat[2,i]<-corr_re_high$p.value
  
  coormat[1,i]<-corr_re_low$estimate
  coormat[2,i]<-corr_re_high$estimate
  
  
  my_comparisons <- list(c("Low", "High"))
  
  p<-ggboxplot(temp,x="LAR group",y="TIDE",notch=F,xlab="LAR risk",ylab="TIDE",order = c("Low","High"),fill = c("lightgreen","red"))+stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+ggtitle(var_name)
  
  
  assign(paste0("p", j),p)
  j<-j+1

}

p131 <- cowplot::plot_grid(p1,p2,p3,p4,ncol = 2,nrow = 2,labels = 'AUTO',rel_heights = c(1,1),label_size = 20,scale = c(0.9,0.9,0.9,0.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/LARrisk_TIDE_all.pdf' ,
                 plot = p131,width = 10,height = 10)





library("ComplexHeatmap")
library("circlize")
library(dplyr)
plotMutiHeatmap=function(up,down,up_break,up_colors,down_break,down_colors,title){
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
      if(down[i, j]>0.1){
      }
    }
  }
  
  p1 <- Heatmap(up, column_title = title
                , rect_gp = gpar(type = "none")
                , show_heatmap_legend = F
                , cluster_rows = F
                , cluster_columns = F, 
                cell_fun = DiagFunc(up = up, down = down) 
  ) 
  col_fun = colorRamp2(down_break, down_colors) 
  lgd <- Legend(title = "Spearman Corralation", 
                col_fun = col_fun, 
                at = c(-0.5,0,0.5), 
                labels = c("-0.5","0","0.5"),  
                direction = "horizontal" 
  )
  col_fun2 = colorRamp2(up_break, up_colors) 
  lgd2 <- Legend(title = "-log10(p-value)", 
                 col_fun = col_fun2, 
                 at = c(0,1,2,3,4,5), 
                 labels = c('0',"1","2","3","4",">5"),  
                 direction = "horizontal"
  )
  
  draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "bottom"
       ,heatmap_legend_side = "bottom", merge_legend = TRUE)
}


up_break=c(0, 5)
down_break=c(-0.5,0,0.5)
up_colors=c("#FFFFFF","#6f9a8d")
down_colors=c("blue",'white',"red")


up=-log10(t(pmat))
down=t(coormat)
title='Spreman correlation between TIDE and LARscore'
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/02_TIDE/TIDE_cortest.pdf',width = 5,height = 10)
plotMutiHeatmap(up,down,up_break,up_colors,down_break,down_colors,title)
dev.off()





###############################
setwd('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/03_r_chop')
GSE181063_chemo<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/03_r_chop/GSE181063_chemo.txt',header = T,check.names = F,sep = '\t')

GSE10846_chemo<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/03_r_chop/GSE10846_chemo.txt',header = T,check.names = F,sep = '\t')



GSE181063_chemo_final<-merge(GSE181063_lar_score_final,GSE181063_chemo,by='sample')
GSE10846_chemo_final<-merge(GSE10846_lar_score_final,GSE10846_chemo,by='sample')


GSE181063_chemo_final$status<-ifelse(GSE181063_chemo_final$status == 0 ,"Alive","Dead")


GSE10846_chemo_final$status<-ifelse(GSE10846_chemo_final$status == 0 ,"Alive","Dead")

##################
my_comparisons <- list(c("Alive", "Dead"))

p2<-ggboxplot(GSE181063_chemo_final,x="status",y="riskscore",notch=F,xlab="",ylab="LARScore",fill=c("slateblue1","orange"),order=c("Alive", "Dead"))+stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+ggtitle("GSE181063")

results<-as.data.frame(table(GSE181063_chemo_final$status,GSE181063_chemo_final$riskgroup))

colnames(results)<-c("Status","LARgroup","Count")

results1<-results
results1$Count[1]<-153/(153+217)
results1$Count[2]<-217/(153+217)
results1$Count[3]<-317/(317+149)
results1$Count[4]<-149/(317+149)
results1$Count<-round(results1$Count,digits=2)
colnames(results1)<-c("Status","LARgroup","Percentage")

p3<-ggplot(results1,aes(x=LARgroup,y=Percentage,fill=Status))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+geom_text(aes(label = results1$Percentage),position=position_stack(vjust =0.5),size = 5)+ggtitle("GSE181063")+theme(legend.position="top")

library(pROC)
a<-rainbow(20)
# Gene ROC
rocobj<- roc(GSE181063_chemo_final$status, GSE181063_chemo_final$riskscore,plot=TRUE, print.auc=TRUE)#感兴趣的基因

#p1<-cowplot::plot_to_gtable(plot(rocobj,col="blue",plot=TRUE, print.auc=TRUE))
rocobj$auc
p1<-ggroc(rocobj)+annotate("text",x=0.5,y=0.2,label="AUC=0.6862")+ggtitle("GSE181063")+theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


##################
my_comparisons <- list(c("Alive", "Dead"))

p5<-ggboxplot(GSE10846_chemo_final,x="status",y="riskscore",notch=F,xlab="",ylab="LARScore",fill=c("slateblue1","orange"),order=c("Alive", "Dead"))+stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+ggtitle("GSE10846")

results2<-as.data.frame(table(GSE10846_chemo_final$status,GSE10846_chemo_final$riskgroup))

colnames(results2)<-c("Status","LARgroup","Count")

results3<-results2
results3$Count[1]<-47/(47+28)
results3$Count[2]<-28/(47+28)
results3$Count[3]<-107/(107+24)
results3$Count[4]<-24/(107+24)
results3$Count<-round(results3$Count,digits=2)
colnames(results3)<-c("Status","LARgroup","Percentage")

p6<-ggplot(results3,aes(x=LARgroup,y=Percentage,fill=Status))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+geom_text(aes(label = results3$Percentage),position=position_stack(vjust =0.5),size = 5)+ggtitle("GSE10846")+theme(legend.position="top")

library(pROC)
a<-rainbow(20)
# Gene ROC
rocobj<- roc(GSE10846_chemo_final$status, GSE10846_chemo_final$riskscore,plot=TRUE, print.auc=TRUE)#感兴趣的基因

#p1<-cowplot::plot_to_gtable(plot(rocobj,col="blue",plot=TRUE, print.auc=TRUE))
rocobj$auc
p4<-ggroc(rocobj)+annotate("text",x=0.5,y=0.2,label="AUC=0.625")+ggtitle("GSE10846")+theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


p7 <- cowplot::plot_grid(p3,p2,p1,p6,p5,p4,ncol = 3,nrow = 2,labels = 'AUTO',rel_heights = c(1,1,1,1,1,1),label_size = 20,scale = c(0.9,.9,.9,0.9,.9,.9))

cowplot::ggsave2(filename = '/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/03_r_chop/RCHOP.pdf' ,
                 plot = p7,width = 15,height = 15)




############################### GDSC #####


GDSC_exp<-read.table('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GDSC/Cell_line_RMA_proc_basalExp.txt',header = T,check.names = F,sep = '\t')
GDSC_exp<-GDSC_exp[!duplicated(GDSC_exp$GENE_SYMBOLS),]
rownames(GDSC_exp)<-GDSC_exp$GENE_SYMBOLS
GDSC_exp$GENE_title<-NULL
GDSC_exp$sample<-rownames(GDSC_exp)
GDSC_exp$GENE_SYMBOLS<-NULL

res_gdsc_degs<-GSE181063_lasso[[1]]


GDSC_exp_use<-GDSC_exp[which(rownames(GDSC_exp) %in% res_gdsc_degs$symbol),]


res_gdsc_degs<-res_gdsc_degs[which(res_gdsc_degs$symbol %in% rownames(GDSC_exp_use)),]


data_deg_ex_gdsc<-GDSC_exp_use
data_deg_ex_gdsc$sample<-NULL
data_deg_ex_gdsc<-log(data_deg_ex_gdsc+1)
#data_deg_ex_gdsc<-as.data.frame(t(data_deg_ex_gdsc))


data_gene_score<-data_deg_ex_gdsc
data_gene_score[,]=0

library(dplyr)
library(psych)
library(tidyverse)


res_gdsc_degs_use<-res_gdsc_degs %>% dplyr::slice(match(rownames(data_gene_score),res_gdsc_degs$symbol))

for (i in 1:1018)
{
  for(j in 1:9)
  {
    data_gene_score[j,i]<-as.numeric(res_gdsc_degs_use[j,2])*data_deg_ex_gdsc[j,i]
  }
}

data_gene_score<-as.data.frame(t(data_gene_score))



for (i in 1:1018)
{
  data_gene_score$MPIrisk_score[i]<-sum(data_gene_score[i,1:9])
}

data_gene_score_final<-data_gene_score
data_gene_score_final$sample<-rownames(data_gene_score_final)
data_gene_score_final[,1:9]<-NULL

write.table(data_gene_score_final,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/04_GDSC/GDSC_LARrisk_score.txt',row.names=T,col.names=T,quote=F,sep="\t")

rownames(data_gene_score_final)<-gsub('DATA.','',rownames(data_gene_score_final))
data_gene_score_final$sample<-rownames(data_gene_score_final)

gdsc_drug_ic<-read.csv('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GDSC/PANCANCER_IC_Mon Apr 26 16_16_19 2021.csv',header = T)
gdsc_drug_ic[,10:13]<-NULL
gdsc_drug_ic[,5:8]<-NULL
colnames(gdsc_drug_ic)[4]<-"sample"


GDSC_wm_score_drug<-merge(gdsc_drug_ic,data_gene_score_final,by='sample')

drug_name<-as.data.frame(table(GDSC_wm_score_drug$Drug.name))
colnames(drug_name)<-c("drug","counts")

library(psych)
DN<-as.character(drug_name$drug)
mydata1 <- data.frame(Drugs=c(NA),Othergenes=c(NA),estimate.rho=c(NA),pvalue=c(NA),FDR=c(NA))
p<-1
for (i in 1:dim(drug_name)[1])
{
  temp<-GDSC_wm_score_drug[which(GDSC_wm_score_drug$Drug.name == DN[i]),]
  GDSC_wm_score_drug$AUC
  GDSC_wm_score_drug$MPIrisk_score
  res1<-corr.test(temp$AUC,temp$MPIrisk_score,method = "spearman",adjust = "fdr")
  mydata1[p,1]=DN[i]
  mydata1[p,2]='MPIrisk_score'
  mydata1[p,3]=res1$r
  mydata1[p,4]=res1$p
  mydata1[p,5]=res1$se
  p=p+1
}
mydata1$Class<-"ns"
for(i in 1:dim(mydata1)[1])
{
  if(mydata1$estimate.rho[i]>0.1 & mydata1$pvalue[i]<0.05)
  {
    mydata1$Class[i]<-"positive"
  }
  else if(mydata1$estimate.rho[i] < -0.1 & mydata1$pvalue[i]<0.05)
  {
    mydata1$Class[i]<-"negtive"
  }
}
table(mydata1$Class)

write.table(mydata1,'/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/04_GDSC/GDSC_LARrisk_scorecorr.txt',row.names=T,col.names=T,quote=F,sep="\t")

drug_wm<-mydata1[which(mydata1$Class != 'ns'),]
drug_wm<-drug_wm[order(drug_wm$estimate.rho),]
drug_wm$`-log10(pvalue)`<--log10(drug_wm$pvalue)


library(tidyverse)


pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/04_GDSC/drug_LARrisk_score.pdf',width = 10,height = 6)
drug_wm %>% 
  ggplot(aes(reorder(Drugs, estimate.rho), estimate.rho)) + 
  geom_col(aes(fill = `-log10(pvalue)`)) + 
  scale_fill_gradient2(low = "blue", 
                       high = "red", 
                       midpoint = 2) + 
  #coord_flip() + 
  labs(x = "")+
  labs(y = "Rs of drug sensitivity and LARrisk_score")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



drug_pathway<-read.csv('//Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/data/GDSC/PANCANCER_ANOVA_Mon Apr 26 17_21_47 2021.csv',header = T)

drug_pathway[,5:22]<-NULL
drug_pathway[,2:3]<-NULL
drug_pathway<-drug_pathway[which(drug_pathway$drug_name %in% drug_wm$Drugs),]

drug_pathway_use<-as.data.frame(table(drug_pathway))
drug_pathway_use<-drug_pathway_use[which(drug_pathway_use$Freq != 0),]
colnames(drug_pathway_use)<-c("Drugs","pathway","feq")

drug_pathway_use_final<-merge(drug_pathway_use,drug_wm,by='Drugs')
drug_pathway_use_final[,3:4]<-NULL
drug_pathway_use_final[,5:7]<-NULL
pathway<-as.data.frame(table(drug_pathway_use_final$pathway))
pathway_use<-as.character(pathway$Var1)



up<-data.frame(matrix(data = 0,nrow=14,ncol=20))
rownames(up)<-pathway_use
colnames(up)<-drug_pathway_use_final$Drugs

down<-data.frame(matrix(data = 0,nrow=14,ncol=20))
rownames(down)<-pathway_use
colnames(down)<-drug_pathway_use_final$Drugs



for(i in 1:20)
{
  for(j in 1:14)
  {
    for(p in 1:20)
    {
      if(rownames(up)[j] == as.character(drug_pathway_use_final[p,2]) & colnames(up)[i] == as.character(drug_pathway_use_final[p,1]))
      {
        up[j,i]<-(-log10(drug_pathway_use_final$pvalue[p]))
        down[j,i]<-drug_pathway_use_final$estimate.rho[p]
      }
    }
  }
}


library("ComplexHeatmap")
library("circlize")
plotMutiHeatmap1=function(up,down,up_break,up_colors,down_break,down_colors,title){
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
      if(down[i, j]>0.1){
      }
    }
  }
  
  p1 <- Heatmap(up, column_title = title
                , rect_gp = gpar(type = "none")
                , show_heatmap_legend = F
                , cluster_rows = T
                , cluster_columns = T, 
                cell_fun = DiagFunc(up = up, down = down) 
  ) 
  col_fun = colorRamp2(down_break, down_colors) 
  lgd <- Legend(title = "Correlation", 
                col_fun = col_fun, 
                at = c(-0.3,0,0.3), 
                labels = c("-0.3","0","0.3"),  
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
down_break=c(-0.3,0,0.3)
up_colors=c("#FFFFFF","#6f9a8d")
down_colors=c("blue",'white',"red")

up=up
down=down
title=''
pdf('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/04_GDSC/drug_pathway.pdf',width = 10,height = 10)
plotMutiHeatmap1(up,down,up_break,up_colors,down_break,down_colors,title)
dev.off()
###############################




###############################
save.image('/Volumes/Work/projects/BJ/small/active/20211129-GDT210901001/runtime_ver1/output/04_LAR/04_LAR.RData')










































































































