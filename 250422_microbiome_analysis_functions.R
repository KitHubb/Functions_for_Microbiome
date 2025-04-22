# Mon Mar 31 15:30:20 2025 ------------------------------
# volcano 보정 및 lfc 영역 설정 변수 추가 
# Mon Apr 21 16:52:37 2025 ------------------------------
  # volcano, mean abundance 순서 수정
library(phyloseq)
library(dplyr)
library(ggtext)
library(ggplot2)
library(glue)
library(ggprism)
library(ggrepel)
library(patchwork)
library(broom)

# source("./~") 이런식으로 사용 가능

#### 1. alpha diversity #### 


alpha_plot <- function(table, x_, y_ = "Shannon", colx_ = x_, ylim = 7, y_kruskal = 6.5,
                       col = c("#E31A1C", "#1F78B4", "#4D4D4D")){
  
  comp_ <- as.vector(unique(table[, x_]))
  Num <- length(comp_)
  
  
  
  
  if (Num >  3) {
    
    plot <- ggplot(table, aes_string(x =x_, y = y_)) +
      geom_violin(aes(fill =!!rlang::sym(x_), 
                      color = !!rlang::sym(x_)),
                  show.legend = FALSE, outlier.shape = NA, alpha=0.6) +
      
      geom_jitter(aes(color = !!rlang::sym(x_)), width = 0.15, alpha=0.3, size=2) +
      
      geom_boxplot(aes(color = !!rlang::sym(x_)),
                   show.legend = FALSE, outlier.shape = NA, width = 0.2 ) + 
      # geom_boxplot(aes_string(fill =colx_, color = colx_),
      #              show.legend = FALSE, outlier.shape = NA, alpha=0.6) +
      # geom_jitter(aes_string(color = colx_), width = 0.15, alpha=0.6, size=2) +
      scale_y_continuous(limits = c(0, ylim)) +
      scale_fill_manual(values =  col) +
      scale_color_manual(values =  col) +
      theme_test() +
      ylab(y_) + 
      theme(text = element_text(size = 10),
            axis.text.x.bottom = element_text(angle = 45,hjust = 1.1,vjust = 0.9),
            legend.title = element_blank(),
            strip.background = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank() ) +
      theme(# aspect.ratio=1,
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = NA, color = NA) )+
      theme(plot.margin=unit(c(0,0,0,0),"points"))+
      
      stat_compare_means(method = "kruskal.test",  tip.length=0.00 , # 0.05
                         size = 3,
                         label.y = y_kruskal, label.x =((Num+1)/2-0.4),
                         aes(label = paste0("kruskal.test, p = ", ..p.format..))
                         # aes(label = paste0("p = ",..p.format..))
      ) +
      geom_pwc(method = "wilcox_test", 
               label = "p = {ifelse(p < 0.001, '<0.001',p)}{p.signif}",label.size = 3,
               tip.length=0.00) 
    
    
  } else if (Num == 3) {
    com_pairs <- list(comp_[1:2], comp_[c(1, 3)], comp_[2:3])
    
    plot <- ggplot(table, aes_string(x =x_, y = y_)) +
      geom_violin(aes(fill =!!rlang::sym(x_), 
                      color = !!rlang::sym(x_)),
                  show.legend = FALSE, outlier.shape = NA, alpha=0.6) +
      
      geom_jitter(aes(color = !!rlang::sym(x_)), width = 0.15, alpha=0.3, size=2) +
      
      geom_boxplot(aes(color = !!rlang::sym(x_)),
                   show.legend = FALSE, outlier.shape = NA, width = 0.2 ) + 
      # geom_boxplot(aes_string(fill =x_, color = x_),
      #              show.legend = FALSE, outlier.shape = NA, alpha=0.6) +
      # geom_jitter(aes_string(color = x_), width = 0.15, alpha=0.6, size=2) +
      scale_y_continuous(limits = c(0, ylim)) +
      scale_fill_manual(values =  col) +
      scale_color_manual(values =  col) +
      theme_test() +
      ylab(y_) + 
      
      theme(text = element_text(size = 10),
            axis.text.x.bottom = element_text(angle = 45,hjust = 1.1,vjust = 0.9),
            legend.title = element_blank(),
            strip.background = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank() ) +
      theme(# aspect.ratio=1,
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = NA, color = NA) )+
      theme(plot.margin=unit(c(0,0,0,0),"points"))+
      
      stat_compare_means(method = "kruskal.test",  tip.length=0.01, size = 3,
                         label.y = y_kruskal, label.x =((Num+1)/2-0.4),
                         # label="p.signif"
                         aes(label = paste0("kruskal.test, p = ", ..p.format..))
      )+
      
      # https://github.com/kassambara/ggpubr/issues/327
      geom_pwc(method = "wilcox_test", 
               label = "p = {ifelse(p < 0.001, '<0.001',p)}{p.signif}",label.size = 3,
               tip.length=0.00) 
    
    # stat_compare_means(method = "wilcox.test", tip.length=0.00 , # 0.05
    #                    size = 3,
    #                    label.x = ((Num+1)/2-0.2),
    #                    aes(label = paste0("p = ",..p.format..)),
    #                    comparisons = com_pairs)
    
  } else if (Num == 2) {
    
    plot <- ggplot(table, aes_string(x =x_, y = y_)) +
      geom_violin(aes(fill =!!rlang::sym(x_), 
                      color = !!rlang::sym(x_)),
                  show.legend = FALSE, outlier.shape = NA, alpha=0.6) +
      
      geom_jitter(aes(color = !!rlang::sym(x_)), width = 0.15, alpha=0.3, size=2) +
      
      geom_boxplot(aes(color = !!rlang::sym(x_)),
                   show.legend = FALSE, outlier.shape = NA, width = 0.2 ) + 
      # geom_boxplot(aes_string(fill =x_, color = x_),
      #              show.legend = FALSE, outlier.shape = NA, alpha=0.6) +
      # geom_jitter(aes_string(color = x_), width = 0.15, alpha=0.6, size=2) +
      scale_y_continuous(limits = c(0, ylim)) +
      scale_fill_manual(values = col) +
      scale_color_manual(values = col) +
      theme_test() +
      ylab(y_) + 
      theme(text = element_text(size = 10),
            axis.text.x.bottom = element_text(angle = 45,hjust = 1.1,vjust = 0.9),
            legend.title = element_blank(),
            strip.background = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank() ) +
      theme(# aspect.ratio=1,
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = NA, color = NA) )+
      theme(plot.margin=unit(c(0,0,0,0),"points"))+
      geom_pwc(method = "wilcox_test", label = "p = {ifelse(p < 0.001, '<0.001',p)}{p.signif}", 
               tip.length=0.00)
    
    # stat_compare_means(method = "wilcox.test", tip.length=0.00 , # 0.05
    #                    size = 3,
    #                    label.x = 1.4,
    #                    aes(label = paste0("p = ",..p.format..)),
    #                    comparisons = list(comp_))
  }
  
  # result_df <- alpha_df(df = table,# alpha_aiversity = y_,
  #                       group = x_)
  
  
  return(plot) #, result = result_df
  
}


alpha_df <- function(df = df,
                     group = group 
                     # path = path
) {
  alpha_aiversity = c("Observed", "Chao1", "FaithPD",  "InvSimpson", "Simpson", "Shannon")
  
  # setting
  options(scipen = 5) # 숫자 요약에서 풀로 보기
  compar <- as.vector(unique(df[, group])) # 비교하고자 하는 그룹
  compar_num <- length(compar) # 비교하고자 하는 그룹 수 
  permut <- combinat::combn(compar, 2) %>% t() # 순열 
  p_val_df <- permut %>% as.data.frame()  # p-value값이 들어갈 data.frame
  num <- nrow(permut)   
  range <- seq(from = 1, to = num, by = 1) # for문을 위한 수열
  combined_df <- data.frame()
  
  for (index in alpha_aiversity) { # 각 alpha diversity index동안 
    
    # normal distribution : https://stackoverflow.com/questions/35953394/calculating-length-of-95-ci-using-dplyr
    result <- df %>%
      dplyr::group_by(!!rlang::sym(group)) %>%
      dplyr::summarize(
        index = index,
        Mean = mean(!!rlang::sym(index)),           # 평균
        N = n(),                                    # 행 개수
        Sd = sd(!!rlang::sym(index))) %>%           # 표준편차
      dplyr::mutate(se = Sd /sqrt(N),
                    lower = Mean - qnorm(0.975)*se, # 95% 신뢰 구간 하한
                    upper = Mean + qnorm(0.975)*se, # 95% 신뢰 구간 상한
                    CI95per = upper-lower           # 95% 신뢰 구간
      )
    
    assign(paste0(index, "_df"), result) 
    
    
    p_val_df$index <-  index
    
    if (compar_num>=3){ # 3군 이상일때 kruskal wallis 분석
      ###  kruskal wallis  
      krus_p <- kruskal.test(df[, index] ~ df[, group], data = df)$p.value
      
      # p_val_df["kruskal wallis—pval", index] <- krus_p
      p_val_df$Kruskal <- krus_p
    }
    
    for (rows in range)  {  # 2 군씩 나누어서 pairwise wilcoxon rank sun test
      
      ### wilcoxon 
      compar_var <- permut[rows, ]
      compar_df <- df %>% filter((!!rlang::sym(group)) %in% compar_var) 
      p.val<- wilcox.test(compar_df[, index] ~ compar_df[, group], data = compar_df)$p.value
      
      p_val_df[p_val_df$V1 %in% compar_var[1]& 
                 p_val_df$V2 %in% compar_var[2]&
                 p_val_df$index %in% index, "Wilcoxon" ]  <- p.val
      
      # Total_result[Total_result$index %in% index, "Stat"] <- "wilcoxon rank sum test"
      # Total_result[Total_result$index %in% index, "P.val"] <- p_val_df
      # 
    }
    p_val_df[, "Wilcoxon.adjust"] <- p.adjust(p_val_df[, "Wilcoxon"], method = "fdr")
    
    combined_df <- rbind(combined_df, p_val_df)
    
    
  }
  
  
  
  ### Save_mean value 
  Total_result <- rbind(Observed_df,
                        Chao1_df,
                        Shannon_df,
                        Simpson_df,
                        InvSimpson_df,
                        FaithPD_df )
  
  # write.csv(Total_result, paste0(path, "alpha_",group , "_mean_and_CI.csv"), col.names = T, row.names = T)
  
  #### Save_statistic value
  
  # print(p_val_df)
  # write.csv(p_val_df, paste0(path, "alpha_",group , "_StatisticTest.csv"), col.names = T, row.names = T)
  
  return( list(Total_result=Total_result, p_val_df=combined_df))
  
}

#### 2. Beta diversity #### 
#' Beta diversity analysis and visualization
#'
#' This function performs beta diversity analysis and visualization using phyloseq data.
#'
#' @param phyloseq A phyloseq object containing microbiome data.
#' @param type A character string indicating the variable for comparison.
#' @param shap A character string indicating the shape variable for plotting.
#' @param seed An integer specifying the random seed for reproducibility.
#' @param plot A character string specifying the ordination method ("PCoA" by default).
#' @param SampleID A character string specifying the sample ID column in sample data.
#' @param type_col A vector of colors for plotting.
#' @param col_inout A character string indicating the position of the legend ("out" by default).
#' @param col_right_left A character string indicating the justification of the legend ("right" by default).
#'
#' @return A list containing the total plot and individual plots for each beta diversity index.
#'
#' @examples
#' beta_plot4(phyloseq_obj, "Treatment", shap = NULL, seed = 42, plot = "PCoA", type_col = c("#1F78B4", "#4D4D4D"))
#'
#' @import ggplot2
#' @import ggtext
#' @import dplyr
#' @import phyloseq
#' @importFrom MASS lda
#' @importFrom gridExtra ggarrange
#' @importFrom ggpubr ggarrange




# ps.T, type = "Skin.site", type_col = c("#FA8334", "#94832F", "#8eb1c1")
# phyloseq = ps;type = "Skin.site"; shap = NULL;type_col = c("#1F78B4", "#4D4D4D") ;
# seed=42 ;plot="PCoA";
# index = "bray"


# 2024.04.01

# phyloseq = ITS1
# type = "Drug_type"
# type_col = col4
# plot="PCoA"
# SampleID = "SampleID"
# col_inout = c("out")
# col_right_left =c("right")
# indices =  c("bray", "jaccard",  "wunifrac", "unifrac")
# seed =42
# shap=NULL

beta_plot <- function(phyloseq, type, shap=NULL, seed=42, plot="PCoA", 
                      SampleID = "SampleID", type_col, 
                      col_inout = c("out"),
                      col_right_left =c("right"), 
                      indices =  c("bray", "jaccard", "unifrac", "wunifrac")) { 
  
  # Function to round p-values
  value <- function(val) {
    if (val > 0.05) {
      val2 = round(val, 3)
    } else if (val <= 0.05 & val > 0.001) {
      val2 = round(val, 3)
    } else if (val <= 0.001) {
      val2 = "<0.001"
    }
    return(val2)
  }
  
  
  # Function to extract variable name
  out_name <- function(v1) {
    deparse(substitute(v1))
  }
  
  
  plots <- list()
  result_list <- list()
  
  for (index in indices){
    
    set.seed(seed)
    x.dist <- phyloseq::distance(phyloseq, method = index)
    dist <- out_name(x.dist)
    meta <- phyloseq %>% sample_data %>% data.frame
    
    # PERMANOVA
    set.seed(seed)
    Perm <- adonis2(as.formula (glue("{dist} ~ {type}")), data=data.frame(sample_data(phyloseq)),
                    permutations=9999, method=index)
    Perm.p <- value(Perm$`Pr(>F)`[1])
    Perm.R2 = round(Perm$R2[1], 3)
    
    # ANOSIM
    set.seed(seed)
    Ans <- anosim(x.dist, data.frame(sample_data(phyloseq))[, type],  permutations = 9999)
    Ans.p <- value(Ans$signif[1])
    Ans.R <- round(Ans$statistic, 3)
    
    # PERMDISP (beta dispersion test)
    set.seed(seed)
    PMD <- anova(betadisper(x.dist, data.frame(sample_data(phyloseq))[, type]),  permutations = 9999)
    PMD.p <- value(PMD$`Pr(>F)`[1])
    PMD.f <- value(PMD$`F value`[1])
    
    
    # arrange result 
    result<- data.frame(
      Statistical.test = c("PERMANOVA", "ANOSIM", "PERMDISP"),
      DF = c(Perm$Df[1], NA, PMD$Df[1]),
      Sum.Sq = c(Perm$SumOfSqs[1], NA, PMD$`Sum Sq`[1]),
      Mean.Sq = c(NA, NA, PMD$`Mean Sq`[1]),
      F.Statist = c(Perm$F[1], NA, PMD$`F value`[1]),
      R.Squared = c(Perm.R2, Ans.R, NA),
      P.value = c(Perm$`Pr(>F)`[1], Ans$signif[1], PMD$`Pr(>F)`[1])
    )    
    result_list[[index]] <- result
    
    
    
    # Ordinate
    set.seed(seed)
    ord <- ordinate(phyloseq, plot, index)
    pcoa_df <- data.frame(meta, ord$vectors[, 1:2])
    PC1 <- round(ord$values["Relative_eig"][1,] * 100, 1)
    PC2 <- round(ord$values["Relative_eig"][2,] * 100, 1)

    Title <- switch(index,bray="Bray-curtis",
                    jaccard="Jaccard",
                    unifrac="Unweighted UniFrac",
                    wunifrac = "Weighted UniFrac")
    
    # 4 .draw beta plot 
    main.plot <- pcoa_df %>%
      ggplot(aes(x = Axis.1, y=Axis.2)) +
      geom_vline(xintercept = 0, colour = "grey80") +
      geom_hline(yintercept = 0, colour = "grey80") +
      geom_point(aes_string(shape = shap, color=type),alpha = 0.5, size=2) +
      stat_ellipse(aes_string(color= type) ) +
      scale_color_manual(values = type_col) +
      labs(y = paste0("PCoA2 (", PC2, "%)"), 
           x = paste0("PCoA1 (", PC1, "%)")
           ) +
  annotate("text", hjust = 0, vjust = -1, x = -Inf, y = -Inf, size = 2,
               label = paste0("    ", Title,
                              '\nPERMANOVA: R2 = ', Perm.R2, ', p-value = ', Perm.p,
                              '\nANOSIM: R = ', Ans.R, ', p-value = ', Ans.p,
                              '\nPERMADISP: F = ', PMD.f, ', p-value = ', PMD.p)) +
      theme_test() +
      theme(plot.caption = element_text(hjust = 0))+
      theme(plot.caption = element_markdown(),
            aspect.ratio = 1,
            legend.title = element_blank(),
            legend.margin = margin(6, 6, 6, 6),
            legend.background = element_rect(fill = NA, color = NA),
            plot.margin = unit(c(0, 0, 0, 0), "points"))
    
                 
    #            label= paste0(
    #              "    ",Title, '
    # PERMANOVA: R2 = ', Perm.R2,  ', p-value=', Perm.p ,
    #              '
    # ANOSIM: R = ', Ans.R,  ', p-value=', Ans.p,
    #              '
    # PERMDISP: F = ', PMD.f,  ', p-value=', PMD.p 
    #            ))
    # 
    #  "**PERMANOVA**  R<sup>2</sup>=",Perm.R2, ", *p*-value=", Perm.p, "<br>",
    #  "**ANOSIM**  R=",    Ans.R,   ", *p*-value=", Ans.p, "<br>",
    # "**ANOVA**  F=",     Anv.f,   ", *p*-value=", Anv.p ))
    
    # ggtext::geom_richtext(label.color = NA, size = 2, fill = NA,
    #                       hjust = 0, vjust =0, x = -Inf, y = -Inf,   # left bottom
    #                       #  hjust = 0, vjust =1, x = -Inf, y = Inf, # left upper 
    #                       label= paste0( "**", Title,"** ",  "<br>",
    #                                      "**PERMANOVA**  R<sup>2</sup>=",Perm.R2, ", *p*-value=", Perm.p, "<br>", 
    #                                      "**ANOSIM**  R=",    Ans.R,   ", *p*-value=", Ans.p, "<br>", 
    #                                      "**ANOVA**  F=",     Anv.f,   ", *p*-value=", Anv.p ) ) 
    
    
    if (col_inout == "in" & col_right_left == "right" ){
      main.plot <- main.plot + theme(
        legend.position = c(1, 1 ),
        legend.justification = c("right", "top"),
        legend.box.just = "right"
      )  
    } else if(col_inout == "in" & col_right_left == "left" ){
      main.plot <- main.plot + theme(
        legend.position = c(0, 1),
        legend.justification = c("left", "top"),
        legend.box.just = "left"
      )  
    }
    
    assign(paste0(index, "_p"), main.plot) 
    plots[[index]] <- main.plot
  }
  
  
  plot_list <- lapply(indices, function(index) {
    get(paste0(index, "_p")) # envir = parent.frame()
  })
  # ggarrange를 사용하여 plot을 배열
  plots[["Total"]] <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)
  # plots[["Total"]] <- ggarrange(bray_p, jaccard_p, unifrac_p, wunifrac_p, ncol = 2, nrow = 2)
  
  out <- list(plots = plots, results = result_list)
  return(out)
}


beta_plot.envfit <- function(phyloseq, type, shap = NULL, seed = 42, plot = "PCoA", 
                             SampleID = "SampleID", type_col, 
                             col_inout = c("out"),
                             col_right_left = c("right"),
                             indices = c("bray", "jaccard", "unifrac", "wunifrac"),
                             envfit = FALSE, PERMANOVA = FALSE,
                             formula = NULL) {
  
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(ggrepel)
  library(glue)
  library(ggpubr)
  library(dplyr)
  
  value <- function(val) {
    if (val > 0.05) {
      val2 = round(val, 3)
    } else if (val <= 0.05 & val > 0.001) {
      val2 = round(val, 3)
    } else if (val <= 0.001) {
      val2 = "<0.001"
    }
    return(val2)
  }
  
  plots <- list()
  result_list <- list()
  permanova_out <- list()
  
  for (index in indices) {
    set.seed(seed)
    x.dist <- phyloseq::distance(phyloseq, method = index)
    meta <- data.frame(sample_data(phyloseq))
    
    # PERMANOVA for pairwise group comparison (used for caption)
    Perm <- adonis2(as.formula(glue("x.dist ~ {type}")), data = meta, permutations = 9999, method = index)
    Perm.p <- value(Perm$`Pr(>F)`[1])
    Perm.R2 <- round(Perm$R2[1], 3)
    
    # ANOSIM
    Ans <- anosim(x.dist, meta[[type]], permutations = 9999)
    Ans.p <- value(Ans$signif[1])
    Ans.R <- round(Ans$statistic, 3)
    
    # PERMDISP (beta dispersion test)
    PMD <- anova(betadisper(x.dist, meta[[type]]), permutations = 9999)
    PMD.p <- value(PMD$`Pr(>F)`[1])
    PMD.f <- value(PMD$`F value`[1])
    
    result_list[[index]] <- data.frame(
      Statistical.test = c("PERMANOVA", "ANOSIM", "PERMDISP"),
      DF = c(Perm$Df[1], NA, PMD$Df[1]),
      Sum.Sq = c(Perm$SumOfSqs[1], NA, PMD$`Sum Sq`[1]),
      Mean.Sq = c(NA, NA, PMD$`Mean Sq`[1]),
      F.Statist = c(Perm$F[1], NA, PMD$`F value`[1]),
      R.Squared = c(Perm.R2, Ans.R, NA),
      P.value = c(Perm$`Pr(>F)`[1], Ans$signif[1], PMD$`Pr(>F)`[1])
    )
    
    ord <- ordinate(phyloseq, plot, index)
    pcoa_df <- data.frame(meta, ord$vectors[, 1:2])
    PC1 <- round(ord$values["Relative_eig"][1,] * 100, 1)
    PC2 <- round(ord$values["Relative_eig"][2,] * 100, 1)
    
    Title <- switch(index,
                    bray = "Bray-curtis",
                    jaccard = "Jaccard",
                    unifrac = "Unweighted UniFrac",
                    wunifrac = "Weighted UniFrac")
    
    main.plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2)) +
      geom_vline(xintercept = 0, colour = "grey80") +
      geom_hline(yintercept = 0, colour = "grey80") +
      geom_point(aes_string(shape = shap, color = type), alpha = 0.5, size = 2) +
      stat_ellipse(aes_string(color = type)) +
      scale_color_manual(values = type_col) +
      labs(y = paste0("PCoA2 (", PC2, "%)" ), x = paste0("PCoA1 (", PC1, "%)")) +
      annotate("text", hjust = 0, vjust = -1, x = -Inf, y = -Inf, size = 2,
               label = paste0("    ", Title,
                              '\n    PERMANOVA: R2 = ', Perm.R2, ', p-value = ', Perm.p,
                              '\n    ANOSIM: R = ', Ans.R, ', p-value = ', Ans.p,
                              '\n    PERMADISP: F = ', PMD.f, ', p-value = ', PMD.p)) +
      theme_test() +
      theme(plot.caption = element_text(hjust = 0))+
      theme(plot.caption = element_markdown(),
            aspect.ratio = 1,
            legend.title = element_blank(),
            legend.margin = margin(6, 6, 6, 6),
            legend.background = element_rect(fill = NA, color = NA),
            plot.margin = unit(c(0, 0, 0, 0), "points"))
    
    #### Envfit #### 
    if (envfit && !is.null(formula)) {
      envfir_matrix <- ord$vectors[, 1:2]
      
      set.seed(seed)
      envfit_result <- vegan::envfit(as.formula(glue("envfir_matrix ~ {formula}")),data = meta,
                                     permutations = 9999)
      
      # Numeric variable
      arrow_df <- as.data.frame(scores(envfit_result, display = "vectors"))
      arrow_df$Variable <- rownames(arrow_df)
      colnames(arrow_df)[1:2] <- c("Axis.1", "Axis.2")
      # Category variable
      centroid_df <- as.data.frame(scores(envfit_result, display = "factors"))
      centroid_df$Variable <- rownames(centroid_df)
      colnames(centroid_df)[1:2] <- c("Axis.1", "Axis.2")
      
      
      # Plot
      main.plot <- main.plot +
        # Line
        geom_segment(data = arrow_df, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
                     arrow = arrow(length = unit(0.2, "cm")), color = "blue", linewidth = 0.4) +
        geom_text_repel(data = arrow_df, aes(x = Axis.1, y = Axis.2, label = Variable), color = "blue", size = 3) +
        # Point
        geom_point(data = centroid_df, aes(x = Axis.1, y = Axis.2), shape = 17, size = 2.5, color = "black") +
        geom_text_repel(data = centroid_df, aes(x = Axis.1, y = Axis.2, label = Variable), size = 3, color = "black", fontface = "italic")
      
      result_list[[index]]$Envfit <- envfit_result
      
    }
    
    
    
    
    ## Plot 정렬
    if (col_inout == "in" & col_right_left == "right") {
      main.plot <- main.plot + theme(legend.position = c(1, 1), legend.justification = c("right", "top"), legend.box.just = "right")
    } else if (col_inout == "in" & col_right_left == "left") {
      main.plot <- main.plot + theme(legend.position = c(0, 1), legend.justification = c("left", "top"), legend.box.just = "left")
    }
    
    plots[[index]] <- main.plot
    
    
    #### Global PERMANOVA by terms ####
    if (PERMANOVA && !is.null(formula)) {
      full_formula <- as.formula(paste0("x.dist ~ ", formula))
      full_perm <- adonis2(full_formula, data = meta, permutations = 9999, method = index, by = "terms")
      perma.df <- data.frame(
        Group = rownames(full_perm)[1:(nrow(full_perm)-1)],
        R2 = full_perm$R2[1:(nrow(full_perm)-1)],
        P.val = full_perm$`Pr(>F)`[1:(nrow(full_perm)-1)]
      ) %>%
        mutate(
          Explained = R2 * 100,
          p.label = ifelse(P.val < 0.001, "p<0.001", paste0("p=", round(P.val, 3))),
          significance = ifelse(P.val < 0.05, "Significant", "Not Significant")
        )
      perma.df.forplot <- perma.df[perma.df$Group != "Residual", ]
      bar <- ggplot(perma.df.forplot, aes(x = reorder(Group, Explained), y = Explained, fill = significance)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = c("Significant" = "coral", "Not Significant" = "grey")) +
        geom_text(aes(label = p.label), hjust = -0.1, size = 4.2) +
        coord_flip() +
        scale_y_continuous(limits = c(0, max(perma.df.forplot$Explained) + 2), expand = c(0, 0)) +
        theme_minimal(base_size = 14) +
        labs(y = expression("PERMANOVA R"^2), x = NULL) +
        theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12),
              panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
              plot.margin = margin(10, 30, 10, 10))
      
      permanova_out[[index]] <- list(plot = bar, table = perma.df)
    }
  }
  
  plots[["Total"]] <- ggarrange(plotlist = plots[indices], ncol = 2, nrow = 2)
  
  return(list(plots = plots, results = result_list, permanova = permanova_out))
}
#### 3. Taxonomy diversity #### 
#####   Abund_cal ##### 
Abund_cal <- function(ps.glom, tax_level, group, path) {
  
  # melt
  melt <- psmelt(ps.glom)
  
  # setting
  Taxonomy <- as.vector(unique(melt[, tax_level]))
  meta <- data.frame(sample_data(ps.glom))
  meta_com <- meta[, group, drop = TRUE] %>% unique
  tax_num <- length(Taxonomy)
  
  # reset
  Total_result <- NULL
  
  ## Total abundance calculation
  Total_result <- melt %>%
    dplyr::group_by(!!rlang::sym(tax_level)) %>%
    dplyr::summarize(
      Total.Mean = mean(Abundance), # 평균
      Total.N    = n(),             # 행 개수
      Total.Sd   = sd(Abundance)    # 표준편차
    )
  
  ## Each group abundance calculation
  for (com in meta_com) {
    # 각 그룹별 계산
    melt.2 <- melt %>% filter(!!rlang::sym(group) == com)
    result <- melt.2 %>%
      dplyr::group_by(!!rlang::sym(tax_level)) %>%
      dplyr::summarize(
        !!paste0(com, ".Mean") := mean(Abundance),  # 평균
        !!paste0(com, ".N")    := n(),              # 행 개수
        !!paste0(com, ".Sd")  := sd(Abundance)      # 표준편차
      ) %>%   
      dplyr::mutate(!!paste0(com, ".se")      := !!rlang::sym(paste0(com, ".Sd"))    / sqrt(!!rlang::sym(paste0(com, ".N"))),           # 표준오차
                    !!paste0(com, ".lower")   := !!rlang::sym(paste0(com, ".Mean"))  - qnorm(0.975) * !!rlang::sym(paste0(com, ".se")), # 95% 신뢰 구간 하한
                    !!paste0(com, ".upper")   := !!rlang::sym(paste0(com, ".Mean"))  + qnorm(0.975) * !!rlang::sym(paste0(com, ".se")), # 95% 신뢰 구간 상한
                    !!paste0(com, ".CI95per") := !!rlang::sym(paste0(com, ".upper")) - !!rlang::sym(paste0(com, ".lower"))              # 95% 신뢰 구간
      )
    
    # abundance 결과 합치기
    Total_result <- bind_cols(Total_result, result[, -1])
  }
  print(Total_result)
  write.csv(Total_result, paste0(path, "Taxa_", group, "_stat_", tax_level, ".csv"), row.names = F)
}

##### 1) tax_process ##### 
# tax_process <- function(.melt., .taxa., .otu) {
#   
#   
#   # 1. select Top taxa ____________________________________________________________ 
#   Order = tapply(.melt.$Abundance, .melt.[, .taxa.], sum) %>% sort(decreasing = T)  # abundance 순으로 정렬 
#   
#   # 만약 Genus 에서 Bacteria Unclassified가 나오면 Other로 처리하기
#   # Top_g <- names(Order[c(1:.top.)])
#   Names <- names(Order[c(1:.top.)])
#   if ("Bacteria Unclassified" %in% Names ) {
#     Names.2 <- names(Order[c(1:(.top.+1))])
#     Names.2 <- Names.2[Names.2 !="Bacteria Unclassified"]
#     Top_g <- Names.2
#   } else (Top_g = Names) 
#   
#   
#   Top_p <- .melt.[.melt.[, .taxa.] %in% Top_g, "Phylum"]%>% unique()  
#   
#   p_tax_table <- .melt.[.melt.[, .taxa.] %in% Top_g, c("Phylum", .taxa.)]  %>% 
#     .[!duplicated(.[ , .taxa.]),] 
#   
#   
#   # 2. define Other _____________________________________________________________
#   .melt.2 <- .melt. # Back up
#   
#   # Not Top15 Phylum <-  Other
#   .melt.2[!.melt.2[, "Phylum"] %in% Top_p, "Phylum"] <- "Others"
#   .melt.2[!.melt.2[, "Phylum"] %in% Top_p, .taxa.] <- "Others"
#   
#   
#   
#   if (.taxa. != "Species"){
#     for ( i in Top_p ) {
#       G <- p_tax_table[p_tax_table[, "Phylum" ] ==  i , .taxa.]
#       .melt.2[.melt.2[, "Phylum" ] == i & !.melt.2[, .taxa.] %in% G, .taxa. ] <- paste0(i, "_Others")
#     }
#     for ( i in Top_p ) {
#       G <- p_tax_table[p_tax_table[, "Phylum"] ==  i , .taxa.]
#       for (g in G){
#         .melt.2[.melt.2[, .taxa.] == g, .taxa.] <- paste0(i,"_", g)
#       }
#     }
#   } else {
#     for ( i in Top_p ) {
#       G <- p_tax_table[p_tax_table[, "Phylum" ] ==  i , .taxa.]
#       .melt.2[.melt.2[, "Phylum" ] == i & !.melt.2[, .taxa.] %in% G, .taxa. ] <- paste0("Others")
#       .melt.2[.melt.2[, "Phylum" ] == i & !.melt.2[, .taxa.] %in% G, "Phylum" ] <- paste0("Others")
#     }
#     
#   }
#   
#   # check
#   .melt.2[.melt.2[, "Phylum"] %in% Top_p ,.taxa.] %>% unique() 
#   
#   
#   # 3. order _____________________________________________________________________
#   .melt.3 <- .melt.2
#   table <- .melt.3[.melt.3[, "Phylum" ] %in% Top_p, c("Abundance", "Phylum", .taxa.)] 
#   
#   
#   # Phylum Order  ############################################################## 수정 부분 
#   p_order <- table %>%  .[,"Phylum" ]%>% unique
#   if("Firmicutes" %in% p_order) {
#     P_levels <- c("Firmicutes", as.vector(p_order[p_order != "Firmicutes"]), "Others")
#   } else {
#     P_levels <- p_order %>% sort()
#   }
#   
#   .melt.3[,"Phylum"]  <- factor(.melt.3[,"Phylum" ],  levels = P_levels)
#   
#   # Genus order 
#   # Genus order  by phylum and Abundance 
#   table.2 <- table %>% dplyr::group_by(Phylum, !!rlang::sym(.taxa.)) %>%
#     dplyr::summarise(sum.Abundance=sum(Abundance), .groups = 'drop') %>%
#     as.data.frame()
#   
#   
#   g_order <- table.2 %>% 
#     dplyr::arrange( -sum.Abundance) %>% 
#     dplyr::arrange(Phylum) %>%
#     .[,.taxa. ]
#   
#   
#   .melt.3[,.taxa.] <- factor(.melt.3[,.taxa.], levels = c(g_order, "Others"))
#   # check
#   .melt.3[,.taxa.]%>% levels()
#   .melt.4 <- .melt.3
#   
#   
#   # 6. palette ___________________________________________________________________
#   table.3 <-table.2 %>% arrange( -sum.Abundance) %>% 
#     arrange(Phylum) %>% 
#     select(Phylum, !!rlang::sym(.taxa.))
#   df <- table.3
#   
#   # Find how many colour categories to create and the number of colours in each
#   categories <- aggregate(as.formula(paste(.taxa., "Phylum", sep="~" )), df, 
#                           function(x) length(unique(x)))
#   
#   
#   # Create a color-coded print list
#   color_list.names  <- categories[, "Phylum"]
#   color_list <- vector("list", length(color_list.names)) 
#   names(color_list) <- color_list.names
#   
#   
#   '%!in%' <- function(x,y)!('%in%'(x,y))
#   
#   
#   P <- categories$Phylum
#   basic_p <-c("Bacteroidetes", "Proteobacteria", "Firmicutes", "Fusobacteria", "Actinobacteria")
#   
#   # main color 
#   
#   Blues_c   <-c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B") 
#   Reds_c    <-c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
#   Purples_c <-c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D")
#   Greens_c  <-c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")
#   YlOrBr_c  <-c("#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506")
#   
#   
#   
#   # for 문 
#   for (i in P) {
#     if (i  == "Firmicutes")     { 
#       Fi_num <- categories[categories[, "Phylum"] %in% c("Firmicutes"), .taxa.]
#       
#       ifelse(Fi_num == 1,              Firmicutes_color <- rev(brewer.pal(9, "Blues")[5]),
#              ifelse(Fi_num == 2,              Firmicutes_color <- rev(brewer.pal(9, "Blues")[c(3, 7)]),
#                     ifelse(Fi_num >=3 & Fi_num <= 9, Firmicutes_color <- rev(brewer.pal(Fi_num, "Blues")),
#                            Firmicutes_color <- colorRampPalette(Blues_c)(Fi_num)))) 
#       color_list$Firmicutes <- Firmicutes_color
#       
#     } else if (any(i == "Actinobacteria"))    { 
#       Ac_num <- categories[categories[, "Phylum"]  %in% c("Actinobacteria"), .taxa.]
#       
#       ifelse(Ac_num == 1,              Actinobacteria_color <- rev(brewer.pal(9, "Reds")[5]),
#              ifelse(Ac_num == 2,              Actinobacteria_color <- rev(brewer.pal(9, "Reds")[c(2, 7)]),
#                     ifelse(Ac_num >=3 & Ac_num <= 9, Actinobacteria_color <- rev(brewer.pal(Ac_num, "Reds")),
#                            Actinobacteria_color <- colorRampPalette(Reds_c)(Ac_num)))) 
#       color_list$Actinobacteria <- Actinobacteria_color
#       
#       
#       
#     } else if (i == "Bacteroidetes")   { 
#       Ba_num <- categories[categories[, "Phylum"]  %in% c( "Bacteroidetes"), .taxa.]
#       
#       ifelse(Ba_num == 1,              Bacteroidetes_color <- rev(brewer.pal(9, "Purples")[5]),
#              ifelse(Ba_num == 2,              Bacteroidetes_color <- rev(brewer.pal(9, "Purples")[c(3, 7)]),
#                     ifelse(Ba_num >=3 & Ba_num <= 9, Bacteroidetes_color <- rev(brewer.pal(Ba_num, "Purples")),
#                            Bacteroidetes_color <- colorRampPalette(Purples_c)(Ba_num))))  
#       color_list$Bacteroidetes <- Bacteroidetes_color
#       
#     } else if (i  == "Proteobacteria") { 
#       Pr_num <- categories[categories[, "Phylum"]  %in% c( "Proteobacteria"), .taxa.]
#       
#       ifelse(Pr_num == 1,              Proteobacteria_color <- rev(brewer.pal(9, "Greens")[5]),
#              ifelse(Pr_num == 2,              Proteobacteria_color <- rev(brewer.pal(9, "Greens")[c(3, 7)]),
#                     ifelse(Pr_num >=3 & Pr_num <= 9, Proteobacteria_color <- rev(brewer.pal(Pr_num, "Greens")),
#                            Proteobacteria_color <- colorRampPalette(Greens_c)(Pr_num))))  
#       color_list$Proteobacteria <- Proteobacteria_color
#       
#       
#     } else if (i  == "Fusobacteria")   { 
#       Fu_num <- categories[categories[, "Phylum"]  %in% c( "Fusobacteria"), .taxa.]
#       
#       
#       ifelse(Fu_num == 1,              Fusobacteria_color <- rev(brewer.pal(9, "YlOrBr")[5]),
#              ifelse(Fu_num == 2,              Fusobacteria_color <- rev(brewer.pal(9, "YlOrBr")[c(3, 7)]),
#                     ifelse(Fu_num >=3 & Fu_num <= 9, Fusobacteria_color <- rev(brewer.pal(Fu_num, "YlOrBr")),
#                            Fusobacteria_color <- colorRampPalette(YlOrBr_c)(Fu_num)))) 
#       color_list$Fusobacteria <- Fusobacteria_color
#       
#     } else if (i %!in% basic_p)          { 
#       other_p <- categories[categories$Phylum %!in% c("Bacteroidetes", "Proteobacteria", "Firmicutes", 
#                                                       "Fusobacteria", "Actinobacteria"), "Phylum"] 
#       
#       
#       # others_colors <- c("BrBG", "RdPu",  "YlOrBr","Greens","Blues","Reds", "Purples",   "Spectral", "PuOr",   "YlOrBr")
#       RdPu_c <-    c( "#FFF7F3","#FA9FB5" , "#DD3497" ,"#7A0177", "purple")
#       YlOrBr_c <-  c( "#FFFFE5", "#FFF7BC",  "#FEE391",  "#FEC44F",  "#FE9929",  "#EC7014", "#CC4C02", "#993404", "#662506")
#       BrBG_c <-    c( "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E" )
#       Blues_c   <-c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B") 
#       Reds_c    <-c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
#       PuOr_c <-    c("#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788")
#       Spectral_c <-c( "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD") # c("red","yellow","springgreen","royalblue")
#       Purples_c <-c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D")
#       Greens_c  <-c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")
#       YlOrBr_c  <-c("#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#993404", "#662506")
#       
#       others_colors2 <- list(RdPu_c,  YlOrBr_c, Greens_c,Reds_c, Blues_c,  PuOr_c, YlOrBr_c, BrBG_c, PuOr_c,  Purples_c,  Spectral_c)
#       # for ( i in 1:length(other_p) ) {
#       #   p <- other_p[i]
#       #   color <-  others_colors[i]
#       #   num <- categories[categories$Phylum == p, .taxa.]
#       #   ifelse(num == 1,              col <- rev(brewer.pal(9, color)[5]),
#       #          ifelse(num == 2,              col <- rev(brewer.pal(9, color)[c(3, 7)]),
#       #                 ifelse(num >=3 & num <= 9,    col <- rev(brewer.pal(num, color)),
#       #                        stop("Number of Other Phylum is more than 9, Use manual way.")                       
#       #                 )))
#       #   name <- paste0(p, "_color")
#       #   assign(name, col)
#       #   color_list[[p]] <- get(name)
#       # }
#       
#       for ( i in 1:length(other_p) ) {
#         p <- other_p[i]
#         color2 <-  others_colors2[i]
#         num <- categories[categories$Phylum == p, .taxa.]
#         ifelse(num == 1, col <- rev(colorRampPalette(color2[[1]])(9)[5]),
#                ifelse(num >= 2, col <- rev(colorRampPalette(color2[[1]])(num)),
#                       col <- rev(colorRampPalette(color2[[1]])(num))
#                ))
#         name <- paste0(p, "_color")
#         assign(name, col)
#         color_list[[p]] <- get(name)
#       }
#       
#     }
#   }
#   name <- names(color_list)
#   `%nin%` = Negate(`%in%`)
#   
#   color_vector <- color_list %>% unlist %>% unname
#   Final_col <- c(color_vector, "#D3D3D3")
#   
#   # result <- list(Final_color=Final_col, data=.melt.4)
#   return(list(Final_color=Final_col, data=.melt.4))
# }
# ##### 2) taxa_plot ##### 
# tax_process_Genus <- function(melt, Genus = "Staphylococcus") {
#   # 1. select Top taxa ____________________________________________________________
#   
#   Sp_df <- melt[melt$Genus %in% Genus, "Species"] %>% unique
#   sp_p <-  melt[melt$Genus %in% Genus, "Phylum"] %>% unique
#   # 2. define Other _____________________________________________________________
#   melt.2 <- melt # Back up
#   
#   # Not Top15 Phylum <-  Other
#   melt.2[!melt.2[, "Phylum"] %in% sp_p, "Phylum"] <- "Others"
#   melt.2[!melt.2[, "Genus"] %in% Genus, "Genus"] <- "Others"
#   melt.2[!melt.2[, "Genus"] %in% Genus, "Species"]  <- "Others"
#   
#   # melt.2$Species <- factor(melt.2$Species , levels=c("Others", Sp_df))
#   
#   n <- length(Sp_df)
#   
#   max_color <- c(brewer.pal(12, "Paired"), 
#                  brewer.pal(12, "Set3"), 
#                  brewer.pal(8, "Accent")) 
#   
#   if(n == 1){color <- rev(brewer.pal(9, "Blues")[5])
#   }else if(n == 2) {
#     color <- rev(brewer.pal(9, "Blues")[c(3, 7)])
#   }else if (n >=3 & n <= 9) {
#     color <- rev(brewer.pal(n, "Blues"))
#   }else { 
#     color <- max_color
#   }
#   # color <- terrain.colors(n)
#   
#   
#   Final_col <- c( "#D3D3D3", color)
#   # re <- list(Final_color=Final_col, data=melt.2)
#   return(list(Final_color=Final_col, data=melt.2))
# }

# taxa_plot <- function(tax_data, level = "Species", x_ = "SampleID", y_  = "Abundance"){
#   plot <- ggplot(data=tax_data$data, aes_string(x=x_, y=y_, fill=level)) +
#     geom_bar(aes(), position="fill", stat="identity") +
#     labs(y = "Relative Abundance", x = "") + 
#     theme(axis.text.x.bottom = element_text(angle = 45,hjust = 1.1,vjust = 0.9),
#           axis.text = element_text(size = 10),
#           legend.text = element_text(size = 10),
#           legend.position = "right",
#           legend.background = element_rect(fill='transparent'),
#           legend.box.background = element_rect(fill='transparent', color=NA),
#           plot.background = element_rect(fill='transparent', color=NA),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank()) +
#     if(level == "Species"){
#       scale_fill_manual(name =level, values = tax_data$Final_color)
#     } else {scale_fill_manual(name =paste0('Phylum - ', level), values = tax_data$Final_color)}
#   
#   return(plot)  
# }



getwd()


taxa_plot <- function(melt, taxa, tax_otu, x_axis, phylum_or = NULL){ # 2024 12 11
  
  #### Function ####
  # F1.process_data = function(melt, tax_level,  tax_otu) {
  #   
  #   
  #   tax_tab <-   melt[, c("OTU", "Phylum", tax_level)] %>% unique
  #   # subset taxa by OTU ID
  #   tax_tab2 <- tax_tab[tax_tab$OTU %in% tax_otu, ]
  #   tax_phylum <- tax_tab2$Phylum %>% unique
  #   tax_index <- tax_tab2[,  tax_level,  drop=T]
  #   # Reconstruct taxa classified as Others
  #   melt.2 <- melt # Back up
  #   melt.2[!melt.2[, "Phylum"] %in% tax_phylum, "Phylum"] <- "Others"
  #   melt.2[!melt.2[, "Phylum"] %in% tax_phylum, tax_level] <- "Others"
  #   
  #   out=list(melt = melt.2,  tax_table=tax_tab2,  Top_p= tax_phylum,  Top_taxa=tax_index)
  #   return(out)
  # }
  
  F1.process_data = function(melt , taxa, tax_otu) {
    
    tax_tab <-   melt[, c("OTU", "Phylum", taxa)] %>% unique
    # tax_tab
    tax_tab2 <- tax_tab[tax_tab$OTU %in% tax_otu, ]
    tax_phylum <- tax_tab2$Phylum %>% unique
    tax_index <- tax_tab2[,  taxa,  drop=T]
    
    # Reconstruct taxa classified as Others
    melt.2 <- melt # Back up
    melt.2[!melt.2[, "Phylum"] %in% tax_phylum, "Phylum"] <- "Other"
    melt.2[!melt.2[, "Phylum"] %in% tax_phylum, taxa] <- "Other"
    
    
    # Genus 와 Phylum정렬
    if (taxa != "Species") {
      for (i in tax_phylum) {
        G <-tax_tab2[tax_tab2[, "Phylum"] == i, taxa]
        
        melt.2[melt.2[, "Phylum"] == i, taxa]
        
        
        
        melt.2[melt.2[, "Phylum"] == i & !melt.2[, taxa] %in% G, taxa] <- paste0(i, "_Other")
      }
      for (i in tax_phylum) {
        G <- tax_tab2[tax_tab2[, "Phylum"] == i, taxa]
        for (g in G) {
          melt.2[melt.2[, taxa] == g, taxa] <- paste0(i, "_", g)
        }
      }
    } else {
      for (i in tax_phylum) {
        G <- tax_tab2[tax_tab2[, "Phylum"] == i, taxa]
        melt.2[melt.2[, "Phylum"] == i & !melt.2[, taxa] %in% G, taxa] <- "Other"
        melt.2[melt.2[, "Phylum"] == i & !melt.2[, taxa] %in% G, "Phylum"] <- "Other"
      }
    }
    return(list(df = melt.2,
                Phylum_list = tax_phylum))
    
  }
  # F2.Order_data
  F2.Order_data = function(processed_data, Top_p, taxa){
    # phylum level 
    table <- processed_data[processed_data[, "Phylum" ] %in% Top_p, c("Abundance", "Phylum", taxa)]
    p_order <- table %>%  .[,"Phylum" ]%>% unique
    
    processed_data[,"Phylum"]  <- factor(processed_data[,"Phylum" ],  levels = c(sort(p_order), "Other"))
    
    # Genus order
    table_2 <- table %>% 
      dplyr::group_by(Phylum, !!rlang::sym(taxa)) %>%
      dplyr::summarise(sum.Abundance=sum(Abundance), .groups = 'drop') %>%
      dplyr::arrange( -sum.Abundance) %>%
      ungroup() %>% 
      as.data.frame()
    
    g_order <- table_2 %>% 
      dplyr::arrange(Phylum) %>%
      select(all_of(taxa))  %>% .[[1]]
    
    processed_data[ ,taxa] <- factor(processed_data[,taxa], levels = c(g_order, "Other"))
    return(list(df = processed_data, summary_df = table_2))
    
  }
  # get_palette_colors
  get_palette_colors = function(palette_name, num_taxa) {
    if (num_taxa == 1) {
      colors <- rev(brewer.pal(9, palette_name)[5])
    } else if (num_taxa == 2) {
      colors <- rev(brewer.pal(9, palette_name)[c(3, 7)])
    } else if (num_taxa >= 3 & num_taxa <= 9) {
      colors <- rev(brewer.pal(num_taxa, palette_name))
    } else {
      color_list <- rev(brewer.pal(9, palette_name))
      colors <- colorRampPalette(color_list)(num_taxa)
    }
    return(colors)
  } 
  # F3.generate_colors
  F3.generate_colors = function(df, taxa) {
    
    ## arrange by abundance and phylum
    table_3 <- df %>%
      dplyr::arrange(-sum.Abundance) %>%
      dplyr::arrange(Phylum) %>%
      dplyr::select(Phylum, !!rlang::sym(taxa))
    
    ## count taxa  (2024.09.24)
    categories <- table_3 %>% 
      group_by(Phylum) %>% 
      summarise(Taxa = n())
    colnames(categories)[2] <- taxa
    # categories <- aggregate(as.formula(paste(taxa, "Phylum", sep = "~")), 
    #                         table_3, 
    #                         function(x) length(unique(x))) 
    
    ## First phylum 
    # P_levels <- table_3$Phylum %>% unique
    
    ## color list 
    color_list.names <- categories[, "Phylum"]
    color_list <- vector("list", length(color_list.names))
    names(color_list) <- color_list.names
    
    ## 
    phylum_color_map <- list(
      # Bacteria
      Actinobacteria = "Reds",
      Actinobacteriota = "Reds",
      Actinomycetota = "Reds",
      
      Firmicutes = "Blues",
      Bacillota = "Blues",
      Firmicutes_A = "Blues",
      Firmicutes_B = "BuGn",
      Firmicutes_C = "PuBu",
      Firmicutes_D = "YlGnBu", 
      
      Bacteroidetes = "Purples",
      Bacteroidota = "Purples",
      
      Proteobacteria = "Greens",
      Pseudomonadota = "Greens",
      
      Fusobacteria = "YlOrBr",
      Fusobacteriota = "YlOrBr",
      
      # Fungi 
      Ascomycota = "RdPu",
      Basidiomycota = "YlOrBr"
      
      
      
      
      
    )
    
    basic_p <- names(phylum_color_map)
    other_colors <- c("BrBG", "Spectral", "BrBG", "PuOr"# "RdPu", "YlOrBr",
    )
    other_p <- categories$Phylum[!categories$Phylum %in% basic_p]
    In_p <- categories$Phylum[categories$Phylum %in% basic_p]
    
    phylum_color_map2 <- phylum_color_map[In_p]
    if (is.na(other_p[1])){
      phylum_color_map2 <- c(phylum_color_map2)
    } else {
      other_color_map <- setNames(other_colors[1:length(other_p)], other_p)
      phylum_color_map2 <- c(phylum_color_map2, other_color_map)  
      
    }
    
    for (phylum in names(phylum_color_map2)) {
      num_taxa <- categories[categories$Phylum == phylum, taxa]
      palette_name <- phylum_color_map2[[phylum]]
      
      if (length(num_taxa) > 0) {
        color_list[[phylum]] <- get_palette_colors(palette_name, as.numeric(num_taxa))
      }
    }
    
    
    color_vector <- unlist(color_list, use.names = FALSE)
    final_color <- c(color_vector, "#D3D3D3")
    return(final_color)
  }
  
  F5.taxa_plot <- function(df, color, taxa, x_axis){
    p <- ggplot(df, aes(x = !!rlang::sym(x_axis), y = Abundance, fill = !!rlang::sym(taxa))) +
      geom_bar(stat = "identity", position="fill") +
      labs(y = "Relative Abundance") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),  # x축 라벨 각도 조정
            plot.title = element_text(hjust = 0.5),  # 제목 가운데 정렬
            legend.position = "right",  # 범례를 하단에 위치
            legend.title = element_blank()) +  # 범례 제목 제거
      scale_fill_manual(values = color)   # 색상 팔레트 변경
    
    print(p)
  }
  
  F6.sampleID_order <- function(df, phylum) {
    sample_order <- df%>% 
      data.frame() %>%
      dplyr::group_by(!!rlang::sym(x_axis)) %>% 
      mutate(Abundance = Abundance / sum(Abundance)) %>%
      filter(Phylum == phylum_or) %>% 
      dplyr::group_by(!!rlang::sym(x_axis)) %>% 
      dplyr::summarise(Abundance = sum(Abundance)) %>% 
      dplyr::arrange(Abundance) %>%
      pull(!!rlang::sym(x_axis)) %>% as.character()
    
    df[, x_axis] <- factor(df[, x_axis], levels = sample_order)
    
    return(df)
  }
  
  
  
  
  #### Run ####
  out2 <- F1.process_data(melt = melt, taxa = taxa, tax_otu = tax_otu)
  
  out3 <- F2.Order_data(processed_data = out2$df,Top_p =  out2$Phylum_list, taxa = taxa)
  color_code <- F3.generate_colors(df = out3$summary_df, taxa=taxa )
  # return(list(out = out3$df, colors = color_code))  
  plot <- F5.taxa_plot(out3$df, color_code, taxa, x_axis)
  
  if (!is.null(phylum_or)){
    out4 <- out3
    out4$df <- F6.sampleID_order(out3$df, phylum_or)
    plot <- F5.taxa_plot(out4$df, color_code, taxa, x_axis)
    
  } 
  return(out = list(plot = plot, 
                    data = out3,
                    color = color_code))
  
  
}
##### 2-2) Species taxa plot 
tax_plot.sp <- function(melt, taxa, tax_otu, x_axis, phylum_or = NULL){ 
  # Tue Apr 22 10:49:40 2025 ------------------------------
  
  
  F1.Process_data <- function(melt, taxa, otu) {
    
    
    # pick taxa
    tax.summary <- melt[melt$OTU %in% otu, 
                        c("OTU", "Phylum", "Genus", "Species")] %>% unique
    
    tax_phylum  <- tax.summary$Phylum %>% unique
    tax_genus   <- tax.summary$Genus %>% unique
    tax_species <- tax.summary[,  taxa,  drop=T]
    
    ## 1. 일단 otu list에 없는 Phylum이나 Genus를 other로 치환 _________________
    melt.2 <- melt # Back up
    
    melt.2[!melt.2$Genus %in% tax_genus, "Genus"] <- "Other"
    melt.2[!melt.2$Genus %in% tax_genus, taxa] <- "Other"
    
    ## 2. Genus 들 중에 species level에서 otu list와 일치하지 않는 친구 -> Genus_other(Species 열에 추가)
    for (i in tax_genus) {
      G <-tax.summary[tax.summary[, "Genus"] == i, taxa]
      
      melt.2[melt.2[, "Genus"] == i & !melt.2[, taxa] %in% G, taxa] <- paste0(i, "_Other")
    }
    
    # check
    # melt.2$Species %>% unique
    
    ## 3.Phylum도 Genus list 에 없는 친구들 Other로 치환하기  __________________
    melt.2[!melt.2$Phylum %in% tax_phylum, "Phylum"] <- "Other"
    melt.2[!melt.2$Genus %in% tax_genus, "Phylum"] <- "Other"
    
    
    
    ## 4. 순서 지정하기 (chr to foctor)  _______________________________________
    tax.summary2 <- melt.2[, c("Phylum", "Genus", taxa)] %>% unique
    
    # Species level
    tax.summary3  <- tax.summary2 %>%
      dplyr::mutate( is_other = Species == paste0(Genus, "_other")) %>%
      dplyr::arrange( Phylum, Genus, desc(is_other), Species) %>%
      dplyr::select(-is_other)  # 최종적으로는 is_other 컬럼 제거
    
    Species_order <- tax.summary3 %>% pull(Species)
    
    # Species level
    melt.3 <- melt.2
    melt.3[ ,taxa] <- factor(melt.3[,taxa], levels = c(Species_order))
    # melt.3[ ,taxa] %>% unique
    
    
    return(out = list(df=melt.3, tax_summary = tax.summary3) )
    
  }
  
  
  
  #############################################################################
  # get_palette_colors
  get_palette_colors  <- function(palette_name, num_taxa) {
    if (num_taxa == 1) {
      colors <- rev(brewer.pal(9, palette_name)[5])
    } else if (num_taxa == 2) {
      colors <- rev(brewer.pal(9, palette_name)[c(3, 7)])
    } else if (num_taxa >= 3 & num_taxa <= 9) {
      colors <- rev(brewer.pal(num_taxa, palette_name))
    } else {
      color_list <- rev(brewer.pal(9, palette_name))
      colors <- colorRampPalette(color_list)(num_taxa)
    }
    return(colors)
  }
  
  
  # F3.generate_colors
  F2.generate_colors = function(tax_summary, taxa) {
    ## count taxa  (2024.09.24)
    categories <- tax_summary[!tax_summary$Phylum %in%"Other", ] %>%
      group_by(Phylum) %>%
      summarise(Taxa = n())
    colnames(categories)[2] <- taxa
    
    ## color list
    color_list.names <- categories[["Phylum"]]
    color_list <- vector("list", length(color_list.names))
    names(color_list) <- color_list.names
    ##
    phylum_color_map <- list(
      # Bacteria
      Actinobacteria = "Reds", Actinobacteriota = "Reds", Actinomycetota = "Reds", 
      Firmicutes = "Blues", Bacillota = "Blues", Firmicutes_A = "Blues", 
      Firmicutes_B = "BuGn", Firmicutes_C = "PuBu", Firmicutes_D = "YlGnBu", 
      Bacteroidetes = "Purples", Bacteroidota = "Purples", Proteobacteria = "Greens", 
      Pseudomonadota = "Greens", Fusobacteria = "YlOrBr", Fusobacteriota = "YlOrBr",
      # Fungi
      Ascomycota = "RdPu",  Basidiomycota = "YlOrBr"
    )
    
    basic_p <- names(phylum_color_map)
    other_colors <- c("BrBG", "Spectral", "BrBG", "PuOr", "RdPu", "YlOrBr")
    other_p <- categories$Phylum[!categories$Phylum %in% basic_p]
    In_p <- categories$Phylum[categories$Phylum %in% basic_p]
    
    phylum_color_map2 <- phylum_color_map[In_p]
    if (is.na(other_p[1])){
      phylum_color_map2 <- c(phylum_color_map2)
    } else {
      other_color_map <- setNames(other_colors[1:length(other_p)], other_p)
      phylum_color_map2 <- c(phylum_color_map2, other_color_map)
      
    }
    
    for (phylum in names(phylum_color_map2)) {
      num_taxa <- categories[categories$Phylum == phylum, taxa]
      palette_name <- phylum_color_map2[[phylum]]
      
      if (length(num_taxa) > 0) {
        color_list[[phylum]] <- get_palette_colors(palette_name, as.numeric(num_taxa))
      }
    }
    
    
    color_vector <- unlist(color_list, use.names = FALSE)
    final_color <- c(color_vector, "#D3D3D3")
    return(final_color)
  }
  
  ##############################################################################
  F3.taxa_plot <- function(df, color, taxa, x_axis){
    p <- ggplot(df, aes(x = !!rlang::sym(x_axis), y = Abundance, fill = !!rlang::sym(taxa))) +
      geom_bar(stat = "identity", position="fill") +
      labs(y = "Relative Abundance") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),  # x축 라벨 각도 조정
            plot.title = element_text(hjust = 0.5),  # 제목 가운데 정렬
            legend.position = "right",  # 범례를 하단에 위치
            legend.title = element_blank()) +  # 범례 제목 제거
      scale_fill_manual(values = color)   # 색상 팔레트 변경
    
    print(p)
  }
  
  sampleID_order <- function(df, phylum) {
    sample_order <- df%>%
      data.frame() %>%
      dplyr::group_by(!!rlang::sym(x_axis)) %>%
      mutate(Abundance = Abundance / sum(Abundance)) %>%
      filter(Phylum == phylum_or) %>%
      dplyr::group_by(!!rlang::sym(x_axis)) %>%
      dplyr::summarise(Abundance = sum(Abundance)) %>%
      dplyr::arrange(Abundance) %>%
      pull(!!rlang::sym(x_axis)) %>% as.character()
    
    df[, x_axis] <- factor(df[, x_axis], levels = sample_order)
    
    return(df)
  }
  
  #### Run ####
  out <- F1.Process_data(melt = melt, taxa = taxa, otu = tax_otu)
  
  color_code <- F2.generate_colors(tax_summary = out$tax_summary, taxa = taxa)
  
  plot <- F3.taxa_plot(out$df, color_code, taxa, x_axis)
  
  
  return(out = list(plot = plot,
                    data = out,
                    color = color_code))
  
  
}



##### 3) sampleID_order ##### 


sampleID_order <- function(tax_data, ph) {
  sample_order <-  tax_data$data %>% 
    data.frame() %>%
    dplyr::group_by(SampleID) %>% 
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    filter(Phylum == ph) %>% 
    dplyr::group_by(SampleID) %>% 
    dplyr::summarise(Abundance = sum(Abundance)) %>% 
    dplyr::arrange(Abundance) %>%
    pull(SampleID) %>% as.character()
  
  tax_data$data$SampleID <- factor(tax_data$data$SampleID, levels = sample_order)
  
  return(tax_data)
}



sampleID_order_g <- function(tax_data, gn = NULL) {
  sample_order <- tax_data$data %>% 
    data.frame() %>%
    dplyr::group_by(SampleID) %>% 
    mutate(Abundance = Abundance / sum(Abundance))
  sample_order <- sample_order %>%
    filter(Genus == gn) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(Abundance = sum(Abundance)) %>%
    dplyr::arrange(Abundance) %>%
    pull(SampleID) %>% as.character()
  
  tax_data$data$SampleID <- factor(tax_data$data$SampleID, levels = sample_order)
  
  return(tax_data)
}





#### 4. compar mean abundance ####


# melt= ps.10.gn.rel.melt;
# type = "Target";
# tax_level = "Genus";
# topn = p10.gn.top;
# col = c("red", "blue", "grey");
# Phylum = c("Show");
# elements = c("mean_bar", "logfold_bar" , "p_val_ano","p_adj_ano")



##### 1. data melt/ define Other ###### ________________________________________
prepare_data = function(melt, tax_level,  tax_otu) {
  
  # make taxonomy table 
  tax_tab <-   melt[, c("OTU", "Phylum", tax_level)] %>% unique
  
  # subset taxa by OTU ID
  tax_tab2 <- tax_tab[tax_tab$OTU %in% tax_otu, ]
  tax_phylum <- tax_tab2$Phylum %>% unique
  tax_index <- tax_tab2[,  tax_level,  drop=T]
  
  # Reconstruct taxa classified as Others
  melt.2 <- melt # Back up
  melt.2[!melt.2[, "Phylum"] %in% tax_phylum, "Phylum"] <- "Others"
  melt.2[!melt.2[, "Phylum"] %in% tax_phylum, tax_level] <- "Others"
  
  
  out=list(melt = melt.2,  tax_table=tax_tab2,  Top_p= tax_phylum,  Top_taxa=tax_index)
  return(out)
}


##### 2. Phylum level + tax name? ##### ________________________________________
process_not_phylum <- function(ps.mlt2, tax_level, Top_p, p_tax_table) {
  for (i in Top_p) {
    G <- p_tax_table[p_tax_table[, "Phylum"] == i, tax_level]
    ps.mlt2[ps.mlt2[, "Phylum"] == i & !ps.mlt2[, tax_level] %in% G, tax_level] <- "Others"
  }
  ps.mlt3 <- ps.mlt2
  table <- ps.mlt3[ps.mlt3[, "Phylum"] %in% Top_p, c("Abundance", tax_level)]
  table.2 <- table %>% group_by(!!rlang::sym(tax_level)) %>% dplyr::summarise(sum.Abundance = sum(Abundance), .groups = 'drop')
  g_order <- table.2 %>%
    dplyr::arrange(-sum.Abundance) %>%
    select(!!rlang::sym(tax_level))
  ps.mlt3[, tax_level] <- factor(ps.mlt3[, tax_level], levels = rev(c(g_order[g_order != "Others"], "Others")))
  return(ps.mlt3)
}

process_show_phylum <- function(ps.mlt2, tax_level, Top_p, p_tax_table) {
  for (i in Top_p) {
    G <- p_tax_table[p_tax_table[, "Phylum"] == i, tax_level]
    ps.mlt2[ps.mlt2[, "Phylum"] == i & !ps.mlt2[, tax_level] %in% G, tax_level] <- paste0(i, "_Others")
  }
  
  for (i in Top_p) {
    G <- p_tax_table[p_tax_table[, "Phylum"] == i, tax_level]
    for (g in G) {
      ps.mlt2[ps.mlt2[, tax_level] == g, tax_level] <- paste0(i, "_", g)
    }
  }
  
  ps.mlt3 <- ps.mlt2
  
  table <- ps.mlt3[ps.mlt3[, "Phylum"] %in% Top_p, c("Abundance", "Phylum", tax_level)]
  table.2 <- table %>% group_by(!!rlang::sym("Phylum"), !!rlang::sym(tax_level)) %>%
    dplyr::summarise(sum.Abundance = sum(Abundance), .groups = 'drop') %>%
    as.data.frame()
  
  table.2[grepl("_Other", table.2[, tax_level]), "Other_taxa"] <- TRUE
  
  g_order <- table.2 %>%
    dplyr::arrange(sum.Abundance) %>%
    dplyr::arrange(desc(-Other_taxa)) %>%
    # dplyr::arrange((!!rlang::sym(tax_level) == paste0(Phylum, "_Others"))) %>%
    dplyr::arrange(desc(Phylum)) %>% .[, tax_level]
  p_order <- table %>% .[, "Phylum"] %>% unique %>% sort
  
  ps.mlt3[, "Phylum"] <- factor(ps.mlt3[, "Phylum"], levels = rev(p_order))
  ps.mlt3[, tax_level] <- factor(ps.mlt3[, tax_level], levels = c("Others", g_order))
  return(ps.mlt3)
}

##### 3. Mean Abundance Calculation  #####  ____________________________________
Mean_Abundance_Calculation <- function(ps.mlt3, tax_level, type){
  df <- plyr::ddply(ps.mlt3, c(tax_level, type), 
                    dplyr::summarise, 
                    Mean=mean(Abundance),
                    Sd=sd(Abundance, na.rm=TRUE), 
                    N=length(Abundance), 
                    Se=Sd/sqrt(N))
  return(df)
}
##### 4. statistical test ##### ________________________________________________
kruskal_wallis_test <- function(data, type, tax_level) {
  
  kruskal_df <- data %>%
    tidyr::nest(data = -tax_level) %>%
    dplyr::mutate(
      test = purrr::map(.x = data, ~ kruskal.test(as.formula(glue("Abundance ~ {type}")), data = .x) %>% tidy))  %>%
    tidyr::unnest(test) %>%
    dplyr::mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
    select(!!rlang::sym(tax_level), p.value, p.adjust) %>%
    as.data.frame()
  
  # 추가적인 처리
  kruskal_df <- kruskal_df %>%
    dplyr::mutate(fil = "f") %>%
    dplyr::mutate(rank = as.numeric(rank(as.factor(kruskal_df[, tax_level]))))
  
  colnames(kruskal_df)[2:3] <- c("p.val", "p.adj")
  return(kruskal_df)
}

wilcoxon_test <- function(data, type, tax_level) {
  
  wilcox_df <- data %>%
    tidyr::nest(data = -tax_level) %>%
    dplyr::mutate(test = purrr::map(.x = data, ~ wilcox.test(as.formula(glue("Abundance ~ {type}")), data = .x) %>% tidy)) %>%
    tidyr::unnest(test) %>%
    dplyr::mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
    dplyr::select(!!rlang::sym(tax_level), p.value, p.adjust) %>%
    as.data.frame()
  
  wilcox_df <- wilcox_df %>%
    dplyr::mutate(fil = "f") %>%
    dplyr::mutate(rank = as.numeric(rank(as.factor(wilcox_df[, tax_level]))))
  
  colnames(wilcox_df)[2:3] <- c("p.val", "p.adj")
  
  return(wilcox_df)
}
##### 5. mean bar ##### ________________________________________________________

# Mon Apr 21 12:58:15 2025 ------------------------------
  #   # 단점은 error bar  중심이 bar 랑 잘 안맞음..근데 이건 고치는 법을 잘 모르겠다;

generate_mean_bar <- function(df, tax_level, type, col) {
  
  
  mean_bar <- ggplot(df,  aes(y = !!rlang::sym("Mean"), 
                              x = !!rlang::sym(tax_level), 
                              fill = !!rlang::sym(type)))  +
    # geom_bar(stat = "identity", position = "dodge") +
    geom_bar(stat = "identity",  position=position_dodge2(reverse = TRUE)) +
    scale_fill_manual(values = col) + 
    ggplot2::coord_flip()
  
  df
  mean_bar2 <- mean_bar +
    geom_errorbar(aes(ymin = Mean - Se, ymax = Mean + Se), 
                  position = position_dodge2(reverse = TRUE, padding = 0.5),
                  # width =1 
    ) +
    theme_classic() 
  mean_bar3 <- mean_bar2 +
    GGally::geom_stripped_cols() + 
    labs(y = "Mean Relative Abundance", x = NULL) + 
    theme(axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth  = 0.5),
          axis.ticks.x = element_line(linewidth  = 0.5),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(margin = ggplot2::margin(r = 0)), 
          axis.text.y = element_text(size = 10, color = "black", margin = ggplot2::margin(b = 6), hjust = 0),
          axis.title.x = ggplot2::element_text(size = 10, color = "black", hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = ggplot2::unit(0.1, "cm"),
          legend.direction = "vertical",
          legend.justification = "left",
          legend.text = ggplot2::element_text(size = 8, face = "bold"),
          legend.box.just = "right", 
          plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")
    ) 
  print(mean_bar3)
}

# generate_mean_bar <- function(df, tax_level, type, col) {
#   
#   
#   mean_bar <- ggplot(df,  aes_string(y = "Mean", x = tax_level, fill = type)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     scale_fill_manual(values = col) + 
#     ggplot2::coord_flip()
#   
#   mean_bar2 <- mean_bar +
#     geom_errorbar(aes(ymin = Mean - Se, ymax = Mean + Se), position = position_dodge(0.9), width = 0.2) +
#     theme_classic() 
#   
#   mean_bar3 <- mean_bar2 +
#     GGally::geom_stripped_cols() + 
#     labs(y = "Mean Relative Abundance (%)", x = NULL) + 
#     theme(axis.ticks.y = element_blank(),
#           axis.line.y = element_blank(),
#           axis.line.x = element_line(linewidth  = 0.5),
#           axis.ticks.x = element_line(linewidth  = 0.5),
#           panel.grid.major.y = element_blank(),
#           panel.grid.major.x = element_blank(),
#           axis.text = element_text(size = 10, color = "black"),
#           axis.text.x = element_text(margin = ggplot2::margin(r = 0)), 
#           axis.text.y = element_text(size = 10, color = "black", margin = ggplot2::margin(b = 6), hjust = 0),
#           axis.title.x = ggplot2::element_text(size = 10, color = "black", hjust = 0.5),
#           legend.title = element_blank(),
#           legend.position = "top",
#           legend.key.size = ggplot2::unit(0.1, "cm"),
#           legend.direction = "vertical",
#           legend.justification = "left",
#           legend.text = ggplot2::element_text(size = 8, face = "bold"),
#           legend.box.just = "right", 
#           plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")
#     ) 
#   print(mean_bar3)
# }

##### 6. mean bar ##### ________________________________________________________
generate_heatmap <- function(melt2, tax_level, type){
  data <- melt2[,c(tax_level, "Abundance",type)  ] %>% 
    dplyr::group_by(!!rlang::sym(tax_level), !!rlang::sym(type)) %>% 
    summarise(Abundance2 = sum(Abundance))# %>% 
  # dplyr::mutate(z_score = scale(Abundance2))
  
  data2 <- data %>%
    group_by(!!rlang::sym(tax_level)) %>%
    mutate(MMS_Abundance  = (Abundance2 - min(Abundance2)) / (max(Abundance2) - min(Abundance2))) %>%
    ungroup()
  
  heat <- ggplot(data2, aes(x = !!rlang::sym(type), y = !!rlang::sym(tax_level), fill = MMS_Abundance )) +
    geom_tile(color = "grey", #"white", 
              lwd = 1.0, linetype = 1) +
    # scale_fill_continuous(breaks = c(0, 1), labels = c("low", "high"), type = "viridis")+
    scale_fill_gradient2 ( low =  "white" ,  high = "red" )  +
    labs(y = NULL, x = NULL) +
    theme_void()+
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(margin = ggplot2::margin(r = 0), angle = 90, vjust = 0.5, hjust = 1), 
          axis.text.y =  element_blank(), 
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = ggplot2::unit(0.2, "cm"),
          legend.direction = "vertical",
          legend.justification = "left",
          legend.text = ggplot2::element_text(size = 8, face = "bold"),
          legend.box.just = "right", 
          plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")) 
  return(heat)
}
##### 7. p-value #####  ________________________________________________________
plot_p_values <- function(df, adj) {
  label_val <- NA
  if(adj == "Yes"){
    label_val <- "p.adj"
  } else if(adj == "No") {
    label_val <- "p.val"
  }
  
  
  p_val.p <- df %>% 
    ggplot() + 
    geom_text(aes(fil, rank,
                  label = ifelse(df[,label_val] < 0.001, "< 0.001", round(df[,label_val], 3))),
              fontface = ifelse(df[,label_val] < 0.05, "bold", "plain")
    )  +
    labs(x = paste(label_val), y = NULL) + 
    scale_y_discrete(position = "right") + 
    theme(axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(c(0,0.2, 0, 0), "cm"),
          legend.position = "non")
  print(p_val.p)
  return(p_val.p)
}

##### 8. log fold bar #####  ___________________________________________________
log2fold_bar <- function(df, tax_level, col){
  # logfold_bar: Log2 fold change
  tax_list <- unique(df[, tax_level])
  df2 <- data.frame(row.names = tax_list, tax = tax_list)
  
  for (i in tax_list) {
    Mean.df <- df[df[, tax_level] %in% i, ]$Mean
    df2[df2$tax %in%i, "log_2_fold_change"] <- log2(Mean.df[1]/Mean.df[2]) 
  }
  
  
  logfold_df <- df2 %>% 
    dplyr::mutate(new_col = case_when(
      is.nan(log_2_fold_change) ~ "0",
      is.infinite(log_2_fold_change) & log_2_fold_change > 0 ~ "> 0",
      is.infinite(log_2_fold_change) & log_2_fold_change < 0 ~ "< 0",
      log_2_fold_change < 0 ~ "< 0",
      log_2_fold_change >= 0 ~ "> 0" ))
  logfold_df$new_col <- factor(logfold_df$new_col, levels =  c("> 0", "< 0", "0"))
  
  
  ggplot(logfold_df, aes(x = tax, y = log_2_fold_change, fill = new_col)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) + 
    GGally::geom_stripped_cols() + 
    labs(y = 'Log2 Fold Change', x =NULL ) + 
    theme_classic() +
    scale_fill_manual(values = col) + # Update color values
    scale_color_manual(values = col)+ 
    theme(axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(size = 0.5),
          axis.ticks.x = element_line(size = 0.5),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 10,color = "black"),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 10, color = "black",margin = margin(b = 6)),
          axis.title.x = element_text(size = 11,color = "black", hjust = 0.5),
          legend.position = "non") + 
    geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed", color = "black") + 
    ggplot2::coord_flip()
}

##### 8. Run function mean plot 2 ##### ________________________________________



mean_plot <- function(melt,
                      tax_level = "Genus",
                      tax_otu ,
                      type = "Group",
                      palet,
                      Phylum = c("Show")){
  
  # if(is.null(palet) ){ palet = c("red", "blue", "grey")}
  
  # 1. data melt
  out_data <- prepare_data(melt, tax_level, tax_otu )
  
  # 2. change label by Top N
  if (Phylum == "Not") {   
    out_data2 <- process_not_phylum(ps.mlt2 = out_data$melt,
                                    tax_level = tax_level, 
                                    Top_p = out_data$Top_p, 
                                    p_tax_table =  out_data$tax_table)
  } else if (Phylum == "Show" ){
    out_data2 <- process_not_phylum(ps.mlt2 = out_data$melt,
                                    tax_level = tax_level, 
                                    Top_p = out_data$Top_p, 
                                    p_tax_table =  out_data$tax_table)
  }
  
  
  # 3. Mean_Abundance_Calculation
  df <- Mean_Abundance_Calculation(out_data2, tax_level, type)
  type_num <- length(unique(df[,type]))
  print("step2")
  
  
  # 4.Graph Generation
  
  if(type_num>2){
    mean_bar <- generate_mean_bar(df = df, tax_level = tax_level, type = type, col = palet)
    kruskal_df <- kruskal_wallis_test(out_data2, type, tax_level)
    p_val_ano <- plot_p_values(df=kruskal_df, adj = "No") 
    p_adj_ano <- plot_p_values(df=kruskal_df, adj = "Yes")
    heatmap <- generate_heatmap(melt2 = out_data2, tax_level = tax_level, type = type)
    statistical_result <- kruskal_df
    
    
  } else if(type_num==2){
    
    mean_bar <- generate_mean_bar(df = df, tax_level = tax_level, type = type, col = palet)
    logfold_bar <- log2fold_bar(df = df, tax_level = tax_level, col = palet)
    wilcoxon_df <- wilcoxon_test(data = out_data2, 
                                 type = type,
                                 tax_level = tax_level)
    p_val_ano <- plot_p_values(df = wilcoxon_df, adj = "No") 
    p_adj_ano <- plot_p_values(df = wilcoxon_df, adj = "Yes")
    statistical_result <- wilcoxon_df
    
  }
  
  
  # 5. Show plot
  
  # weights <- c(mean_bar = 1.5, logfold_bar = 0.7, p_val_ano = 0.5, p_adj_ano = 0.5, heatmap = 1.2)
  # selected_weights <- weights[names(selected_plots)]
  
  if(type_num>2){
    
    plot <- mean_bar + heatmap + p_val_ano + p_adj_ano +
      patchwork::plot_layout(ncol = 4, nrow = 1, widths = c(1.5, 0.4, 0.4,0.4))
    
  } else if(type_num==2){
    plot <- mean_bar + logfold_bar + p_val_ano + p_adj_ano +
      patchwork::plot_layout(ncol = 4, nrow = 1, widths = c(1.5, 0.7, 0.4,0.4))
    
  }
  
  
  
  #   if (all(elements %in% c("mean_bar","p_val_ano", "p_adj_ano" ))){
  #     plot <- get(elements[1]) + get(elements[2]) + get(elements[3]) +
  #       patchwork::plot_layout(ncol = 3, nrow = 1, widths = c(1.3, 0.7, 0.5))
  #   } else if (all(elements %in% c("mean_bar","logfold_bar","p_adj_ano" ))){
  #     plot <- get(elements[1]) + get(elements[2]) + get(elements[3]) +
  #       patchwork::plot_layout(ncol = 3, nrow = 1, widths = c(1.3, 0.7, 0.5))
  #   } else if (all(elements %in% c("mean_bar","p_adj_ano" ))){    
  #     plot <- get(elements[1]) + get(elements[2]) +
  #       patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(1.5,  0.5))
  #   } else if (all(elements %in% c("mean_bar","p_val_ano" ))){    
  #     plot <- get(elements[1]) + get(elements[2]) +
  #       patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(1.5,  0.5))
  #   } else if (all(elements %in% c("mean_bar","logfold_bar","p_val_ano" ))){
  #     plot <- get(elements[1]) + get(elements[2]) + get(elements[3]) +
  #       patchwork::plot_layout(ncol = 3, nrow = 1, widths = c(1.3, 0.7, 0.5))
  #   } else if (all(elements %in% c("mean_bar", "logfold_bar" , "p_val_ano","p_adj_ano"))){
  #     plot <- mean_bar + logfold_bar + p_val_ano + p_adj_ano +
  #       patchwork::plot_layout(ncol = 4, nrow = 1, widths = c(1.5, 0.7, 0.4,0.4))
  #   }
  #   print(plot)
  
  
  out =  list(plot = plot, 
              statistical_result = statistical_result, 
              df=df )
  return(out)
  
}



#### 5. volcano plot #### 

volcano_plot <- function(ph_glom , tax_rank, comparision, compar1, compar2, x_lim = 6, y_lim=2.5, 
                         mycolors =c("blue", "red", "grey"), lfc = 1 ) {
  
  sample_d <- sample_data(ph_glom) %>% data.frame()
  if ( !c(compar1) %in% sample_d[, comparision] ) { cat("Compar1 is not in Comparision column")
  } else if ( !c(compar1) %in% sample_d[, comparision] ) { cat("Compar2 is not in Comparision column") } 
  
  tax_d <- tax_table(ph_glom)  %>% as.data.frame() 
  otu_d <- otu_table(ph_glom)  %>% as.data.frame() 
  table_ot <- merge(tax_d, otu_d, by = "row.names") %>%
    column_to_rownames("Row.names")
  
  otu_table <- table_ot[, c(tax_rank, colnames(otu_d))]
  
  
  summarize_otu_table <- function(otu_table) {
    summarized_table <- aggregate(as.formula(glue(". ~ ", tax_rank)), data = otu_table, 
                                  FUN = function(x) sum(x))
    return(summarized_table)
  }
  
  # 함수 적용
  result <- summarize_otu_table(otu_table)
  result <- result %>% column_to_rownames(tax_rank)
  
  
  
  table <- result %>% as.matrix %>% t()
  
  p_value <- apply(table, 2,
                   function(x){wilcox.test(as.numeric(x[sample_d[,comparision] == compar1]),
                                           as.numeric(x[sample_d[,comparision] == compar2]),
                                           correct = F)$p.value})
  # print(p_value)
  
  
  table.2 <- cbind(table %>% t() %>% data.frame(), data.frame(p_value))
  table.2$q_value <- p.adjust(table.2$p_value, "BH")
  
  table_compar1 <- table[sample_d[,comparision] == compar1, ]
  table_compar2 <- table[sample_d[,comparision] == compar2, ]
  
  mean.compar1 <- colMeans(data.matrix(table_compar1))
  mean.compar2 <- colMeans(data.matrix(table_compar2))
  
  table.2$log2fold_1_2 <- log2(mean.compar1/mean.compar2)
  
  # Nan, Inf 제거
  # table.3 <- table.2[!is.infinite(rowSums(table.2)),]
  table.3 <- table.2[!is.nan(rowSums(table.2)),]
  table.3$threshold <- table.3$p_value < 0.05
  table.3$BH_threshold <- table.3$q_value < 0.05
  table.3$fold_threshold <- abs(table.3$log2fold_1_2) > lfc
  
  # sort 
  table.3_ord <- table.3[order(table.3$p_value), ] 
  table.3_ord$fold_change <- table.3_ord$log2fold_1_2 > 0
  table.3_ord$Species <- rownames(table.3_ord)
  
  # p-value cutoff
  table.3_ord$DA <- NA
  table.3_ord[table.3_ord$log2fold_1_2 > 0 & 
                table.3_ord$fold_threshold == T &  
                table.3_ord$threshold == T, "DA"] <- compar1
  table.3_ord[table.3_ord$log2fold_1_2 < 0 & 
                table.3_ord$fold_threshold == T&
                table.3_ord$threshold == T, "DA"] <- compar2
  
  table.3_ord$DA[is.na(table.3_ord$DA)] <- "Not Sig"
  # table.3_ord$DA %>% table()
  
  # q-value cutoff
  table.3_ord$BH_DA <- NA
  table.3_ord[table.3_ord$log2fold_1_2 > 0 & 
                table.3_ord$fold_threshold == T &  
                table.3_ord$BH_threshold == T, "BH_DA"] <- compar1
  table.3_ord[table.3_ord$log2fold_1_2 < 0 & 
                table.3_ord$fold_threshold == T&
                table.3_ord$BH_threshold == T, "BH_DA"] <- compar2
  
  table.3_ord$BH_DA[is.na(table.3_ord$BH_DA)] <- "Not Sig"
  
  ## Change point color (https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html)
  
  names(mycolors) <- c(compar2, compar1, "Not Sig")
  
  p1 <- ggplot(table.3_ord, aes(x=log2fold_1_2, y=-log10(p_value)) )  +
    geom_point(aes(colour=DA, alpha = fold_threshold)) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_vline(xintercept = 1, linetype="dotted") +
    geom_vline(xintercept = -1, linetype="dotted") +
    geom_text_repel(aes(label = ifelse(fold_threshold == T &  threshold == T, Species, "")),
                    max.overlaps =15 ) +
    ylim(0, y_lim) +
    xlim(-x_lim, x_lim) +
    scale_alpha_manual(values=c(0.3, 1)) +
    scale_shape_manual(values=c(16, 15)) +
    scale_color_manual(values=mycolors) +
    xlab("Log2 fold change") + 
    ylab("-Log10 p-value") +
    theme_bw() + 
    theme(legend.position = "right",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title = element_text(size=10))  
  
  p2 <- ggplot(table.3_ord, aes(x=log2fold_1_2, y=-log10(q_value)) )  +
    geom_point(aes(colour=BH_DA, alpha = fold_threshold)) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_vline(xintercept = 1, linetype="dotted") +
    geom_vline(xintercept = -1, linetype="dotted") +
    geom_text_repel(aes(label = ifelse(fold_threshold == T &  BH_threshold == T, Species, "")),
                    max.overlaps =15 ) +
    ylim(0, y_lim) +
    xlim(-x_lim, x_lim) +
    scale_alpha_manual(values=c(0.3, 1)) +
    scale_shape_manual(values=c(16, 15)) +
    scale_color_manual(values=mycolors) +
    xlab("Log2 fold change") + 
    ylab("-Log10 p-value (FDR)") +
    theme_bw() + 
    theme(legend.position = "right",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.title = element_text(size=10))  
  
  return(out = list(Original=p1, 
                    BH_FDR_correction=p2)
  )
  
  
  
}

##### blast_arrange ##### 
blast_arrange <- function( blast_output, ASV) {
  
  table <- data.frame(
    `ASV` = blast_output[, ASV],
    `Before` = blast_output$sscinames,
    `Genus`  = NA,
    `Species` = NA,
    `Genus_Species` = NA
  )
  
  
  Clear_Species <- function(data) {
    Sp <- strsplit(data, " ")[[1]]          # 문자열 분할
    Genus <- Sp[1]
    Species <- Sp[2]
    G_S   <- paste(Sp[1], Sp[2]) # , sep = " ") # 첫 번째 및 두 번째 단어 출력
    
    return(c(Genus, Species, G_S ))
  }
  
  for (i in blast_output$sscinames ) {
    table[table$Before == i, "Genus"] <- Clear_Species(i)[1]
    table[table$Before == i, "Species"] <-Clear_Species(i)[2]
    table[table$Before == i, "Genus_Species"] <- Clear_Species(i)[3]
  }
  return(table)
}

# ETC
library(glue)
seq_to_fas <- function(seq, path) {
  fas <- data.frame(
    ID = seq,
    Sequence =  seq)
  
  
  # save as fas
  Xfasta <- character(nrow(fas) * 2)
  Xfasta[c(TRUE, FALSE)] <- paste0(">", fas$ID)
  Xfasta[c(FALSE, TRUE)] <- fas$Sequence
  
  print(glue("save .fas file in ",path ))
  
  writeLines(Xfasta, path)
  
}


#### other function #### 
'%!in%' <- function(x,y)!('%in%'(x,y))


prevalence_filt <- function(ps, threshold){
  prevelancedf = apply(X = otu_table(ps),
                       MARGIN = 1,
                       FUN = function(x){sum(x > 0)})
  
  # Add taxonomy and total read counts to this data.frame
  prevelancedf = data.frame(Prevalence = prevelancedf,
                            TotalAbundance = taxa_sums(ps),
                            tax_table(ps))
  
  
  prevelancedf1 = subset(prevelancedf, Phylum %in% get_taxa_unique(ps, taxonomic.rank = "Phylum"))
  
  prevalenceThreshold = threshold * nsamples(ps)
  keepTaxa = rownames(prevelancedf1)[(prevelancedf1$Prevalence >= prevalenceThreshold)]
  
  ps = prune_taxa(keepTaxa, ps)
  return(ps)
  
} 

ps_tax_update <- function(physeq_obj) {
  # Tue Apr  8 15:00:37 2025 ------------------------------
  
  require(phyloseq)
  require(dplyr)
  require(stringr)
  
  # Extract taxonomy and sequence data
  tax_df <- tax_table(physeq_obj) %>% data.frame()
  seq_df <- refseq(physeq_obj) %>%
    as.data.frame() %>%
    dplyr::rename(Sequence = x) %>%
    mutate(length = str_length(Sequence))
  
  # Combine taxonomy and sequence information
  merged_df <- cbind(tax_df, seq_df)
  
  # Add abundance and relative abundance
  merged_df$Abundance <- taxa_sums(physeq_obj)
  merged_df$RelativeAbundance <- round(merged_df$Abundance * 100 / sum(merged_df$Abundance), 4)
  
  # Determine digit length based on number of taxa
  num_taxa <- ntaxa(physeq_obj)
  digit_length <- nchar(as.character(num_taxa))
  asv_format <- paste0("ASV_%0", digit_length, "d")
  
  # Sort and assign ASV IDs
  merged_df <- merged_df %>%
    arrange(-Abundance) %>%
    mutate(ASV = sprintf(asv_format, 1:nrow(.)))
  
  # Update tax_table in phyloseq object
  updated_physeq <- physeq_obj
  tax_table(updated_physeq) <- tax_table(as.matrix(merged_df))
  
  return(list(
    updated_phyloseq = updated_physeq,
    tax_table_df = tax_table(updated_physeq) %>% data.frame()
  ))
}
