library(tidyverse)
library(ggrepel)
library(viridis)
library(ComplexHeatmap)
library(VennDiagram)

# Up e Downs

fn_ORA_gene_diff_topN <- function(x, lfc = 1, p_valor = 0.05, N = 30) {
  
  x %>%
    drop_na() %>%
    dplyr::mutate(
      abs_logFC = abs(logFC) # Calcula o valor absoluto do logFC
    ) %>%
    dplyr::mutate(Diff = case_when(
      logFC < -lfc & P.Value < p_valor ~ "Down",
      logFC > lfc & P.Value < p_valor ~ "Up",
      TRUE ~ "Ns"
    )) %>%
    dplyr::filter(Diff != "Ns") %>% 
    dplyr::arrange(desc(abs_logFC)) %>% # Ordena em ordem decrescente pelo valor absoluto do logFC
    head(N) %>% # Seleciona os 30 primeiros genes
    pull(Symbol)
}

# Volcano plot

volcanoplot <- function(x,plot_title ="", lfc = 1, p_valor = 0.05){
  
  x_trans <- select(x, logFC, P.Value, Symbol, UNIPROT) %>% 
    mutate(log_p_valor = -log10(P.Value)) %>% 
    dplyr::mutate(Diff = case_when(
      logFC < -lfc & P.Value < p_valor ~ "Down",
      logFC > lfc & P.Value < p_valor ~ "Up",
      TRUE ~ "Ns"
    )) %>% 
    mutate(Diff = factor(Diff, c("Up","Ns","Down")))
  
  top20 <- x_trans %>% 
    filter(Diff != "Ns") %>% 
    mutate(labs = case_when(
      is.na(Symbol) ~ UNIPROT,
      TRUE ~ Symbol
    )) %>% 
    slice_max(log_p_valor, n = 20,with_ties = T)
  
  x_trans %>% 
    ggplot(aes(x = logFC, y = log_p_valor, col = Diff)) +
    geom_point(alpha = .75) +
    geom_hline(yintercept = 1.30103,
               linetype = "dashed") +
    geom_vline(xintercept = c(-1,1),
               linetype = "dashed") +
    geom_label_repel(data = top20,
                     aes(label = labs),
                     show.legend = FALSE,
                     max.overlaps = 20) +
    scale_color_manual(values = c(
      "Up" = "firebrick4",
      "Ns" = "black",
      "Down" = "dodgerblue")) +
    labs(x = "log2FC",
         y = "-log10(p-value)",
         col = "",
         title = plot_title) +
    theme_bw() +
    theme(text = element_text(family = "serif"),
          plot.title = element_text(size = 22,
                                    face = "bold"),
          axis.title = element_text(size = 18,
                                    face = "bold"),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16,
                                      face = "bold"))
}

de_tibble <- list("AAxNT" = AAxNT_tibble,
                  "AA_CKDxAA" = AA_CKDxAA_tibble,
                  "AA_CKD_vs_CKD" = AA_CKD_vs_CKD_tibble)


nms_comps <- names(de_tibble)

volcano_list <- map2(de_tibble,nms_comps,
                     ~volcanoplot(.x,plot_title = .y))


walk2(volcano_list, nms_comps,
      .f = \(x,y) ggsave(paste0("Output/pub_plots/volcano_", y,".png"),
                         plot = x,
                         scale = 3,
                         bg = "white",
                         dpi = 600))

# Heatmap

total_diff

hm_data_total <-
  raw_data %>% 
  dplyr::select(Symbol = PG.Genes,
                ends_with(c("_R1","_R2","_R3","_R4"))) %>% 
  dplyr::filter(Symbol %in% total_diff)

#hm1

hm1_mat <-
  hm_data_total %>% 
  select(Symbol, starts_with(c("AA_R","non"))) %>% 
  dplyr::filter(Symbol %in% fn_ORA_gene_diff(AAxNT_tibble)) %>% 
  mutate(across(where(is.numeric),
                \(x) log10(x+1))) %>% 
  mutate(across(where(is.numeric),
                \(x) (x-mean(x))/sd(x)
                )
         ) %>%
  # pivot_longer(-Symbol,
  #              names_to = "Groups",
  #              values_to = "z_score") %>% 
  column_to_rownames("Symbol") %>% 
  as.matrix()

# Calculando a média por linha
row_means <- hm_data_total %>% 
  select(Symbol, starts_with(c("AA_R","non"))) %>% 
  dplyr::filter(Symbol %in% fn_ORA_gene_diff(AAxNT_tibble)) %>% 
  mutate(across(where(is.numeric),
                \(x) log10(x+1))) %>% 
  rowwise() %>% 
  mutate(mean_value = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  ungroup() %>% 
  pull(mean_value)

# Criando a anotação lateral direita com um ponto por linha
right_annotation <- rowAnnotation(
  Average = anno_points(row_means,
                      size = unit(2, "mm"),
                      gp = gpar(col = "black"))
)

tiff(filename = "Output/pub_plots/AAxNT_hm.tif",
     bg = "white",
     compression = "lzw",
     width = 3000,
     height = 3000,
     res = 600)
Heatmap(hm1_mat,
        name = "Z-score",
        border_gp = gpar(lwd = 2),
        cluster_columns = F,
        rect_gp = gpar(lwd = 1),
        right_annotation = right_annotation)
dev.off()
# hm2

paper_proteins <- c("Chil3","Lcn2","Spp1","Mmp8",
                    "Orm1","Orm2","Ltf","Apcs",
                    "Ambp","Atf6","Pon3","Ces2a",
                    "Ces2e","Ces3a")

hm2_mat <-
  hm_data_total %>% 
  select(Symbol, starts_with(c("AA_R","AA_CKD"))) %>% 
  #dplyr::filter(Symbol %in% fn_ORA_gene_diff_topN(AA_CKDxAA_tibble,N = 20)) %>% 
  dplyr::filter(Symbol %in% paper_proteins) %>% 
  mutate(across(where(is.numeric),
                \(x) log10(x+1))) %>% 
  mutate(across(where(is.numeric),
                \(x) (x-mean(x))/sd(x)
                )
  ) %>% 
  # pivot_longer(-Symbol,
  #              names_to = "Groups",
  #              values_to = "z_score") %>% 
  column_to_rownames("Symbol") %>% 
  as.matrix()



# Calculando a média por linha
row_means <- hm_data_total %>% 
  select(Symbol, starts_with(c("AA_R","AA_CKD"))) %>% 
  #dplyr::filter(Symbol %in% fn_ORA_gene_diff_topN(AA_CKDxAA_tibble, N = 20)) %>% 
  dplyr::filter(Symbol %in% paper_proteins) %>% 
  mutate(across(where(is.numeric), \(x) log10(x + 1))) %>% 
  rowwise() %>% 
  mutate(mean_value = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  ungroup() %>% 
  pull(mean_value)

# Criando a anotação lateral direita com um ponto por linha
right_annotation <- rowAnnotation(
  Average = anno_points(row_means,
                      size = unit(2, "mm"),
                      gp = gpar(col = "black"))
)

tiff(filename = "Output/pub_plots/AA_CKDxAA_hm_n20.tif",
     bg = "white",
     compression = "lzw",
     width = 4000,
     height = 12000,
     res = 600)
Heatmap(hm2_mat,
        name = "Z-score",
        border_gp = gpar(lwd = 2),
        cluster_columns = F,
        rect_gp = gpar(lwd = 1),
        right_annotation = right_annotation)
dev.off()

# hm 3

hm3_mat <-
  hm_data_total %>% 
  select(Symbol, starts_with(c("CKD_R","AA_CKD"))) %>% 
  dplyr::filter(Symbol %in% fn_ORA_gene_diff(AA_CKD_vs_CKD_tibble)) %>% 
  mutate(across(where(is.numeric),
                \(x) log10(x+1))) %>% 
  mutate(across(where(is.numeric),
                \(x) (x-mean(x))/sd(x)
                )
  ) %>% 
  # pivot_longer(-Symbol,
  #              names_to = "Groups",
  #              values_to = "z_score") %>% 
  column_to_rownames("Symbol") %>% 
  as.matrix()

# Calculando a média por linha
row_means <- hm_data_total %>% 
  select(Symbol, starts_with(c("CKD_R","AA_CKD"))) %>% 
  dplyr::filter(Symbol %in% fn_ORA_gene_diff(AA_CKD_vs_CKD_tibble)) %>% 
  mutate(across(where(is.numeric),
                \(x) log10(x+1))) %>% 
  rowwise() %>% 
  mutate(mean_value = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>% 
  ungroup() %>% 
  pull(mean_value)

# Criando a anotação lateral direita com um ponto por linha
right_annotation <- rowAnnotation(
  Average = anno_points(row_means,
                      size = unit(2, "mm"),
                      gp = gpar(col = "black"))
)

tiff(filename = "Output/pub_plots/AA_CKD_vs_CKD_hm.tif",
     bg = "white",
     compression = "lzw",
     width = 4000,
     height = 4000,
     res = 600)
Heatmap(hm3_mat,
        name = "Z-score",
        border_gp = gpar(lwd = 2),
        cluster_columns = F,
        rect_gp = gpar(lwd = 1),
        right_annotation = right_annotation)
dev.off()

# Barplot vias

ora_list <- list("AAxNT" = ora_result_AAxNT,
                 "AA_CKDxAA" = ora_result_AA_CKDxAA,
                 "AA_CKD_vs_CKD" = ora_result_AA_CKD_vs_CKD)

fn_barplot <- function(x, plot_title = "") {
  x@result %>%
    filter(pvalue < 0.05) %>%
    mutate(ID = str_split_fixed(ID, "_", 2)[,2]) %>% 
    mutate(ID = str_replace_all(ID,"_"," ")) %>% 
    arrange(pvalue) %>%  # Ordena primeiro por Count (descendente) e depois por pvalue (ascendente)
    mutate(ID = factor(ID, levels = rev(unique(ID)))) %>%  # Garante a ordenação
    ggplot(aes(x = Count, y = ID, fill = -log10(pvalue))) +
    geom_col(color = "grey50") +
    scale_fill_viridis(option = "F") +
    labs(title = plot_title,
         y = '') +
    theme_bw()
}


ORA_plot_list <- map2(ora_list, nms_comps,
                      ~fn_barplot(.x,plot_title = .y))
ORA_plot_list$AA_CKDxAA
# fn_barplot(ora_list$AAxNT)
# fn_barplot(ora_list$AA_CKDxAA)
# fn_barplot(ora_list$AA_CKD_vs_CKD)
paste0("Output/pub_plots/barplot_",
       nms_comps,
       ".png")

walk2(ORA_plot_list, nms_comps,
      .f = \(x,y) ggsave(paste0("Output/pub_plots/barplot_", y,".png"),
                         plot = x,
                         width = 12,
                         height = 8,
                         bg = "white",
                         dpi = 600),
      .progress = T)


# Venn diagram
deglist <- list()


deglist$AA_CKD_vs_CKD <-
  fn_ORA_gene_diff(AA_CKD_vs_CKD_tibble)

deglist$AA_CKDxAA <-
  fn_ORA_gene_diff(AA_CKDxAA_tibble)

deglist$AAxNT <-
  fn_ORA_gene_diff(AAxNT_tibble)

venn.diagram(deglist,
             filename = "Output/pub_plots/venn_diffprot.png",
             imagetype="png" ,
             height = 1500 , 
             width = 1500 , 
             resolution = 600,
             main = "Differentially expressed proteins",
             main.cex = 1,
             compression = "lzw",
             lwd = 2,
             col = c("firebrick","forestgreen","royalblue"),
             #label.col = c(rep("black",3),"firebrick",rep("black",3)),
             cex = .65,
             fontface = "bold",
             fontfamily = "serif",
             output = F,
             margin = 0.15,
             category.names = c("AA CKD vs. CKD",
                                "AA CKD vs. AA",
                                "AA vs. NT"),
             cat.col = c("firebrick","forestgreen","royalblue"),
             cat.dist = c(0.085,0.085,0.085),
             cat.cex = .5
             )


oralist <- list()

fn_filter_ORA <- function(x) {
  x@result %>% 
    filter(pvalue < 0.05) %>%
    #mutate(ID = str_split_fixed(ID, "_", 2)[,2]) %>% 
    #mutate(ID = str_replace_all(ID,"_"," ")) %>% 
    pull(ID)
}

oralist$AA_CKD_vs_CKD <-
  fn_filter_ORA(ora_result_AA_CKD_vs_CKD)

oralist$AA_CKDxAA <-
  fn_filter_ORA(ora_result_AA_CKDxAA)

oralist$AAxNT <-
  fn_filter_ORA(ora_result_AAxNT)



venn.diagram(x = oralist,
             filename = "Output/pub_plots/venn_oralist.png",
             imagetype="png" ,
             height = 1500 , 
             width = 1500 , 
             resolution = 600,
             main = "Enriched terms",
             main.cex = 1,
             compression = "lzw",
             lwd = 2,
             col = c("firebrick","forestgreen","royalblue"),
             #label.col = c(rep("black",3),"firebrick",rep("black",3)),
             cex = .65,
             fontface = "bold",
             fontfamily = "serif",
             output = F,
             margin = 0.15,
             euler.d = F,
             scaled = F,
             category.names = c("AA CKD vs. CKD",
                                "AA CKD vs. AA",
                                "AA vs. NT"),
             cat.col = c("firebrick","forestgreen","royalblue"),
             cat.dist = c(0.085,0.085,0.085),
             cat.cex = .5
)

# Table diff prot

feature_data %>% 
  filter(PG.Genes %in% fn_ORA_gene_diff(AA_CKD_vs_CKD_tibble)) %>% 
  write_csv(file = "Output/pub_plots/AA_CKD_vs_CKD.csv")

feature_data %>% 
  filter(PG.Genes %in% fn_ORA_gene_diff(AA_CKDxAA_tibble)) %>% 
  write_csv(file = "Output/pub_plots/AA_CKD_vs_AA.csv")

feature_data %>% 
  filter(PG.Genes %in% fn_ORA_gene_diff(AAxNT_tibble)) %>% 
  write_csv(file = "Output/pub_plots/AAxNT.csv")
