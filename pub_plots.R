library(ggrepel)
library(viridis)
# Up e Downs



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
         col = "Significance",
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

fn_get_diff_mat


# Barplot vias

ora_list <- list("AAxNT" = ora_result_AAxNT,
                 "AA_CKDxAA" = ora_result_AA_CKDxAA,
                 "AA_CKD_vs_CKD" = ora_result_AA_CKD_vs_CKD)

fn_barplot <- function(x, plot_title = "") {
  x@result %>%
    filter(pvalue < 0.05) %>%
    arrange((Count), pvalue) %>%  # Ordena primeiro por Count (descendente) e depois por pvalue (ascendente)
    mutate(ID = factor(ID, levels = unique(ID))) %>%  # Garante a ordenação
    ggplot(aes(x = Count, y = ID, fill = -log10(pvalue))) +
    geom_col(color = "grey50") +
    scale_fill_viridis(option = "F") +
    labs(title = plot_title,
         y = '') +
    theme_bw()
}


ORA_plot_list <- map2(ora_list,nms_comps,
                      ~fn_barplot(.x,plot_title = .y))

# fn_barplot(ora_list$AAxNT)
# fn_barplot(ora_list$AA_CKDxAA)
# fn_barplot(ora_list$AA_CKD_vs_CKD)

walk2(ORA_plot_list, nms_comps,
      .f = \(x,y) ggsave(paste0("Output/pub_plots/barplot_", y,".png"),
                         plot = x,
                         scale = 3,
                         bg = "white",
                         dpi = 600))
