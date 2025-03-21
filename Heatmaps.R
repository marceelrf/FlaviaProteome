fn_ORA_gene_diff(AA_CKD_vs_CKD_tibble)
fn_ORA_gene_diff(AA_CKDxAA_tibble)
fn_ORA_gene_diff(AAxNT_tibble)

total_diff <-
  unique(c(fn_ORA_gene_diff(AAxNT_tibble),
       fn_ORA_gene_diff(AA_CKDxAA_tibble),
       fn_ORA_gene_diff(AA_CKD_vs_CKD_tibble)))


fn_filter_diffs <- function(x){
  x %>%
    dplyr::filter(Symbol %in% total_diff) %>% 
    dplyr::select(logFC,Symbol,P.Value,adj.P.Val) %>%
    dplyr::distinct() %>% 
    tidyr::drop_na()
}

agg_AAxNT_diff <- fn_filter_diffs(AAxNT_tibble) %>% 
  dplyr::filter(P.Value < 0.05) %>% 
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(logFC = mean(logFC))

agg_AA_CKD_vs_CKD_diff <- fn_filter_diffs(AA_CKD_vs_CKD_tibble) %>% 
  dplyr::filter(P.Value < 0.05) %>% 
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(logFC = mean(logFC))
agg_AA_CKDxAA_diff <- fn_filter_diffs(AA_CKDxAA_tibble) %>% 
  dplyr::filter(P.Value < 0.05) %>% 
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(logFC = mean(logFC))


agg_list <- list(AAxNT = agg_AAxNT_diff,
                 AA_CKD_x_CKD = agg_AA_CKD_vs_CKD_diff,
                 AA_CKDxAA = agg_AA_CKDxAA_diff)

bind_rows(agg_list,.id = "Comp") %>% 
  complete(Comp, Symbol) %>% 
  mutate(Symbol = factor(Symbol)) %>%
  mutate(Symbol = fct_rev(Symbol)) %>%
  ggplot(aes(x = Comp, y = Symbol, fill = logFC)) +
  geom_tile(col = "black") +
  scale_fill_gradient2(low = "dodgerblue",
                       mid = "white",
                       high = "firebrick4",
                       midpoint = 0,
                       breaks = c(-5,-2.5,-1,0,1,2.5,5)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  labs(title = "Differential expressed proteins",
       subtitle = "Filter by p-value < 0.05 per group",
       x = "", y = "") +
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(face = "bold",
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        axis.text.y = element_text(size = 5),
        plot.title = element_text(face = "bold",size = 22),
        panel.border = element_rect(fill = NA,linewidth = 2))

ggsave(filename = "Output/heatmap.png",dpi = 450,
       height = 18,width = 6)
