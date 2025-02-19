library(ComplexHeatmap)
# Functions ---------------------------------------------------------------

dotplot_mrf <- function(data, ncat = 10, empates = F,
                        titulo = "" ){ 
  data@result %>%
    mutate(Description = fct_reorder(Description, Count)) %>%  # Ordenar pelo número de genes
    slice_max(Count, n = ncat, with_ties = empates) %>% 
    ggplot(aes(x = Count, y = Description, size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_viridis_b(option = "C",direction = -1) +  # Escala de cores perceptual
    scale_size(range = c(3, 10)) +  # Ajustar o tamanho dos pontos
    labs(
      x = "Número de Genes",
      y = "Termos Enriquecidos",
      color = "-log10(p.adjust)",
      size = "Contagem",
      title = titulo
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right")
  }

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

heatmap_enrich <- function(data) {
  data@result %>%
    pivot_longer(cols = starts_with("Gene_"), 
                 names_to = "Gene", values_to = "Expression") %>%
    mutate(Gene = fct_reorder(Gene, Expression, .fun = mean),
           Description = fct_reorder(Description, -p.adjust)) %>%
    ggplot(aes(x = Gene, y = Description, fill = Expression)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "B") +  # Escala de cores perceptual
    labs(x = "Genes", y = "Termos Enriquecidos", fill = "Expressão") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



# dotplot -----------------------------------------------------------------

dotplot_mrf(ora_result_AA_CKDxAA,
            empates = T,
            titulo = "AA+CKD vs AA")
dotplot_mrf(ora_result_AAxNT,
            empates = T,
            titulo = "AA vs NT")
dotplot_mrf(ora_result_AA_CKD_vs_CKD,
            empates = T,
            titulo = "AA+CKD vs CKD")


# heatmap -----------------------------------------------------------------

ora_result_AAxNT@result %>% 
  separate_rows(geneID) %>% 
  ggplot(aes(x = geneID, y = Description)) +
  geom_tile(col = "black")

ora_result_AA_CKD_vs_CKD@result %>% 
  separate_rows(geneID) %>% 
  select(Description,geneID) %>% 
  group_by(Description,geneID) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = geneID,values_from = n,values_fill = 0) %>% 
  column_to_rownames("Description") %>% 
  as.matrix() %>% 
  Heatmap(clustering_distance_rows = "binary",
          clustering_distance_columns = "binary")
