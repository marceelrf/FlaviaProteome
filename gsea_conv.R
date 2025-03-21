mapping <- read_tsv(file = "Data/idmapping_AA_CKD_vs_CKD.tsv")


AAxNT_tibble <- AAxNT %>% 
  separate_rows(UNIPROT,sep = ";") %>% 
  left_join(mapping, by = join_by(UNIPROT == From)) %>% 
  select(-SYMBOL,-ENTREZID) %>% 
  rename(Symbol = To)

AA_CKD_vs_CKD_tibble <- AA_CKD_vs_CKD %>% 
  separate_rows(UNIPROT,sep = ";") %>%
  left_join(mapping, by = join_by(UNIPROT == From)) %>% 
  select(-SYMBOL,-ENTREZID) %>% 
  rename(Symbol = To)

AA_CKDxAA_tibble <- AA_CKDxAA %>% 
  separate_rows(UNIPROT,sep = ";") %>%
  left_join(mapping, by = join_by(UNIPROT == From)) %>% 
  select(-SYMBOL,-ENTREZID) %>% 
  rename(Symbol = To)



# ORA ---------------------------------------------------------------------

fn_ORA_gene_diff <- function(x, lfc = 1, p_valor = 0.05) {

  x %>% 
    drop_na() %>% 
    dplyr::mutate(Diff = case_when(
      logFC < -lfc & P.Value < p_valor ~ "Down",
      logFC > lfc & P.Value < p_valor ~ "Up",
      TRUE ~ "Ns"
    )) %>%
    dplyr::filter(Diff != "Ns") %>% 
    pull(Symbol)
}

# fn_ORA_gene_diff(AA_CKDxAA_tibble)

ora_result_AAxNT <- enricher(fn_ORA_gene_diff(AAxNT_tibble),
                             TERM2GENE = msig_data,
                             minGSSize = 1,
                             maxGSSize = 500,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "BH")

ora_result_AA_CKDxAA <- enricher(fn_ORA_gene_diff(AA_CKDxAA_tibble),
                                 TERM2GENE = msig_data,
                                 minGSSize = 1,
                                 maxGSSize = 500,
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1,
                                 pAdjustMethod = "BH")
ora_result_AA_CKD_vs_CKD <- enricher(fn_ORA_gene_diff(AA_CKD_vs_CKD_tibble),
                                     TERM2GENE = msig_data,
                                     minGSSize = 1,
                                     maxGSSize = 500,
                                     pvalueCutoff = 1,
                                     qvalueCutoff = 1,
                                     pAdjustMethod = "BH")

# 
# ora_result_AA_CKDxAA@result %>% 
#   mutate(qscore = -log(p.adjust, base=10)) %>% 
#   barplot(x="qscore")
