AA_CKD_vs_CKD_GL <- AA_CKD_vs_CKD %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(stat = logFC*-log10(adj.P.Val)) %>% 
  arrange(desc(stat)) %>% 
  select(SYMBOL,stat) %>% 
  deframe()
AAxNT_GL <- AAxNT %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(stat = logFC*-log10(adj.P.Val)) %>% 
  arrange(desc(stat)) %>% 
  select(SYMBOL,stat) %>% 
  deframe()
AA_CKDxAA_GL <- AA_CKDxAA %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(stat = logFC*-log10(adj.P.Val)) %>% 
  arrange(desc(stat)) %>% 
  select(SYMBOL,stat) %>% 
  deframe()


AA_CKD_vs_CKD_gsea_result <- 
  GSEA(AA_CKD_vs_CKD_GL, TERM2GENE = msig_data)
AAxNT_gsea_result <- 
  GSEA(AAxNT_GL, TERM2GENE = msig_data)
AA_CKDxAA_gsea_result <- 
  GSEA(AA_CKDxAA_GL, TERM2GENE = msig_data)


