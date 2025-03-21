library(fgsea)
AA_CKD_vs_CKD_GL <- AA_CKD_vs_CKD_tibble %>% 
  filter(!is.na(Symbol)) %>% 
  mutate(stat = logFC*-log10(adj.P.Val)) %>% 
  group_by(Symbol) %>% 
  summarise(stat = mean(stat)) %>% 
  drop_na() %>%
  arrange(desc(stat)) %>%
  select(Symbol,stat) %>% 
  deframe()
AAxNT_GL <- AAxNT_tibble %>% 
  filter(!is.na(Symbol)) %>% 
  mutate(stat = logFC*-log10(adj.P.Val)) %>% 
  group_by(Symbol) %>% 
  summarise(stat = mean(stat)) %>% 
  drop_na() %>%
  arrange(desc(stat)) %>%
  select(Symbol,stat) %>% 
  deframe()
AA_CKDxAA_GL <- AA_CKDxAA_tibble %>% 
  filter(!is.na(Symbol)) %>% 
  mutate(stat = logFC*-log10(adj.P.Val)) %>% 
  group_by(Symbol) %>% 
  summarise(stat = mean(stat)) %>% 
  drop_na() %>%
  arrange(desc(stat)) %>%
  select(Symbol,stat) %>% 
  deframe()


gs_list <-
  msig_data %>% 
  group_split(gs_name)

# gsea_AA_CKDxAA <- clusterProfiler::GSEA(geneList = AA_CKDxAA_GL,
#                       TERM2GENE = msig_data,
#                       minGSSize = 1,
#                       maxGSSize = 500,
#                       pvalueCutoff = .05,
#                       pAdjustMethod = "BH")


?fgsea()
gsea_AAxNT <- 
  fgsea(gs_list,AAxNT_GL,minSize=1)

gsea_AA_CKD_vs_CKD <- 
  fgsea(gs_list,sort(AA_CKDxAA_GL + rnorm(length(AA_CKDxAA_GL), mean = 0, sd = 1e-5),decreasing = T))
