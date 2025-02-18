mapping <- read_tsv(file = "Data/idmapping_AA_CKD_vs_CKD.tsv")


AAxNT_tibble <- AAxNT %>% 
  separate_rows(UNIPROT,sep = ";") %>% 
  left_join(mapping, by = join_by(UNIPROT == From)) %>% 
  select(-SYMBOL,-ENTREZID) %>% 
  rename(Symbol = To)

AA_CKD_vs_CKD_tibble <- AA_CKD_vs_CKD %>% 
  left_join(mapping, by = join_by(UNIPROT == From)) %>% 
  select(-SYMBOL,-ENTREZID) %>% 
  rename(Symbol = To)

AA_CKDxAA_tibble <- AA_CKDxAA %>% 
  left_join(mapping, by = join_by(UNIPROT == From)) %>% 
  select(-SYMBOL,-ENTREZID) %>% 
  rename(Symbol = To)
