# AA CKD vs CKD -----------------------------------------------------------

teste_ora_1 <- readxl::read_xlsx(path = "Output/ora/ora_result_AA_CKD_vs_CKD.xlsx")
teste_ora_2 <- readxl::read_xlsx(path = "Output2/enrichment/AA_CKDxCKD.xlsx")


length(teste_ora_1$Description)
length(teste_ora_2$Description)

teste_ora_1$Description[teste_ora_1$Description %in% teste_ora_2$Description]


teste_ora_1 %>% 
  separate_rows(geneID,sep = "\\/") %>% 
  pull(geneID) %>% 
  unique()

teste_ora_2 %>% 
  separate_rows(geneID,sep = "\\/") %>% 
  pull(geneID) %>% 
  unique()


# AA vs NT ----------------------------------------------------------------


teste_ora_1 <- readxl::read_xlsx(path = "Output/ora/ora_result_AAxNT.xlsx")
teste_ora_2 <- readxl::read_xlsx(path = "Output2/enrichment/AAxNT.xlsx")


length(teste_ora_1$Description)
length(teste_ora_2$Description)

teste_ora_1$Description[teste_ora_1$Description %in% teste_ora_2$Description]


teste_ora_1 %>% 
  separate_rows(geneID,sep = "\\/") %>% 
  pull(geneID) %>% 
  unique()

teste_ora_2 %>% 
  separate_rows(geneID,sep = "\\/") %>% 
  pull(geneID) %>% 
  unique()


# AA CKD vs AA ------------------------------------------------------------

teste_ora_1 <- readxl::read_xlsx(path = "Output/ora/ora_result_AA_CKDxAA.xlsx")
teste_ora_2 <- readxl::read_xlsx(path = "Output2/enrichment/AA_CKDxAA.xlsx")


length(teste_ora_1$Description)
length(teste_ora_2$Description)

teste_ora_1$Description[teste_ora_1$Description %in% teste_ora_2$Description]


teste_ora_1 %>% 
  separate_rows(geneID,sep = "\\/") %>% 
  pull(geneID) %>% 
  unique()

teste_ora_2 %>% 
  separate_rows(geneID,sep = "\\/") %>% 
  pull(geneID) %>% 
  unique()

