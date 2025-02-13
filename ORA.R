library(msigdbr)

msigdbr::msigdbr_collections()
mmu_c2 <- msigdbr::msigdbr(species = "Mus musculus", category = c("C2"))
mmu_c5 <- msigdbr::msigdbr(species = "Mus musculus", category = c("C5"))

sigs <-
  mmu_c2 %>% 
  bind_rows(mmu_c5) %>% 
  filter(str_detect(gs_description,"renal|Inflamassome|phosphate|phosphorus"))

sigs %>% 
  group_by(gs_description) %>% 
  count()

msig_data <- sigs %>% select(gs_name, gene_symbol)

genes_interesse <- AnnotationDbi::select(x = org.Mm.eg.db,
                                         keys = AAxNT_diff,
                                         columns = "SYMBOL",
                                         keytype = "UNIPROT") %>% 
  drop_na() %>% 
  pull(SYMBOL)

# Rodar ORA com clusterProfiler
ora_result <- enricher(genes_interesse, TERM2GENE = msig_data)
