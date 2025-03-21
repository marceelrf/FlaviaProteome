library(msigdbr)
library(enrichplot)

msigdbr::msigdbr_collections()
mmu_c2 <- msigdbr::msigdbr(species = "Mus musculus", category = c("C2"))
mmu_c5 <- msigdbr::msigdbr(species = "Mus musculus", category = c("C5"))

sigs <-
  mmu_c2 %>% 
  bind_rows(mmu_c5) %>% 
  filter(str_detect(gs_description,"renal|Inflamassome|phosphate|phosphorus|calcif"))


sigs %>% 
  group_by(gs_description) %>% 
  count()

sigs %>% 
  filter(str_detect(gs_description,"calcif")) %>% 
  View()

msig_data <- sigs %>% select(gs_name, gene_symbol)

genes_interesse_AAxNT <- AnnotationDbi::select(x = org.Mm.eg.db,
                                         keys = AAxNT_diff,
                                         columns = "SYMBOL",
                                         keytype = "UNIPROT") %>% 
  drop_na() %>% 
  pull(SYMBOL)

genes_interesse_AA_CKDxAA <- AnnotationDbi::select(x = org.Mm.eg.db,
                                               keys = AA_CKDxAA_diff,
                                               columns = "SYMBOL",
                                               keytype = "UNIPROT") %>% 
  drop_na() %>% 
  pull(SYMBOL)

genes_interesse_AA_CKD_vs_CKD <- AnnotationDbi::select(x = org.Mm.eg.db,
                                               keys = AA_CKD_vs_CKD_diff,
                                               columns = "SYMBOL",
                                               keytype = "UNIPROT") %>% 
  drop_na() %>% 
  pull(SYMBOL)

# Rodar ORA com clusterProfiler
ora_result_AAxNT <- enricher(genes_interesse_AAxNT, TERM2GENE = msig_data)
ora_result_AA_CKDxAA <- enricher(genes_interesse_AA_CKDxAA, TERM2GENE = msig_data)
ora_result_AA_CKD_vs_CKD <- enricher(genes_interesse_AA_CKD_vs_CKD, TERM2GENE = msig_data)

annot_gene <- AnnotationDbi::select(x = org.Mm.eg.db,
                      keys = AA_CKDxAA_diff,
                      columns = c("SYMBOL","ENTREZID","UNIPROT"),
                      keytype = "UNIPROT")
AA_CKDxAA <- AA_CKDxAA %>% 
  left_join(annot_gene)

annot_gene <- AnnotationDbi::select(x = org.Mm.eg.db,
                                    keys = AAxNT_diff,
                                    columns = c("SYMBOL","ENTREZID","UNIPROT"),
                                    keytype = "UNIPROT")
AAxNT <- AAxNT %>% 
  left_join(annot_gene)

annot_gene <- AnnotationDbi::select(x = org.Mm.eg.db,
                                    keys = AA_CKD_vs_CKD_diff,
                                    columns = c("SYMBOL","ENTREZID","UNIPROT"),
                                    keytype = "UNIPROT")
AA_CKD_vs_CKD <- AA_CKD_vs_CKD %>% 
  left_join(annot_gene)

AA_CKD_vs_CKD_GL <- AA_CKD_vs_CKD %>% 
  filter(!is.na(SYMBOL)) %>% 
  select(SYMBOL,logFC) %>% 
  deframe()
AAxNT_GL <- AAxNT %>% 
  filter(!is.na(SYMBOL)) %>% 
  select(SYMBOL,logFC) %>% 
  deframe()
AA_CKDxAA_GL <- AA_CKDxAA %>% 
  filter(!is.na(SYMBOL)) %>% 
  select(SYMBOL,logFC) %>% 
  deframe()

cnetplot(x = ora_result_AAxNT,foldChange = AAxNT_GL)
cnetplot(x = ora_result_AA_CKDxAA,foldChange = AA_CKDxAA_GL)
cnetplot(x = ora_result_AA_CKD_vs_CKD,foldChange = AA_CKD_vs_CKD_GL)
