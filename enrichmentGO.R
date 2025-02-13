library(clusterProfiler)
library(org.Mm.eg.db)



# AA x NT (Asfotase Alfa x Non-Treated) -----------------------------------
# coef = 2

AAxNT <- topTable(fit2,
         coef = 2,#o contraste avaliado
         adjust.method = "BH",#Benjamin hocheback - Falsos positivos corrigidos
         number = Inf #numero de linhas para serem exibidas
)

AAxNT$UNIPROT <- rownames(AAxNT)

AAxNT %>% 
  mutate(Diff = case_when(
    logFC < -1 & P.Value < 0.05 ~ "Down",
    logFC > 1 & P.Value < 0.05 ~ "Up",
    TRUE ~ "Ns"
  )) %>% 
  count(Diff)

AAxNT_diff <- 
  AAxNT %>% 
  mutate(Diff = case_when(
    logFC < -1 & P.Value < 0.05 ~ "Down",
    logFC > 1 & P.Value < 0.05 ~ "Up",
    TRUE ~ "Ns"
  )) %>%
  filter(Diff != "Ns") %>% 
  pull(UNIPROT)

annot_gene <- AnnotationDbi::select(x = org.Mm.eg.db,
                                    keys = AAxNT_diff,
                                    columns = c("SYMBOL","ENTREZID","UNIPROT"),
                                    keytype = "UNIPROT")

AAxNT_diff_go <- enrichGO(gene = AAxNT_diff,
         OrgDb = org.Mm.eg.db,
         keyType = "UNIPROT",
         ont = "ALL",
         pvalueCutoff = 1,
         qvalueCutoff = 1,
         minGSSize = 1)

# AA+CKD x AA -------------------------------------------------------------
# coef = 4


AA_CKDxAA <- topTable(fit2,
                  coef = 4,#o contraste avaliado
                  adjust.method = "BH",#Benjamin hocheback - Falsos positivos corrigidos
                  number = Inf #numero de linhas para serem exibidas
)

AA_CKDxAA$UNIPROT <- rownames(AA_CKDxAA)

AA_CKDxAA_diff <-
  AA_CKDxAA %>% 
  mutate(Diff = case_when(
    logFC < -1 & P.Value < 0.05 ~ "Down",
    logFC > 1 & P.Value < 0.05 ~ "Up",
    TRUE ~ "Ns"
  ))  %>%
  filter(Diff != "Ns") %>% 
  pull(UNIPROT)

annot_gene <- AnnotationDbi::select(x = org.Mm.eg.db,
                                    keys = AA_CKDxAA_diff,
                                    columns = c("SYMBOL","ENTREZID","UNIPROT"),
                                    keytype = "UNIPROT")

AA_CKDxAA_diff_go <- enrichGO(gene = AA_CKDxAA_diff,
                          OrgDb = org.Mm.eg.db,
                          keyType = "UNIPROT",
                          ont = "ALL",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          minGSSize = 1)

# AA+CKD x CKD ------------------------------------------------------------
# coef = 5

AA_CKD_vs_CKD <- topTable(fit2,
                      coef = 5,#o contraste avaliado
                      adjust.method = "BH",#Benjamin hocheback - Falsos positivos corrigidos
                      number = Inf #numero de linhas para serem exibidas
)

AA_CKD_vs_CKD$UNIPROT <- rownames(AA_CKD_vs_CKD)

AA_CKD_vs_CKD_diff <- AA_CKD_vs_CKD %>% 
  mutate(Diff = case_when(
    logFC < -1 & P.Value < 0.05 ~ "Down",
    logFC > 1 & P.Value < 0.05 ~ "Up",
    TRUE ~ "Ns"
  )) %>%
  filter(Diff != "Ns") %>% 
  pull(UNIPROT)

annot_gene <- AnnotationDbi::select(x = org.Mm.eg.db,
                                    keys = AA_CKD_vs_CKD_diff,
                                    columns = c("SYMBOL","ENTREZID","UNIPROT"),
                                    keytype = "UNIPROT")

AA_CKD_vs_CKD_diff_go <- enrichGO(gene = AA_CKD_vs_CKD_diff,
                              OrgDb = org.Mm.eg.db,
                              keyType = "UNIPROT",
                              ont = "ALL",
                              pvalueCutoff = 1,
                              qvalueCutoff = 1,
                              minGSSize = 1)