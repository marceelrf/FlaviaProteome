library(limma)
library(Biobase)

fData(expset)

limma::plotDensities(expset,legend = F)
boxplot(exprs(expset)[,1:16],
        main="Boxplot normalized Intensities")
expset_log2 <- expset


Biobase::exprs(expset_log2) <- log2(Biobase::exprs(expset))
limma::plotDensities(expset_log2,legend = F)

boxplot(exprs(expset_log2)[,1:16],
        main="Boxplot normalized Intensities")

eset_log_norm <- expset_log2
exprs(eset_log_norm) <- normalizeBetweenArrays(exprs(expset_log2),
                                               method = "quantile")
limma::plotDensities(eset_log_norm,legend = F)

boxplot(exprs(eset_log_norm)[,1:16],
        main="Boxplot normalized Intensities")
# DEP ---------------------------------------------------------------------
pData(eset_log_norm)

(design <- model.matrix(~ 0 + factor(pData(eset_log_norm)$Condition)))

colnames(design) <- c("AA","AA_CKD","CKD","nt")

fit <- lmFit(exprs(eset_log_norm),
             design)

# cm <- makeContrasts(b1 = "mv_CKD + vsm_CKD - mv_normal - vsm_normal", #Doença
#                     b2 = "mv_CKD + mv_normal - vsm_CKD - vsm_normal", #Fonte
#                     b3 = "mv_normal - vsm_normal", #Fonte - comparando normal
#                     b4 = "vsm_normal - vsm_CKD", #Doença - comparando nas VSM
#                     levels = design)
# Define the contrasts
cm <- makeContrasts(
  AA_vs_CKD = AA - CKD,
  AA_vs_nt = AA - nt,
  CKD_vs_nt = CKD - nt,
  AA_CKD_vs_AA = AA_CKD - AA,
  AA_CKD_vs_CKD = AA_CKD - CKD,
  AA_CKD_vs_nt = AA_CKD - nt,
  levels = design
)

comparison_names <- c(
  "AA_vs_CKD",
  "AA_vs_nt",
  "CKD_vs_nt",
  "AA_CKD_vs_AA",
  "AA_CKD_vs_CKD",
  "AA_CKD_vs_nt"
)

cm
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

topTable(fit2,
         coef = 1,#o primeiro contraste avaliado
         adjust.method = "BH",#Benjamin hocheback - Falsos positivos corrigidos
         number = Inf #numero de linhas para serem exibidas
)



for(i in 1:6){
  tmp_tab <- NULL
  tmp_tab <- topTable(fit2,
           coef = 1,#o primeiro contraste avaliado
           adjust.method = "BH",#Benjamin hocheback - Falsos positivos corrigidos
           number = Inf #numero de linhas para serem exibidas
  )
  
  tmp_tab %>% 
    rownames_to_column("PG.ProteinAccessions") %>% 
    mutate(Sig = case_when(
      logFC > 1 & P.Value < 0.05 ~ "Up",
      logFC < -1 & P.Value < 0.05 ~ "Down",
      TRUE ~ "NotSig"
    )) %>% 
    left_join(feature_data) %>% 
    writexl::write_xlsx(
      path = paste0("Output/",
                    comparison_names[i],".xlsx"))
  
}
