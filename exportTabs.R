library(tidyverse)

writexl::write_xlsx(x = ora_result_AA_CKD_vs_CKD@result,
                    path = "Output/ora/ora_result_AA_CKD_vs_CKD.xlsx")
writexl::write_xlsx(x = ora_result_AA_CKDxAA@result,
                    path = "Output/ora/ora_result_AA_CKDxAA.xlsx")
writexl::write_xlsx(x = ora_result_AAxNT@result,
                    path = "Output/ora/ora_result_AAxNT.xlsx")

#create networks
fn_create_net  <- function(x) {
  x@result %>%
    filter(pvalue < 0.05) %>% 
    select(ID,geneID) %>% 
    separate_rows(geneID,sep = "\\/")
}

fn_create_net(ora_result_AA_CKD_vs_CKD) %>% 
  write_tsv(file = "Output/ora/net_AA_CKD_vs_CKD.txt")
fn_create_net(ora_result_AA_CKDxAA) %>% 
  write_tsv(file = "Output/ora/net_AA_CKDxAA.txt")
fn_create_net(ora_result_AAxNT) %>% 
  write_tsv(file = "Output/ora/net_AAxNT.txt")
