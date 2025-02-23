library(biomaRt)

ensembl <- useMart("ENSEMBL_MART_SNP")

ensembl <- useDataset("hsapiens_snp", mart = ensembl)

listDatasets(ensembl) |> View()
listFilters(ensembl) |> View()


snp_ids <- clipr::read_clip()

snp_info <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "allele","minor_allele"),
                  filters = "snp_filter", values = snp_ids, mart = ensembl)

snp_info |>
  dplyr::filter(chr_name %in% 1:22)
