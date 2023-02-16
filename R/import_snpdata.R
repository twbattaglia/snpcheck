#' Import data for SNP checkk
#'
#' @param input_folder A folder with VCF files from SNP data
#' @return Dataframe of formatted SNP data
#' @export
import_snpdata = function(input_folder) {

  # Load libraries
  library(tidyverse)
  library(progress)

  # Import files
  snp.list = list.files(input_folder, recursive = T, full.names = T, pattern = ".vcf")
  names(snp.list) = list.files(input_folder, pattern = ".vcf", recursive = T) %>%
    basename() %>%
    str_remove("_OpenArrayCalls.vcf") %>%
    str_remove("_.*")

  # Error checking (TODO)
  message("Importing SNP data...")

  import.res = snp.list %>%
    map_dfr(.id = "sampleId", .f = function(x){
      data.table::fread(x) %>%
        separate(col = contains("_"), into = c("GT", "GQ"), sep = "\\:", remove = F) %>%
        mutate(Genotype = case_when(
          GT == "0/0" ~ "A",
          GT == "1/0" ~ "B",
          GT == "0/1" ~ "B",
          GT == "1/1" ~ "C",
          GT == "1/2" ~ "D",
          GT == "2/1" ~ "D",
        )) %>%
        mutate(SNP = paste0("chr", `#CHROM`, ":", POS)) %>%
        select(SNP, Genotype)
    }) %>%
    mutate(Genotype = replace_na(Genotype, "-")) %>%
    pivot_wider(id_cols = sampleId, names_from = SNP, values_from = Genotype) %>%
    filter(!grepl("NTC", sampleId)) %>%
    mutate(CODE = str_extract(sampleId, "VE[:digit:]+|AKP[0-9]+|K[0-9]+")) %>%
    rowwise() %>%
    mutate(STUDY = stringr::str_split_fixed(sampleId, "-", 3)[1]) %>%
    mutate(temp = stringr::str_split_fixed(sampleId, "-", 3)[2]) %>%
    mutate(SAMPLE = str_match(pattern = regex("^[^0-9]*([0-9]+)(.*)"), string = temp)[2]) %>%
    select(-temp) %>%
    mutate(SAMPLE = str_pad(SAMPLE, 3, pad = "0")) %>%
    relocate(sampleId, STUDY, SAMPLE, CODE) %>%
    rename(Study_ID = sampleId)

  message("Finished importing SNP data...")

  return(import.res)
}


