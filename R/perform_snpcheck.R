#' Perform SNP check
#'
#' @param snp_table A number
#' @param study Study identifier to include
#' @param sample Sample identifier to include
#' @param usable_minimum Minimum usable SNP's to include in analysis
#' @param identical_minimum Minimum identical SNP's to determine similarity
#' @param distance_max Maximun distance of changing SNP loci
#' @return Dataframe of pairwise comparisons
#' @export
perform_snpcheck = function(snp_table, study = NULL, sample = NULL, usable_minimum = 18, identical_minimum = 18, distance_max = 4, check_blood = F) {

  # Load libraries
  library(tidyverse)
  library(progress)

  # Error checking (TODO)
  message("Importing global SNP data...")
  snps = snp_table

  # Perform sample-sample or sample-blood comparisons
  if(!check_blood){
    # Filter to keep a selected set of samples that have the same ID
    if(!is.null(sample)) {
      snps = snps %>%
        filter(STUDY %in% study) %>%
        filter(SAMPLE %in% sample)
      message(paste0("Found ", nrow(snps), " samples..."))
    } else {
      stop("Please enter a sample identifier (e.g 066, 079, 111)")
    }
  } else{
    snps = snps %>%
      filter(grepl("BL", Study_ID) | grepl("PBMC", Study_ID) | Study_ID == sample)
    message(paste0("Found ", nrow(snps), " samples..."))
  }

  # Remove unneeded columns
  snps = snps %>%
    mutate(uniqueId = paste(Study_ID, VENUMBER, sep = "_")) %>%
    relocate(uniqueId) %>%
    group_by(uniqueId) %>%
    mutate(row = row_number()) %>%
    mutate(uniqueId = paste0(uniqueId, "_" ,row)) %>%
    select(-row) %>%
    select(-STUDY:-VENUMBER, -Study_ID)

  # get combinations
  message("Getting pairwise combinations...")
  if(!check_blood){
    combinations = combn(snps$uniqueId, 2) %>%
      t() %>%
      as.data.frame()
  } else{
    combinations = combn(snps$uniqueId, 2) %>%
      t() %>%
      as.data.frame() %>%
      filter((grepl(sample, V1) | (grepl(sample, V2))))
   }

  # Iterate over each pairwise combination
  message(paste0("Computing ", nrow(combinations), " combinations"))
  all_comparisons <- list()
  all_tables <- list()
  pb <- progress_bar$new(total = nrow(combinations))
  for(i in 1:nrow(combinations)){
    pb$tick()

    # Get reduced table
    snps.sub = snps %>%
      filter(uniqueId %in% combinations[i,]) %>%
      pivot_longer(cols = -c("uniqueId")) %>%
      pivot_wider(names_from = uniqueId, values_from = value) %>%
      rename(locus = name)

    # Sample names
    sample1.name = colnames(snps.sub)[2]
    sample2.name = colnames(snps.sub)[3]

    # Compute total skipped (needs to be chr)
    total_skipped = 0

    # Compute total usable sites
    total_usable = snps.sub %>%
      rowwise() %>%
      filter(across(everything(), ~ .x != "-")) %>%
      nrow()

    # Computer total failed sites
    total_failed = 26 - total_usable

    # (A-A=0 , A-B=1 , A-C=2 , B-B=0 , C-B=1)
    snps.sub.compare = snps.sub %>%
      filter(across(everything(), ~ .x != "-")) %>%
      mutate(Dist = case_when(
        .[[2]] == "A" & .[[3]] == "A" ~ 0,
        .[[2]] == "B" & .[[3]] == "B" ~ 0,
        .[[2]] == "C" & .[[3]] == "C" ~ 0,

        .[[2]] == "A" & .[[3]] == "B" ~ 1,
        .[[2]] == "B" & .[[3]] == "A" ~ 1,

        .[[2]] == "A" & .[[3]] == "C" ~ 2,
        .[[2]] == "C" & .[[3]] == "A" ~ 2,

        .[[2]] == "C" & .[[3]] == "B" ~ 1,
        .[[2]] == "B" & .[[3]] == "C" ~ 1
      )) %>%
      mutate(Match = case_when(
        .[[2]] == "A" & .[[3]] == "A" ~ 1,
        .[[2]] == "C" & .[[3]] == "C" ~ 1,
        .[[2]] == "B" & .[[3]] == "B" ~ 1
      )) %>%
      mutate(Match = replace_na(Match, 0))

    # Compute identical
    total_identical = sum(snps.sub.compare$Match)
    percent_identical = total_identical / total_usable

    # Compute different
    total_different = total_usable - total_identical
    percent_different = total_different / total_usable

    # Compute total distance
    total_distance =  sum(snps.sub.compare$Dist)

    # Compute output
    if(total_usable < usable_minimum){
      FINAL = 'FAILED'
      FINAL_REASON = paste0("Not enough usable positions: ", total_usable, " < ", usable_minimum)
    } else if (total_distance > distance_max) {
      FINAL = "FAILED"
      FINAL_REASON = paste0("Distance too large: ", total_distance, " > ", distance_max)
    } else if (total_identical < identical_minimum) {
      FINAL = "FAILED"
      FINAL_REASON = paste0("Not enough identical calls: ", total_identical, " < ", identical_minimum)
    } else {
      FINAL = "PASS"
      FINAL_REASON = ""
    }

    # Output into dataframe
    res = tibble("Sample_1" = sample1.name,
                 "Sample_2" = sample2.name,
                 "Total_skipped" = total_skipped,
                 "Total_usable" = total_usable,
                 "Total_failed" = total_failed,
                 "Total_identical" = total_identical,
                 "Percent_identical" = percent_identical,
                 "Total_different" = total_different,
                 "Percent_different" = percent_different,
                 "Total_distance" = total_distance,
                 "DATA" = list(tibble(snps.sub.compare)),
                 "FINAL" = FINAL,
                 "Remarks" = FINAL_REASON)
    all_comparisons[[i]] = res
    all_tables[[i]] = snps.sub.compare
  }

  # Merge all results together
  message("Finished...")
  output = do.call(rbind, all_comparisons)
  return(output)
}


