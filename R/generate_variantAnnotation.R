generate_variantAnnotation <- function(rowData, infoMetrics, variantLocation, variantConsequence, exacAlleleFreq){


  annotatedVCF_temp <- rowData %>% as.data.frame() %>% tibble::rowid_to_column("QUERYID") %>% # make a queryID column from rd row index
    dplyr::select(ExACname, chromosome = seqnames, start, end, REF, ALT, QUERYID) %>% # extract relevant columns from rd
    cbind(infoMetrics) %>%  # merge rd/rowData with infoMetrics
    dplyr::left_join(variantLocation, by = "QUERYID") %>% # merge variantLocation containing all variant location annotation
    dplyr::left_join(variantConsequence, by = "QUERYID") %>% # merge variantConsequence containing coding variant consequence
    dplyr::left_join(exacAlleleFreq, by = c("ExACname.x" = "ExACname")) # merge exacAlleleFreq by ExACname columns

  annotatedVCF_temp %>% head()

}
