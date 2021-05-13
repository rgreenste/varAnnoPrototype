generate_variantAnnotation <- function(rowData, infoMetrics, variantLocation, variantConsequence, exacAlleleFreq){

  # check the class of input file
  if(!class(rowData) == "GRanges" )
    stop(" Error:'rowData' should be of class 'GRanges'. Please read in a vcf file with the importExpandedVCF function and extract the rowData with 'extract_rowData()' before proceeding")

  # check the class of input file
  if(!class(infoMetrics) == "data.frame" )
    stop("Error: 'infoMetrics' should be of class 'data.frame'. Please generate 'infoMetrics' via the 'extract_infoMetrics()' function before proceeding.")

  # check lengths match - should be as they come from same ExpandedVCF object with no subsetting of rows
  if(!(nrow(infoMetrics) == length(rowData)))
    stop("Error: 'infoMetrics' and 'rowData' are of different lengths. Check that both were generated from the same vcf file with no subsetting.")

  # check the class of input file
  if(!"data.frame" %in% class(variantLocation))
    stop("Error: 'variantLocation' should be of class 'data.frame'. Please generate 'variantLocation' via the 'extract_variantLocation()' function before proceeding.")

  # check the class of input file
  if(!"data.frame" %in% class(variantConsequence))
    stop("Error: 'variantConsequence' should be of class 'data.frame'. Please generate 'variantConsequence' via the 'extract_variantConsequence()' function before proceeding.")

  # check the class of input file
  if(!"data.frame" %in% class(exacAlleleFreq))
    stop("Error: 'exacAlleleFreq' should be of class 'data.frame'. Please generate 'exacAlleleFreq' via the 'extract_ExACalleleFreq()' function before proceeding.")

  annotatedVCF_temp <- rowData %>% as.data.frame() %>% tibble::rowid_to_column("QUERYID") %>% # make a queryID column from rd row index
    dplyr::select(ExACname, chromosome = seqnames, start, end, REF, ALT, QUERYID) %>% # extract relevant columns from rd
    cbind(infoMetrics) %>%  # merge rd/rowData with infoMetrics
    dplyr::left_join(variantLocation, by = "QUERYID") %>% # merge variantLocation containing all variant location annotation
    dplyr::left_join(variantConsequence, by = "QUERYID") %>% # merge variantConsequence containing coding variant consequence
    dplyr::left_join(exacAlleleFreq, by = c("ExACname.x" = "ExACname")) # merge exacAlleleFreq by ExACname columns

  annotatedVCF_temp %>% head()

}
