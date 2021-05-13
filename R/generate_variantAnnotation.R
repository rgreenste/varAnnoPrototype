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

  ## perform checks on the above joins before proceeding ##

  # make sure only coding changes have a consequence
  if(!(all(table(annotatedVCF_temp$CONSEQUENCE, annotatedVCF_temp$LOCATION) %>% colSums() > 0) %in% c(1,rep(0,6))))
    warning("Error: Consequences detected for variants that are noncoding. Annotation discrepancy between 'variantLocation' and 'variantConsequence'. Proceed with caution as annoations may be incorrect.")

  # check that GENEIDs from variantLocation and variantConsequence the same?
  if(!mean(annotatedVCF_temp$GENEID.x == annotatedVCF_temp$GENEID.y, na.rm = T) == 1)
    warning("Error: Mismatch beween GENEID columns detected. Annotation discrepancy between 'variantLocation' and 'variantConsequence'. Proceed with caution as annoations may be incorrect.")

   # check that ExAC names from rowData and exacAlleleFreq match
  if(!mean(annotatedVCF_temp$ExACname.x == annotatedVCF_temp$ExACname.y, na.rm = T) == 1)
    warning("Error: Mismatch beween ExACname columns detected. Annotation discrepancy between 'rowData' and 'exacAlleleFreq'. Proceed with caution as annoations may be incorrect.")

   # check that length of output is the same as rd/rowData
  if(!nrow(annotatedVCF_temp) == length(rowData))
    warning("Error: 'rowData' input and function output have different lengths. Some variants have been dropped from original file. Proceed with caution as annotations may be incorrect.")

  # check that all QUERYID are represented
  if(!all(annotatedVCF_temp$QUERYID %in% seq(1:length(rowData))))
    warning("Error: Some variants have been dropped from original file. Proceed with caution as annotations may be incorrect.")

  # select final columns and rename
  annotatedVCF<-annotatedVCF_temp %>%
    dplyr::select(-QUERYID, -GENEID.y, -ExACname.y, # drop redundant columns and QUERYID which is just an index value
                  ExACname = ExACname.x, # rename column for clarity
                  GENEID_Entrez = GENEID.x) # rename column for clarity

  annotatedVCF
}
