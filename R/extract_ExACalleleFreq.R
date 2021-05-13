extract_ExACalleleFreq <- function(rowData){

  # check the class of input file
  if(!class(rowData) == "GRanges" )
    stop("'rowData' should be of class 'GRanges'. Please read in a vcf file with the importExpandedVCF function and extract the rowData with 'extract_rowData()' before proceeding")

  # check for correct genome
  if(!GenomeInfoDb::genome(rowData)[1] == "hg19")
    stop("for this release only 'hg19' is currently supported")

  # using a POST call perform a bulk query of ExAC database using the variant names
  res<- httr::POST("http://exac.hms.harvard.edu/rest/bulk/variant/variant", body = jsonlite::toJSON(rowData$ExACname))

  # convert the data from JSON format to list
  dat<-httr::content(res, "text") %>% jsonlite::fromJSON()

  # check that the returned list dat is same length as rowData
  if(!length(dat) == length(rowData))
    stop("Error: Query and result do not have the same length")

  # extract the ExAC allele frequencies from the dat list
  exac_allele_freq<-sapply(dat, function(x){
    x$allele_freq
  }) %>%
    unlist(use.names = T) %>%
    tibble::as_tibble(rownames = NA) %>%
    dplyr::mutate(ExACname = rownames(.), ExAC_allele_freq = value) %>%
    dplyr::select(ExACname, ExAC_allele_freq)

  exac_allele_freq
}
