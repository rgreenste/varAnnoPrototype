#' Extract Allele Frequencies from ExAC API
#'
#' Query the ExAC API for variants based on their ExAC name and return allele frequencies of variants from that database
#'
#' @param rowData GRanges object containing genomic coordinates for each variant and Reference (REF) and Alternate (ALT) alleles
#'
#' @return a tibble object containing allele frequencies and ExAC variant names for variants present within the ExAC database
#' @export
#'
#' @examples \dontrun{extract_ExACalleleFreq(rowData)}
#'
#' @importFrom rlang .data
extract_ExACalleleFreq <- function(rowData){

  # check the class of input file
  if(!class(rowData) == "GRanges" )
    stop("Error: 'rowData' should be of class 'GRanges'. Please read in a vcf file with the importExpandedVCF function and extract the rowData with 'extract_rowData()' before proceeding")

  # check for correct genome
  if(!GenomeInfoDb::genome(rowData)[1] == "hg19")
    stop("Error: For this release only 'hg19' is currently supported")

  message("Querying the ExAC API for information on variants in the vcf file")

  # using a POST call perform a bulk query of ExAC database using the variant names
  res<- httr::POST("http://exac.hms.harvard.edu/rest/bulk/variant/variant", body = jsonlite::toJSON(rowData$ExACname))

  message("Extracting content from POST request to ExAC API")

  # convert the data from JSON format to list
  dat<-httr::content(res, "text") %>% jsonlite::fromJSON()

  # check that the returned list dat is same length as rowData
  if(!length(dat) == length(rowData))
    stop("Error: Query and result do not have the same length")

  message("Extracting allele frequencies from ExAC API content")

  # extract the ExAC allele frequencies from the dat list
  exac_allele_freq<-sapply(dat, function(x){
    x$allele_freq
  }) %>%
    unlist(use.names = T) %>%
    tibble::as_tibble(rownames = NA) %>% # convert to tibble
    tibble::rownames_to_column(var = "ExACname") %>% # make the rownames into a new column
    dplyr::mutate(ExAC_allele_freq = .data$value) %>% # rename column
    dplyr::select(.data$ExACname, .data$ExAC_allele_freq) # keep desired columns

  exac_allele_freq
}
