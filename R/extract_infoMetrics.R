#' Extract Various Metrics from INFO Field of vcf File
#'
#' Extract metrics from the INFO field of a vcf file. Returns a data.frame containing one row for each variant and columns containing information about the type of variant, read depth, number of reads supporting the variant allele, number of reads supporting the reference allele, and the percentage of reads supporting the variant
#'
#' @param expandedVCFfile object of ExpandedVCF class
#'
#' @return dataframe object containing the requested information for the INFO field of the vcf file
#' @export
#'
#' @examples \dontrun{extract_infoMetrics(expandedVCFfile)}
#'
#' @importFrom rlang .data
extract_infoMetrics <- function(expandedVCFfile) {

  # check the class of input file
  if(!class(expandedVCFfile) == "ExpandedVCF" )
    stop("Error: 'expandedVCFfile' should be of class 'ExpandedVCF'. Please read in a vcf file with the importExpandedVCF function")

  # check for appropriate info fields
  if(!all(c("TYPE", "DP", "AO", "RO") %in% (VariantAnnotation::info(expandedVCFfile) %>% as.data.frame() %>% names())))
    stop("Error: Not all requested INFO columns are present - vcf file cannot be annotated as requested")

  # extract various metrics from INFO field of vcf file
  VariantAnnotation::info(expandedVCFfile) %>% as.data.frame() %>%
    dplyr::mutate(variant_type = .data$TYPE, # type of variant
                  total_read_depth = .data$DP, #depth of sequencing coverage at site of variation
                  variant_read_count = .data$AO, #number of reads supporting the variant
                  reference_read_count = .data$RO, #number of of reads supporting the reference
                  percent_variant_reads = round(100*.data$AO/.data$DP,2)) %>% #percentage of reads supporting the variant
    dplyr::select(.data$variant_type, .data$total_read_depth, .data$variant_read_count, .data$reference_read_count, .data$percent_variant_reads)

}

