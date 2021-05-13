#' Import vcf file and generate ExpandedVCF object
#'
#'
#' @param vcfFile character
#' @param genome character
#'
#' @return ExpandedVCF
#' @export
#'
#' @examples importExpandedVCF("myVariants.vcf", "hg19")
importExpandedVCF <- function(vcfFile, genome = "hg19") {

  # check file format of vcf
  if(!is.character(vcfFile))
    stop("Error: 'vcfFile' should be character")
  if(!(stringr::str_sub(vcfFile, start = -3, end = -1) %>% stringr::str_to_lower()) == "vcf")
    stop("Error:'vcfFile' must be a vcf file")

  # check name of genome, must be one of possible human genome assemblies
  if(!genome %in% c("hg18", "hg19", "hg38"))
    stop("Error: 'genome' should be one of 'hg18', 'hg19', 'hg38'")

  # for this release only hg19 is currently supported
  if(!genome == "hg19")
    stop("Error: For this release only 'hg19' is currently supported")

  # check file exists then read in the vcf file and expand
  if(file.exists(vcfFile)) {
    VariantAnnotation::readVcf(vcfFile, genome) %>% VariantAnnotation::expand()
  } else {
    print("Error: The file does not exist. Please use a valid file.")
  }
}
