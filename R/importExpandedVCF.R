importExpandedVCF <- function(vcfFile, genome = "hg19") {

  # check dependencies
  if(!require(dplyr))
    stop("the 'dplyr' package must be be installed first")
  if(!require(VariantAnnotation))
    stop("the 'variantAnnotation' package must be be installed first")

  # check file format of vcf
  if(!is.character(vcfFile))
    stop("'vcfFile' should be character")
  if(!(stringr::str_sub(vcfFile, start = -3, end = -1) %>% stringr::str_to_lower()) == "vcf")
    stop("'vcfFile' must be a vcf file")

  # check file exists then read in the vcf file and expand
  if(file.exists(vcfFile)) {
    VariantAnnotation::readVcf(vcfFile, genome) %>% VariantAnnotation::expand()
  } else {
    print("the file does not exist")
  }
}
