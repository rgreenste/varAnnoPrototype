importExpandedVCF <- function(vcfFile, genome = "hg19") {

  # check dependencies
  if(!require(dplyr))
    stop("the 'dplyr' package must be be installed first")
  if(!require(VariantAnnotation))
    stop("the 'variantAnnotation' package must be be installed first")

  VariantAnnotation::readVcf(vcfFile, genome) %>% VariantAnnotation::expand()
}
