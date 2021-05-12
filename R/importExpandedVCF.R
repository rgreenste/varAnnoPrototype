importExpandedVCF <- function(vcfFile, genome = "hg19") {
  VariantAnnotation::expand(VariantAnnotation::readVcf(vcfFile, genome))
}
