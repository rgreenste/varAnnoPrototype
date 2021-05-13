extract_variantLocation <- function (rowData) {

   # check for correct genome
  if(!genome(rd)[1] == "hg19")
    stop("for this release only 'hg19' is currently supported")

  # annotate variants based on txdb
  allvar <- VariantAnnotation::locateVariants(rd, # rowRanges from vcf renamed to have UCSC style chromosomes
                                              TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, #UCSC txdb object
                                              AllVariants()) # what variants to annotate - all of them
  allvar
}
