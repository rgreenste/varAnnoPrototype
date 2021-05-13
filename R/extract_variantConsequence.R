extract_variantConsequence <- function (rowData, expandedVCFfile){

  # predict coding variants using rowData, txdb, Hspaiens BSGenome object, and alt accessor of vcf
  coding <- VariantAnnotation::predictCoding(rowData, # rowRanges from vcf renamed to have UCSC style chromosomes
                                             TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, #UCSC txdb object containing feature annotations
                                             seqSource = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, # BSGenome object containing hg19 human genome sequence
                                             VariantAnnotation::alt(expandedVCFfile)) # alt accessor of EpandedVCF object used to extract the alternate allele for each variant
  coding
  }
