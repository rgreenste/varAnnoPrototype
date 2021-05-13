extract_variantConsequence <- function (rowData, expandedVCFfile){

  # check the class of input file
  if(!class(rowData) == "GRanges" )
    stop("'rowData' should be of class 'GRanges'. Please read in a vcf file with the importExpandedVCF function and extract the rowData with 'extract_rowData()' before proceeding")

  # check for correct genome
  if(!GenomeInfoDb::genome(rowData)[1] == "hg19")
    stop("for this release only 'hg19' is currently supported")

  # check the class of input file
  if(!class(expandedVCFfile) == "ExpandedVCF" )
    stop("'expandedVCFfile' should be of class 'ExpandedVCF'. Please read in a vcf file with the importExpandedVCF function")

  # predict coding variants using rowData, txdb, Hspaiens BSGenome object, and alt accessor of vcf
  coding <- VariantAnnotation::predictCoding(rowData, # rowRanges from vcf renamed to have UCSC style chromosomes
                                             TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, #UCSC txdb object containing feature annotations
                                             seqSource = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, # BSGenome object containing hg19 human genome sequence
                                             VariantAnnotation::alt(expandedVCFfile)) # alt accessor of EpandedVCF object used to extract the alternate allele for each variant

  # explicitly set the levels of consequence by order of severity - frameshift/nonsense are similar without more specific info
  coding$CONSEQUENCE <- factor(coding$CONSEQUENCE, levels = c("frameshift", "nonsense", "nonsynonymous", "synonymous"))

  coding
   }
