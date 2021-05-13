extract_rowData <- function(expandedVCFfile) {

  # check the class of input file
  if(!class(expandedVCFfile) == "ExpandedVCF" )
    stop("'expandedVCFfile' should be of class 'ExpandedVCF'. Please read in a vcf file with the importExpandedVCF function")

  # extract the rowRanges from the VCF file
  rd<-rowRanges(expandedVCFfile)

  # check that only main chromosomes have variants
  if(!all(names(which(table(GenomeInfoDb::seqnames(rd)) > 0)) %in% c(seq(1:22), "X", "Y")))
    warning('vcf file contains variants on nonstandard chromsomes - these will not be considered in this analysis')

  # keep only the standard chromosomes
  rd<-GenomeInfoDb::keepSeqlevels(rd, c(seq(1:22), "X", "Y"))

  # convert to UCSC styles
  newStyle <- GenomeInfoDb::mapSeqlevels(seqlevels(rd), "UCSC")
  rd <- GenomeInfoDb::renameSeqlevels(rd, newStyle)

  rd
}
