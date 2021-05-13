extract_rowData <- function(expandedVCFfile) {
  # check the class of input file
  if(!class(expandedVCFfile) == "ExpandedVCF" )
    stop("'expandedVCFfile' should be of class 'ExpandedVCF'. Please read in a vcf file with the importExpandedVCF function")

  # extract the rowRanges from the VCF file
  rd<-rowRanges(expandedVCFfile)
  rd
}
