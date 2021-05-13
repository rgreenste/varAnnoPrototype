extract_rowData <- function(expandedVCFfile) {
  # extract the rowRanges from the VCF file
  rd<-rowRanges(expandedVCFfile)
  rd
}
