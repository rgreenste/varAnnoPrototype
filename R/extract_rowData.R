#' Extract rowData from ExpandedVCF file
#'
#' Extract the rowData from ExpandedVCF as GRanges object and convert seqlevels/chromsome names to UCSC style
#'
#' @param expandedVCFfile object of ExpandedVCF class
#'
#' @return GRanges object containing genomic coordinates for each variant and in the metadata columns Reference (REF) and Alternate (ALT) alleles, as well as quality scores (QUAL) and filter status
#' @export
#'
#' @examples \dontrun{extract_rowData(expandedVCFfile)}
extract_rowData <- function(expandedVCFfile) {

  # check the class of input file
  if(!class(expandedVCFfile) == "ExpandedVCF" )
    stop("'expandedVCFfile' should be of class 'ExpandedVCF'. Please read in a vcf file with the importExpandedVCF function")

  # extract the rowRanges from the VCF file
  rd<-SummarizedExperiment::rowRanges(expandedVCFfile)

  # check that only main chromosomes have variants
  if(!all(names(which(table(GenomeInfoDb::seqnames(rd)) > 0)) %in% c(seq(1:22), "X", "Y")))
    warning('vcf file contains variants on nonstandard chromsomes - these will not be considered in this analysis')

  # keep only the standard chromosomes
  rd<-GenomeInfoDb::keepSeqlevels(rd, c(seq(1:22), "X", "Y"))

  # convert to UCSC styles
  newStyle <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(rd), "UCSC")
  rd <- GenomeInfoDb::renameSeqlevels(rd, newStyle)

  rd
}
