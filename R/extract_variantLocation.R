#' Annotate Variants Based on Genomic Feature Location
#'
#' @param rowData GRanges object containing genomic coordinates for each variant and Reference (REF) and Alternate (ALT) alleles
#'
#' @return tibble object containing information about the location of each variant with respect to genomic features and if located within a gene additionally the Gene ID number for that annotation
#' @export
#'
#' @examples \dontrun{extract_variantLocation(rowData)}
extract_variantLocation <- function (rowData) {

   # check for correct genome
  if(!genome(rd)[1] == "hg19")
    stop("for this release only 'hg19' is currently supported")

  # annotate variants based on txdb
  allvar <- VariantAnnotation::locateVariants(rd, # rowRanges from vcf renamed to have UCSC style chromosomes
                                              TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, #UCSC txdb object
                                              AllVariants()) # what variants to annotate - all of them

  # explicitly reorder the levels of the LOCATION factor by order of assumed 'severity' if mutated
  allvar$LOCATION <- factor(allvar$LOCATION, levels = c("coding", "spliceSite", "fiveUTR", "promoter", "threeUTR", "intron", "intergenic"))

  # group by QUERYID (the row number from the vcf file / rd (rowData))
  # keep rows based on location type in order of worst "consequence"
  # keep one entry per gene (based on smallest TXID) if found in multiple transcripts
  allvar_collapsed <- allvar %>% data.frame() %>% dplyr::mutate(TXID_noNA = tidyr::replace_na(TXID, "*")) %>% # convert NAs in TXID to *
    dplyr::group_by(QUERYID) %>% dplyr::slice_min(LOCATION) %>% # keep the smallest value (most "severe") of LOCATION factor per QUERYID
    dplyr::slice_min(TXID_noNA) %>% # keep the smallest transcript based on TXID_noNA if variant found in more than on TXID
    dplyr::select(LOCATION, QUERYID, GENEID) # keep columns of interest for simplicity

  # make sure there is only one row per QUERYID
  if(!table(allvar_collapsed$QUERYID) %>% max() == 1)
    stop("There is more than one row per variant. Subsequent joins will fail. Additional redundancies must be evaluated")

  allvar_collapsed
}
