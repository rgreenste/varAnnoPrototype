#' Annotate Coding Variants by Biological Consequence
#'
#' Extracts consequences for coding variants, returning only the most deleterious consequence per variant if there are multiple consequences.
#'
#' @param rowData GRanges object containing genomic coordinates for each variant and Reference (REF) and Alternate (ALT) alleles
#'
#' @return a tibble object containing information about the consequence of each variant in the CONSEQUENCE column along with additional information of interest about each variant including the reference allele codon and amino acid (REFCODON, REFAA), the alternative allele codon and amino acid (VARCODON, VARAA) and the variant name in the ExAC style
#' @export
#'
#' @examples \dontrun{extract_variantConsequence(rowData)}
#'
#' @importFrom rlang .data
extract_variantConsequence <- function (rowData){

  # check the class of input file
  if(!class(rowData) == "GRanges" )
    stop("Error: 'rowData' should be of class 'GRanges'. Please read in a vcf file with the importExpandedVCF function and extract the rowData with 'extract_rowData()' before proceeding")

  # check for correct genome
  if(!GenomeInfoDb::genome(rowData)[1] == "hg19")
    stop("Error: For this release only 'hg19' is currently supported")

  # check for required metadata column containing ALT allele
  if(!"ALT" %in% names(GenomicRanges::mcols(rowData)))
    stop("Error: alternative allele not included in this GRanges object. Generate rowData from ExpandedVCF file with 'extract_rowData().")

  message("Predicting the effect of coding variants")

  # predict coding variants using rowData, txdb, Hspaiens BSGenome object, and alt accessor of vcf
  coding <- VariantAnnotation::predictCoding(rowData, # rowRanges from vcf renamed to have UCSC style chromosomes
                                             TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, #UCSC txdb object containing feature annotations
                                             seqSource = BSgenome.Hsapiens.UCSC.hg19::Hsapiens, # BSGenome object containing hg19 human genome sequence
                                             rowData$ALT) # alternative allele from the rowData

  # explicitly set the levels of consequence by order of severity - frameshift/nonsense are similar without more specific info
  coding$CONSEQUENCE <- factor(coding$CONSEQUENCE, levels = c("frameshift", "nonsense", "nonsynonymous", "synonymous"))

  message("Extracting most deleterious consequence per variant")

  # group by QUERYID (the row number from the vcf file / rd (rowData))
  # keep one row based on worst coding variant CONSEQUENCE, one row per gene
  # generate the ExAC name based on convetion
  coding_collapsed <- coding %>% data.frame() %>% dplyr::group_by(.data$QUERYID) %>% # group by QUERYID
    dplyr::slice_min(.data$CONSEQUENCE) %>% # keep one row per variant based on worst CONSEQUENCE, smallest level of factor
    dplyr::slice_min(.data$TXID) %>% # keep one entry per gene (based on smallest TXID) if found in multiple transcripts
    dplyr::mutate(ExACname = paste(stringr::str_remove(.data$seqnames, "chr"), .data$start, .data$REF, .data$ALT, sep="-")) %>% # add a variant name in the ExAC style chrom#-position-REF-ALT
    dplyr:: select(.data$QUERYID, .data$GENEID, .data$CONSEQUENCE, .data$REFCODON, .data$VARCODON, .data$REFAA, .data$VARAA, .data$ExACname)  # keep columns of interest for simplicity

  # make sure there is only one row per QUERYID
  if(!table(coding_collapsed$QUERYID) %>% max() == 1)
    stop("Error: There is more than one row per variant. Subsequent joins will fail. Additional redundancies must be evaluated")

  coding_collapsed

   }
