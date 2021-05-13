#' Generate Final Variant Annotation Dataframe for Export
#'
#'
#'
#' @param rowData GRanges object containing genomic coordinates for each variant and Reference (REF) and Alternate (ALT) alleles
#' @param infoMetrics dataframe object containing the requested information for the INFO field of the vcf file including type of variant, depth of sequencing coverage at site of variation, number of reads supporting the variant, number of of reads supporting the reference and percentage of reads supporting the variant
#' @param variantLocation tibble object containing information about the location of each variant with respect to genomic features and if located within a gene additionally the Gene ID number for that annotation
#' @param variantConsequence a tibble object containing information about the consequence of each variant in the CONSEQUENCE column along with additional information of interest about each variant including the reference allele codon and amino acid (REFCODON, REFAA), the alternative allele codon and amino acid (VARCODON, VARAA) and the variant name in the ExAC style
#' @param exacAlleleFreq tibble object containing allele frequencies and ExAC variant names for variants present within the ExAC database
#'
#' @return Dataframe object containing annotation information about variants from an original vcf file. The following columns are included:
#' \describe{
##'  \item{\code{ExACname}}{Variant name in the ExAC naming style chromsome-position-referenceAllele-alternateAllele}
##'  \item{\code{chromosome}}{Chromosome location of variant}
##'  \item{\code{start}}{Starting position on chromsome of variant}
##'  \item{\code{end}}{Ending position on chromsome of variant}
##'  \item{\code{REF}}{Reference allele at this postion. Maybe be one or more nucleotides}
##'  \item{\code{ALT}}{Alternate allele at this postion. Maybe be one or more nucleotides}
##'  \item{\code{variant_type}}{Type of variant. May be one of: \describe{
##'  \item{\code{snp}}{single nucleotide polymorphism}
##'  \item{\code{mnp}}{multiple nucleotide polymorphism}
##'  \item{\code{ins}}{insertion}
##'  \item{\code{del}}{deletion}
##'  \item{\code{complex}}{complex events including composite insertion and substitution events}}
##' }
##'  \item{\code{total_read_depth}}{depth of sequencing coverage at site of variation}
##'  \item{\code{variant_read_count}}{number of reads supporting the variant}
##'  \item{\code{reference_read_count}}{number of of reads supporting the reference}
##'  \item{\code{percent_variant_reads}}{percentage of reads supporting the variant}
##'  \item{\code{LOCATION}}{Variant location with respect to genomic features. May be one of: \describe{
##'  \item{\code{coding}}{Variant is located within a coding sequence that can be translated into a polypeptide}
##'  \item{\code{spliceSite}}{Variant is located within a splice site. This is defined as within the first or last two nucleotides of an intron}
##'  \item{\code{promoter}}{Variant is located within the promoter region of a gene/transcript}
##'  \item{\code{fiveUTR}}{Variant is located with the 5' untranslated region (UTR) of a transcript}
##'  \item{\code{threeUTR}}{Variant is located with the 3' untranslated region (UTR) of a transcript}
##'  \item{\code{intron}}{Variant is located within an intron}
##'  \item{\code{intergenic}}{complex events including composite insertion and substitution events}}
##' }
##'  \item{\code{GENEID_Entrez}}{Entrez Gene Identification numbers for variants annotated to a gene by being located within a transcript or promoter}
##'  \item{\code{CONSEQUENCE}}{Biological consequence for variants located within coding regions. May be one of: \describe{
##'  \item{\code{frameshift}}{Variant causes a frameshift in the resulting transcript. This leads to an offset in the codons with respect to the three nucleotide frame and will most likely change all subsequent amino acid sequences.}
##'  \item{\code{nonsense}}{Variant causes change of codon to stop codon. This leads to a truncated polypeptide sequence.}
##'  \item{\code{nonsynonymous}}{Variant causes a change of codon to that of a different amino acid residue. This may or may not have an effect of the function of the protein.}
##'  \item{\code{synonymous}}{Variant causes a change of codon but that codon is degenerate/redundant to the codon of the reference allele. The translated amino acid will be unchanged.}}
##' }
##'  \item{\code{REFCODON}}{Codon present with the reference allele}
##'  \item{\code{VARCODON}}{Codon present with the variant allele}
##'  \item{\code{REFAA}}{Amino acid encoded by the codon containing the reference allele}
##'  \item{\code{VARAA}}{Amino acid encoded by the codon containing the variant allele}
##'  \item{\code{ExAC_allele_freq}}{Allele frequency for the allele extracted from the ExAC API. Not all variants are present}
##' }
##'
#' @export
#'
#' @examples \dontrun{generate_variantAnnotation(rowData,
#' infoMetrics, variantLocation, variantConsequence, exacAlleleFreq)}
#'
#' @importFrom rlang .data
generate_variantAnnotation <- function(rowData, infoMetrics, variantLocation, variantConsequence, exacAlleleFreq){

  # check the class of input file
  if(!class(rowData) == "GRanges" )
    stop("Error:'rowData' should be of class 'GRanges'. Please read in a vcf file with the importExpandedVCF function and extract the rowData with 'extract_rowData()' before proceeding")

  # check the class of input file
  if(!class(infoMetrics) == "data.frame" )
    stop("Error: 'infoMetrics' should be of class 'data.frame'. Please generate 'infoMetrics' via the 'extract_infoMetrics()' function before proceeding.")

  # check lengths match - should be as they come from same ExpandedVCF object with no subsetting of rows
  if(!(nrow(infoMetrics) == length(rowData)))
    stop("Error: 'infoMetrics' and 'rowData' are of different lengths. Check that both were generated from the same vcf file with no subsetting.")

  # check the class of input file
  if(!"data.frame" %in% class(variantLocation))
    stop("Error: 'variantLocation' should be of class 'data.frame'. Please generate 'variantLocation' via the 'extract_variantLocation()' function before proceeding.")

  # check the class of input file
  if(!"data.frame" %in% class(variantConsequence))
    stop("Error: 'variantConsequence' should be of class 'data.frame'. Please generate 'variantConsequence' via the 'extract_variantConsequence()' function before proceeding.")

  # check the class of input file
  if(!"data.frame" %in% class(exacAlleleFreq))
    stop("Error: 'exacAlleleFreq' should be of class 'data.frame'. Please generate 'exacAlleleFreq' via the 'extract_ExACalleleFreq()' function before proceeding.")

  annotatedVCF_temp <- rowData %>% as.data.frame() %>% tibble::rowid_to_column("QUERYID") %>% # make a queryID column from rd row index
    dplyr::select(.data$ExACname, chromosome = .data$seqnames, .data$start, .data$end, .data$REF, .data$ALT, .data$QUERYID) %>% # extract relevant columns from rd
    cbind(infoMetrics) %>%  # merge rd/rowData with infoMetrics
    dplyr::left_join(variantLocation, by = "QUERYID") %>% # merge variantLocation containing all variant location annotation
    dplyr::left_join(variantConsequence, by = "QUERYID") %>% # merge variantConsequence containing coding variant consequence
    dplyr::left_join(exacAlleleFreq, by = c("ExACname.x" = "ExACname")) # merge exacAlleleFreq by ExACname columns

  ## perform checks on the above joins before proceeding ##

  # make sure only coding changes have a consequence
  if(!(all(table(annotatedVCF_temp$CONSEQUENCE, annotatedVCF_temp$LOCATION) %>% colSums() > 0) %in% c(1,rep(0,6))))
    warning("Warning: Consequences detected for variants that are noncoding. Annotation discrepancy between 'variantLocation' and 'variantConsequence'. Proceed with caution as annoations may be incorrect.")

  # check that GENEIDs from variantLocation and variantConsequence the same?
  if(!mean(annotatedVCF_temp$GENEID.x == annotatedVCF_temp$GENEID.y, na.rm = T) == 1)
    warning("Warning: Mismatch beween GENEID columns detected. Annotation discrepancy between 'variantLocation' and 'variantConsequence'. Proceed with caution as annoations may be incorrect.")

   # check that ExAC names from rowData and exacAlleleFreq match
  if(!mean(annotatedVCF_temp$ExACname.x == annotatedVCF_temp$ExACname.y, na.rm = T) == 1)
    warning("Warning: Mismatch beween ExACname columns detected. Annotation discrepancy between 'rowData' and 'exacAlleleFreq'. Proceed with caution as annoations may be incorrect.")

   # check that length of output is the same as rd/rowData
  if(!nrow(annotatedVCF_temp) == length(rowData))
    warning("Warning: 'rowData' input and function output have different lengths. Some variants have been dropped from original file. Proceed with caution as annotations may be incorrect.")

  # check that all QUERYID are represented
  if(!all(annotatedVCF_temp$QUERYID %in% seq(1:length(rowData))))
    warning("Warning: Some variants have been dropped from original file. Proceed with caution as annotations may be incorrect.")

  # select final columns and rename
  annotatedVCF<-annotatedVCF_temp %>%
    dplyr::select(-.data$QUERYID, -.data$GENEID.y, -.data$ExACname.y, # drop redundant columns and QUERYID which is just an index value
                  ExACname = .data$ExACname.x, # rename column for clarity
                  GENEID_Entrez = .data$GENEID.x) # rename column for clarity

  annotatedVCF
}
