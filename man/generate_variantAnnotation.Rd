% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_variantAnnotation.R
\name{generate_variantAnnotation}
\alias{generate_variantAnnotation}
\title{Generate Final Variant Annotation Dataframe for Export}
\usage{
generate_variantAnnotation(
  rowData,
  infoMetrics,
  variantLocation,
  variantConsequence,
  exacAlleleFreq
)
}
\arguments{
\item{rowData}{GRanges object containing genomic coordinates for each variant and Reference (REF) and Alternate (ALT) alleles}

\item{infoMetrics}{dataframe object containing the requested information for the INFO field of the vcf file including type of variant, depth of sequencing coverage at site of variation, number of reads supporting the variant, number of of reads supporting the reference and percentage of reads supporting the variant}

\item{variantLocation}{tibble object containing information about the location of each variant with respect to genomic features and if located within a gene additionally the Gene ID number for that annotation}

\item{variantConsequence}{a tibble object containing information about the consequence of each variant in the CONSEQUENCE column along with additional information of interest about each variant including the reference allele codon and amino acid (REFCODON, REFAA), the alternative allele codon and amino acid (VARCODON, VARAA) and the variant name in the ExAC style}

\item{exacAlleleFreq}{tibble object containing allele frequencies and ExAC variant names for variants present within the ExAC database}
}
\value{
Dataframe object containing annotation information about variants from an original vcf file. The following columns are included:
\describe{
\item{\code{ExACname}}{Variant name in the ExAC naming style chromsome-position-referenceAllele-alternateAllele}
\item{\code{chromosome}}{Chromosome location of variant}
\item{\code{start}}{Starting position on chromsome of variant}
\item{\code{end}}{Ending position on chromsome of variant}
\item{\code{REF}}{Reference allele at this postion. Maybe be one or more nucleotides}
\item{\code{ALT}}{Alternate allele at this postion. Maybe be one or more nucleotides}
\item{\code{variant_type}}{Type of variant. May be one of: \describe{
\item{\code{snp}}{single nucleotide polymorphism}
\item{\code{mnp}}{multiple nucleotide polymorphism}
\item{\code{ins}}{insertion}
\item{\code{del}}{deletion}
\item{\code{complex}}{complex events including composite insertion and substitution events}}
}
\item{\code{total_read_depth}}{depth of sequencing coverage at site of variation}
\item{\code{variant_read_count}}{number of reads supporting the variant}
\item{\code{reference_read_count}}{number of of reads supporting the reference}
\item{\code{percent_variant_reads}}{percentage of reads supporting the variant}
\item{\code{LOCATION}}{Variant location with respect to genomic features. May be one of: \describe{
\item{\code{coding}}{Variant is located within a coding sequence that can be translated into a polypeptide}
\item{\code{spliceSite}}{Variant is located within a splice site. This is defined as within the first or last two nucleotides of an intron}
\item{\code{promoter}}{Variant is located within the promoter region of a gene/transcript}
\item{\code{fiveUTR}}{Variant is located with the 5' untranslated region (UTR) of a transcript}
\item{\code{threeUTR}}{Variant is located with the 3' untranslated region (UTR) of a transcript}
\item{\code{intron}}{Variant is located within an intron}
\item{\code{intergenic}}{complex events including composite insertion and substitution events}}
}
\item{\code{GENEID_Entrez}}{Entrez Gene Identification numbers for variants annotated to a gene by being located within a transcript or promoter}
\item{\code{CONSEQUENCE}}{Biological consequence for variants located within coding regions. May be one of: \describe{
\item{\code{frameshift}}{Variant causes a frameshift in the resulting transcript. This leads to an offset in the codons with respect to the three nucleotide frame and will most likely change all subsequent amino acid sequences.}
\item{\code{nonsense}}{Variant causes change of codon to stop codon. This leads to a truncated polypeptide sequence.}
\item{\code{nonsynonymous}}{Variant causes a change of codon to that of a different amino acid residue. This may or may not have an effect of the function of the protein.}
\item{\code{synonymous}}{Variant causes a change of codon but that codon is degenerate/redundant to the codon of the reference allele. The translated amino acid will be unchanged.}}
}
\item{\code{REFCODON}}{Codon present with the reference allele}
\item{\code{VARCODON}}{Codon present with the variant allele}
\item{\code{REFAA}}{Amino acid encoded by the codon containing the reference allele}
\item{\code{VARAA}}{Amino acid encoded by the codon containing the variant allele}
\item{\code{ExAC_allele_freq}}{Allele frequency for the allele extracted from the ExAC API. Not all variants are present}
}
}
\description{
Generate Final Variant Annotation Dataframe for Export
}
\examples{
\dontrun{generate_variantAnnotation(rowData,
infoMetrics, variantLocation, variantConsequence, exacAlleleFreq)}

}
