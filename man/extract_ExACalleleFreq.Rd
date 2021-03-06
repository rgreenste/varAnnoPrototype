% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_ExACalleleFreq.R
\name{extract_ExACalleleFreq}
\alias{extract_ExACalleleFreq}
\title{Extract Allele Frequencies from ExAC API}
\usage{
extract_ExACalleleFreq(rowData)
}
\arguments{
\item{rowData}{GRanges object containing genomic coordinates for each variant and Reference (REF) and Alternate (ALT) alleles}
}
\value{
a tibble object containing allele frequencies and ExAC variant names for variants present within the ExAC database
}
\description{
Query the ExAC API for variants based on their ExAC name and return allele frequencies of variants from that database
}
\examples{
\dontrun{extract_ExACalleleFreq(rowData)}

}
