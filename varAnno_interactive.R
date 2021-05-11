# An interactive/hard coded R script to build variant annotation tool for one vcf file

# load dependencies
library(tidyverse)
library(VariantAnnotation)


# read in the vcf file to a ExpandedVCF object - one variant per line
vcf <- readVcf("inst/extdata/Challenge_data_(1).vcf", "hg19") %>% expand()



