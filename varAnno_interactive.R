# An interactive/hard coded R script to build variant annotation tool for one vcf file

# load dependencies
library(tidyverse)
library(VariantAnnotation)


# read in the vcf file to a ExpandedVCF object - one variant per line
vcf <- readVcf("inst/extdata/Challenge_data_(1).vcf", "hg19") %>% expand()

################################################################################
## examine the parameters of this object ##

# look at the object
vcf

# what coordinates are represented in this file
rowRanges(vcf)

# how many variants are represented here?
length(rowRanges(vcf))

# what info is in the vcf header?
info(header(vcf)) %>% as.data.frame()

# what info is in the genotype/FORMAT field?
geno(header(vcf))
################################################################################

