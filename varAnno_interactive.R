# An interactive/hard coded R script to build variant annotation tool for one vcf file

# load dependencies
library(tidyverse)
library(VariantAnnotation)


# read in the vcf file to a ExpandedVCF object - one variant per line
vcf <- readVcf("inst/extdata/Challenge_data_(1).vcf", "hg19") %>% expand()

################################################################################
#### examine the parameters of this object ####
################################################################################

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
#### Extract metrics from INFO field of vcf file ####
################################################################################

# these are the info metrics I want to extract
info_metrics <- info(vcf) %>% as.data.frame() %>%
  dplyr::mutate(variant_type = TYPE,
                total_read_depth = DP, #depth of sequencing coverage at site of variation
                variant_read_count = AO, #number of reads supporting the variant
                reference_read_count = RO, #number of of reads supporting the reference
                percent_variant_reads = 100*AO/DP) %>% #percentage of reads supporting the variant
  dplyr::select(variant_type, total_read_depth, variant_read_count, reference_read_count, percent_variant_reads)


