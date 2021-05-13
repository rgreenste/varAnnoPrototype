extract_infoMetrics <- function(expandedVCFfile) {

  # extract various metrics from INFO field of vcf file
  VariantAnnotation::info(expandedVCFfile) %>% as.data.frame() %>%
    dplyr::mutate(variant_type = TYPE, # type of variant
                  total_read_depth = DP, #depth of sequencing coverage at site of variation
                  variant_read_count = AO, #number of reads supporting the variant
                  reference_read_count = RO, #number of of reads supporting the reference
                  percent_variant_reads = 100*AO/DP) %>% #percentage of reads supporting the variant
    dplyr::select(variant_type, total_read_depth, variant_read_count, reference_read_count, percent_variant_reads)

}

