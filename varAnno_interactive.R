# An interactive/hard coded R script to build variant annotation tool for one vcf file

# load dependencies
library(tidyverse)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(httr)
library(jsonlite)

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


################################################################################
#### Adjust chromosome names in VCF to match TxDb #### 
################################################################################

# need to convert chromosome names from "1" to "chr1" style
# check which chromosomes have variants 
which(table(seqnames(rowRanges(vcf))) > 0)

# if only main chroms then proceed
rd<-rowRanges(vcf)
# keep only the standard chromosomes
rd<-keepSeqlevels(rd, c(seq(1:22), "X", "Y"))
# convert to UCSC styles
newStyle <- mapSeqlevels(seqlevels(rd), "UCSC")
rd <- renameSeqlevels(rd, newStyle)

# load the txdb object containing hg19 feature data
txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
txdb<-keepSeqlevels(txdb, seqlevels(rd))

# look at the shared chromsomes between rowData and txdb
intersect(seqlevels(rd) , seqlevels(txdb))

# make sure the seqlengths of the chromosomes in each are the same
stopifnot(mean(seqlengths(rd) == seqlengths(txdb)) == 1)

################################################################################
#### Annotate all variants to genomic feature/location #### 
################################################################################

# annotate variants based on txdb
allvar <- locateVariants(rd, txdb, AllVariants())

# look at it
head(allvar)

# what locations are represented
table(allvar$LOCATION)

# need to collapse allvar to keep most severe annotation per variant
length(allvar) == length(rd) 

# what is order of levels of LOCATION factor?
levels(allvar$LOCATION)

# explicitly reorder the levels of the LOCATION factor by order of assumed 'severity' if mutated
allvar$LOCATION <- factor(allvar$LOCATION, levels = c("coding", "spliceSite", "fiveUTR", "promoter", "threeUTR", "intron", "intergenic"))

# what is order of levels of LOCATION factor?
levels(allvar$LOCATION)

# group by QUERYID (the row number from the vcf file / rd (rowData))
# keep rows based on location type in order of worst "consequence" 
# keep one entry per gene (based on smallest TXID) if found in multiple transcripts
allvar_collapsed <- allvar %>% data.frame() %>% mutate(TXID_noNA = replace_na(TXID, "*")) %>% # convert NAs in TXID to *
  group_by(QUERYID) %>% slice_min(LOCATION) %>% # keep the smallest value (most "severe") of LOCATION factor per QUERYID
  slice_min(TXID_noNA) %>% # keep the smallest transcript based on TXID_noNA if variant found in more than on TXID
  dplyr::select(LOCATION, QUERYID, GENEID) # keep columns of interest for simplicity

# make sure there is only one row per QUERYID
stopifnot(table(allvar_collapsed$QUERYID) %>% max() == 1)

# look at it
allvar_collapsed %>% head()


################################################################################
#### determine CONSEQUENCE of coding variants #### 
################################################################################

# predict coding variants using rowData, txdb, Hspaiens BSGenome object, and alt accessor of vcf
coding <- predictCoding(rd, txdb, seqSource=Hsapiens, alt(vcf) )

# how many of each kind?
table(coding$CONSEQUENCE)

# explicitly set the levels of consequence - frameshift/nonsense are similar without more specific info
coding$CONSEQUENCE <- factor(coding$CONSEQUENCE, levels = c("frameshift", "nonsense", "nonsynonymous", "synonymous"))

# what is order of levels of CONSEQUENCE factor?
levels(coding$CONSEQUENCE)

# need to collapse coding to keep most severe annotation per variant
length(coding) == length(rd) 

# group by QUERYID (the row number from the vcf file / rd (rowData))
# keep one row based on worst coding variant CONSEQUENCE, one row per gene
# generate the ExAC name based on convetion
coding_collapsed <- coding %>% data.frame() %>% group_by(QUERYID) %>% # group by QUERYID 
  slice_min(CONSEQUENCE) %>% # keep one row per variant based on worst CONSEQUENCE, smallest level of factor
  slice_min(TXID) %>% # keep one entry per gene (based on smallest TXID) if found in multiple transcripts
  mutate(ExACname = paste(str_remove(seqnames, "chr"), start, REF, ALT, sep="-")) %>% # add a variant name in the ExAC style chrom#-position-REF-ALT
  dplyr:: select(QUERYID, GENEID, CONSEQUENCE, REFCODON, VARCODON, REFAA, VARAA, ExACname)  # keep columns of interest for simplicity

# make sure there is only one row per QUERYID
stopifnot(table(coding_collapsed$QUERYID) %>% max() == 1)

# look at it
coding_collapsed %>% head()


################################################################################
#### generate ExAC IDs for all variants in vcf/rd #### 
################################################################################

# make a new metadata column on the rowData in the ExAC variant naming style
mcols(rd)$ExACname <- paste(str_remove(seqnames(rd), "chr"), start(rd), rd$REF, rd$ALT, sep="-")



################################################################################
#### perform bulk query to ExAC database for all variants by ExAC ID #### 
################################################################################

# establish the ExAC base url for database queries
baseURL<-"http://exac.hms.harvard.edu/"

# using a POST call perform a bulk query of ExAC database using the variant names
res<- POST(paste0(baseURL, "/rest/bulk/variant/variant"), body = toJSON(rd$ExACname))

# convert the data from JSON format to list
dat<-content(res, "text") %>% fromJSON()             

# extract the ExAC allele frequencies from the dat list
exac_allele_freq<-sapply(dat, function(x){
  x$allele_freq
}) %>% 
  unlist(use.names = T) %>% 
  as_tibble(rownames = NA) %>% 
  dplyr::mutate(ExACname = rownames(.), ExAC_allele_freq = value) %>% 
  dplyr::select(ExACname, ExAC_allele_freq)


################################################################################
#### Merge Existing Data Objects #### 
################################################################################

# check lengths match - should be as they come from same ExpandedVCF object with no subsetting of rows
stopifnot(nrow(info_metrics) == length(rd))

annotatedVCF_test <- rd %>% as.data.frame() %>% tibble::rowid_to_column("QUERYID") %>% # make a queryID column from rd row index
  dplyr::select(ExACname, chromosome = seqnames, start, end, REF, ALT, QUERYID) %>% # extract relevant columns from rd
  cbind(info_metrics) %>%  # merge rd/rowData with info_metrics
  left_join(allvar_collapsed, by = "QUERYID") %>% # merge allvar_collapsed containing all variant location annotation
  left_join(coding_collapsed, by = "QUERYID") %>% # merge coding_collapsed containing coding variant consequence
  left_join(exac_allele_freq, by = c("ExACname.x" = "ExACname")) # merge exac_allele_freq by ExACname columns

# look at it
annotatedVCF_test %>% head(n=15)

# do some checking on the above joins
table(annotatedVCF_test$CONSEQUENCE, annotatedVCF_test$LOCATION) # make sure only coding changes have a consequence
stopifnot((table(annotatedVCF_test$CONSEQUENCE, annotatedVCF_test$LOCATION) %>% colSums() > 0) == c(1,rep(0,6)))
stopifnot(mean(annotatedVCF_test$GENEID.x == annotatedVCF_test$GENEID.y, na.rm = T) == 1) # are the GENEIDs from allvar and coding the same?
stopifnot(mean(annotatedVCF_test$ExACname.x == annotatedVCF_test$ExACname.y, na.rm = T) == 1) # are the ExAC names from different tables the same?
stopifnot(nrow(annotatedVCF_test) == length(rd)) # check that length of output is the same as rd/rowData
stopifnot(annotatedVCF_test$QUERYID == seq(1:length(rd))) # check that all QUERYID are represented

# select final columns and rename
annotatedVCF<-annotatedVCF_test %>% 
  dplyr::select(-QUERYID, -GENEID.y, -ExACname.y, # drop redundant columns and QUERYID which is just an index value
                                    ExACname = ExACname.x, # rename column for clarity
                                    GENEID_Entrez = GENEID.x) # rename column for clarity
# look at it
annotatedVCF %>% head()


# plan
# extract metrics from info section of vcf - done
# annotate all variants to region/feature - done
# annotate coding variants for consequence - done
# generate ExAC styles names for all variants - done
# query ExAC database for all variants - done
# extract allele_frequency for variants in ExAC database - done
# giant merge at the end - done
# incremental merge at each step?
# export to csv
