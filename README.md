
<!-- README.md is generated from README.Rmd. Please edit that file -->

# varAnnoPrototype - a prototype variant annotation tool

<!-- badges: start -->
<!-- badges: end -->

vcf files are complex and can be difficult to interpret.
varAnnoPrototype is a prototype variant annotation tool that extracts
some relevant information from vcf files in addition to annotating the
variants contained within the file to the location and type of genomic
feature in which they are located as well as their biological
consequence. This package contains a set of functions to read in a vcf
file, extract the relevant information, perform the genomic annotations,
and output a dataframe containing the result which can be written to
file, viewed in external software, and shared with others.

## Installation

The development version of this package can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rgreenste/varAnnoPrototype")
```

## Example Usage

Here is an example of how you could use this package to annotate a vcf
file. Please not that at this time only vcf files for human data created
with the h19 reference are supported. Additional support for other
genome builds including coordinate liftover could be included in future
updates.

### Load the Library

After installing the package load the library with the following R code:

``` r
library(varAnnoPrototype)
## basic example code
```

### Read in a vcf File

A sample vcf file is contained within the extdata of this package.
Extract the filepath and load the file with the following code using the
`importExpandedVCF()` function:

``` r
# extract the path to the vcf file
vcfFile<-system.file("extdata", "Challenge_data_(1).vcf", package="varAnnoPrototype")

# import the vcf file with package function importExpandedVCF()
vcfTest<-importExpandedVCF(vcfFile, genome = "hg19")
#> Reading in the vcf file

# take a quick look at the result
vcfTest
#> class: ExpandedVCF 
#> dim: 7569 2 
#> rowRanges(vcf):
#>   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
#> info(vcf):
#>   DataFrame with 43 columns: NS, DP, DPB, AC, AN, AF, RO, AO, PRO, PAO, QR, ...
#> info(header(vcf)):
#>            Number Type    Description                                          
#>    NS      1      Integer Number of samples with data                          
#>    DP      1      Integer Total read depth at the locus                        
#>    DPB     1      Float   Total read depth per bp at the locus; bases in rea...
#>    AC      A      Integer Total number of alternate alleles in called genotypes
#>    AN      1      Integer Total number of alleles in called genotypes          
#>    AF      A      Float   Estimated allele frequency in the range (0,1]        
#>    RO      1      Integer Reference allele observation count, with partial o...
#>    AO      A      Integer Alternate allele observations, with partial observ...
#>    PRO     1      Float   Reference allele observation count, with partial o...
#>    PAO     A      Float   Alternate allele observations, with partial observ...
#>    QR      1      Integer Reference allele quality sum in phred                
#>    QA      A      Integer Alternate allele quality sum in phred                
#>    PQR     1      Float   Reference allele quality sum in phred for partial ...
#>    PQA     A      Float   Alternate allele quality sum in phred for partial ...
#>    SRF     1      Integer Number of reference observations on the forward st...
#>    SRR     1      Integer Number of reference observations on the reverse st...
#>    SAF     A      Integer Number of alternate observations on the forward st...
#>    SAR     A      Integer Number of alternate observations on the reverse st...
#>    SRP     1      Float   Strand balance probability for the reference allel...
#>    SAP     A      Float   Strand balance probability for the alternate allel...
#>    AB      A      Float   Allele balance at heterozygous sites: a number bet...
#>    ABP     A      Float   Allele balance probability at heterozygous sites: ...
#>    RUN     A      Integer Run length: the number of consecutive repeats of t...
#>    RPP     A      Float   Read Placement Probability: Phred-scaled upper-bou...
#>    RPPR    1      Float   Read Placement Probability for reference observati...
#>    RPL     A      Float   Reads Placed Left: number of reads supporting the ...
#>    RPR     A      Float   Reads Placed Right: number of reads supporting the...
#>    EPP     A      Float   End Placement Probability: Phred-scaled upper-boun...
#>    EPPR    1      Float   End Placement Probability for reference observatio...
#>    DPRA    A      Float   Alternate allele depth ratio.  Ratio between depth...
#>    ODDS    1      Float   The log odds ratio of the best genotype combinatio...
#>    GTI     1      Integer Number of genotyping iterations required to reach ...
#>    TYPE    A      String  The type of allele, either snp, mnp, ins, del, or ...
#>    CIGAR   A      String  The extended CIGAR representation of each alternat...
#>    NUMALT  1      Integer Number of unique non-reference alleles in called g...
#>    MEANALT A      Float   Mean number of unique non-reference allele observa...
#>    LEN     A      Integer allele length                                        
#>    MQM     A      Float   Mean mapping quality of observed alternate alleles   
#>    MQMR    1      Float   Mean mapping quality of observed reference alleles   
#>    PAIRED  A      Float   Proportion of observed alternate alleles which are...
#>    PAIREDR 1      Float   Proportion of observed reference alleles which are...
#>    MIN_DP  1      Integer Minimum depth in gVCF output block.                  
#>    END     1      Integer Last position (inclusive) in gVCF output record.     
#> geno(vcf):
#>   List of length 10: GT, GQ, GL, DP, DPR, RO, QR, AO, QA, MIN_DP
#> geno(header(vcf)):
#>           Number Type    Description                                           
#>    GT     1      String  Genotype                                              
#>    GQ     1      Float   Genotype Quality, the Phred-scaled marginal (or unc...
#>    GL     G      Float   Genotype Likelihood, log10-scaled likelihoods of th...
#>    DP     1      Integer Read Depth                                            
#>    DPR    A      Integer Number of observation for each allele                 
#>    RO     1      Integer Reference allele observation count                    
#>    QR     1      Integer Sum of quality of the reference observations          
#>    AO     A      Integer Alternate allele observation count                    
#>    QA     A      Integer Sum of quality of the alternate observations          
#>    MIN_DP 1      Integer Minimum depth in gVCF output block.
```

### Extract rowData from ExpandedVCF File

Next extract the rowData from the the ExpandedVCF object using the
`extract_rowData()` function. rowData contains the genomic coordinates
of the variants along with the Reference (REF) and Alternate (ALT)
alleles. A column contain the variant name in the ExAC notation style
(chromosome-position-REF-ALT) is created from the information in
rowData:

``` r
# extract info from the ExpandedVCF object
rd<-extract_rowData(vcfTest)

# take a quick look at the result
rd
#> GRanges object with 7569 ranges and 6 metadata columns:
#>          seqnames              ranges strand | paramRangeID            REF
#>             <Rle>           <IRanges>  <Rle> |     <factor> <DNAStringSet>
#>      [1]     chr1              931393      * |           NA              G
#>      [2]     chr1              935222      * |           NA              C
#>      [3]     chr1             1277533      * |           NA              T
#>      [4]     chr1             1284490      * |           NA              G
#>      [5]     chr1             1571850      * |           NA              G
#>      ...      ...                 ...    ... .          ...            ...
#>   [7565]     chrX           155011829      * |           NA              T
#>   [7566]     chrX           155011926      * |           NA              T
#>   [7567]     chrX           155233098      * |           NA              T
#>   [7568]     chrX           155239499      * |           NA              G
#>   [7569]     chrX 155239824-155239827      * |           NA           ACAA
#>                     ALT        QUAL      FILTER              ExACname
#>          <DNAStringSet>   <numeric> <character>           <character>
#>      [1]              T 2.17938e-13           .          1-931393-G-T
#>      [2]              A 1.68667e+04           .          1-935222-C-A
#>      [3]              C 2.81686e+04           .         1-1277533-T-C
#>      [4]              A 6.30056e+03           .         1-1284490-G-A
#>      [5]              A 0.00000e+00           .         1-1571850-G-A
#>      ...            ...         ...         ...                   ...
#>   [7565]              G     53198.2           .       X-155011829-T-G
#>   [7566]              C    100027.0           .       X-155011926-T-C
#>   [7567]              C     13266.6           .       X-155233098-T-C
#>   [7568]              A         0.0           .       X-155239499-G-A
#>   [7569]           GCAG     16743.5           . X-155239824-ACAA-GCAG
#>   -------
#>   seqinfo: 24 sequences from hg19 genome
```

### Extract Metrics from INFO Field of vcf File

Next extract some metrics from the INFO field of the vcf file via the
ExpandedVCF object using the `extract_infoMetrics()` function. Metrics
include type of variant, depth of sequencing coverage at site of
variation, number of reads supporting the variant, number of of reads
supporting the reference and percentage of reads supporting the variant.
Variant types include `snp` or single nucleotide polymorphism, `mnp` or
multinucleotide polymorphism, `ins` or insertion, `del` or deletion, or
`complex` which includes events such as composite insertion and
substitutions.

``` r
# extract info from the ExpandedVCF object
infoMetrics<-extract_infoMetrics(vcfTest)

# take a quick look at the result
head(infoMetrics)
#>   variant_type total_read_depth variant_read_count reference_read_count
#> 1          snp             4124                 95                 4029
#> 2          snp             1134                652                  480
#> 3          snp              786                786                    0
#> 4          snp              228                228                    0
#> 5          snp             4055                 94                 3961
#> 6          snp             3456                 26                 3430
#>   percent_variant_reads
#> 1                  2.30
#> 2                 57.50
#> 3                100.00
#> 4                100.00
#> 5                  2.32
#> 6                  0.75
```

### Annotate Variants Based on Genomic Feature Location

Next extract Location information about each variant with respect to
genomic features using the `extract_variantLocation()` function. Variant
positions are queried against a file of genomic features called
`TxDb.Hsapiens.UCSC.hg19.knownGene` from the R/Bioconductor package of
the same name. Information about this package can be found \[here\]
(<https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html>).
Information is returned about the location of each variant with respect
to genomic features and, if located within a gene, additionally the Gene
ID number for that annotation.

Possible genomic feature locations include `coding`, `spliceSite`,
`promoter`, `fiveUTR`, `threeUTR`, `intron`, and `intergenic`.

``` r
# extract the variant location with respect to genomic features
varLoc<-extract_variantLocation(rd)
#> Determining the location of variants with respect to genomic features
#> 'select()' returned many:1 mapping between keys and columns
#> 'select()' returned many:1 mapping between keys and columns
#> 'select()' returned many:1 mapping between keys and columns
#> 'select()' returned many:1 mapping between keys and columns
#> 'select()' returned many:1 mapping between keys and columns
#> 'select()' returned many:1 mapping between keys and columns
#> Extracting the feature annotation likely to have the most deleterious consequence if mutated

# take a quick look at the result
varLoc
#> # A tibble: 7,559 x 3
#> # Groups:   QUERYID [7,559]
#>    LOCATION   QUERYID GENEID
#>    <fct>        <int> <chr> 
#>  1 intergenic       1 <NA>  
#>  2 coding           2 57801 
#>  3 coding           3 1855  
#>  4 fiveUTR          4 1855  
#>  5 intron           5 984   
#>  6 intron           6 984   
#>  7 intron           7 984   
#>  8 coding           8 984   
#>  9 coding           9 984   
#> 10 coding          10 728642
#> # … with 7,549 more rows
```

### Annotate Coding Variants by Biological Consequence

Next extract information about the biological consequence of variants
within coding regions of the genome using the
`extract_variantConsequence()` function. Coding variant positions with
respect to genes are determined by querying against the
`TxDb.Hsapiens.UCSC.hg19.knownGene` genomic feature file. The
surrounding genomic sequence and nucleotide consequence of the variant
is determined by querying the hg19 human genome reference sequence found
the BioStrings object `BSgenome.Hsapiens.UCSC.hg19` from the
R/Bioconductor package of the same name. Information about this package
can be found \[here\]
(<https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html>).
Information is returned about the consequence of each variant in the
CONSEQUENCE column. Possible consequences include `frameshift`,
`nonsense`, `nonsynonymous`, and `synonymous`.

Additional information of interest about each variant is also returned,
including the reference allele codon and amino acid (REFCODON, REFAA),
the alternative allele codon and amino acid (VARCODON, VARAA) and the
variant name in the ExAC style.

``` r
# extract the variant consequence with respect to genomic features
varConseq<-extract_variantConsequence(rd)
#> Predicting the effect of coding variants
#> Extracting most deleterious consequence per variant

# take a quick look at the result
varConseq
#> # A tibble: 3,262 x 8
#> # Groups:   QUERYID [3,262]
#>    QUERYID GENEID CONSEQUENCE   REFCODON VARCODON REFAA VARAA ExACname     
#>      <int> <chr>  <fct>         <chr>    <chr>    <chr> <chr> <chr>        
#>  1       2 57801  nonsynonymous AGG      AGT      R     S     1-935222-C-A 
#>  2       3 1855   synonymous    CCA      CCG      P     P     1-1277533-T-C
#>  3       8 984    nonsynonymous GGG      AGG      G     R     1-1575784-C-T
#>  4       9 984    synonymous    GCG      GCA      A     A     1-1577180-C-T
#>  5      10 728642 synonymous    AAA      AAG      K     K     1-1635004-T-C
#>  6      11 728642 synonymous    CGG      CGT      R     R     1-1635749-C-A
#>  7      12 728642 synonymous    TAC      TAT      Y     Y     1-1636044-G-A
#>  8      13 728642 nonsynonymous GGG      AGG      G     R     1-1638994-C-T
#>  9      15 984    synonymous    GAA      GAG      E     E     1-1647814-T-C
#> 10      16 984    synonymous    CGA      CGG      R     R     1-1647871-T-C
#> # … with 3,252 more rows
```

### Extract Allele Frequencies from ExAC API

Next perform a bulk query to ExAC API for all variants based on their
ExAC name using the `extract_ExACalleleFreq()` function. Information
containing allele frequencies and ExAC variant names for variants
present within the ExAC database are returned.

``` r
# perform bulk query to ExAC database for all variants by ExAC ID and return allele frequencies
exacAF<-extract_ExACalleleFreq(rd)
#> Querying the ExAC API for information on variants in the vcf file
#> Extracting content from POST request to ExAC API
#> No encoding supplied: defaulting to UTF-8.
#> Extracting allele frequencies from ExAC API content

# take a quick look at the result
head(exacAF)
#> # A tibble: 6 x 2
#>   ExACname        ExAC_allele_freq
#>   <chr>                      <dbl>
#> 1 2-142567910-T-C          0.234  
#> 2 17-21207835-G-A          0.362  
#> 3 10-97990583-A-G          0.638  
#> 4 11-65664346-T-C          0.508  
#> 5 18-31325884-T-C          0.0198 
#> 6 14-38061572-G-A          0.00515
```

### Generate Final Variant Annotation for Export

Next a series of joins of the previously generated datasets are
performed to yield a final variant annotation file for export using the
`generate_variantAnnotation()` function.

``` r
# generate final annotation dataframe
varAnno<-generate_variantAnnotation(rd, infoMetrics, varLoc, varConseq, exacAF)

# take a quick look at the result
head(varAnno)
#>        ExACname chromosome   start     end REF ALT variant_type
#> 1  1-931393-G-T       chr1  931393  931393   G   T          snp
#> 2  1-935222-C-A       chr1  935222  935222   C   A          snp
#> 3 1-1277533-T-C       chr1 1277533 1277533   T   C          snp
#> 4 1-1284490-G-A       chr1 1284490 1284490   G   A          snp
#> 5 1-1571850-G-A       chr1 1571850 1571850   G   A          snp
#> 6 1-1572579-A-G       chr1 1572579 1572579   A   G          snp
#>   total_read_depth variant_read_count reference_read_count
#> 1             4124                 95                 4029
#> 2             1134                652                  480
#> 3              786                786                    0
#> 4              228                228                    0
#> 5             4055                 94                 3961
#> 6             3456                 26                 3430
#>   percent_variant_reads   LOCATION GENEID_Entrez   CONSEQUENCE REFCODON
#> 1                  2.30 intergenic          <NA>          <NA>     <NA>
#> 2                 57.50     coding         57801 nonsynonymous      AGG
#> 3                100.00     coding          1855    synonymous      CCA
#> 4                100.00    fiveUTR          1855          <NA>     <NA>
#> 5                  2.32     intron           984          <NA>     <NA>
#> 6                  0.75     intron           984          <NA>     <NA>
#>   VARCODON REFAA VARAA ExAC_allele_freq
#> 1     <NA>  <NA>  <NA>               NA
#> 2      AGT     R     S     6.611108e-01
#> 3      CCG     P     P     9.975711e-01
#> 4     <NA>  <NA>  <NA>     8.972994e-01
#> 5     <NA>  <NA>  <NA>     1.113759e-05
#> 6     <NA>  <NA>  <NA>     2.176891e-03
```

### Write to File

The `varAnno` dataframe can be written to file for viewing in external
software and sharing with others.

``` r
# export results to csv file
write.csv(varAnno, file = "annotated_challenge_data_vcf.csv", row.names = F)
```

## Get Help

Help information an documentation for functions in this package can be
accessed by typing a ‘?’ before the function name. An example is
included below:

``` r
# get help for a function
?generate_variantAnnotation
```
