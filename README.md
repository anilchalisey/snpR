The snpR package is designed to identify trait-associated SNPs from the
GWAS catalog, and use this list of SNPs to:

  - Prune the list to generate independent SNPs - defined as those with
    low LD (r2 \< 0.1) to each other and beyond 500kb from each other.
  - For each independent (index) SNP determine all SNPs in LD (r2
    defined by user) from 1K Genomes Project Phase 3 CEU data.
  - Generate a BED file of the index SNP locus - the genomic region
    spanning the distance between the most-5’ and most-3’ SNPs in LD at
    that threshold.
  - Determine enrichment of the overlap between the SNPs and genomic
    regions of interest.

The website for snpR may be found
[here](https://anilchalisey.github.io/snpR/).

# Installation

## Command line tools

The workflow used by snpR requires the installation of PLINK. PLINK
binaries exist for both windows and linux-based systems. If PLINK is not
found in the executable path, then the necessary binary will be
automatically downloaded and used.

### R package dependencies

There are a number of dependencies to the snpR package detailed in the
`DESCRIPTION` file. These will be installed automatically.

### Install snpR

The package is easily installed as follows:

``` r
devtools::install_github("anilchalisey/snpR", build_vignettes = TRUE)
```

## Required data

In order to perform analysis, the package requires PLINK binary genotype
formatted 1000 Genomes Project data which is available from the [Broad
Institute at this
link](http://www.broadinstitute.org/mpg/depict/depict_download/1kg/1000_genomes_project_phase3_CEU.tar.gz).
This may be done within R using the following commands (WARNING - this
file is very
large):

``` r
download.file("http://www.broadinstitute.org/mpg/depict/depict_download/1kg/1000_genomes_project_phase3_CEU.tar.gz", "plink_data.tar.gz", mode = "wb")
untar("plink_data.tar.gz", exdir = "plink_data")
```

Also required is precomputed LD r2 boundaries for each 1000 genomes
project phase 3 SNP. I have provided this in a Rdata format with SNPs
binned by MAF. The file is too large to be included with this package
but may be downloaded from
[here](https://www.dropbox.com/s/s9roomzkq98zvhb/collection.rda?dl=0).

## Usage

The commands below demonstrate the usage of this package. In this
example, we are interested in GWAS hits identified from studies on
Alzheimer’s disease.

``` r
# Generate a list of SNPS
gwas <- get_GWAS(catalog = "download", filter.criteria = "alzheimer")

# Prune the list of SNPs (N.B. If plink is installed, then specify the path
# to the plink executable.  Otherwise, it will be automatically downloaded and 
# installed. It is assumed, in this example, that the genotype data was downloaded
# as suggested above
pruned_snps <- prune_SNPs(plink = NULL, snps = gwas, 
                          genotypeData = "plink_data")

# The downloaded PLINK binary is located at "plink/plink.exe" (on a Windows machine)
# If using Linux or MAC OS then set this appropriately
plink <- "plink/plink.exe"

# Identify all SNPs in LD (r2 > 0.8) to the independent SNPs
ld_snps <- ld_SNPs(plink = plink, genotypeData = "plink_data", 
                   independent_snps = pruned_snps, r2 = 0.8)

# Create a BED file of the independent SNP loci
loci_snps <- ld_SNPs2BED(ld_snps, "alzheimer_snp_loci.bed")

# Create a set of matching SNP loci - in this example, we only use 100 permutations.
# In practice, at least 1000 should be used
matching_loci <- matching_SNPs(snpdatabase = "collection.rda", 
                               nperm = 100, ld.snps = ld_snps)

# Finally we run the enrichment analysis.  For this we use a BED file containing
# regions generated from a ChIP-seq analysis.
bed <- system.file("extdata", "sample_peaks.narrowPeak", package = "snpR")
snp_enrichment <- run_enrichment(reg = bed, ld.snps = ld_snps, matched.list = matching_loci)
```

Alternatively, the entire pathway may be run using a single call to the
function `snp_enrichment()`.

``` r

result <- snp_enrichment(catalog = "gwas_hg38.txt.gz", filter = "alzheimer", 
                         plink = "plink/plink.exe", nperm = 100, genotypeData = "plink_data", 
                         r2 = 0.8, outdir = "results", 
                         snpdatabase = "collection.rda", 
                         reg = system.file("extdata", "sample_peaks.narrowPeak", package = "snpR"))
```

## About snpR

The snpR package has been developed in the Chris O’Callaghan Group at
the Centre for Cellular and Molecular Biology, University of Oxford by
Anil Chalisey and Chris O’Callaghan.
