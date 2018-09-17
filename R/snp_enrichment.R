#' Run the SNP enrichment pipeline
#' 
#' @description The snpR package is designed to identify trait-associated SNPs from the GWAS catalog, 
#' and use this list of SNPs to: 
#' 
#' (1) Prune the list to generate independent SNPs - defined as those with low LD (r2 < 0.1) to each 
#' other and beyond 500kb from each other; 
#' 
#' (2) For each independent (index) SNP determine all SNPs in LD (r2 defined by user) from 1K 
#' Genomes Project Phase 3 CEU data;
#' 
#' (3) Generate a BED file of the index SNP locus - the genomic region spanning the distance between 
#' the most-5' and most-3' SNPs in LD at that threshold; 
#' 
#' (4) Determine enrichment of the overlap between the SNPs and genomic regions of interest.
#' 
#' @section snpdatabase:
#' The snpdatabase is an rdata object containing all snps in the 1000G 
#' Phase3 along with the MAF, number of genes in the locus boundary 
#' (defined by r2 0.8) and number of ld buddies, split by MAF bin.  It may be
#' downloaded from \url{https://www.dropbox.com/s/s9roomzkq98zvhb/collection.rda?dl=0}.
#' 
#' @param catalog character; either a file path to a previously downloaded
#' version of the GWAS catalog or "download", to download the latest version.
#' [Default = "download"]
#' @param filter character; filtering criteria - list of terms, separated
#' bu '|' as a single character string, e.g. "diabetes | ischaemic | stroke". If no
#' filter is chosen, then all SNPs are outputted [Default = NA]
#' @param plink Path to PLINK (v1.9)
#' @param nperm Number of permutations to perform (i.e. sets of matching SNPs to create)
#' @param genotypeData Path to the directory containing the PLINK binary genotype formatted 1000 
#' Genomes Project data which may be downloaded as a tar.gz file from this link: 
#' \url{http://www.broadinstitute.org/mpg/depict/depict_download/1kg/1000_genomes_project_phase3_CEU.tar.gz}.
#' @param r2 numeric value 0-1 indicating the r2 value to calculate SNPs in LD
#' @param outdir output directory in which to save results
#' @param snpdatabase character string specifying path to the snp database collection.
#' @param reg either a character vector to a BED file or a GRanges object 
#' in which to check for snp enrichment
#'
#' @export

snp_enrichment <- function(catalog = "download", filter = NULL, 
                           plink = NULL, nperm = 100,
                           genotypeData = "plink_data", 
                           r2 = 0.8, outdir = "snp_result", snpdatabase = NULL,
                           reg = NULL) {
  
  gwas <- get_GWAS(catalog = catalog, filter.criteria = filter)
  pruned_snps <- prune_SNPs(plink = plink, snps = gwas, 
                            genotypeData = genotypeData)
  
  if (is.null(plink)) plink <- "plink/plink.exe"
  
  ld_snps <- ld_SNPs(plink = plink, genotypeData = genotypeData, 
                     independent_snps = pruned_snps, r2 = r2)
  
  loci_snps <- ld_SNPs2BED(ld_snps, "snp_loci.bed")
  
  matching_loci <- matching_SNPs(snpdatabase = snpdatabase, 
                                 nperm = nperm, ld.snps = ld_snps)
  
  snp_enrichment <- run_enrichment(reg = reg, ld.snps = ld_snps, matched.list = matching_loci)
  
  tomove <- c(list.files(pattern = "*snp_loci.bed"), "matched_SNPs_GR.rda", "gwas_hg38.txt.gz", 
              "overlapping_snps.txt", list.files(pattern = "plink"), "snpenrichment.png")
  
  lapply(tomove, function(x) file.move(from = x, to = file.path(outdir, x)))
  
  return(list(ld_snps, loci_snps, snp_enrichment))
}
  