#' Function to generate independent SNPs
#'
#' Given a list of SNPs, this function will utilise PLINK to generate independent
#' SNPs (r2 < 0.2 and at least 500kb apart).  It utilises the
#'
#' @param plink Path to PLINK (v1.9)
#' @param snps character string to file containing output from \code{get_GWAS} or
#' the name of the data.table object created by \code{getGWAS}
#' @param genotypeData Path to the directory containing the PLINK binary genotype formatted 1000 
#' Genomes Project data which may be downloaded as a tar.gz file from this link: 
#' \url{http://www.broadinstitute.org/mpg/depict/depict_download/1kg/1000_genomes_project_phase3_CEU.tar.gz}.
#'
#' @return
#' data.table object
#'
#' @importFrom readr write_tsv
#' @importFrom data.table fread
#'
#' @export

prune_SNPs <- function(plink = NULL, snps, genotypeData) {
  
  # Check that plink command works
  tryCatch({
    null <- system(command = paste(plink, "--help"), intern = TRUE)
  }, error = function(e) plink <<- install_plink())
  
  if (!is.character(snps)) {
    readr::write_tsv(snps, path = "tmp_snps.txt", col_names = TRUE)
    snps <- "tmp_snps.txt"
  }
  
  gtdf <- file.path(genotypeData, 
                    "ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes")
  code.clumping <- sprintf(
    "%s --bfile %s --clump-p1 5e-8 --clump-kb 500 --clump-r2 0.1 --clump %s --out clumpedSNPs",
    plink, gtdf, snps)
  system(code.clumping)
  snps.clumped <- data.table::fread("clumpedSNPs.clumped",
                                    header = T, stringsAsFactors = F)
  unlink(snps)
  unlink(list.files(pattern = "clumpedSNPs*"))
  return(snps.clumped)
}
