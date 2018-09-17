#' Create BED file of SNP loci
#'
#' @param ld_snps list object created by \code{ld_SNPs}
#' @param file character string specifying name of output BED file
#'
#' @return
#' a data.frame object
#'
#' @importFrom readr write_tsv
#' @importFrom data.table fwrite
#' @importFrom parallel detectCores
#'
#' @export

ld_SNPs2BED <- function(ld_snps, file = NA) {
  x <- ld_snps[[1]]
  results.bed <- x[, c(2, 6, 7, 1)]
  results.bed$chr <- paste0("chr", results.bed$chr)
  results.bed$locus_upstream_boundary <-
    gsub("^.*:", "", results.bed$locus_upstream_boundary)
  results.bed$locus_downstream_boundary <-
    gsub("^.*:", "", results.bed$locus_downstream_boundary)
  threads <- parallel::detectCores() - 1
  data.table::fwrite(x = results.bed, append = FALSE, file = file,
                     col.names = FALSE, nThread = threads, row.names = FALSE,
                     sep = "\t")
  return(results.bed)
}
