#' Download and parse GWAS catalog data
#'
#' @description This function will download the latest version of the GWAS
#' catalog from UCSC with hg38 co-ordinates, retain only those
#' SNPs with a genome-wide association p-value less than 5x10-8, and then
#' filter this further to create a data.frame containing SNPs determined
#' by the user-specified filtering criteria.
#'
#' @param catalog character; either a file path to a previously downloaded
#' version of the GWAS catalog or "download", to download the latest version.
#' [Default = "download"]
#' @param filter.criteria character; filtering criteria - list of terms, separated
#' bu '|' as a single character string, e.g. "diabetes | ischaemic | stroke". If no
#' filter is chosen, then all SNPs are outputted [Default = NA]
#'
#' @return
#' a data.frame of the gwas catalog pruned
#'
#' @importFrom dplyr select arrange filter distinct mutate
#' @importFrom utils download.file
#' @importFrom readr read_tsv
#'
#' @export

get_GWAS <- function(catalog = "download", filter.criteria = NA) {
  
  chrom <- chromEnd <- SNP <- marker <- P <- trait <- NULL

  if (catalog == "download") {
    utils::download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gwasCatalog.txt.gz",
                  "gwas_hg38.txt.gz")
    catalog <- "gwas_hg38.txt.gz"
  } else {
    if (!file.exists(catalog)) stop("The specified file does not exist.")
  }

  gwas.hg38 <- readr::read_tsv(catalog, col_names = FALSE, skip = 1)
  if (ncol(gwas.hg38) != 23) stop("The wrong number of columns are present in the
                                file - are you sure this is the correct file?")

  colnames(gwas.hg38) <- c("bin", "chrom", "chromStart", "chromEnd", "SNP",
                           "pubMedID", "author", "pubDate", "journal", "title",
                           "trait", "initSample", "replSample", "region",
                           "genes", "riskAllele", "riskAlFreq", "P",
                           "pValueDesc", "orOrBeta", "ci95", "platform",
                           "cnv")

  gwas.hg38.pruned <- gwas.hg38 %>%
    dplyr::mutate(marker = paste(gsub("chr", "", chrom), chromEnd, sep = ":")) %>%
    dplyr::select(SNP, marker, P, trait) %>%
    dplyr::arrange(P) %>%
    dplyr::filter(P < 5e-8)

  if (!is.na(filter.criteria)) {
    filter.criteria <- gsub (" \\| ", "\\|", filter.criteria)
    filter.criteria <- trimws(filter.criteria)
    gwas.hg38.pruned <- gwas.hg38.pruned[grepl(filter.criteria,
                                               gwas.hg38.pruned$trait,
                                               ignore.case = TRUE), ]
  }
  return(gwas.hg38.pruned)
}
