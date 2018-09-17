#' Identify SNPs in LD
#'
#' This function identifies all SNPs in LD with a list of independent_snps
#' generated using the \code{pruneSNPs()} function
#'
#' @inheritParams prune_SNPs
#' @param independent_snps data.table object created by \code{prune_SNPs}
#' @param r2 numeric value 0-1 indicating the r2 value to calculate SNPs in LD
#'
#' @return
#' a data.table object
#'
#' @export
#' 
#' @importFrom readr write_tsv
#' @importFrom data.table fread setnames :=
#' @importFrom dplyr arrange mutate select summarise group_by

ld_SNPs <- function(plink, genotypeData, independent_snps, r2 = 0.8) {
  
  CHR_A <- BP_A <- CHR_B <- BP_B <- R2 <- SNP_B <- comb <- indexSNP <- ldPos <- NULL
  
  if (!is.character(independent_snps)) {
    readr::write_tsv(independent_snps, path = "tmp_snps.txt", col_names = TRUE)
    snps <- "tmp_snps.txt"
  }
  out <- paste0("ldpartners_", gsub("\\.", "", r2))
  gtdf <- file.path(genotypeData, 
                    "ALL.chr_merged.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes") 
  code.ld <- sprintf(
    "%s --bfile %s --r2 --ld-snp-list tmp_snps.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 %s --out %s",
    plink, gtdf, r2, out)
  system(code.ld)
  ld.partners <- data.table::fread(paste0(out, ".ld"), header = T,
                                   stringsAsFactors = FALSE)
  ld.partners.r2 <- ld.partners
  ld.partners.r2[, c("indexPos", "ldPos") := {
    list(paste(CHR_A, BP_A, sep = ":"),
         paste(CHR_B, BP_B, sep = ":"))
  }
  ]
  ld.partners.r2 <- ld.partners.r2[, c("SNP_A", "indexPos", "SNP_B", "ldPos", "R2")]
  oldnames <- colnames(ld.partners.r2)
  newnames <- c("indexSNP", "indexPos", "ldSNP", "ldPos", "r2")
  data.table::setnames(ld.partners.r2, old = oldnames, new = newnames)

  clump <- independent_snps$SP2
  clump <- lapply(clump, function(x) {
    x <- unlist(strsplit(x, ","))
    x <- gsub("\\(1\\)", "", x)
  })
  ld <- ld.partners[, 6:7]
  clump.ld <- lapply(clump, function(x) ld[ld$SNP_B %in% x, ])
  clump.ld <- lapply(clump.ld, function(x) {
    x %>% dplyr::arrange(desc(R2)) %>%
      dplyr::mutate(comb = paste0(SNP_B, " (", round(R2, 2), ")")) %>%
      dplyr::select(comb)
  })
  clump.ld <- lapply(clump.ld, function(x) {
    x <- do.call(c, x)
    x <- paste(x, collapse = ", ")
  })
  independent_snps$SP2_ld <- unlist(clump.ld)
  results <- ld.partners.r2 %>% dplyr::group_by(indexSNP) %>%
    dplyr::summarise(locus_up = min(ldPos), locus_down = max(ldPos))

  results <- merge(independent_snps, results, by.x = "SNP", by.y = "indexSNP",
                   keep.all = T)
  results.snps <- results[ , c(1, 2, 4, 5, 13, 14, 15)]
  oldnames <- colnames(results.snps)
  newnames <- c('snp_name', 'chr', 'pos','pvalue',
                'plink_ld_partners', 'locus_upstream_boundary',
                'locus_downstream_boundary')
  data.table::setnames(results.snps, oldnames, newnames)
  unlink(c("tmp_snps.txt", list.files(pattern = "ldpartners_")))
  return(list(results.snps, ld.partners.r2))
}
