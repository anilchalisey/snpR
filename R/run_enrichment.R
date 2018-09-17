#' Determine SNP enrichment within genomic loci
#' 
#' This function determines whether a set of SNPs are enriched, compared to 
#' matched SNPs, within a defined set of genomic regions.
#'
#' @param reg either a character vector to a BED file or a GRanges object
#' @param ld.snps list object created by \code{ld_SNPs()}
#' @param matched.list list of matched SNPs created by \code{matching_SNPs()}
#'
#' @export
#' 
#' @importFrom dplyr mutate select 
#' @importFrom GenomicRanges makeGRangesFromDataFrame countOverlaps
#' @import grDevices
#' @import graphics
#' @importFrom rtracklayer import
#' @importFrom utils write.table
#' @importFrom stats dnorm sd
#' @importFrom tools file_ext
#' @importFrom IRanges subsetByOverlaps

run_enrichment <- function(reg = NULL, ld.snps, matched.list) {
  
  ldPos <- end <- chr <- start <- ldSNP <- NULL
  
  # Determine the number of overlaps for the actual SNPs and then the matched 
  # SNPs
  if (is.character(reg)) {
    if (tools::file_ext(reg) == "narrowPeak") peaks <- peak2Granges(reg)
    if (tools::file_ext(reg) == "broadPeak") peaks <- peak2Granges(reg, broad = TRUE)
    if (tools::file_ext(reg) == "bed|BED") peaks <- rtracklayer::import(reg, format = "BED")
  } else {
    peaks <- reg
  }
  
  # actual overlap
  ld.partners.r2.GR <- ld.snps[[2]] %>% 
    dplyr::mutate(chr = paste0("chr", gsub(":.*$", "", ldPos)), 
           end = as.numeric(gsub("^.*:", "", ldPos)), 
           start = (end - 1)) %>%
    dplyr::select(chr, start, end, snp = ldSNP)
  ld.partners.r2.GR <- ld.partners.r2.GR[stats::complete.cases(ld.partners.r2.GR), ]
  ld.partners.r2.GR <- GenomicRanges::makeGRangesFromDataFrame(ld.partners.r2.GR, 
                                                               keep.extra.columns = T)
  oversnps <- sum(GenomicRanges::countOverlaps(ld.partners.r2.GR, peaks))
  
  # overlap of snp sets
  oversnps.matched <- lapply(matched.list, function(x) {
    sum(GenomicRanges::countOverlaps(x, peaks))
  })
  
  # Calculate p-value of overlap
  nperm <- length(matched.list)
  r <- sum(oversnps.matched > oversnps)
  p.value <- (r+1)/(nperm + 1) #r+1/number of permutations + 1
  
  # to get a p-value of 0.05 we would have had to observe...
  r.ns <- ceiling(0.05*(nperm + 1) - 1)
  o.ns <- unlist(oversnps.matched)[order(unlist(oversnps.matched), 
                                         decreasing = T)][r.ns]
  
  # Create a histogram
  try({
    graphics::plot.new()
    grDevices::pdf(file = "snpenrichment.png", 
                   width = 11.7, height = 8.3, 
                   pointsize = 9)
    graphics::par(mar=c(5.1, 4.1, 5.1, 2.1))
    g <- graphics::hist(unlist(oversnps.matched), breaks = round(oversnps, -1), 
              col = "lightgray", density = 10, xlab = "", 
              main = "",
              xlim = c(0, max(unlist(oversnps.matched), oversnps)+1), 
              axes = FALSE,
              ylim = c(0, max(round(table(unlist(oversnps.matched))), -1)))
    graphics::title(paste0("Distribution of number of overlaps for matched sets of SNPs\n (",
                 nperm, " permutations)"), line = +2.5)
    graphics::axis(side = 2, pos = 0.2)
    xfit <- seq(min(unlist(oversnps.matched)), 
                max(unlist(oversnps.matched)), 
                length = 100)
    yfit <- stats::dnorm(xfit, mean = mean(unlist(oversnps.matched)),
                  sd = stats::sd(unlist(oversnps.matched))) 
    yfit <- yfit*diff(g$mids[1:2])*length(unlist(oversnps.matched))
    graphics::lines(xfit, yfit, col="black", lwd=2)
    graphics::abline(v = oversnps, col = "darkgreen", lwd = 3)
    graphics::abline(v = o.ns, col = "red", lwd = 3)
    graphics::mtext(side = 1, text = "value for p-value = 0.05", font = 2, col = "red", 
          at = o.ns)
    graphics::mtext(side = 1, text = paste0(o.ns, " overlaps"), font = 2, col = "red",
          line = 1, at = o.ns)
    graphics::mtext(side = 3, font = 2, col = "darkgreen", at = oversnps, 
          text = paste0("observed = ", oversnps))
    graphics::mtext(side = 3, font = 2, col = "darkgreen", at = oversnps, line = 1, 
          paste0("p-value = ", format.pval(p.value, digits = 2)))
    grDevices::dev.off()
    grDevices::graphics.off()
  })
  
  # 5) Output the overlapping SNPs and annotate
  oversnps.list <- IRanges::subsetByOverlaps(ld.partners.r2.GR, peaks)
  names(oversnps.list) <- oversnps.list$snp
  
  utils::write.table(oversnps.list, file = "overlapping_snps.txt", 
              col.names = T, row.names = F, sep = "\t", quote = F)
  
  cat("\n\nAnalysis complete.")
  if (p.value < 0.05) {
    cat("\n\nAt a significance threshold of 0.05, there was enrichment for the
        inputted SNPs.\n")
  } else {
    cat("\n\nAt a significance threshold of 0.05, there was no enrichment for the
        inputted SNPs.\n")
  }
  cat("Number of overlaps = ", oversnps, "\n")
  cat("Fold-enrichment = ", oversnps/mean(unlist(oversnps.matched)), "\n")
  cat("p-value = ", p.value, "\n")
  
  unlink("tmp.txt")
}