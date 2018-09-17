#' @title Create a set of matching SNP loci
#'
#' @description This function creates a set of matching SNP loci using data from 
#' the 1000 Genomes Phase 3 Project.  The SNPs in each set are matched for MAF, 
#' number of genes in the locus and LD 'buddies' at r2 0.8.
#' 
#' @details The snpdatabase is an rdata object containing all snps in the 1000G 
#' Phase3 along with the MAF, number of genes in the locus boundary 
#' (defined by r2 0.8) and number of ld buddies, split by MAF bin.  It may be
#' downloaded from \url{https://www.dropbox.com/s/s9roomzkq98zvhb/collection.rda?dl=0}.
#'  
#' @param snpdatabase character string specifying path to the snp database collection.
#' @param nperm Number of permutations to perform (i.e. sets of matching SNPs to create)
#' @param ld.snps list object created by \code{ld_snps}
#' 
#' @importFrom data.table := data.table rbindlist setkey set setattr
#' @importFrom dplyr select mutate
#' @importFrom stats complete.cases
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @return List of GRanges whose length is equal to nperm
#' @export

matching_SNPs <- function(snpdatabase = NULL, nperm = 1000, ld.snps) {
  
  collection <- ldSNP <- ldPos <- snp_maf <- gene_count <- friends_ld08 <- NULL
  HGNC_nearest_gene_snpsnap <- dist_nearest_gene_snpsnap <- gene_count <- NULL
  rowid <- J <- snpID <- end <- chr <- start <- rsID <- NULL
  
  cat("\n\nGenerating matched SNPs for permutation testing.  SNPs will be 
      matched for MAF, number of genes in SNP locus, and LD 'buddies' at 
      r2 0.8.\n")
  cat("Number of permutations =", nperm)
  cat("\nLoading SNP database - this may take a while.  If this takes longer than 2-3 
      minutes, you may not have sufficient RAM to proceed.\n")
  load(snpdatabase)
  cat("\nDatabase loaded.\n")
  
  ld.partners.r2.dt <- data.table::data.table(ld.snps[[2]])
  ld.partners.r2.dt <- lapply(collection, function(x) {
    tmp <- merge(ld.partners.r2.dt, x, by.x = "ldPos", 
                 by.y = "snpID", keep.x = T) %>%
      dplyr::select(ldSNP, ldPos, snp_maf, gene_count, friends_ld08, 
                    nearest_gene = HGNC_nearest_gene_snpsnap, 
                    dist_nearest_gene = dist_nearest_gene_snpsnap)
  })
  
  cat("\nCreating sets of matched SNPs for each permutation - this step may be 
      time-consuming.\n")
  
  matched <- vector("list", 49)
  
  ################## TO DO: rewrite this as a foreach loop
  for (i in seq_along(ld.partners.r2.dt)) {
    x <- ld.partners.r2.dt[[i]]
    y <- list()
    for (j in seq_len(nrow(x))) {
      gc <- as.numeric(x[j, gene_count])
      ld <- as.numeric(x[j, friends_ld08])
      s1 <- collection[[i]][as.numeric(gene_count) < 1.3*gc & 
                              as.numeric(gene_count) > 0.7*gc]
      s2 <- s1[as.numeric(friends_ld08) < 1.3*ld & 
                 as.numeric(friends_ld08) > 0.7*ld]
      y[[j]] <- s2[sample(nrow(s2), nperm, replace = T), ]
    }
    matched[[i]] <- y
    names(matched)[i] <- names(collection)[i]
  }
  
  matched2 <- unlist(matched, recursive = F)
  # This bit where we take the ith row from each dataframe in matched2 and rbind
  # potentially takes forever!  I have finally come up with a data.table solution
  # that is about 1000x faster!!  
  
  # convert to data.table
  invisible(lapply(matched2, data.table::setattr, name = "class", 
                   value = c("data.table", "data.frame")))
  # make one big table
  bigdata <- data.table::rbindlist(matched2)
  # generate an index - the number of dataframes will be the number of 
  # permutations
  index <- as.character(seq_len(nperm))
  bigdata[, `:=`(rowid, index)]
  # set a key based on the row index
  data.table::setkey(bigdata, rowid)
  # split on this
  matched.list <- lapply(index, function(i, j, x) x[i = J(i)], x = bigdata)
  # Drop the `row id` column
  invisible(lapply(matched.list, function(x) data.table::set(x, j = "rowid", value = NULL)))
  # Convert back to data.frame
  invisible(lapply(matched.list, data.table::setattr, name = 'class', 
                   value = c('data.frame')))
  
  # Convert the matched SNPs into a GRanges object list
  matched.list <- 
    lapply(matched.list, function(x) {
      a <- x %>% dplyr::mutate(chr = paste0("chr", gsub(":.*$", "", snpID)), 
                        end = as.numeric(gsub("^.*:", "", snpID)), 
                        start = (end - 1)) %>%
        dplyr::select(chr, start, end, snp = rsID, 
                      nearest_gene = HGNC_nearest_gene_snpsnap, 
                      dist_nearest_gene = dist_nearest_gene_snpsnap)
      a <- a[stats::complete.cases(x), ]
      return(a)
    })
  
  matched.list.GR <- lapply(matched.list, function(x) {
    a <- as.data.frame(x)
    GenomicRanges::makeGRangesFromDataFrame(a,  keep.extra.columns = TRUE, 
                                            seqnames.field = "chr", start.field = "start",
                                            end.field = "end")
  })
  save(matched.list.GR, file = paste0("matched_SNPs_GR.rda"))
  return(matched.list.GR)
}