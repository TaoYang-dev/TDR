#' Add cushion for sequencing data.
#'
#' Add a cushion for sequencing data which reduces the discrepancy of read counts that
#' has no biological meaning but greatly increases the discordance in ranks of the data.
#'
#' @param dat a matrix that has two columns, each of them represent a replicate of data,
#' i.g., sequencing counts for genomic loci.
#' @param cush an integer, the differece in sequencing counts that are deemed the same.
#' For example, if a cush of 2 is chosen, the read count of 4 at a loci of replicate 1 and
#' the read count of 6 at the same loci of replicate 2 are deemed the same. The 6 in replicate 2
#' is corrected to 4.
#' @return This function returns a matrix of two colunms of data, i.g., seqencing counts. Each
#' column is a vector a counts after correction to eliminate the discordance between two
#' replicates in ranks.
#' @export
#' @examples
#' data(Chipseq_TF)
#' cush <- 2
#' seqcush(Chipseq_TF, cush)

seqcush = function(dat, cush) {
  if (cush < 0 && !is.numeric(cush))
    stop("invalid cushion : alpha\n")
  else{
  idx = which(abs(dat[,1]-dat[,2])<=cush)
  dat[idx,1] = dat[idx,2]
  }
  return(dat)
}
