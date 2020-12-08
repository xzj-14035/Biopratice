#' @title  read.fasta
#' @description This function reads the fasta files, extract the sequence from
#' it and calculate the proportion of each kind of base in sequence.
#' @param fa_files A character vector in which all elements are the
#'     absolute path of fasta files.
#' @export
read.fasta <- function(fa_files){
  fa_sum <- vector()
  for (f in fa_files){
    f_temp <- readLines(f)
    fa_sum <- c(fa_sum, f_temp)
  }
  label.index <- grep('^>', fa_sum)
  seq.index <- label.index + 1
  seq.list <- fa_sum[seq.index]
  names(seq.list) <- sub('^>(\\w+)\\|?', '\\1', fa_sum[label.index])
  seq.len <- vector()
  result <- array(0, dim = c(length(label.index), 17))
  colnames(result) <- c('A', 'T', 'C', 'G', 'U', "R", 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', '-')
  rownames(result) <- names(seq.list)
  for (i in 1:length(seq.list)){
    s <- seq.list[i]
    s.list <- unlist(strsplit(toupper(s), split = ''))
    s.len <- length(s.list)
    names(s.len) <- names(s)
    seq.len <- c(seq.len, s.len)
    for (b in s.list){
      result[names(s), b] <- result[names(s), b] + 1
    }
    result[names(s), ] <- round(result[names(s), ]/s.len, 2)
  }
  structure(
    list(seq.list = seq.list,
         result = result,
         Labels = sub('>', '', fa_sum[label.index]),
         len = seq.len
    ),
    class = "readfasta"
  )
}
