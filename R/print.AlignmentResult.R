#' @title print.AlignmentResult
#' @description This function prints the object returned by function
#' global_align.
#' @param x The object returned by function global_align
#' @param ... Other parameter
#' @method print AlignmentResult
#' @export
print.AlignmentResult <- function(x, ...){
  cat("Sequence X:\t")
  cat(x$Sequences[1])
  cat('\n')
  cat("Sequence Y:\t")
  cat(x$Sequences[2])
  cat('\n')
  cat('Scoring System:\t')
  cat(x$Scoring_system[1])
  cat('\t')
  cat('for match;')
  cat('\t')
  cat(x$Scoring_system[2])
  cat('\t')
  cat('for mismatch;')
  cat('\t')
  cat(x$Scoring_system[3])
  cat('\t')
  cat('for gap;')
  cat('\t\n\n')
  cat('Dynamic programing matrix:\n')
  print(x$Score_matrix, ...)
  cat('\n')
  cat('Alignment:\n')
  for (i in 1:length(x$Aln_result)){
    temp_result <- x$Aln_result[[i]]
    seq1 <- unlist(strsplit(temp_result[1], split = ''))
    seq2 <- unlist(strsplit(temp_result[2], split = ''))
    bars <- ''
    for ( i in 1:length(seq1)){
      if ((seq1[i]==seq2[i])&(seq1[i]!='-')&(seq2[i]!='-')){
        bars <- paste(bars, '|', sep = '')
      }
      else{
        bars <- paste(bars, ' ', sep = '')
      }
    }
    bars <- paste(bars, '\n', sep = '')
    cat(temp_result[1])
    cat('\n')
    cat(bars)
    cat(temp_result[2])
    cat('\n\n')
  }
  cat('Optimum alignment mat:\t')
  cat(x$Score_matrix[dim(x$Score_matrix)[1], dim(x$Score_matrix)[2]])
  cat('\n')
}
