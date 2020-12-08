#' @title print.readfasta
#' @description This function prints the object returned by function read.fasta.
#' @param x The object returned by function read.fasta
#' @param ... Other parameter pass to print
#' @method print readfasta
#' @export
print.readfasta <- function(x, ...){
  cat(length(x$seq.list))
  cat(' DNA sequence in binary format stored in a list.\n\n')
  if (length(x$seq.list) == 1){
    cat('Sequence length:')
    cat(x$len)
    cat('\n\n')
    cat('Label:\n')
    cat(x$Labels)
    cat('\n\n')
  }
  else{
    max_len <- max(x$len)
    min_len <- min(x$len)
    mean_len <- round(mean(x$len), 2)
    cat('Mean sequence length: ')
    cat(max_len)
    cat('\n\t')
    cat('Shortest sequence: ')
    cat(min_len)
    cat('\n\t')
    cat('Longest sequence: ')
    cat(max_len)
    cat('\n\n')
    cat('Label:\n')
    for (label in x$Labels){
      Index <- which(x$Labels == label)
      if (Index >= 4){
        break
      }
      else{
        cat(label)
        cat('\n')
      }
    }
  }
  cat('\nBase composition:\n')
  if (dim(x$result)[1] >= 4){
    print(x$result[1:4, ])
  }
  else{
    print(x$result)
  }
}
