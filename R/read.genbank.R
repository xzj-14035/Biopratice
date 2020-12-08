#' @title  read.genbank
#' @description This function reads the genbank files, extracts the sequence
#' from it, calculates the proportion of each kind of base in sequence
#' and generates a fasta file which contains the sequence in it.
#' @param gb_files A character vector in which all elements are the
#' absolute path of genbank files
#' @export
read.genbank <- function(gb_files){
  fa_files <- vector()
  for (gb_file in gb_files){
    gb <- readLines(gb_file)
    acc <- gsub("ACCESSION\\s+", '', gb[grep("ACCESSION", gb)])
    acc_file <- paste(acc, ".fasta", sep = '')
    accession <- paste(">", acc, sep = '')
    seq <- gb[(grep("ORIGIN", gb)+1):(grep('//', gb)-1)]
    seq <- gsub('\\d+', '', seq)
    seq <- gsub('\\s+', '', seq)
    seq <- paste(seq, collapse = '')
    seq.len <- length(unlist(strsplit(seq, split = '')))
    fasta <- c(accession, seq)
    gb_path <- unlist(strsplit(gb_file, split = '/'))
    fa_path <- paste('/', paste(gb_path[2:(length(gb_path)-1)], collapse = '/'), sep = '')
    writeLines(fasta, paste(fa_path, acc_file, sep = '/'))
    fa_files <- c(fa_files, paste(fa_path, acc_file, sep = '/'))
  }
  x <- read.fasta(fa_files)
  print(x)
}
