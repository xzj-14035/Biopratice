#' @title  global_aln
#' @description This function reads the fasta files, extract the sequence from
#'     it and calculate the proportion of each kind of base in sequence.
#' @param X A sequence to be aligned
#' @param Y A sequence to be aligned
#' @param Scoring_system A vector containing the score of match, mismatch and
#' inserting a gap
#' @export
global_aln <- function(X, Y, Scoring_system){
  Scoring_system <- Scoring_system
  names(Scoring_system) <- c('match', 'mismatch', 'gap')
  seq1 <- X
  seq2 <- Y
  seq1 <- unlist(strsplit(seq1, split = ''))
  seq2 <- unlist(strsplit(seq2, split = ''))
  seq1.len <- length(seq1)
  seq2.len <- length(seq2)
  score_matrix <- matrix(0, nrow = seq1.len + 1, ncol = seq2.len + 1)
  for (Row in 1:(seq1.len+1)){
    for (Col in 1:(seq2.len+1)){
      if (Row == 1){
        if (Col != 1){score_matrix[Row, Col] <- score_matrix[Row, Col-1] + Scoring_system['gap']}
      }
      else{
        if (Col == 1){score_matrix[Row, Col] <- score_matrix[Row-1, Col] + Scoring_system['gap'];next;}
        else{
          insert1 <- score_matrix[Row, Col] <- score_matrix[Row-1, Col] + Scoring_system['gap']
          insert2 <- score_matrix[Row, Col] <- score_matrix[Row, Col-1] + Scoring_system['gap']
          if (seq1[Row-1] == seq2[Col-1]){mis_or_match <- score_matrix[Row-1, Col-1] + Scoring_system['match']}
          else{mis_or_match <- score_matrix[Row-1, Col-1] + Scoring_system['mismatch']}
          score_matrix[Row, Col] <- max(insert1, insert2, mis_or_match)
        }
      }
    }
  }
  trace_back <- function(temp_seq1='', temp_seq2='', Row=seq1.len+1, Col=seq2.len+1, aln_result = list()){
    if (Row == 1&Col == 1){
      result <- list(c(temp_seq1, temp_seq2))
      aln_result <- append(aln_result, result)
      aln_result
    }
    else{
      if(Row>1&Col>1){
        if ((seq1[Row-1] == seq2[Col-1])&(score_matrix[Row-1, Col-1] + Scoring_system['match']==score_matrix[Row, Col])){
          Seq1 <- paste(seq1[Row-1], temp_seq1, sep = '')
          Seq2 <- paste(seq2[Col-1], temp_seq2, sep = '')
          aln_result1 <- trace_back(Seq1, Seq2, Row = Row-1, Col = Col-1)
          aln_result <- append(aln_result, aln_result1)
        }
        if ((seq1[Row-1] != seq2[Col-1])&(score_matrix[Row-1, Col-1] + Scoring_system['mismatch']==score_matrix[Row, Col])){
          Seq1 <- paste(seq1[Row-1], temp_seq1, sep = '')
          Seq2 <- paste(seq2[Col-1], temp_seq2, sep = '')
          aln_result2 <- trace_back(Seq1, Seq2, Row = Row-1, Col = Col-1)
          aln_result <- append(aln_result, aln_result2)
        }
        if ((score_matrix[Row-1, Col] + Scoring_system['gap']==score_matrix[Row, Col])){
          Seq1 <- paste(seq1[Row-1], temp_seq1, sep = '')
          Seq2 <- paste('-', temp_seq2, sep = '')
          aln_result3 <- trace_back(Seq1, Seq2, Row = Row-1, Col = Col)
          aln_result <- append(aln_result, aln_result3)
        }
        if ((score_matrix[Row, Col-1] + Scoring_system['gap']==score_matrix[Row, Col])){
          Seq1 <- paste('-', temp_seq1, sep = '')
          Seq2 <- paste(seq2[Col-1], temp_seq2, sep = '')
          aln_result4<- trace_back(Seq1, Seq2, Row = Row, Col = Col-1)
          aln_result <- append(aln_result, aln_result4)
        }
      }
      else{
        if(Row>1){
          Seq1 <- paste(seq1[Row-1], temp_seq1, sep = '')
          Seq2 <- paste('-', temp_seq2, sep = '')
          aln_result5 <- trace_back(Seq1, Seq2, Row = Row-1, Col = Col)
          aln_result <- append(aln_result, aln_result5)
        }
        if(Col>1){
          Seq1 <- paste('-', temp_seq1, sep = '')
          Seq2 <- paste(seq2[Col-1], temp_seq2, sep = '')
          aln_result6 <- trace_back(Seq1, Seq2, Row = Row, Col = Col-1)
          aln_result <- append(aln_result, aln_result6)
        }
      }
    }
    aln_result
  }
  aln_result <- trace_back()
  structure(
    list(
      Sequences = c(X, Y),
      Scoring_system = Scoring_system,
      Score_matrix = score_matrix,
      Aln_result = aln_result
    ),
    class = 'AlignmentResult'
  )
}
