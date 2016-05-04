#' This function filter out sequences that has undefined nucleotides.
#'
#' For each forward and reverse read the function looks for unvalid nucleotides.
#' If either the forward or reverse read has a unvalid nucleotide, both of them
#' are considered unvalid reads, hence rendering that row in the dataframe
#' unused.
#' @param Forward.df: (Dataframe) Forward dataframe. This dataframes must have a
#' column named 'Sequence' where the sequence should be.
#' @param Reverse.df: (Dataframe) Reverse dataframe. This dataframes must have a
#' column named 'Sequence' where the sequence should be.
#' @param validRows:  (Bool[])    Array with bools telling which sequences do you
#' want to keep (usually all TRUE). Each one corresponde to a row in the dataframes in the
#' same order. If the value is TRUE it means that we keep it, if the filter pass. If the value is
#' false, it means that we discard it, regardless of the result of the filter.
#' @return (Bool[]) Array with bools, with the non valid rows set to FALSE.
filterByNucleotides <- function(forward.df, reverse.df, validRows) {

  # Get the sequences into two arrays
  candidateForwardSequences  <- forward.df$"Sequence"
  candidateReverseSequences  <- reverse.df$"Sequence"

  # Look for invalid nucleotides, if you get a FALSE, it means that sequence DO NOT HAVE an invalid nucleotide
  unvalidForward <- grepl("[^atcgATCG]", candidateForwardSequences)
  unvalidReverse <- grepl("[^atcgATCG]", candidateReverseSequences)

  # Make the OR operation of both vector, if there is a TRUE one of them then must be deleted (hence the OR)
  result <- unvalidForward + unvalidReverse

  # We keep those which are FALSE, so we need to negate the result
  result <- !result

  # Make an AND operation with the validRow argument.
  result <- result * validRows

  # Return result, but converted to logical values (TRUE, FALSE)
  return(result > 0)
}

#' This filter out sequences which have bad quality readings.
#'
#' For each forward and reverse read, the function looks for sequence qualities
#' and decide if the quality is good enought to be taken into account. If not,
#' then is label as a bad sequence and won't be processed.
#' @param Forward.df: (Dataframe) Forward dataframe. This dataframes must have a
#' column named 'Quality' where the sequence should be.
#' @param Reverse.df: (Dataframe) Reverse dataframe. This dataframes must have a
#' column named 'Quality' where the sequence should be.
#' @param minimum:    (Int)       This is the minimum quality that we accept for
#' every nucleotide. For example, if we have a sequence with nucleotides which have quality
#' 50-50-50-50-10, and we set the minimum to 30, the whole sequence will be a bad sequence.
#' The minimum is set to 0 by default.
#' @param average:    (Int)       This is what the average score of the quality of
#' sequence should be (or greater). For example, if we have a sequence with nucleotides which have
#' quality 70-70-70, the average would be 70. If set the average to 70 or less the sequence will
#' pass. If we set the average to 71 the sequence will not pass. The average is set to 0 by default.
#' @param validRows:  (Bool[])    Array with bools telling which sequences do you
#' want to keep (usually all TRUE). Each one corresponde to a row in the dataframes in the
#' same order. If the value is TRUE it means that we keep it, if the filter pass. If the value is
#' false, it means that we discard it, regardless of the result of the filter.
#' @return (Bool[]) Array with bools, with the non valid rows set to FALSE.
#'
filterByQuality <- function(forward.df, reverse.df, validRows, minimum = 0, average = 0) {

  # Get the quality for each sequence
  candidateForwardQuality  <- forward.df$"Quality"
  candidateReverseQuality  <- reverse.df$"Quality"

  # For each sequence
  for(i in 1:length(candidateForwardQuality)){

    # If is a valid row, find out the qualities
    if(validRows[i]==TRUE){

      # Transform the quality string into an array of char
      qualityArrayForward <- paste( strsplit(as.character(candidateForwardQuality[i]),"")[[1]] , collapse=" ")
      qualityArrayReverse <- paste( strsplit(as.character(candidateReverseQuality[i]),"")[[1]] , collapse=" ")

      # Transform the array of char into an array of integer
      qualityArrayForwardNumeric <- strtoi(charToRaw(qualityArrayForward), base=16L)
      qualityArrayReverseNumeric <- strtoi(charToRaw(qualityArrayReverse), base=16L)

      # Find the average in each
      averageForward <- sum(qualityArrayForwardNumeric)/length(qualityArrayForwardNumeric)
      averageReverse <- sum(qualityArrayReverseNumeric)/length(qualityArrayReverseNumeric)

      # Find which one of those are above the minimum
      aboveMinimumForward <- (qualityArrayForwardNumeric >= minimum)
      aboveMinimumReverse <- (qualityArrayReverseNumeric >= minimum)

      # Find out if we have any that is not above the minimum, if these variable are TRUE, all of them are above minimum
      aboveMinimumForward <- sum(aboveMinimumForward)==length(aboveMinimumForward)
      aboveMinimumReverse <- sum(aboveMinimumReverse)==length(aboveMinimumReverse)

      # Find out if the average is above the threshold
      aboveAverageForward <- (averageForward >= average)
      aboveAverageReverse <- (averageReverse >= average)

      # If both above average and above minimum are TRUE, do nothin. Otherwise, set the row as invalid.
      if(!(aboveMinimumForward && aboveMinimumReverse && aboveAverageForward && aboveAverageReverse)){
        validRows[i] <- FALSE
      }
    }
  }
  return (validRows)
}
