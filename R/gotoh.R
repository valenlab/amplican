#' This function uses the Gotoh algorithm to align two nucleotides sequences.
#'
#' \code{gRCPP} This is the main aligning function. The function allocate
#'              memory for six matrices of size N+1 x M+1 , where N and M
#'              are the length of the two sequences that you want to align. The
#'              function return a string representation with verbose
#'              information about the alignment, indels and mismatches,
#'              position of those, length, etc... This information is latter
#'              to be processed in R.
#' @param pattern (string) See subject.
#' @param subject (string) These two strings are the string representation of
#'                          the sequences that we want to align. Passed as
#'                          constant references (don't duplicate and don't
#'                          modify). Both sequences are nucleotide sequences
#'                          and we restrict the alphabet to {a,t,g,c,A,T,G,C}.
#' @param scoringMatrix (string) Select a scoring matrix for the nucleotides.
#'                               The default matrix is NUC44. You can choose
#'                               from the following matrices:
#'                               - NUC44 (+5 for match, -4 for miss)
#' @param gapOpening (int) The score penalty for opening a new gap. Default
#'                         value is 50.
#' @param gapExtension (int) The score penalty for extending the gap. Default
#'                           value is 0.
#' @param gapEnding (bool) Set to TRUE if you want that the alignment score
#'                         be affected by the gaps at the end of the alignment,
#'                         FALSE for otherwise. (FALSE default).
#' @param farIndels (bool) If a gap start at the first positions, or ends at
#'                         the last position, then is a far away indel. In
#'                         some cases, you might want to consider that a
#'                         deletion/insertion, and in some other you don't.
#'                         Set to TRUE if you want that the start or ending
#'                         gap counts as deletion/insertion. FALSE othewise.
#'                         (TRUE default).
#' @return (string) The function return a string with all the information. The
#'                  string is divided in 5 parts. These parts are:
#'                  - Verbose alignment result.
#'                  - List of events from the alignment perspective.
#'                  - List of events from the subject perspective.
#'                  - Pattern alignment string.
#'                  - Pattern subject string.
#'
#'                  You can split each part using
#'                  the substring '++++' (four plus symbols).
#'
#'                  The parts mean the following:
#'
#'                  -- Verbose alignment result --
#'
#'                  In here you will find a verbose and comprehensive summary
#'                  of the alignment. With all the variables logged,
#'                  the comparison between the pattern and the subject, and the
#'                  list of insertions, deletions, and mismatches.
#'                  Example:
#'
#'                   ****************
#'                   SEQUENCES
#'                   ----------------
#'                   Pattern: aaaaa
#'                   Subject: aaddttaaaa
#'                   ****************
#'                   SCORING INFO
#'                   ----------------
#'                   Score matrix:       NUC44
#'                   Gap opening:        50
#'                   Gap extension:      0
#'                   Gap ending penalty: 1
#'                   ****************
#'                   RESULTS
#'                   ----------------
#'                   Alignment Length: 10
#'                   Score: 16
#'
#'                   -----aaaaa    10
#'                        .||||
#'                   aaddttaaaa
#'                   ****************
#'                   ALIGNMENT INFO
#'                   ----------------
#'                   Total Insertions: 0
#'                   Total Deletions:  1
#'                   Total Mismatches: 1
#'                   ----------------
#'                   INSERTIONS:
#'                   ----------------
#'                   --Insertions respect patterns coordinates--
#'                   ----------------
#'                   DELETIONS:
#'                   ----------------
#'                   In Pattern
#'                   START: 1
#'                   END:   5
#'                   --Deletions respect subject coordinates--
#'                   In Pattern
#'                   START: 1
#'                   END:   5
#'                   ----------------
#'                   MISMATCHES:
#'                   ----------------
#'                   START:    5
#'                   ORIGINAL: t
#'                   MUTATED:  a
#'                   --Mismatches respect subject coordinates--
#'                   START:    5
#'                   ORIGINAL: t
#'                   MUTATED:  a
#'                   ****************
#'                   -- List of events from the alignment perspective --
#'                   The next string is compact representation of the
#'                   insertions, deletions, and mismatches. The representation
#'                   follows the following format:
#'                   (at) <total N insertions> (at)
#'                   <insertion 1, start position>,<insertion 1, end position>
#'                   *
#'                   ...
#'                   *
#'                   <insertion N, start position>,<insertion N, end position>
#'                   !
#'                   (at) <total M deletions> (at)
#'                   <deletion 1, start position>,<deletion 1, end position>
#'                   *
#'                   ...
#'                   *
#'                   <deletion M, start position>,<deletion M, end position>
#'                   !
#'                   (at) <total P mismatches> (at)
#'                   <mismatch 1, position>,
#'                   <character of the original nucleotide>,
#'                   <character of the new nucleotide>,
#'                   *
#'                   ...
#'                   *
#'                   <mismatch P, position>,
#'                   <character of the original nucleotide>,
#'                   <character of the new nucleotide>,
#'
#'                   This is an example:
#'                   (at)0(at)!(at)1(at)1,5*!(at)1(at)5,t,a*
#'                   -- List of events from the subject perspective --
#'                   The next string is compact representation of the
#'                   insertions, deletions, and mismatches. The twist here is
#'                   that the information is given from the subject coordinates
#'                   instead of the alignment coordinates.
#'
#'                   For example, if we have the alignment
#'
#'                   12345678901234567801 | Alignment coordinates
#'                   --------------------
#'                   AAAATTTTAAAA----AAAA | Pattern
#'                   AAAA----AAAACCCCAAAA | Subject
#'                   --------------------
#'                   1234----567890123457 | Subject coordinates
#'
#'                   You can see that there is a deletion in the pattern. This
#'                   deletion goes from 13 to 16, both included.
#'
#'                   However, this deletion have different coordinates, if we
#'                   use the subject coordinates. In this case, the coordinates
#'                   tells that there is a deletion in the pattern that goes
#'                   from 9 to 12 both included.
#'
#'                   The format in which this information is portrayed is the
#'                   same as in the previous list of events.
#'                   -- Pattern alignment string --
#'
#'                   A string with the characters that represent the resulted
#'                   alignment for the pattern. In our example:
#'                   -----aaaaa
#'                   -- Pattern subject string --
#'
#'                   A string with the characters that represent the resulted
#'                   alignment for the subject. In our example:
#'
#'                   aaddttaaaa
#' #@examples
#' #gRCPP('aaaaa', 'aaddttaaaa')
#' #gRCPP('aaaaa', 'aaddttaaaa', 'NUC44', 30, 5)
#' #gRCPP('aaaaa', 'aaddttaaaa', 'NUC44', 30, 5, TRUE, TRUE)
#' @useDynLib amplican, .registration = TRUE
#' @importFrom Rcpp sourceCpp
gRCPP <- function(pattern,
                  subject,
                  scoringMatrix = "NUC44",
                  gapOpening = 50,
                  gapExtension = 0,
                  gapEnding = FALSE,
                  farIndels = TRUE) {

    result <- gotohRCPP(pattern,
                        subject,
                        scoringMatrix,
                        gapOpening,
                        gapExtension,
                        gapEnding,
                        farIndels)
    return(result)
}
