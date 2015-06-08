if(FALSE){source("http://bioconductor.org/biocLite.R")
biocLite("DECIPHER")}

### R code from vignette source 'vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.Rnw'

###################################################
### code chunk number 1: ArtOfAlignmentInR.Rnw:52-54
###################################################
options(continue=" ")
options(width=60)


###################################################
### code chunk number 2: expr0
###################################################
library(DECIPHER)

n_points <- 10
N0 <- ceiling(2^seq(5, 13, length.out=n_points))
N1 <- ceiling(2^seq(5, 12, length.out=n_points))
N2 <- ceiling(2^seq(5, 13, length.out=n_points))
N3 <- ceiling(2^seq(5, 16, length.out=n_points))
timings0 <- setNames(rep(0, length(N0)), N0)
timings1 <- setNames(rep(0, length(N1)), N1)
timings2 <- setNames(rep(0, length(N2)), N2)
timings3 <- setNames(rep(0, length(N3)), N3)

print("We are here: 1")
cat(n_points)
print("We are here: 2")
print(n_points)

for (i in seq_len(length(N0))) {
  for (j in 0:3) {
    N <- eval(parse(text=paste("N", j, sep="")))
    
    
    # simulate sequences with 15% distance
    string1 <- DNAStringSet(paste(sample(DNA_ALPHABET[1:4], N[i], replace = TRUE), collapse = ""))
    string2 <- replaceAt(string1,
                         at=IRanges(sample(N[i], ceiling(N[i]/5)), width=1),
                         sample(c(DNA_ALPHABET[1:4], ""), ceiling(N[i]/5), replace = TRUE))
    
    #print("Sequences to align:")
    #print(string1)
    #print(string2)
    
    # align the sequences using two methods
    if (j==0) {
      timings0[i] <- system.time(pairwiseAlignment(string1, string2))[["user.self"]]
    } else if (j==1) {
      timings1[i] <- system.time(AlignProfiles(string1, string2, restrict=-Inf, anchor=NA, processors=1))[["user.self"]]
    } else if (j==2) {
      timings2[i] <- system.time(AlignProfiles(string1, string2, anchor=NA, processors=1))[["user.self"]]
    } else { # j == 3
      timings3[i] <- system.time(AlignProfiles(string1, string2, processors=1))[["user.self"]]
    }
  }
}

c0 <- lm(timings0 ~ N0 + I(N0^2))
c1 <- lm(timings1 ~ N1 + I(N1^2))
c2 <- lm(timings2 ~ N2)
c3 <- lm(timings3 ~ N3)

N <- seq(1, 46340, length.out=1000) # prediction range
plot(N0, timings0,
     xlab = "Sequence length (nucleotides)",
     ylab = "Elapsed Time (sec.)",
     main = "",
     ylim=c(range(timings0, timings1, timings2, timings3)),
     xlim=c(0, max(N3)))
points(N, predict(c0,
                  data.frame(N0 = N)),
       type="l", lty=3)
points(N1, timings1,
       col="blue", pch=0)
points(N, predict(c1,
                  data.frame(N1 = N)),
       type="l", lty=3, col="blue")
points(N2, timings2,
       col="red", pch=5)
points(N, predict(c2,
                  data.frame(N2 = N)),
       type="l", lty=3, col="red")
N <- seq(1, max(N3), length.out=1000) # prediction range
points(N3, timings3,
       col="green", pch=2)
points(N, predict(c3,
                  data.frame(N3 = N)),
       type="l", lty=3, col="green")
legend("bottomright",
       c("pairwiseAlignment",
         "AlignProfiles (unrestricted, unanchored)",
         "AlignProfiles (restricted, unanchored)",
         "AlignProfiles (restricted, anchored)"),
       pch=c(1, 0, 5, 2), col=c("black", "blue", "red", "green"), lty=3)


###################################################
### code chunk number 3: expr1
###################################################
library(DECIPHER)

# OR find the example sequence file used in this tutorial:
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")

dna <- readDNAStringSet(fas)
dna # the unaligned sequences

print (dna)


###################################################
### code chunk number 4: expr2 (eval = FALSE)
###################################################
## AA <- AlignTranslation(dna, asAAStringSet=TRUE) # align the translation
## BrowseSequences(AA, highlight=1) # view the alignment
## 
## DNA <- AlignTranslation(dna) # align the translation then reverse translate
## DNA <- AlignSeqs(dna) # align the sequences directly without translation


###################################################
### code chunk number 5: expr3 (eval = FALSE)
###################################################
## # form guide tree using pairwiseAlignment
## l <- length(dna)
## d <- matrix(0, nrow=l, ncol=l)
## pBar <- txtProgressBar(style=3)
## for (j in 2:l) {
## 	d[j, 1:(j - 1)] <- pairwiseAlignment(rep(dna[j], j - 1),
## 		dna[1:(j - 1)],
## 		scoreOnly=TRUE)
## 	setTxtProgressBar(pBar, j/l)
## }
## close(pBar)
## 
## # rescale the distance scores from 0 to 1
## m <- max(d)
## d[lower.tri(d)] <- d[lower.tri(d)] - m
## m <- min(d)
## d[lower.tri(d)] <- d[lower.tri(d)]/m
## 
## # form a guide tree from the distance matrix
## gT <- IdClusters(d, cutoff=seq(0, 1, 0.01))
## 
## # use the guide tree as input for alignment
## DNA <- AlignSeqs(dna, guideTree=gT) # align directly
## DNA <- AlignTranslation(dna, guideTree=gT) # align by translation


###################################################
### code chunk number 6: expr4 (eval = FALSE)
###################################################
## # form guide tree using inexact clustering
## gT <- IdClusters(myXStringSet=dna, method="inexact", cutoff=seq(0.05, 0.9, 0.05))
## 
## # use the guide tree as input for alignment
## DNA <- AlignSeqs(dna, guideTree=gT) # align directly
## DNA <- AlignTranslation(dna, guideTree=gT) # align by translation


###################################################
### code chunk number 7: sessinfo
###################################################
toLatex(sessionInfo(), locale=FALSE)


