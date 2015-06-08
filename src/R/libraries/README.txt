RCPP doesn't seem to work well with .h files.

What I'm doing now is, the "nw.h" file that is in /src/CPP is copy into /src/tutorial or /src/R (here) and renamed to "nw.cpp".

You also need to change the DEBUG flag to false.

You also need to change the header of nwRCPP and take out the const and reference arguments. The new header is in the code ready to be copypasted.

All of these changes are done in the cpp-developing branch of git

The final mission is to get Rcpp to work with nw.h with the constant headers as in soruceRcpp("/../CPP/nw.h"). But until then, follow this instructions.

IMPORTANT!!!

THE CPP GOES INTO R AND NOT THE OTHER WAY AROUND, IF YOU NEED TO CHANGE SOMETHING DO IT IN THE CPP-DEVELOPING BRANCH AND NOT IN THE R-DEVELOPING BRANCH

TRANSPASERS WILL BE SHOOT ON SIGHT!!
