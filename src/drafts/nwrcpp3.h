/*Rcpp script for aligning two sequences*/

#include <iostream>
#include <string>
#include <algorithm> //reverse
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>

//#include <Rcpp.h>


using namespace std;

int max(int A, int B, int C, char* direction)

void  dpm_init(int** grid, char** traceback, int lengthPattern, int lengthSubject, int gapPenalty)


int nw_align(int** grid, char** traceback, string sequencePattern,
             string sequenceSubject, string& sequencePatternAlignment,
             string& sequenceSubjectAlignment,int gapPenalty)

int nw(string sequencePattern,string sequenceSubject ,string& sequencePatternAlignment, string& sequenceSubjectAlignment)
