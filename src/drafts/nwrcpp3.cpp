/*Rcpp script for aligning two sequences*/

#include <iostream>
#include <string>
#include <algorithm> //reverse
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>

//#include <Rcpp.h>


using namespace std;

/**
	
	Return the maximum of three integers. The character indicates
	the direction of the traceback matrix propagation, like this:
	 
	-----------
	| B  |  A |
	-----------
	| C  |  X |
	-----------

    @param int A,B,C: integer from where we get the maximum
    @return 
*/
int max(int A, int B, int C, char* direction){
	int  max = -999999 ;

    if(A>=B && A>=C){
		max = A;
        *direction = '|';
	}
    else if(B>C){
		max = B ;
        *direction = '\\' ;
	}
	
    else{
		max = C ;
        *direction = '-' ;
	}

    return  max ;   
}

/**

	Initialize a matrix representing the integesr grid

    @param 
    @return 
*/
void  dpm_init(int** grid, char** traceback, int lengthPattern, int lengthSubject, int gapPenalty){
	
	//The start is set to 0
	grid[0][0] = 0 ;
    traceback[0][0] = 'x' ;

	//Initialize the first row (Pattern)
    for(int i=1; i< lengthPattern; i++){

		grid[0][i] = -i * gapPenalty;
        traceback[0][i] = '-';
	}
	traceback[0][lengthPattern] = '-';

	//Initialize the first column (Subject)
	for(int i=1; i< lengthSubject; i++){

		grid[i][0] =  -i * gapPenalty ;
        traceback[i][0] =  '|' ;
	}
	traceback[lengthSubject][0] =  '|' ;

}


/**
    

    @param 
    @return 
*/
int nw_align(int** grid, char** traceback, string sequencePattern,
             string sequenceSubject, string& sequencePatternAlignment,
             string& sequenceSubjectAlignment,int gapPenalty){
	
	int k = 0;
	int x = 0;
	int y = 0;
	
	int up = 0;
	int diagonal = 0;
	int left = 0;
	
    char direction = 'x';
    char nucleotide = 'x';
    
    int i=0;
    int j=0;

	/* The scoring matrix goes as follow:
	---------------------
	|   | A | T | C | G |
	---------------------
	| A | 2 |-1 |-1 |-1 |
	---------------------
	| T |-1 | 2 |-1 |-1 |
	---------------------
	| C |-1 |-1 | 2 |-1 |
	---------------------
 	| G |-1 |-1 |-1 | 2 |
	*/

    const int  scoreMatrix[4][4] = {{ 2, -1, -1, -1 },    
									{-1,  2, -1, -1 },
									{-1, -1,  2, -1 },
									{-1, -1, -1,  2 }};

    int  lengthPattern = sequencePattern.length();
    int  lengthSubject = sequenceSubject.length();

	//For each of the subject nucleotides
    for(i=1; i<=lengthSubject; i++){

		//For each of the pattern nucleotides
		for(j=1; j <= lengthPattern; j++){

			//Get the nucleotide of the Pattern
			nucleotide = sequencePattern[j-1];

			switch(nucleotide){
				case 'A': x = 0; break;
                case 'T': x = 1; break;
                case 'C': x = 2; break;
                case 'G': x = 3; break;
			}

			//Get the nucleotide of the subject
			nucleotide = sequenceSubject[i-1];

			switch(nucleotide){
				case 'A': y = 0; break;
				case 'T': y = 1; break;
				case 'C': y = 2; break;
				case 'G': y = 3; break;
			}

			//Find the score from the three directions
			up       = grid[i-1][j] - gapPenalty;
            diagonal = grid[i-1][j-1] + scoreMatrix[x][y];
            left     = grid[i][j-1] - gapPenalty;

			//Get the biggest one and write it into the grid
            grid[i][j] = max(up,diagonal,left,&direction);

			//Write the direction into the traceback matrix
            traceback[i][j] = direction;
		}
	}

    //Now, lets go backwards and write the alignment result
    i--;
    j--;

	while(i>0 || j>0){
		switch(traceback[i][j]){
			case '|':
				sequencePatternAlignment += '-'; 
                sequenceSubjectAlignment += sequenceSubject[i-1]; 
                i-- ;
                break;

			case '\\':
				sequencePatternAlignment += sequencePattern[j-1]; 
                sequenceSubjectAlignment += sequenceSubject[i-1]; 
                i--;
                j--;
                break;

			case '-' :
				sequencePatternAlignment += sequencePattern[j-1]; 
                sequenceSubjectAlignment += '-'; 
                j--;
                break;
		}
		k++ ;
	}

	//The final result is the traceback path, but reversed
	reverse(sequencePatternAlignment.begin(), sequencePatternAlignment.end());
	reverse(sequenceSubjectAlignment.begin(), sequenceSubjectAlignment.end());

	return  0 ;
}

/**
    

    @param 
    @return 
*/

// [[Rcpp::export]]
string nw(string sequencePattern,string sequenceSubject ,string& sequencePatternAlignment, string& sequenceSubjectAlignment){
	
	int  gapPenalty = 2 ;

    int lengthPattern = sequencePattern.length();
    int lengthSubject = sequenceSubject.length();

    //Create the grid matrix
    int** grid = new int* [lengthSubject+1];
    for(int i = 0; i <= lengthSubject; i++){
		//cout<<i<<endl;		
		grid[i] = new int[lengthPattern+1];
	}
	
    //Create the traceback matrix
    char** traceback = new char* [lengthSubject+1];
    for(int i = 0; i <= lengthSubject; i++){
		traceback[i] = new char[lengthPattern+1];
	}

    //Initialize traceback and grid
    dpm_init(grid, traceback, lengthPattern, lengthSubject, gapPenalty);

    //Do the alignment
    nw_align(grid, traceback, sequencePattern, sequenceSubject, sequencePatternAlignment, sequenceSubjectAlignment, gapPenalty);
	
	//Delete datastructures
	for(int i=0; i<=lengthSubject; i++){
		//cout<<i<<endl;
		delete[] grid[i];
	}
    delete[] grid;
    
	for(int i=0; i<=lengthSubject; i++){
		delete[] traceback[i];
	}
    delete[] traceback;

    return  sequencePatternAlignment+"@"+sequenceSubjectAlignment ;
}
