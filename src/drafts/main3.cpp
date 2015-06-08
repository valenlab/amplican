/*Rcpp script for aligning two sequences*/

#include <iostream>
#include <string>
#include <algorithm> //reverse
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <ctime>     //random stuff
#include <time.h> 

#include "nwrcpp3.cpp"

using namespace std;


int  main( int argc, char ** argv ){

		srand (time(NULL));

        //Sequences candidate to align
        string  sequencePattern="ACGTGATGACGT";
        string  sequenceSubject="GATGA";

        //Results of the alignment
        string  sequencePatternAlignment;
        string  sequenceSubjectAlignment;

        string results = "";

		cout<<"Sequence: "<<sequencePattern<<endl;
		cout<<"Mutated: "<<sequenceSubject<<endl;
						
		results = nw(sequencePattern, sequenceSubject, sequencePatternAlignment, sequenceSubjectAlignment);

		cout<<results;

        // Print alignment
        //print_al(sequencePatternAlignment, sequenceSubjectAlignment ) ;        

        return  0 ;
}
