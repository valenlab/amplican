/**
 * @file
 * @author  Rafael Nozal Cañadas <rafaelc@cbu.uib.no>
 * @version 0.1
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * 	Rcpp script for aligning several sequences with a single subject                
 *	sequence, one by one.
 *
 *	The script allocate memory for one matrix only. After an alignment,
 *	the subject sequence remain in place an a new sequence is given.
 *
 *  The results are given in a string witht the following format:
 *
 *	SUBJECT SEQUENCE@SEQUENCE 1 
 *	Indels:
 *	Start@End@<Pattern,Subject>
 *	Mismatches:
 *	Possition@Subject Base@Pattern Base
 *	SUBJECT SEQUENCE@SEQUENCE 2
 *	...
 *
 *	The class Alignment_Result is also implemented.
 *
 *  Compiled and builded with:
 *  g++ -std=c++11 -Wall -c "%f"
 *  g++ -std=c++11 -Wall -o "%e" "%f"
 */

/** @todo : Spelling! */

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <list> 
#include <algorithm> //reverse() when reverse the traceback result
#include <limits>    // in limits<int>:min()
#include <typeinfo>  //in printMatrix<T>, string or int
#include <ctype.h>   // characters to uppercase

// #include <Rcpp.h>  //If you are not working with R, comment this line
					 //otherwise you have to install Rcpp.h and you need admin
					 //privelege to do that.

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]						 

using namespace std;

//Some constant declarations
const int MINUS_INF = std::numeric_limits<int>::min()/2;
const int DEBUG = false;
string PAUSE = "x";


class Mutation_Info{

	private:
		const unsigned int start;
		const char original;
		const char mutated;

	public:

		/**
		Constructor of the class, set atttributes to s,o and m.

		Use example:
			Mutation_Info example(20,C,G)
			cout<<example.to_string()<<endl;
		*/
		Mutation_Info(const int s, const char o, const char m):start(s),original(o),mutated(m){}

		/**
			This function return an string representation of the object.

			The string returned has the following format:
			START:    <start>
			ORIGINAL: <original>
			MUTATED:  <mutated>
			
			@return string with the specified format.
		*/
		string to_string() const{

			string toReturn="START:    " + std::to_string(start) + "\n"+
							"ORIGINAL: " + original              + "\n"+
							"MUTATED:  " + mutated               + "\n";

			return toReturn;
		};

		/**

			This function return a string representation of the mutationl object
			in a compact format, like this:
			<start>,<original>,<mutated>

		*/
		string to_string_lite() const{

			return std::to_string(start) +","+ original +","+ mutated;

		};

		/**
			Get the 'start' attribute of the class.

			@return const unsigned int with the start
		*/
		inline const unsigned int get_start(){return start;};

		/**
			Get the 'original' attribute of the class.

			@return const char with the original nucleotide
		*/
		inline const unsigned int get_original(){return original;};

		/**
			Get the 'mutated' attribute of the class.

			@return const char with the mutation nucleotide
		*/
		inline const bool get_pattern(){return mutated;};

};


/**
Class that represent Indel information.

Two sequences that are aligned can have several insertion or deletions. In order
to track them we have a list of Indel_Info in each alignment.

The Indel_Info specify where the insertion or deletion start and ends. Also
refers if the indel is in the pattern or the subject. For example:

         0123456789012345678
         -------------------
Pattern  AATATAAC----GCACGTG 
Subject  AA--TAACGGAAGCACGTG 

The Pattern have an insertion from 2 to 3.
The Pattern have a deletion from 8 to 11.

The Subject have the oposite respect the pattern:

The Subject have a deletion from 2 to 3.
The Subject have an insertion from 8 to 11.


@invariant : end must be greater or equal than start
@invariant : the length of the indel would be end-start+1

*/
class Indel_Info{
	
	private:
		const unsigned int start;
		const unsigned int end;
		const bool pattern;
	public:
	
		/**
		Constructor of the class, set atttributes to s,e and p.

		Use example:
			Indel_Info example(20,30,false)
			cout<<example.to_string()<<endl;
		*/
		Indel_Info(const int s, const int e, const bool p):start(s),end(e),pattern(p){}

		/**
			This function return an string representation of the object.

			The string returned has the following format:
			START: <start>
			END:   <end>
			<"In Pattern","In Subject">

			@return string with the specified format.
		*/
		string to_string() const{

			string toReturn="START: "+std::to_string(start)+"\n"+
							"END:   "+std::to_string(end)  +"\n";

			if(pattern==true){ toReturn = "In Pattern\n" + toReturn;}
			else{              toReturn = "In Subject\n" + toReturn;}

			return toReturn;
		};

		/**

			This function return a string representation of the indel object in
			a compact format, like this:
			<start>,<end>

		*/
		string to_string_lite() const{

			return std::to_string(start)+","+ std::to_string(end);

		}

		/**
			Get the 'start' attribute of the class.

			@return const unsigned int with the start
		*/
		inline const unsigned int get_start(){return start;};

		/**
			Get the 'end' attribute of the class.

			@return const unsigned int with the end
		*/
		inline const unsigned int get_end(){return end;};

		/**
			Get the 'pattern' attribute of the class.

			@return const bool with the pattern
		*/
		inline const bool get_pattern(){return pattern;};	

};


/**

This class represent various information about an alignment. The members of the
class are described here:

	@see Indel_Info, in here you will find information about how the list of
					  indels is actually implemented.
					  
	@todo Is there actually a standard way to describe members?

    const string pattern: One of the sequences that are aligned.
    const string subject: The other sequence aligned.

    const string scoreMatrix: Score matrix used to make this alignment.
    const int gapOpening:     Score given for opening gaps.
    const int gapExtension:   Score given for extending gaps.

    const int alignmentLength: How long are both alignments counting gaps. 
    const int score:           Final score of the alignment.

	const string patternAlignment: Alignment of the pattern.
	const string markup:           String that compare visually both alignments.
	const string subjectAlignment: Alignment of the subject.

	const int totalInsertions: How many Insertion are in the alignments.
							   An Insertion is an insertion in the pattern or a
							   deletion in the subject.
	const int totalDeletions:  How many Deletions are in the alignments.
							   A Deletion is a deletion in the pattern or an
							   insertion in the subject.
	const int totalMismatches: How many mismatches are in between the two
							   sequences.

	const std::list<Indel_Info> insertions: List with the insertions information
	const std::list<Indel_Info> deletions:  List with the deletions information.
  

	@invariant: length(patternAlignment) = length(subjectAignment) =
	            length(markup) = alignmentLength

	            length(pattern)<= alignmentLength
	            length(subject)<= alignmentLength

	            size(insertion) = totalInsertions
	            size(deletions) = totalDeletions

	            totalMismatches <= min(length(pattern),length(subject))

*/
/* I'm implementing this class so we can have an independent RCPP functions from
everything else. The idea is to have a to_string method here so we can make
an alignment function align() that return a Alignment_Results. And the RCPP
function will just be alignRCPP(){return align().to_string()} so the core
align function is independent.

Maybe in the future we decide to develop this library a little bit more and make
many RCPP dependencies.


	@todo Better container than just a fixed array? Is O(k), can't think of a
	      better k in this case.

	Originally I thought of making this arrays for insertions and deletions

	const Indel_Info* insertions: 
	const Indel_Info* deletions:

	I realized that we don't know indel in advance and until we finish the grid
	matrix. So we need a mutable container. Afterward we would need to copy the
	list into the array and allocate memory again; hence using the original list
	and save the time for the second allocation.

*/
class Alignment_Result{
	
	private:
		const string pattern;
		const string subject;

		const string scoreMatrix;
		const int gapOpening;
		const int gapExtension;
		const bool gapEnding;

		const int alignmentLength;
		const int score;

		const string patternAlignment;
		const string markup;
		const string subjectAlignment;

		const int totalInsertions;
		const int totalDeletions;
		const int totalMismatches;

		const std::list<Indel_Info> insertions;
		const std::list<Indel_Info> deletions;

		const std::list<Indel_Info> insertionsPattern;
		const std::list<Indel_Info> deletionsSubject;
		
		const std::list<Mutation_Info> mutations;
		const std::list<Mutation_Info> mutationsRelative;

	public:
		/**
			Constructor of the class. Every element is initialize as const
			since the result is fixed.
		*/
		//This is ugly :-(
		Alignment_Result(const string pattern, const string subject,
						 const string scoreMatrix, const int gapOpening,
						 const int gapExtension, const bool gapEnding,
						 const int alignmentLength,
						 const int score, const string patternAlignment,
						 const string markup, const string subjectAlignment,
						 const int totalInsertions,	const int totalDeletions,
						 const int totalMismatches,
						 const std::list<Indel_Info> insertions,
						 const std::list<Indel_Info> deletions,
						 const std::list<Indel_Info> insertionsPattern,
						 const std::list<Indel_Info> deletionsSubject,
						 const std::list<Mutation_Info> mutations,
						 const std::list<Mutation_Info> mutationsRelative
						 ):

						 pattern(pattern),subject(subject),
						 scoreMatrix(scoreMatrix),gapOpening(gapOpening),
						 gapExtension(gapExtension), gapEnding(gapEnding),
						 alignmentLength(alignmentLength),score(score),
						 patternAlignment(patternAlignment), markup(markup),
						 subjectAlignment(subjectAlignment),
						 totalInsertions(totalInsertions),
						 totalDeletions(totalDeletions),
						 totalMismatches(totalMismatches),insertions(insertions)
						 ,deletions(deletions),
						 insertionsPattern(insertionsPattern),
						 deletionsSubject(deletionsSubject),
						 mutations(mutations),
						 mutationsRelative(mutationsRelative){}

		/**
			This function return an string representation of the object.

			The string returned has the following format:

			*********************
			SEQUENCES INFORMATION
			*********************
			SCORING INFORMATION
			*********************
			RESULTS
			*********************
			INDELS AND MUTATIONS
			*********************

			If the instance of the object is initialize to NULL it return a
			warning string message

			@return string with the specified format.
		*/
		string to_string(){

			string toReturn =  "";
			
			if(this==NULL){

				toReturn = "WARNING: Alignment_Result init to NULL";

			}
			else{

				int totalChunks = 0;
				int lastChunkLength = 0;
				string patternPiece = "";
				string markupPiece  = "";
				string subjectPiece = "";
				string indelInfoPiece = "";

				toReturn = toReturn +
				"****************\n"+
				"SEQUENCES \n"+
				"----------------\n"+
				"Pattern: "+pattern + "\n"+
				"Subject: "+subject + "\n"+
				"****************\n"+
				"SCORING INFO \n"+
				"----------------\n"+
				"Score matrix:       "+scoreMatrix+"\n"+
				"Gap opening:        "+std::to_string(gapOpening)+"\n"+
				"Gap extension:      "+std::to_string(gapExtension)+"\n"+
				"Gap ending penalty: "+std::to_string(gapEnding)+"\n"+
				"****************\n"+
				"RESULTS \n"+
				"----------------\n"+
				"Alignment Length: "+std::to_string(alignmentLength)+"\n"+
				"Score: "+std::to_string(score) + "\n\n";

				//Now lets separate the aligns in chunks of lenght 80
				totalChunks = alignmentLength/80;
				lastChunkLength = alignmentLength%80;

				//Print all the chunks of length 80
				for(int i=0; i<totalChunks; i++){
					patternPiece = patternAlignment.substr(80*i,80);
					markupPiece  = markup.substr(80*i,80);
					subjectPiece = subjectAlignment.substr(80*i,80);

					toReturn = toReturn +
					patternPiece + "    " + std::to_string(80*(i+1)) + "\n"+
					markupPiece + "    \n"+
					subjectPiece + "    \n\n";
				}

				//Print the last chunk that didn't reach to 80
				if(lastChunkLength>0){
					patternPiece = patternAlignment.substr(totalChunks*80,lastChunkLength);
					markupPiece  = markup.substr(totalChunks*80,lastChunkLength);
					subjectPiece = subjectAlignment.substr(totalChunks*80,lastChunkLength);

					toReturn = toReturn +
					patternPiece + "    " + std::to_string(alignmentLength) + "\n"+
					markupPiece + "    \n"+
					subjectPiece + "    \n";
				}


				toReturn = toReturn +
				"****************\n"+
				"ALIGNMENT INFO \n"+
				"----------------\n"+
				"Total Insertions: "+std::to_string(totalInsertions)+"\n"+
				"Total Deletions:  "+std::to_string(totalDeletions)  +"\n"+
				"Total Mismatches: "+std::to_string(totalMismatches) +"\n"+
				"----------------\n"+
				"INSERTIONS: \n"+
				"----------------\n";

				//Show the list of all insertions
				for (list<Indel_Info>::const_iterator it = insertions.begin();
				     it!=insertions.end(); it++){

					indelInfoPiece += (*it).to_string() + "\n";					
				}

				toReturn = toReturn + indelInfoPiece + "\n" +

				"--Insertions respect patterns coordinates-- \n\n";

				indelInfoPiece = "";

				//Show the list of all deletions from the subject coordinates.
				for (list<Indel_Info>::const_iterator it = insertionsPattern.begin();
				     it!=insertionsPattern.end(); it++){

					indelInfoPiece += (*it).to_string() + "\n";					
				}

				toReturn += indelInfoPiece + "\n" +

				
				"----------------\n"+
				"DELETIONS: \n"+
				"----------------\n";
				
				indelInfoPiece = "";

				//Show the list of all deletions
				for (list<Indel_Info>::const_iterator it = deletions.begin();
				     it!=deletions.end(); it++){

					indelInfoPiece += (*it).to_string() + "\n";					
				}

				toReturn = toReturn + indelInfoPiece + "\n" +

				"--Deletions respect subject coordinates-- \n\n";

				indelInfoPiece = "";
				
				//Show the list of all deletions from the subject coordinates.
				for (list<Indel_Info>::const_iterator it = deletionsSubject.begin();
				     it!=deletionsSubject.end(); it++){

					indelInfoPiece += (*it).to_string() + "\n";					
				}

				toReturn += indelInfoPiece + "\n" +
			
				"----------------\n"+
				"MISMATCHES: \n"+
				"----------------\n";

				indelInfoPiece = "";

				//Show the list of all missmatches
				for (list<Mutation_Info>::const_iterator it = mutations.begin();
				     it!=mutations.end(); it++){

					indelInfoPiece += (*it).to_string() + "\n";					
				}

				toReturn = toReturn + indelInfoPiece + "\n" +

				"--Deletions respect subject coordinates-- \n\n";

				indelInfoPiece = "";

				//Show the list of all missmatches
				for (list<Mutation_Info>::const_iterator it = mutationsRelative.begin();
				     it!=mutationsRelative.end(); it++){

					indelInfoPiece += (*it).to_string() + "\n";					
				}

				toReturn += indelInfoPiece + "\n" +
				
				"\n"+
				"****************\n";
			}
			
			return toReturn;
		}
		
		/**
			This function return an string representation of the object. In this
			case the information given is the minimum that ggplot2 needs to draw
			the information about the alignments

			The string returned have this format:

			<number of insertions>@<insertion 0 start>,<insertion 0 end>*
			<insertion 1 start>,<insertion 1 end>*...,<insertion N end>*!
			<number of deletions>@<deletion 0 start>,<deletion 0 end>*
			<deletion 1 start>,<deletion 1 end>*...,<deletion M end>*!
			<number of mutations>@<mutation 0 start>,<mutation 0 original>,
			<mutation 0 mutated>*<mutation 1 start>,<mutation 1 original>,
			<mutation 1 mutated>*...,<mutation X start>,<mutation X original>,
			<mutation X mutated>*

			 @see ggplot2
			 @todo info to ggplot2 specifications
			 @see ampliCan.R
			 @todo link to specific lines of code in the R script
			 @todo CIGAR 

		*/
		string to_ggplot2_string(){

			string toReturn = "";

			int totalInsertions = insertions.size();
			int totalDeletions = deletions.size();

			//Add the insertions info
			toReturn += "@" + std::to_string(totalInsertions) + "@";

			for (list<Indel_Info>::const_iterator it = insertions.begin();
				it!=insertions.end(); it++){

				toReturn += (*it).to_string_lite() + "*";
			}

			//Make a character that separate the insertions from deletions
			toReturn += "!";

			//Add the deletions info
			toReturn += "@" + std::to_string(totalDeletions) + "@";

			for (list<Indel_Info>::const_iterator it = deletions.begin();
				it!=deletions.end(); it++){

				toReturn += (*it).to_string_lite() + "*";
			}

			//Make a character that separate the deletions from mutations
			toReturn += "!";

			//Add the mutations info
			toReturn += "@" + std::to_string(totalMismatches) + "@";

			for (list<Mutation_Info>::const_iterator it = mutations.begin();
				it!=mutations.end(); it++){

				toReturn += (*it).to_string_lite() + "*";
			}

			 return toReturn;

		}



		string to_subject_coordiantes_string(){

			string toReturn = "";

			int totalInsertions = insertions.size();
			int totalDeletions = deletions.size();

			//Add the insertions info
			toReturn += "@" + std::to_string(totalInsertions) + "@";

			for (list<Indel_Info>::const_iterator it = insertions.begin();
				it!=insertions.end(); it++){

				toReturn += (*it).to_string_lite() + "*";
			}

			//Make a character that separate the insertions from deletions
			toReturn += "!";

			//Add the deletions info
			toReturn += "@" + std::to_string(totalDeletions) + "@";

			for (list<Indel_Info>::const_iterator it = deletionsSubject.begin();
				it!=deletionsSubject.end(); it++){

				toReturn += (*it).to_string_lite() + "*";
			}

			//Make a character that separate the deletions from mutations
			toReturn += "!";

			//Add the mutations info
			toReturn += "@" + std::to_string(totalMismatches) + "@";

			for (list<Mutation_Info>::const_iterator it = mutationsRelative.begin();
				it!=mutationsRelative.end(); it++){
					
				toReturn += (*it).to_string_lite() + "*";
			}

			 return toReturn;

		}


		/**
			Return only the insertions information in a string format:
			Start@End@<Pattern,Subject>
		*/
		string insertions_to_string(){return "";}

		/**
			Return only the insertions information in a string format:
			Start@End@<Pattern,Subject>
		*/		
		string delitions_to_string(){return "";}

		string get_pattern_alignment(){ return patternAlignment ;}
		string get_subject_alignment(){ return subjectAlignment ;}

};

/**

	Print a 2D matrix with generic variable type T, that has a pattern and a
	subject related to it. Also print the pattern on top and the subject on the
	left side.

	@precond
		The matrix is size N x M
		The pattern must be of length M-1
		The subject must be of length N-1

	@todo Can T be converted to const T? function outside class is not allowed
	      to have printMatrix() const {} declaration.

    @param T** matrix: The matrix we want to print on standard output
    @param string sequencePattern: String with the pattern sequence
	@param string sequenceSubject: String with the subject sequence
	
    @return VOID
*/
template <class T>
void printMatrix(T** matrix, const string &sequencePattern,
                  const string &sequenceSubject){

	//Get the length of each sequence
	int lengthPattern = sequencePattern.length();
    int lengthSubject = sequenceSubject.length();

	//Print the first row which is the pattern on top
    cout << "    ";
    for(int i=0; i<lengthPattern;i++){
		cout << sequencePattern[i] << " ";
	}
    cout << "\n";

	//Print the rest of the lines starting with the subject nucleotide
	for(int i=0;i<=lengthSubject;i++){
		
		//If it is the first character, show the subject nucleotide
		if(i>0){
			cout <<sequenceSubject[i-1]<< " ";
		}
		//Otherwise leave a nice spacing from the nucleotide
		else{
			cout<<"  ";
		}

		//Continue with the rest of the row
		for(int j = 0; j<=lengthPattern; j++){

			//If we are printing numbers
			if (typeid(T) == typeid(int)){

				if(matrix[i][j] <= MINUS_INF){   //If the number is very low
												 //is the same as -infinity
					cout<<"-∞ ";
				}
				else{
					//Positives numbers carry a "+" to compensate aligning
					//with negatives numbers that carry a "-" symbol.
					if(matrix[i][j]>=0){
						cout<<"+";
					}
					cout<<matrix[i][j]<< " ";
				}
			}
			//If we print something else
			else{
				cout<<matrix[i][j]<< " ";
			}
		}
		cout<<"\n";
	}
}

/**

	Return the maximum of every integer that is on a row of a 2D matrix.

*/
inline int maxRow(int** matrix, const int &columnIndex, const int &totalRows,
				  unsigned int &maximumIndex, const int &offset = 0){


	int maximum = matrix[offset][columnIndex];
	maximumIndex = offset;
	
	//For each of the elements in the row
    for(int i=offset+1; i<=totalRows; i++){ //+1, maximum is [i][0], start in [i][1]
		if(maximum < matrix[i][columnIndex]){

			maximum = matrix[i][columnIndex];
			maximumIndex = i;

		}		
	}

	return maximum;

}

/**

	Return the maximum of every integer that is on a column of a 2D matrix.

*/
inline int maxColumn(int** matrix, const int &rowIndex, const int &totalColumns,
				  unsigned int &maximumIndex, const int &offset = 0){


	int maximum = matrix[rowIndex][offset];
	maximumIndex = offset;
	
	//For each of the elements in the row
    for(int i=offset+1; i<=totalColumns; i++){ //+1, maximum is [i][0], start in [i][1]

		if(maximum < matrix[rowIndex][i]){

			maximum = matrix[rowIndex][i];
			maximumIndex = i;

		}		
	}

	return maximum;

}

/**
	
	Return the maximum of three integers. The character indicates
	the direction of the traceback matrix propagation, like this:
	 
	-----------
	| B  |  A |
	-----------
	| C  |  X |
	-----------

    @param int A,B,C: Integer from where we get the maximum, which will
					  max(A,B,C).

	@param direction: A character that represent the direction that we shall
					  follow from X in order to get to the maximum. '-' for C
					  '\' for B and '|' for A. WILL BE MODIFY.
					  
    @return int max: A copy of either, A, B or C.
*/
inline int max(int A, int B, int C, char &direction){
	int  max = -999999 ;

	if(B>=A && B>=C){
		max = B ;
        direction = '\\' ;
	}
    else if(A>=C){
		max = A;
        direction = '|';
	}
    else{
		max = C ;
        direction = '-' ;
	}

    return  max ;   
}

/**
	
	Return the maximum of three integers. The character indicates the direction
	of the traceback matrix propagation. A '\' character means that we continue
	in the traceback grid matrix. A 'V' character means that we should go to the
	vertical traceback matrix. A 'H' chracter means that we should go to the
	horizontal traceback matrix.
	 
    @param int grid,vertical,horizontal: Integers from where we get the maximum,
										 which will max(G,V,H).

	@param direction: A character that represent the direction that we shall
					  follow in order to get to the maximum. 'H' for horizontal
					  '\' for grid and 'V' for vertical. WILL BE MODIFY.
					  
    @return int max: A copy of either, grid, vertical or horizontal.
*/
inline int max3D(int grid, int vertical, int horizontal, char &direction){
	int  max = -999999 ;

    if(grid>=vertical && grid>=horizontal){
		max = grid;
        direction = '\\';
	}
    else if(vertical>=horizontal){
		max = vertical ;
        direction = 'H' ;
	}
    else{
		max = horizontal ;
        direction = 'V' ;
	}

    return  max ;   
}


/**

	Return the maximum of two integer. The maximum represent which path is
	better to take in order to make the alignment.

	If the maximum is 'open' it means that we need to go to the grid traceback
	and check where does it take us. We label direction with a 'G'.

	If the maximum is 'extend' it means that we need to keep going in the
	vertical traceback. We label direction with a '|'

	@param int open, extend: Integer from where we get the maximum. The maximum
							 will be max(open,extend)

	@param direction: 		 A character that represent which direction shall
							 we follow; 'G' for open, and '|' for extend.
							 WILL BE MODIFY.
*/
inline int maxHorizontal(const int &open, const int &extend, char &direction){

	int  max = -999999 ;

    if(open>=extend){
		max = open;
        direction = 'G';
	}
    else{
		max = extend;
        direction = '-';
	}
	
    return  max ;   
}

/**

	Return the maximum of two integer. The maximum represent which path is
	better to take in order to make the alignment.

	If the maximum is 'open' it means that we need to go to the grid traceback
	and check where does it take us. We label direction with a 'G'.

	If the maximum is 'extend' it means that we need to keep going in the
	horizontal traceback. We label direction with a '-'

	@param int open, extend: Integer from where we get the maximum. The maximum
							 will be max(open,extend)

	@param direction: 		 A character that represent which direction shall
							 we follow; 'G' for open, and '-' for extend.
							 WILL BE MODIFY.
*/
inline int maxVertical(const int &open, const int &extend, char &direction){

	int  max = -999999;

    if(open>=extend){
		max = open;
        direction = 'G';
	}
    else{
		max = extend;
        direction = '|';
	}
	
    return  max ;   
}


/**

	Initialize SIX matrices.

	One represent the integer grid where the algorithm propagates and find
	the maximum alignment score.

	It is initialize as:
	grid[0][0] = 0
	grid[i][0] = -infinity (i!=0)
	grid[0][j] = -infinity (j!=0)

	The other two integer matrices represent the gap extension algorithm.
	In this case one resolve horizontally while the other does vertically.

	They are initialize as:
	vertical[0][j] = - (gapOpening + gapExtension*j)
	vertical[i][0] = - infinity (i!=0)

	horizontal[i][0] = - (gapOpening + gapExtension*i)
	horizontal[0][j] = - infinity (j!=0)

	The other ones represent the direction of the propagatation. Encoded
	in chars with an ASCII style.

	They are initialize as:
	tracebackGrid[0][0] = x (representing the finish point)
	tracebackGrid[0][j] = H (representing "go to horizontal matrix") (j!=0)
	tracebackGrid[i][0] = V (representing "go to vertical   matrix") (i!=0)

	tracebackVertical[0][j] = G (representing "go to grid") 
	tracebackVertical[i][0] = | (representing "go up     ") (i!=0)

	tracebackHorizontal[i][0] = G (representing "go to grid")
	tracebackHorizontal[0][j] = - (representing "go to left") (j!=0)

	Here is an example for every matrix initialize with:
	Gap opening = -10
	Gap extension = -2

	GRID:
	    
	-----------------------------------------
	|   |   | A | T | C | G | G | G | A | G |
	-----------------------------------------
	|   | 0 |-in|-in|-in|-in|-in|-in|-in|-in|
	-----------------------------------------
	| A |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------
	| T |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------
	| C |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------
 	| G |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------

	VERTICAL:
	    
	-----------------------------------------
	|   |   | A | T | C | G | G | G | A | G |
	-----------------------------------------
	|   |-10|-12|-14|-16|-18|-20|-22|-24|-26|
	-----------------------------------------
	| A |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------
	| T |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------
	| C |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------
 	| G |-in|   |   |   |   |   |   |   |   |
	-----------------------------------------

	HORIZONTAL:
	    
	-----------------------------------------
	|   |   | A | T | C | G | G | G | A | G |
	-----------------------------------------
	|   |-10|-in|-in|-in|-in|-in|-in|-in|-in|
	-----------------------------------------
	| A |-12|   |   |   |   |   |   |   |   |
	-----------------------------------------
	| T |-14|   |   |   |   |   |   |   |   |
	-----------------------------------------
	| C |-16|   |   |   |   |   |   |   |   |
	-----------------------------------------
 	| G |-18|   |   |   |   |   |   |   |   |
	-----------------------------------------

	TRACEBACK GRID:

	-----------------------------------------
	|   |   | A | T | C | G | G | G | A | G |
	-----------------------------------------
	|   | x | H | H | H | H | H | H | H | H |
	-----------------------------------------
	| A | V |   |   |   |   |   |   |   |   |
	-----------------------------------------
	| T | V |   |   |   |   |   |   |   |   |
	-----------------------------------------
	| C | V |   |   |   |   |   |   |   |   |
	-----------------------------------------
 	| G | V |   |   |   |   |   |   |   |   |
	-----------------------------------------

	TRACEBACK VERTICAL:

	-----------------------------------------
	|   |   | A | T | C | G | G | G | A | G |
	-----------------------------------------
	|   | G | G | G | G | G | G | G | G | G |
	-----------------------------------------
	| A | | |   |   |   |   |   |   |   |   |
	-----------------------------------------
	| T | | |   |   |   |   |   |   |   |   |
	-----------------------------------------
	| C | | |   |   |   |   |   |   |   |   |
	-----------------------------------------
 	| G | | |   |   |   |   |   |   |   |   |
	-----------------------------------------

	TRACEBACK HORIZONTAL:

	-----------------------------------------
	|   |   | A | T | C | G | G | G | A | G |
	-----------------------------------------
	|   | G | - | - | - | - | - | - | - | - |
	-----------------------------------------
	| A | G |   |   |   |   |   |   |   |   |
	-----------------------------------------
	| T | G |   |   |   |   |   |   |   |   |
	-----------------------------------------
	| C | G |   |   |   |   |   |   |   |   |
	-----------------------------------------
 	| G | G |   |   |   |   |   |   |   |   |
	-----------------------------------------	


	NOTE: Since integer are a finite set in C++, we are using the upper or lower
	limit to flag the infinity value. As std::numeric_limits<int>::max(). This
	gives a limit to the maximum score of the matrix of ~32K if we use 16 bit
	representation of signed integers.

	Since we are expecting sequences of maximum 500 nucleotides, and the NUC44
	matrix maximum score is 5; we got a maximum score possible of two equal
	sequence as 5 * 500 = 2500. That gives ample margin to operate and not to
	run into an OVERFLOW situation (which would be around 6500 nucleotides)

	However, for the lower limit std::numeric_limits<int>::min(), we need to
	make an slight ajustment. The minimum will OVERFLOW as soon as we find the
	score of min()-gapOpening. For that reason, the -infinity is going to be
	represented as anything bellow min()/2. We will initialize -infinity to that
	number. Again, with a 16 bit integer representation it will set around the
	number - ~16K. This means that, the worse case scenario is an alignment of
	similarity 0% that has length subject+pattern, ie:
	
		----AAAAAAA
		TTTT-------

	We are expecting sequences of 500 nucleotides.
	
	- Worse case where all nucleotides are align but difference gives us a score
	  of 500 * (-4) = -2500 (again, no problem here).

	- Worse case where all nucleotides are UNALIGN give us a score of:
		-(Opening + 500*Extension)*2 . This equation must never be smaller than
		the minimum set for - infinity. In the case where we set up for -16K:
		2*O + 1000*E > -16K => O + 500*E > -8K. This give us a margin of around
		gap extension ~= 16. This is something worth discussing for minimal
		memory comsumption if 16bits representation of integer is used.
		

	NOTE2: The char** matrices are implemented with chars, only for the user
	convinience. The horizontal and vertical are possible to implement with a
	bool since they only have two states (go <up,left> or go to grid). The grid
	traceback needs THREE states that can be accomplished with only a double
	bool struct. After testing, it would be great to optimize that.
	

    @param int** grid, horizontal, vertical: The three scoring matrices that
											 were PREVIOUSLY CREATED which we
											 are going to initialize properly.

	@param char** tracebackGrid, tracebackHorizontal, tracebackVertical:

											 The three traceback matrices that
											 were PREVIOUSLY CREATED which we
											 are going to initialize properly.

	@param int lengthPattern, lengthSubject:

											 The length of the two sequences.

	@param int gapOpening, gapExtension:

											 The scoring for the gaps.

	@param bool gapEnding:

											 Set to TRUE if you want that the
											 score is affected by the gaps
											 at the end of the alignment. FALSE
											 for otherwise. (FALSE default)

    @return VOID
*/
void  initialize(int** grid, int** vertical, int** horizontal,
				 char** tracebackGrid, char** tracebackHorizontal,
				 char** tracebackVertical, const string &sequencePattern,
				 const string &sequenceSubject, const int &gapOpening,
				 const int &gapExtension, const bool &gapEnding){

	int i = 0; //Reserve memory for index only once

    const int lengthPattern = sequencePattern.length();
    const int lengthSubject = sequenceSubject.length();
	
	//Initialize grid
	grid[0][0] = 0 ;

	//We have two diferent inicializations, depending if we count the ending
	//of the alignments as 0 or whatever gap penalty.
	for(i=1; i<=lengthPattern; i++){
		//grid[0][i] = -(gapOpening + gapExtension * i);
		grid[0][i] = 0;
	}
	for(i=1; i<=lengthSubject; i++){
		//grid[i][0] = -(gapOpening + gapExtension * i);
		grid[i][0] = 0;
	}
	
	if(DEBUG==true){
		cout<<"Grid initialized"<<endl;
	}

	//Initialize horizontal
	for(i=1; i<=lengthPattern; i++){
		horizontal[0][i] = -(gapOpening + gapExtension * i);
	}
    for(i=0; i<=lengthSubject; i++){
		horizontal[i][0] = MINUS_INF;
	}

	if(DEBUG==true){
		cout<<"Horizontal initialized"<<endl;
	}
	
	//Initialize vertical
    for(i=1; i<=lengthSubject; i++){
		vertical[i][0] = -(gapOpening + gapExtension * i);		
	}
	for(i=0; i<=lengthPattern; i++){
		vertical[0][i] = MINUS_INF;		
	}

	if(DEBUG==true){
		cout<<"Vertical initialized"<<endl;
	}
	
	//Initialize tracebackGrid
	tracebackGrid[0][0] = 'x';

    for(i=1; i<=lengthPattern; i++){
		tracebackGrid[0][i] = '-';
	}
	for(i=1; i<=lengthSubject; i++){
		tracebackGrid[i][0] = '|';
	}

	//Initialize tracebackVertical
	for(i=1; i<=lengthSubject; i++){
		tracebackVertical[i][0] = '|';
	}
    for(i=0; i<=lengthPattern; i++){
		tracebackVertical[0][i] = 'G';
	}

	//Initialize tracebackHorizontal
    for(i=1; i<=lengthPattern; i++){
		tracebackHorizontal[0][i] = '-';
	}
	for(i=0; i<=lengthSubject; i++){
		tracebackHorizontal[i][0] = 'G';
	}

	if(DEBUG==true){

		cout<<"DEBUG in initialize"<<endl;
		cout<<"Ending of initialize "<<endl;
		cout<<"Grid matrix"<<endl;
		printMatrix(grid, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Vertical matrix"<<endl;
		printMatrix(vertical, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Horizontal matrix"<<endl; 
		printMatrix(horizontal, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Traceback Grid matrix"<<endl;
		printMatrix(tracebackGrid, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Traceback Vertical matrix"<<endl;
		printMatrix(tracebackVertical, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Traceback Horizontal matrix"<<endl;
		printMatrix(tracebackHorizontal, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"END DEBUG in initialize"<<endl<<endl;
	}

}


/**
    This function allow you to align two sequences, provided that the grid
    matrices and the traceback matrices are inicialized.

    @param 
    @return 
*/
int align(int** grid, int** vertical, int** horizontal, char** tracebackGrid,
		  char** tracebackVertical, char** tracebackHorizontal,
		  const string &sequencePattern, const string &sequenceSubject,
		  const string &matrix, const int &gapOpening,
		  const int &gapExtension, const bool &gapEnding,
		  const bool &farIndels, Alignment_Result* &myResult){

	if(DEBUG==true){

		cout<<"Function Align called with parameters:"<<endl;
		cout<<"G Opening: "<<gapOpening<<endl;
		cout<<"G Extensi: "<<gapExtension<<endl;
		cout<<"G Ending : "<<gapEnding<<endl;

	}

	//Find out coordinates in scoring matrix
	int x = 0;
	int y = 0;

	//Variables for finding max scores
	int up = 0;
	int diagonal = 0;
	int left = 0;
	int newGap = 0;
	int extendGap = 0;

	//Auxiliary variables to store nucleotides and direction chars
    char direction = 'x';
    char nucleotide = 'x';

    //Indexes for the loops
    unsigned int i=0;
    unsigned int j=0;

    //Scoring variables
    unsigned int maximumIndex=0;
    int maxScore = 0;
    
    //The result of the alignment is store in here
    string patternResult = "";
    string subjectResult = "";
    string markup = "";

	//The list of insertions and deletions
    std::list<Indel_Info> insertions;
    std::list<Indel_Info> deletions;
    std::list<Indel_Info> insertionsPattern;
    std::list<Indel_Info> deletionsSubject;
    int indelEnd = 0;
    int indelStart = 0;

	int indelSubjectStart = 0;
	int indelSubjectEnd = 0;
	int indelPatternStart = 0;
	int indelPatternEnd = 0;
	int subjectIndex = 0;
	int patternIndex = 0;

	//The list of mutations
	std::list<Mutation_Info> mutations;
	std::list<Mutation_Info> mutationsRelative;
	//int  mutationStart = 0;
	//char mutationOriginal = 'Z';
	//char mutationMutated  = 'Y';

    //Auxiliary stuff
    bool gapFound = 0;
    char currentNucleotide  = 'x';

	//Scoring matrix structure (init to NUC44 by default)
	int  scoreMatrix[4][4] = {{ 5, -4, -4, -4 },    
							  {-4,  5, -4, -4 },
							  {-4, -4,  5, -4 },
							  {-4, -4, -4,  5}};

	//Get the length of the sequences
    unsigned int  lengthPattern = sequencePattern.length();
    unsigned int  lengthSubject = sequenceSubject.length();

	//Set up the scoring matrix
	/**@todo set up other matrices*/
	if( matrix!="NUC44"){
		cerr<<"Oh no! we don't have more matrices :-("<<endl;
	}

	//For each of the subject nucleotides
    for(i=1; i<=lengthSubject; i++){

		//For each of the pattern nucleotides
		for(j=1; j <= lengthPattern; j++){

			//Get the nucleotide of the Pattern
			nucleotide = sequencePattern[j-1];

			switch(nucleotide){
				case 'A': x = 0; break;
				case 'a': x = 0; break;
                case 'T': x = 1; break;
                case 't': x = 1; break;
                case 'C': x = 2; break;
                case 'c': x = 2; break;
                case 'G': x = 3; break;
                case 'g': x = 3; break;
			}

			//Get the nucleotide of the subject
			nucleotide = sequenceSubject[i-1];

			switch(nucleotide){
				case 'A': y = 0; break;
				case 'a': y = 0; break;
				case 'T': y = 1; break;
				case 't': y = 1; break;
				case 'C': y = 2; break;
				case 'c': y = 2; break;
				case 'G': y = 3; break;
				case 'g': y = 3; break;
			}

			//Find the scores in the vertical and horizontal matrices for i,j
			//Finding for the horizontal
			newGap = grid[i][j-1] - (gapOpening + gapExtension);
			extendGap = horizontal[i][j-1] - gapExtension;
			
			horizontal[i][j] = maxHorizontal(newGap,extendGap,direction);
			tracebackHorizontal[i][j] = direction;

			//Finding for the vertical
			newGap = grid[i-1][j] - (gapOpening + gapExtension);
			extendGap = vertical[i-1][j] - gapExtension;
			
			vertical[i][j] = maxVertical(newGap,extendGap,direction);
			tracebackVertical[i][j] = direction;
			
			//Find the score from the three directions
			up       = horizontal[i-1][j-1];
			diagonal = grid[i-1][j-1];
			left     = vertical[i-1][j-1];

			//Get the biggest one and write it into the grid
            grid[i][j] = max3D(diagonal,up,left,direction) + scoreMatrix[x][y];
            
			//Write the direction into the traceback matrix
            tracebackGrid[i][j] = direction;
		}
	}

    //Now, lets go backwards and write the alignment result

    /*
		In here we can actually find out the insertions and deletions.
		Every time we come back from the vertical we have a deletion (pattern
		respect	the subject), and everytime we combe back from the horizontal
		we have an insertion (pattern respect the subject).
    */
    i--;
    j--;
	direction = 'x';
	maxScore = grid[i][j];

	/*
		For the gap ending penalty:

		EXPERIMENTAL

		If we DO HAVE ending penalty we don't need to touch anything.

		Otherwise, what we do is to find the maximum score in the last column
		xor row, depending of if the bigger sequence is the pattern or the
		subject. Once we find the maximum we adjust the i xor j index so we
		start tracebacking from there. Also adjust the result string accordingly
	*/
	if(gapEnding == false){

		if(DEBUG == true){

			cout<<"Enter gap ending" << endl;
			cout<<"I: "<< i << endl;
			cout<<"J: "<< j << endl;

		}

		//Find out the maximum score in the last row and the last column
		//The maximum score is right now set to the right-bottom corner

		//Get the maximum of the right column
		int maximumColumn = maxRow(grid,j,i,maximumIndex,0);
		unsigned int rowIndex = maximumIndex;
		
		//Get the maximum of the bottom row
		int maximumRow = maxColumn(grid,i,j,maximumIndex,0);
		unsigned int columnIndex = maximumIndex;

		//One of those three is the maximum.
		//If the max Score is the bigger, do nothing, otherwise
		if(maxScore < maximumColumn || maxScore < maximumRow){

			if(maximumColumn>=maximumRow){

				//Adjust the alignments and the traceback start
				while(i!=rowIndex){
					
					patternResult += '-'; 
					subjectResult += sequenceSubject[i-1];
					markup += ' ';  
					i--;
					
				}

			}
			else{

				//Adjust the alignments and the traceback start
				while(j!=columnIndex){
					
					patternResult += sequencePattern[j-1]; 
					subjectResult += '-';
					markup += ' ';  
					j--;
					
				}

			}
		}


		if(DEBUG == true){

			cout<<"Final START after gap ending" << endl;
			cout<<"I: "<< i << endl;
			cout<<"J: "<< j << endl;

		}

		


	
	}

	if(DEBUG==true){

		cout<<"DEBUG Align"<<endl;
		cout<<"Start traceback in: "<<maximumIndex<<endl;
		cout<<"Maximum score found: "<<maxScore<<endl;
		cout<<"Grid matrix"<<endl;
		printMatrix(grid, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Horizontal matrix"<<endl; 
		printMatrix(horizontal, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Vertical matrix"<<endl;
		printMatrix(vertical, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Traceback Grid matrix"<<endl;
		printMatrix(tracebackGrid, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Traceback Horizontal matrix"<<endl; 
		printMatrix(tracebackHorizontal, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"Traceback Vertical matrix"<<endl;
		printMatrix(tracebackVertical, sequencePattern, sequenceSubject);
		cout<<endl;
		cout<<"I: "<<i<<endl;
		cout<<"J: "<<j<<endl;
		cout<<"END DEBUG Align"<<endl<<endl;
		cout<<"Making the traceback"<<endl;
	}

	while(i>0 || j>0){

		if(DEBUG == true && true){

			cout<<"Switching..."<<endl;
			cout<<"Candidate: "<<tracebackGrid[i][j]<<endl;
			cout<<"I: "<<i<<endl;
			cout<<"J: "<<j<<endl;
			cout<<"Current Pattern"<<endl;
			cout<<patternResult<<endl;
			cout<<markup<<endl;
			cout<<subjectResult<<endl;

		}
		
		switch(tracebackGrid[i][j]){

			//The easiest case, just go on diagonal.
			case '\\':
				patternResult += sequencePattern[j-1]; 
                subjectResult += sequenceSubject[i-1];
				if(toupper(sequencePattern[j-1]) == toupper(sequenceSubject[i-1])){
					markup += '|';
				}
				else{
					markup += '.';
				}

				if(DEBUG == true && true){

					cout<<"Case Diagonal found"<<endl;
					cout<<patternResult<<endl;
					cout<<markup<<endl;
					cout<<subjectResult<<endl;
					cin.ignore().get();

				}
				
                i--;
                j--;
                break;

			//If the traceback tell us to go vertical
			case 'V':

				/*
				 If the traceback tell us to go vertical it means that the
				 score came back from the DIAGONAL of the horizontal matrix

				 So we need to mark one last pair, then follow the vertical
				 path until it tell us to come back to grid.
				*/
								
				//First, lets mark the alignment with the insertion
				patternResult += sequencePattern[j-1]; 
                subjectResult += sequenceSubject[i-1];
                if(toupper(sequencePattern[j-1]) == toupper(sequenceSubject[i-1])){
					markup += '|';
				}
				else{
					markup += '.';
				}

				if(DEBUG == true && false){

					cout<<"Case Horizontal found"<<endl;
					cout<<patternResult<<endl;
					cout<<markup<<endl;
					cout<<subjectResult<<endl;

				}

				
                //indelEnd = i;
                i--;
                j--;




                /*
                
				Now we need to continue UP in the VERTICAL TRACEBACK MATRIX
				until we find the reference that transport us back to the
				TRACEBACK GRID.

				*/
                direction = tracebackVertical[i][j];
                
				while(i>0 && direction!='G'){
					subjectResult += sequenceSubject[i-1]; 
					patternResult += '-';
					markup += ' ';  
					i--;
					direction = tracebackVertical[i][j];
				}

				/*
					Since we are going back to grid, it means that the cell to
					where we are going back, was one cell TO THE LEFT from the
					current G marker. We need to adjust and go one more step
					to the LEFT
				*/

				subjectResult += sequenceSubject[i-1]; 
				patternResult += '-';
				markup += ' ';  
				i--;

				break;
				
			//If the traceback tell us to go horizontal
			case 'H':

				/*
				 If the traceback tell us to go horizontal it means that the
				 score came back from the DIAGONAL of the horizontal matrix

				 So we need to mark one last pair, then follow the horizontal
				 path until it tell us to come back to grid.
				*/
				
				//First, lets mark the alignment with the insertion
				patternResult += sequencePattern[j-1]; 
                subjectResult += sequenceSubject[i-1];
                if(toupper(sequencePattern[j-1]) == toupper(sequenceSubject[i-1])){
					markup += '|';
				}
				else{
					markup += '.';
				}

				if(DEBUG == true && false){

					cout<<"Case Horizontal found"<<endl;
					cout<<patternResult<<endl;
					cout<<markup<<endl;
					cout<<subjectResult<<endl;

				}

				
                //indelEnd = j;
                i--;
                j--;

                /*
                
				Now we need to continue LEFT in the HORIZONTAL TRACEBACK MATRIX
				until we find the reference that transport us back to the
				TRACEBACK GRID.

				*/
                direction = tracebackHorizontal[i][j];
                
				while(j>0 && direction!='G'){
					patternResult += sequencePattern[j-1]; 
					subjectResult += '-';
					markup += ' ';  
					j--;
					direction = tracebackHorizontal[i][j];

					if(DEBUG == true && false){

						cout<<"Case Horizontal follow"<<endl;
						cout<<patternResult<<endl;
						cout<<markup<<endl;
						cout<<subjectResult<<endl;

					}
					
				}

				/*
					Since we are going back to grid, it means that the cell to
					where we are going back, was one cell TO THE LEFT from the
					current G marker. We need to adjust and go one more step
					to the LEFT
				*/

				patternResult += sequencePattern[j-1]; 
				subjectResult += '-';
				markup += ' ';  
				j--;

				break;


			//If we reach the top of the traceback grid (NOT THE HORIZONTAL!)
			case '-':

				while(j>0){
					patternResult += sequencePattern[j-1]; 
					subjectResult += '-';
					markup += ' ';  
					j--;
				}
			
				break;

			//If we reach the left of the traceback grid (NOT THE VERTICAL!)
			case '|':

				while(i>0){
					subjectResult += sequenceSubject[i-1]; 
					patternResult += '-';
					markup += ' ';  
					i--;
				}

				break;				

		}
		
	}

	//The final result is the traceback path, but reversed
	reverse(patternResult.begin(), patternResult.end());
	reverse(markup.begin(), markup.end());
	reverse(subjectResult.begin(), subjectResult.end());

	//Find the insertions and deletions once the traceback is finish.
	/*
		I was trying to find the indels as the tracebacks. You can find
		the length of the indel. And you can find the position from the
		rear end. But since you don't know the length of the alignment
		until the traceback is finish, you can't do 1-n to find the inverse
		offset as you go.
	*/

	//Find the insertions (look in subject)
	gapFound = false;
	indelStart = 0;
	indelEnd = 0;
	indelPatternStart = 0;
	indelPatternEnd = 0;
	patternIndex = 0;
	for(i=0; i<subjectResult.size(); i++){

		currentNucleotide = subjectResult[i];

		if(patternResult[i] != '-'){

			patternIndex = patternIndex + 1;
			
		}

		//We found the beggining of a new gap
		if(gapFound == false && currentNucleotide == '-'){

			indelStart = i+1;
			gapFound = true;
			indelPatternStart = patternIndex;


		}
		else{
			//We found the end of the current gap
			if(gapFound == true && currentNucleotide != '-'){

				gapFound = false;

				//Check that we didn't found the first gap in case of an alignment
				//starting with a gap (ie: ----AA CCCCAA, no indel in there)
				if(indelStart != 1){

					indelEnd = i;
					indelPatternEnd = patternIndex-1;

					Indel_Info newInsertion(indelStart,indelEnd,false);
					Indel_Info newInsertionPattern(indelPatternStart,
												   indelPatternEnd,false);
												   
					insertions.push_back(newInsertion);
					insertionsPattern.push_back(newInsertionPattern);

					
				}
			}
		}
	}
	
	//If we count the far right/left indels,we need to close the right insertion
	//if there is one and the pattern alignment ends in gaps.
	if(farIndels==true && gapFound==true){
		indelEnd = subjectResult.size();
		indelPatternEnd = patternIndex;

		Indel_Info newInsertion(indelStart,indelEnd,true);
		Indel_Info newInsertionPattern(indelPatternStart,indelPatternEnd,true);
		
		insertions.push_back(newInsertion);
		insertionsPattern.push_back(newInsertionPattern);
	}	

	//Find the deletions (look in pattern)
	gapFound = false;
	indelStart = 0;
	indelEnd = 0;
	indelSubjectStart = 0;
	indelSubjectEnd = 0;
	subjectIndex = 0;
	for(i=0; i<patternResult.size(); i++){

		currentNucleotide = patternResult[i];

		/*

		 For the subject reference
		 Forward the index of the subject if it is not a gap
		 For example:
			 01234567890123456789012345678
		     aaaaaacccctt-------ttcccaaaaa
		     012345678901       2345678901

			 The first 'T' has index 10 in both the alignment and the subject.

			 The third 'T' has index 12 in the subject, but 19 in the alignment.

		  So we only count +1 if we don't have a gap.
			
		*/
		if(subjectResult[i] != '-'){

			subjectIndex = subjectIndex + 1;
			
		}

		//We found the beggining of a new gap
		if(gapFound == false && currentNucleotide == '-'){

			indelStart = i+1;
			gapFound = true;
			indelSubjectStart = subjectIndex;
		}
		else{
			//We found the end of the current gap
			if(gapFound == true && currentNucleotide != '-'){

				gapFound = false;

				/*

					If we want the start and end of the alignment to count as
					an indel, we need to activate 'farIndels' in this function.

					If farIndels == FALSE

					Check that we didn't found the first gap in case of an
					alignment starting with a gap

					(ie: ----AA CCCCAA, no indel in there)

					If farIndels == TRUE

					(ie: ----AA CCCCAA, valid indel)
					
				*/
				if(indelStart != 1 || farIndels==true){

					indelEnd = i;
					indelSubjectEnd = subjectIndex-1;

					Indel_Info newDeletion(indelStart,indelEnd,true);
					Indel_Info newDeletionSubject(indelSubjectStart,
												  indelSubjectEnd,true); 
					deletions.push_back(newDeletion);
					deletionsSubject.push_back(newDeletionSubject);
				}
			}
		}
	}
	//If we count the far right/left indels, we need to close the right deletion
	//if there is one and the pattern alignment ends in gaps.
	if(farIndels==true && gapFound==true){
		indelEnd = patternResult.size();
		indelSubjectEnd = subjectIndex;

		Indel_Info newDeletion(indelStart,indelEnd,true);
		Indel_Info newDeletionSubject(indelSubjectStart,indelSubjectEnd,true);
		
		deletions.push_back(newDeletion);
		deletionsSubject.push_back(newDeletionSubject);
	}

	//Find the missmatches (look in both)
	//We can optimize this by looking in the insertion / deletions; but for now
	//we are prioritazing readability
	subjectIndex = 0;

	for(i=0; i<patternResult.size(); i++){

		//currentNucleotide  = patternResult[i];


		//First lets check out that we are not in a gap so we can uppercase
		// By definition, there is no mutation in a gap
		if(patternResult[i] != '-' && subjectResult[i] !='-'){

			if(toupper(patternResult[i]) != toupper(subjectResult[i])){

				Mutation_Info newMutation(i,subjectResult[i],patternResult[i]);
				mutations.push_back(newMutation);

				Mutation_Info newMutationRelative(subjectIndex,subjectResult[i],patternResult[i]);
				mutationsRelative.push_back(newMutationRelative);
				
			}

		}

		/*

		 For the subject reference
		 Forward the index of the subject if it is not a gap
		 For example:
			 01234567890123456789012345678
		     aaaaaacccctt-------ttcccaaaaa
		     012345678901       2345678901

			 The first 'T' has index 10 in both the alignment and the subject.

			 The third 'T' has index 12 in the subject, but 19 in the alignment.

		  So we only count +1 if we don't have a gap.
			
		*/
		if(subjectResult[i] != '-'){

			subjectIndex = subjectIndex + 1;
			
		}

	}


	

	//Construct a new alignment object
	myResult = new Alignment_Result(sequencePattern, sequenceSubject, matrix,
									gapOpening, gapExtension, gapEnding,
									patternResult.length(), maxScore,
									patternResult, markup, subjectResult,
									insertions.size(), deletions.size(),
									mutations.size(), insertions, deletions,
									insertionsPattern, deletionsSubject,
									mutations, mutationsRelative);
	

	if(DEBUG==true){
		cout<<"Final result"<<endl;
		cout<<patternResult<<endl;
		cout<<markup<<endl;
		cout<<subjectResult<<endl;
		cout<<"Alignment_Result status"<<endl;
		cout<<myResult->to_string()<<endl;
	}

	return  0 ;
}

/**
	This is the main aligning fucntion. The function allocate memory for six
    matrices of size N+1 x M+1 , where N and M are the length of the two
    sequences that you want to align. The function return a string
    representation with verbose information about the alignmnent, indels and
    missmaches, possition of those, length, etc... This information is latter
    to be processed in R.

	@see documentation.pdf: Accompany this code there is a verbose explanation
							on how this algortihm works. The function runs in
							time O(N²) instead of the classic O(N³) where you
							check for the whole line and column for alignment
							extension score.

    @see nwManyRCPP: This other function align X sequences one by one with a
					 fixed subject sequence. Instead of allocating X times NxM
					 matrices it only allocate one of max(N)xM and overwrite it
					 as the patter sequence is swapped.



    @param string sequencePattern, sequence Subject:

			These two strings are the string representation of the sequences
			that we want to align. Passed as constant references (don't
			duplicate and don't modify).
    

	@param string matrix: Select a scoring matrix for the nucleotides. The
						  default matrix is NUC44. You can choose from the
						  following matrices:
						  - NUC44 (+5 for match, -4 for miss)
						  -
						  -

	@param int gapOpening, gapExtension: The score penalty for the gaps, the
										 default values are 50 for opening
										 and 0 for extension.

	@param bool gapEnding: Set to TRUE if you want that the score is affected by
						   the gaps at the end of the alignment. FALSE for
						   otherwise. (FALSE default)										 

	@param bool debug: Set to true if you want to print verbose information
					   about the matrices in the standard output.

    @return int with 0 if everything goes according to plan.
    
*/
int gotoh(const string &sequencePattern,const string &sequenceSubject,
			  const string &matrix, const int &gapOpening,
			  const int &gapExtension, const bool &gapEnding,
			  const bool &farIndels, Alignment_Result* &myResult){
	
    int lengthPattern = sequencePattern.length();
    int lengthSubject = sequenceSubject.length();

    //Create the grid, vertical and horizontal matrix
    int** grid       = new int* [lengthSubject+1];
    int** vertical   = new int* [lengthSubject+1];
    int** horizontal = new int* [lengthSubject+1];
    for(int i = 0; i <= lengthSubject; i++){
		grid[i]       = new int[lengthPattern+1];
		vertical[i]   = new int[lengthPattern+1];
		horizontal[i] = new int[lengthPattern+1];
	}

	//Create the traceback matrices
    char** tracebackGrid       = new char* [lengthSubject+1];
    char** tracebackVertical   = new char* [lengthSubject+1];
    char** tracebackHorizontal = new char* [lengthSubject+1];
    for(int i = 0; i <= lengthSubject; i++){
		tracebackGrid[i]       = new char[lengthPattern+1];
		tracebackVertical[i]   = new char[lengthPattern+1];
		tracebackHorizontal[i] = new char[lengthPattern+1];
	}

	if(DEBUG==true){
		cout<<"Going to initialize"<<endl;

	}

    //Initialize tracebacks and score matrices
	initialize(grid, vertical, horizontal, tracebackGrid, tracebackHorizontal,
			   tracebackVertical, sequencePattern, sequenceSubject, gapOpening,
			   gapExtension, gapEnding);


	//Show the matrices if we are debugging
	if(DEBUG == true){
		
		cout<<"NW DEBUG"<<endl;
		cout<<"Grid matrix"<<endl;
		printMatrix(grid, sequencePattern, sequenceSubject);
		cout<<"Vertical matrix"<<endl;
		printMatrix(vertical, sequencePattern, sequenceSubject);
		cout<<"Horizontal matrix"<<endl;
		printMatrix(horizontal, sequencePattern, sequenceSubject);

		cout<<"Traceback Grid matrix"<<endl;
		printMatrix(tracebackGrid, sequencePattern, sequenceSubject);
		cout<<"Traceback Vertical matrix"<<endl;
		printMatrix(tracebackVertical, sequencePattern, sequenceSubject);
		cout<<"Traceback Horizontal matrix"<<endl;
		printMatrix(tracebackHorizontal, sequencePattern, sequenceSubject);
		cout<<"END NW DEBUG"<<endl<<endl;

	}
    
    //Do the alignment
	align(grid, vertical, horizontal, tracebackGrid, tracebackVertical,
	      tracebackHorizontal, sequencePattern, sequenceSubject, matrix,
		  gapOpening, gapExtension, gapEnding, farIndels, myResult);

    //Delete all the things!
    for(int i = 0; i <= lengthSubject; i++){
		delete[] grid[i];
		delete[] vertical[i];
		delete[] horizontal[i];
		delete[] tracebackGrid[i];
		delete[] tracebackVertical[i];
		delete[] tracebackHorizontal[i];
	}
	delete[] grid;
	delete[] vertical;
	delete[] horizontal;
	delete[] tracebackGrid;
	delete[] tracebackVertical;
	delete[] tracebackHorizontal;	
	
    return 0;
}


/**
	This is the main aligning fucntion. The function allocate memory for six
    matrices of size N+1 x M+1 , where N and M are the length of the two
    sequences that you want to align. The function return a string
    representation with verbose information about the alignmnent, indels and
    missmaches, possition of those, length, etc... This information is latter
    to be processed in R.

	@see documentation.pdf: Accompany this code there is a verbose explanation
							on how this algortihm works. The function runs in
							time O(N²) instead of the classic O(N³) where you
							check for the whole line and column for alignment
							extension score.

    @see nwManyRCPP: This other function align X sequences one by one with a
					 fixed subject sequence. Instead of allocating X times NxM
					 matrices it only allocate one of max(N)xM and overwrite it
					 as the patter sequence is swapped.

	@see nw:		 Same function, but the information is returned in an a
					 alignment_Result class object instead of a String.

    @param string sequencePattern, sequence Subject:

			These two strings are the string representation of the sequences
			that we want to align. Passed as constant references (don't
			duplicate and don't modify).
    
    @param string sequencePatternAlignment, sequenceSubjectAlignment:

			(THIS WILL GO INTO THE ALIGNMENT CLASS, NOT NECESSARY)
	
			In here, the alignment of the string representation of two sequences
			functions will be stored.

			For example if we give the sequences:

			ATTAGTGAGGTT
			CAGTGCTT

			These two variables will be modify into

			ATTAGTGAGG--TT
			--CAGTGC----TT

	@param string matrix: Select a scoring matrix for the nucleotides. The
						  default matrix is NUC44. You can choose from the
						  following matrices:
						  - NUC44 (+5 for match, -4 for miss)
						  -
						  -

	@param int gapOpening, gapExtension: The score penalty for the gaps, the
										 default values are 50 for opening
										 and 0 for extension.

	@param bool gapEnding: Set to TRUE if you want that the score is affected by
						   the gaps at the end of the alignment. FALSE for
						   otherwise. (FALSE default)										 

	@param bool debug: Set to true if you want to print verbose information
					   about the matrices in the standard output.

    @return string with the following format:

		SUBJECT SEQUENCE ALIGNMENT@PATTERN SEQUENCE ALIGNMENT
		Indels:
		Start@End@<Pattern,Subject>
		Missmatches:
		Possition@Subject Base@Pattern Base

	@todo : Fix the const and default in Rcpp (see bellow)
    
*/
/*

For whatever reason CPP don't like the const and &  declaration, neither the
default value call. While fixing that, use this header instead for testing

string nwRCPP(string sequencePattern, string sequenceSubject,
			  string matrix, int gapOpening, int gapExtension,
			  bool gapEnding, bool farIndels){


*/


// [[Rcpp::export]]
string nwRCPP(const string &sequencePattern,const string &sequenceSubject,
			  const string &matrix="NUC44", const int &gapOpening=50,
			  const int &gapExtension=0, const bool &gapEnding=true,
			  const bool &farIndels=true){

	int status = 0;
	Alignment_Result* myResult = NULL;
	string result = "";

	//Call the interface function
	status = gotoh(sequencePattern, sequenceSubject, matrix, gapOpening,
				gapExtension, gapEnding, farIndels, myResult);

	//If everything is fine return the results in a comprehensive string form
	if (status == 0){

		result = myResult->to_string();
		result += "++++\n";
		result += myResult->to_ggplot2_string() + "\n";
		result += "++++\n";
		result += myResult->to_subject_coordiantes_string() + "\n";
		result += "++++\n";
		result += myResult->get_pattern_alignment() + "\n";
		result += "++++\n";
		result += myResult->get_subject_alignment() + "\n";
		
	}
	//Otherwise return an string error
	else{
		result = "Oh no! :-(";
	}

	return result;
}
