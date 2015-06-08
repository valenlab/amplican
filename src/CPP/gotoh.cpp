/* Main for testing the library*/

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>

#include "gotoh2.h"

using namespace std;

int  main( int argc, char ** argv ){

        //Sequences candidate to align
        string  sequencePattern=
        //"GTGCAGGATGATGCGGTCAGCGTGGAGAGCGGCACACACCCCGAGCGTCCCGACACTCCCACCAACACGGCGAGCGCTCCGGGCAGTGGAAGTCCAAGAAGTGCAAATACTCCTTCAAGTGTGAGATCGGAAGAGCACACGTCTGAACTC";
		//"TTTCCCTACACGACGCTCTTCCGATCTGTGCAGGATGATGCGGTCAGCGTGGAGAGCGGCACACACCCCGAGCGTCCCGACACTCCCACCAACACGGCGAGCGCTCCGGGCAGTGGAAGTCCAAGAAGTGCAAATACTCCTTCAAGTGTG";
        //"GGGAAGTGGAAGTCCAAGAAGTGCAAATACTCCTTCAAGTGTG";
        //"GATCTGTGCAGGATGATGCGGTCAGCGTGGAGAGCGGCACACACCCCGAGCGTCCCGACACTCCCACCAACACGGCGAGCGCTCCGGGCAGGAAGAGCTGGGGGAAGGGGAAGTGGAAGTCCAAGAAGTGCAAATACTCCTTCAAGTGTG";
        //"AGTGTGTCAGGACTGGTTTCATGGCAGGTGAGGTTGCAAAATTTTTAAATGACTGTTAAAAGCTAAATCTTGTACTTATGAATGTGTTTTTTGTGAATGTAATAATATATGTTTTTTCTTTTTGTGTAGTTGTGTTGGTGTGGAGGAGGA";
        //"CCTTCAAGTGTGAGATC";
        //"ccttcaagtgtg";
        //"GTGCAGGATGATGCGGTCAGCGTGGCGAGCGGCACACACCCCGAGCGTCCCGACACTCCCACCAACACGGCGAGCGCTCCGGGCAGGAAGAGCTGGGGGAAGGGGAAGTGGAAGTCCAAGAAGTGCAAATACTCCTTCAAGTGTGAGATC";
        //"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		//"AATATATATATACGCTCTCTAGCTAGGGACTCTCTCTCATTGGGCCTCCCCTCTTTACACGGACTATTTCGGAGCTCTACGCGCTACTAGCGACCTCTCGACTCGACTCTAGCAGCTAGCATCAGCGACACTCTCTACGACGACCTAGCTCGACTG";
        //"AATATATATAACATCG";
        "AATATCATCG";
        //"AGACCAATGG";
        //"AAAAAAAAAAAAAAACCCCTT";

        //"TGATTGAGTGTGTCAGGACTGG";
        //"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga";

        //"GCCTGTGTCGTCTCCCTTATGATGTCACCCGATTCATGATTGAGTGTGTTGGTTCTGAGTTCTGGGTTCGTGGCAGTGCAAAATTTTAAATTTTTTATTAAATGCTAAAACCTGAATCTTGTACGTATGAGTTTTTTTTGTGTGAATGTA";
        //"TTTTAGGAGTGCTTGTTTGACAGGTGAGGATGTTGAAAAATTTTTTAACTGTTAATAGATAATTAATGTTGTTCTTATGAATGTGTTTTTTGTGAATGTAATAATATATGTTTTTTTCTTTTTGTGTAGTTGTGTTGGTGTGGAGGAGGA";
        //"GCCTGTGTCGTCTCCCTTATGATGTCACCCGATTCATGATTGAGTGTGTCAGGACTGGTTTCATGGCAGGTGAGGTTGCAAAATTTTTAAATGACTGTTAAAAGCTAAATCTTGTACTTCTTATGAATGTGTTTTTTGTGAATGTAATAA";
		//"GTGTCAGGACTGGTTTCATGGCAGGTGAGGTTGCAAAATTTTTAAATGACTGTTAAAAGCTAAATCTTGTACTTCTTATGAATGTGTTTTTTGTGAATGTAATAATATATGTTTTTTCTTTTTGTGTAGTTGTGTTGGTGTGGAGGAGGA";
		//"GTGTCAGGACTGGTTT";
        
        string  sequenceSubject=
        //"gtgcaggatgatgcggtcagcgtggagagcggcacacaccCCGagcgtcccgacactcccaccaacacagcgagcgctccgggcaggaagagctgggggaaggggaagtggaagtccaagaagtgcaaatactccttcaagtgtg";
        //"gtgcaggatgatgcggtcagcgtggagagcggcacacaccCCGagcgtcccgacactcccaccaacacagcgagcgctccgggcaggaagagctgggggaaggggaagtggaagtccaagaagtgcaaatactccttcaagtgtg";
        //"gaagtggaagtccaagaagtgcaaatactccttcaagtgtg";
        //"gtgcaggatgatgcggtcagcgtggagagcggcacacaccCCGagcgtcccgacactcccaccaacacagcgagcgctccgggcaggaagagctgggggaaggggaagtggaagtccaagaagtgcaaatactccttcaagtgtg";
        //"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga";
        //"ccttcaagtgtg";
        //"CCTTCAAGTGTGAGATC";
        //"gtgcaggatgatgcggtcagcgtggagagcggcacacaccCCGagcgtcccgacactcccaccaacacagcgagcgctccgggcaggaagagctgggggaaggggaagtggaagtccaagaagtgcaaatactccttcaagtgtg";
        //"AAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAA";
        //"AATATATATATACGCCGGGGATATAGAAGATCCACGCTCTTCGACTCGCTCTCTAGCTAGGGACTCTCTCTCATTGGGCCTCCCCTCTTTACACGGACTATTTCGGAGCTCTAGCTAGCATCAGCGACACTCTCTACGACGACCTAGCTCGACTG";
        //"AATATCATCG";
        "AATATATATAACATCG";
        //"GACCAATGG";
        //"AAAAAAAAAAAAAAACCCCTTTTTTTTTTTTTTT";

        //"tgattgagtgtGATgtttgtcaggactgg";
        //"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga";
        
        //"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga";
        //"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga";
		//"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga";
		//"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtcaggactggtttcatggcaggtgaggttgcaaaatttttaaatgactgttaaaagctaaatcttgtacttcttatgaatgtgttttttgtgaatgtaataatatatgttttttctttttgtgtagttgtgttggtgtggaggagga";
		//"gcctgtgtcgtctcccttatgatgtcacccgattcatgattgagtgtGATgtttgtca";

        string results = "";

		//Show the sequences to the user
		cout<<"Sequence: "<<sequencePattern<<endl;
		cout<<"Mutated:  "<<sequenceSubject<<endl;

		//Resolve the alignment with the RCPP version. This version returns a
		//string with the alignment information.
		results = nwRCPP(sequencePattern, sequenceSubject,"NUC44",30,0,false,true);

		//Show the alignment to the user
		cout<<results;


        return  0 ;
}
