#include <iostream>
#include <string>
#include <algorithm> //reverse
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>

#include <Rcpp.h>


// [[Rcpp::export]]
void myAdding(int a, int b, int& c){

	c = a+b;

}
