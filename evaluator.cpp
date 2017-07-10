/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 * 
 *  Serial polynomial evaluation algorithm function implementations goes here
 * 
 */
#include <stdio.h>
#include "evaluator.h"

double poly_evaluator(const double x, const int n, const double* constants){
    //Implementation
    double rlt = constants[n-1];
    for(int i=n-2; i>=0; i--){
    	rlt = rlt * x + constants[i];
    }
    return rlt;
}
