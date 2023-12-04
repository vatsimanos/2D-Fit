#ifndef _TESTOUT_H_
#define _TESTOUT_H_

#include <iostream>
#include "declare.h"

void printTRate_matrix(TRate_matrix matrix, char* text = NULL);
void printVecDarray(VecD* vec, int x, char* text = NULL);
void printTUntermatrix_vector(std::vector<TUntermatrix> &unter, char* text = NULL);
void printTMultiSteady(TMultiSteady &multi, char* text = NULL);
void printTMultiSteadyP(double* multi, char* text = NULL);

#endif /* _TESTOUT_H_ */