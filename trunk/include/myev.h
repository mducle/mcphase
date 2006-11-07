#include<vector.h>

#ifndef MYEV
#define MYEV

#define VERYSMALL 1e-04

// subs to be able to check and directly diagonalize hermitean
// matrizes, inverse a nearly singular matrix
void myPrintComplexMatrix(FILE * file,ComplexMatrix & M);
void myPrintMatrix(FILE * file,Matrix & M);
void myPrintVector(FILE * file,Vector & M);
void myEigenValuesHermitean (ComplexMatrix & M,Vector & lambda,int & sort,int & maxiter);
void myEigenSystemHermitean (ComplexMatrix & M,Vector & lambda,ComplexMatrix & l,int & sort,int & maxiter);
void myEigenSystemHermiteanGeneral (ComplexMatrix& a, ComplexMatrix& b, Vector & e, ComplexMatrix & T, int & sort, int & maxiter);

#endif
