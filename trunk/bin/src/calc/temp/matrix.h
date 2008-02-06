// $Id: matrix.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: matrix.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

#include "stdinc.h"
#include "stdfunc.h"
#include "strings.h"

class Matrix
{protected:

    MDATA ** pT;
 unsigned nRows,nCols;

    int *indx;

     int  Alloc(const unsigned iRows, const unsigned iCols, MDATA v=0);
     int  Alloc(const Matrix &N);
    void  Dealloc(void);
 MDATA ** ReallocColData(const unsigned nR, const unsigned nC);
     int  CopyData(const Matrix &N);
     int  CopyData(MDATA **N);
     int  LUDcmp(double *d);
     int  LUBakSub(double *b);     
public:
    Matrix(void);
    Matrix(const unsigned iRows, const unsigned iCols, const MDATA v=0);
    Matrix(MDATA **p, const unsigned iRows, const unsigned iCols);
    Matrix::Matrix(double *p, const unsigned iRows, const unsigned iCols);

    Matrix(const Matrix &N);
    virtual ~Matrix(void);

    const MDATA  * operator [](const unsigned i); 
    MDATA & operator ()(const unsigned iRow, const unsigned iCol);

     Matrix &  operator=(const Matrix &N);

     friend Matrix  operator+(Matrix &N1,Matrix &N2);
     friend Matrix  operator-(Matrix &N1,Matrix &N2);
     friend Matrix  operator*(Matrix &N1,Matrix &N2);
     
     friend Matrix  operator*(const MDATA N, Matrix &N2);

    inline   unsigned GetRows(void)       const {return nRows;}
    inline   unsigned GetCols(void)       const {return nCols;} 

    void Print(String &S);
    const MDATA * At(const unsigned iCol); 
      MDATA At(const unsigned iCol,const unsigned iRow) const ; 

     int ReadRow(FILE *InStream, const unsigned iCols, const unsigned nR);
     int ReadRows(LineString &LS, const unsigned iCols, const unsigned nR);
     
      MDATA Det(void);
    int Inv(void);
    int Transp(void);
};

// *************************************************************
                              
#endif //CDATA_H
