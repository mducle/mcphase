// $Id: cdata.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: cdata.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.2  1999/02/22 13:39:05  herbie
// Export format changed to 15.8g
//

#ifndef CDATA_H
#define CDATA_H

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

#include "stdinc.h"
#include "stdfunc.h"
#include "formula.h"
#include "strings.h"

#define  PR_ROW  0x0001U
#define  PR_COL  0x0002U
#define PR_DATA  0x0004U

#define MAX_ORDER 10

class CData
{
 friend class CRData;
 protected:
    MDATA * pCol;
        int nSteps;

     MDATA mMin,mMax;
       int iMin,iMax;
     
    void PrintObj(const char *s, FILE *f=stderr);   

public:
    CData(void);
    CData(const int iSteps, const MDATA v=0);
    CData(const CData &N);
    virtual ~CData(void);

    MDATA &  operator [](const int i);
    CData &  operator=(const CData &N);

    CData  operator+(const MDATA v);
    CData  operator+(const CData &N);

    CData  operator-(const MDATA v);
    CData  operator-(const CData &N);

    CData  operator*(const MDATA v);
    CData  operator*(const CData &N);

    CData  operator/(const MDATA v);
    CData  operator/(const CData &N);

    CData  operator-();

    friend CData  operator+(const MDATA v, CData &V);
    friend CData  operator-(const MDATA v, CData &V);
    friend CData  operator*(const MDATA v, CData &V);
    friend CData  operator/(const MDATA v, CData &V);
    friend MDATA  InProd(const CData &C1, const CData &C2);
    friend MDATA  Abs(const CData &V);
	      
    inline    int GetSteps(void)    const {return nSteps;}
    inline  MDATA GetMin(void)      const {return mMin;}
    inline  MDATA GetMax(void)      const {return mMax;}
    inline    int GetIMax(void)     const {return iMax;}
    inline    int GetIMin(void)     const {return iMin;}

    MDATA At(const int iNo) const;
    void Set(const MDATA NewP, const int iNo, const int bNew=TRUE);

    //int Resize(const int nS);
    int Insert(const MDATA NewP, const int bNew=TRUE);
    int Remove(const int iNo, const int bNew=TRUE);
    int NewRange(const int iStart, int iEnd,const int bNew=TRUE);
    int DelRange(const int iStart, const int iEnd, const int bNew=TRUE);
    int Copy(CData *TDSrc,const int iStart,const int iEnd);
   
    int Reduce(const int iR);

     int NewMinMax(int iStart=0, int iEnd=-1);
     int SortData(void);
     int BasicCalc(const char cOp, const MDATA Val);
     int BasicCalc(const MDATA Val, const char cOp);
     int BasicCalc(const char cOp, const CData &V);
     int BasicCalc(const char cOp, const CData &V1,const CData &V2);
     int Merge(CData *TD);
     int ApplyFunc(MDATA (*fp)(double), double x);
     int CalculateDl_l(const double lfL0);
     int CalculateGap(const double lfCDiam, const double lfCGap);
    void Print(void);
     int SPrintStatistic(char * szB, const unsigned Flag=UINT_MAX);
    void SPrintStatistic(String &S, const unsigned Flag=UINT_MAX);

 };
// **************************************************
class CRData
{protected:

    CData ** pT;
         int nCols;
         int iCx,iCy,iCz;

     int  Alloc(const int iCols, const int iSteps, MDATA v=0);
     int  Alloc(const CRData &N);
    void  Dealloc(void);
 CData ** ReallocColData(const int nC, const int nR);
     int  CopyData(const CRData &N);

public:
    CRData(void);
    CRData(const int iCols, const int iSteps, const MDATA v=0);

    CRData(const CRData &N);
    virtual ~CRData(void);

    CData &  operator [](const int i); //column 1st index, row 2nd index
    CRData &  operator=(const CRData &N);

    inline   int GetCols(void)       const {return nCols;}
    inline   int GetColX(void)       const {return iCx;}
    inline   int GetColY(void)       const {return iCy;}
    inline   int GetColZ(void)       const {return iCz;}

    int GetRows(const int iCol) const; 
    int GetSteps(void) const;
    
    void NewMinMax(int iStart=0, int iEnd=-1);
    MDATA GetColMax(int iCol) const;
    MDATA GetColMin(int iCol) const;

    int GetColIMax(int iCol) const;
    int GetColIMin(int iCol) const;

    void GetMinMax(MPoint &PMin, MPoint & PMax) const ;
    MDATA GetMinX(void) const {return GetColMin(iCx);}
    MDATA GetMinY(void) const {return GetColMin(iCy);}
    MDATA GetMaxX(void) const {return GetColMax(iCx);}
    MDATA GetMaxY(void) const {return GetColMax(iCy);}

    int GetIMinX(void) const {return GetColIMin(iCx);}
    int GetIMinY(void) const {return GetColIMin(iCy);}
    int GetIMaxX(void) const {return GetColIMax(iCx);}
    int GetIMaxY(void) const {return GetColIMax(iCy);}

    MDATA GetPointX(const int iRow) const {return At(GetColX(),iRow);}
    MDATA GetPointY(const int iRow) const {return At(GetColY(),iRow);}
    MDATA GetPointZ(const int iRow) const {return At(GetColZ(),iRow);}

    void SetPointX(const MDATA v, const int iRow) {(*this)[GetColX()][iRow]=v;}
    void SetPointY(const MDATA v, const int iRow) {(*this)[GetColY()][iRow]=v;}
    void SetPointZ(const MDATA v, const int iRow) {(*this)[GetColZ()][iRow]=v;}

     int SetColX(const int iZ);
     int SetColY(const int iZ);
     int SetColZ(const int iZ);
     int AssignXYZ(const int ixyz, const int iCol);
     int ExchgCol(const int iFrom, const int iTo);
     int Remove(const int iNo, const int bNew=TRUE);
     int InsertCol(const int iCol);
      
    void Print(void);
    CData * At(const int iCol); 
      MDATA At(const int iCol,const int iRow) const {return pT[iCol]->At(iRow);} 
        int SortData(const int iC);
        int ApplyFunc(MDATA (*fp)(double), const int iCSrc, const int iCDest);
	
    void SPrintStatistic(char * szB, const unsigned Flag=UINT_MAX);
    void SPrintStatistic(String & S, const unsigned Flag=UINT_MAX);

     int ExportData(FILE *OutSt, const enum EndLine LineT=NOEnd,
                    const int iCS=0, const int iCE=-1, const int iRS=0,
                    const int iRE=-1);
     int ReadCols(FILE *InStream, const int iSteps, const int nC,
                 const int iX=0, const  int iY=1, const  int iZ=2);

     int ReadCols(LineString &LS, const int iSteps, const int nC,
                  const int iX=0, const  int iY=1, const  int iZ=2);

     int Spline(const int iX=-1, const int iY=-1, const int iZ=-1);
//    int SplineIntpol(const double lfX,double &lfY);
     int Calculate(Formula &F);     
     int LinIntpol(MDATA Xp, MDATA &Yp, const int iX=-1, const int iY=-1);
     int MathOper(CRData *TD,const  char cSy, const int iX=-1, const int iY=-1);
     int Deriv( const int nP, const int iX=-1, const int iY=-1);
     int Integrate(const int iX=-1, const int iY=-1);
     int NewColRange(const int iCStart, const int iCEnd);
     int NewRowRange(const int iRStart, const int iREnd);
     int SelVal(const int iC, const double lfS, const double lfE);
     int DelRange(const int iCs, const int iCe, const int iRs, const int iRe);
     int DelVal(const int iC, const double lfS, const double lfE);
     int LookTable(CRData *Table, const int iInsert,
                                  const int iX=-1, const int iY=-1);
     int ZeroShift(MDATA xz, const int iX=-1, const int iY=-1);
     int Reduce(const int iR);
     int MergeRows(CRData &C);
     int PolyFit(Polynom &P, const int iOrder ,const int iX=-1, const int iY=-1);
                      
};

// *************************************************************
                              
#endif //CDATA_H
