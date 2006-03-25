//File:spline.h
// $Id: spline.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: spline.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef SPLINE_H
#define SPLINE_H 1

#include "cdata.h"


// **********************************************************************
class SplineTab
  {//public:

   double * plfX,*plfY;
   double * plfD2;
        int nPoints;
       char szHead[83];

 public:
          SplineTab(double *pX, double *pY, const int nP);
   inline SplineTab() {plfX=plfY=plfD2=0; nPoints=0;}
          SplineTab(CRData *TD);

   ~SplineTab();
   inline int GetLength() {return nPoints;}
   
      void GetHeader(char *szH);
       int ReadTable(const char *szFileName);
       int SortData(void);
      void PrintTable(int bPrD2=0);
       int SaveTable(const char * szFileName);
       int Spline(float fD0=1.E30, float fDn=1.E30);
       int Intpol(const double lfX,double &lfY);
  };
#endif //SPLINE.H
