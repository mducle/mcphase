//File: thefunc.h
// $Id: thefunc.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: thefunc.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#include <limits.h>

#ifndef THEFUNC_H
#define THEFUNC_H 1

#include "strings.h"
#include "cdata.h"

 // Header of new Data File:
 // [FileHeader]
 //Idn=1.0:THE (orig.meas.data)
 //Cell=N
 //p_Cell=10
 //p_Heater=20
 //H_Field=0.0
 //KMon=0.3097
 //DMon=-0.2664
 //HeliPot=4.2
 //CalConst=926.77
 //CalFile=NULL
 //SampleSensor=GaAs9257
 //HeaterSensor=9243
 //L0=1.23
 //SampleID=test 3 col
 //K0=0
 //;
 //StartTime=Fri Mar 01 15:42:19 1996
 //EndTime=Fri Mar 01 15:52:22 1996
 //NoofCols=#
 //DataPoints=6
 //ColID1=Tsample[K]
 //ColID2=C[pF]
 //ColId3=Tsample[V]
 //...
 //ColId#=LastV[??]
 //@EOH
 // N Data Points C1 C2 C3 ... C# separeted by blancs
// **********************************************************************

// ---For old versions should not be used
#define  FTYPE_ILLEGAL -1
#define    FTYPE_ASCII  0
#define FTYPE_MEMBRANE  1
#define      FTYPE_OLD  2
#define      FTYPE_NEW  3
#define     FTYPE_SOME  4
#define      FTYPE_RCP  6
#define    FTYPE_MULTI  7
// ----

typedef struct the_file_header
  { char szHead[82];
    char szText[82];
     int nSteps;
   float fL0;
   float fTSlope;
    char szPCell[20],szPHeat[20];
    char cCell;
    char cInsert;
   float fMField;
   float fKMon,fDMon,fCalConst,fHeliPot;
     //UMon=K*I+D (Powersupply)
  double lfK0;
    char szCalFile[MAXPATH+1];
    char szSSens[20]; char szHSens[20];
    char szSTime[30]; char szETime[30];
     int nCols;
  } THE_HEAD;


// *********************************************************

#define CC_ZERO_FILE 0x0001u
#define CC_CAL1_FILE 0x0002u
#define CC_CAL2_FILE 0x0004u
#define CC_CELL_DIAM 0x0008u
#define  CC_CELL_GAP 0x0010u
#define CC_CELL_PIVO 0x0020u
#define   CC_CELL_C0 0x0040u
#define CC_CELL_ZTEM 0x0080u
#define   CC_CELL_L0  0x0100u

typedef struct cap_cell
   {unsigned uFlag;
	char szZero[MAXPATH+1];
	char szCal[MAXPATH+1];
	char szCal2[MAXPATH+1];
      double fCDiam,fCGap;
      double fCPivot,fZeroTemp;
      double l0, K0;
   }CAP_CELL;

class DataFile;

// *******************************************************
int GetCellPars(CAP_CELL *CC, LineString &PP, LineString &TH);
int SPrintCC(CAP_CELL &CC, char *szB, const unsigned W=UINT_MAX);
int CalcGap(DataFile *DF, LineString &FP, char * szA=NULL);
int Calibrate(DataFile *DF, LineString &PP, char *szA=NULL);
int u2th(CRData *CR, const char * szF, const double lfH);
// *******************************************************

// *************************************************************
// class TiltPlate
// *************************************************************
#define EPS0 0.8854188E-11 //A sec/V m
#define FUNC_MICRO 0
#define FUNC_COMPENS 1
class TiltPlate
{protected:
 double ro; //Outer plate radius
 double ri; //Inner plate radius
 double b;  //Distance between centre of capacitor and pivot
 double ds; //Thickness of Sapphire whasher
 double k0; //pivot distance at T=300K

 double T;
 double C;

 int iErr;

 CRData *DLAg; //Therm. exp. of Ag
 CRData *DLSaph; //Therm. exp. of sapphire

 int RegulaFalsi(const double x0, const double x1, const double Eps,
				  const double Delta, double &x,
                                  int iFuncType=FUNC_MICRO);


 public:
 TiltPlate(const double Nro, const double Nri, const double Nb,
			  const double Nd, const double nk,
			  CRData *Ag, CRData *Saph=0);

  
 inline int IsValid(void){return !iErr;}
 
   void SetT_C(const double NT,const double C_pf){T=NT;C=C_pf*1E-12;}
                                  //Input C in pF formula requires Farad
   void CalcGap(CRData *Meas,int iFuncType=FUNC_MICRO);
 double DistFunc(const double x);
// new cell should be totaly compensated
 double CompDistFunc(const double x);
// Sapphire not relevant vor C(T)
 double SelectFunction(double x, int iFuncType);


};

// ****************************************************************
#endif
//THEFUNC_H
