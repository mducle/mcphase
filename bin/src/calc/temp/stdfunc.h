// File: stdfunc.h
// $Id: stdfunc.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: stdfunc.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.2  1999/02/19 09:32:07  herbie
// Added GetArraysFrom...
//
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifndef STDFUNC_H
#define STDFUNC_H
                    
#include "stdinc.h"
#include "strings.h"

#define LINUX

#define CHECK_POINTER_EXIT(P,V)\
        if(!P){fprintf(stderr,"Out of memory in %s at line %d\n",\
                       __FILE__,__LINE__);\
                exit(V);\
               }
        
#define CHECK_POINTER_RETURN(P,V)\
        if(!P){fprintf(stderr,"Out of memory in %s at line %d\n",\
                       __FILE__,__LINE__);\
	       return(V);\
               }

#define CHECK_POINTER_RET(P)\
        if(!P){fprintf(stderr,"Out of memory in %s at line %d\n",\
                       __FILE__,__LINE__);\
	       return;\
               }

#define CHECK_FILE_POINTER_RETURN(P,NAME,V)\
        if(!P){fprintf(stderr,"Error opening file %s\n",NAME);\
                return(V);\
               }
        
#define CHECK_FILE_POINTER_EXIT(P,NAME,V)\
        if(!P){fprintf(stderr,"Error opening file %s\n",NAME);\
                exit(V);\
               }

#define CHECK_FILE_POINTER_RET(P,NAME)\
        if(!P){fprintf(stderr,"Error opening file %s\n",NAME);\
                return;\
               }

#define CHECK_INDEX_RETURN(I,IMAX,V)\
        if(I<0 || I>=IMAX)\
	   {fprintf(stderr,"%s->%s line %d\nIndex %d out of range (0 ...%d)\n",\
	            __PRETTY_FUNCTION__,__FILE__,__LINE__,I,IMAX-1);\
	    return (V);\
	   }

#define CHECK_INDEX_RET(I,IMAX)\
        if(I<0 || I>=IMAX)\
	  {fprintf(stderr,"%s->%s line %d\nIndex %d out of range (0 ...%d)\n",\
	                  __PRETTY_FUNCTION__,__FILE__,__LINE__,I,IMAX-1);\
	   return;\
	  }


#define PRINT_DEBUG(FORMAT,ARGS...)\
	{fprintf(stderr, FORMAT, ##ARGS);\
         fprintf(stderr,"in %s->%s line %d\n",\
	                __PRETTY_FUNCTION__,__FILE__,__LINE__);}

#define DEBUG_INFO\
        fprintf(stderr,"in %s->%s line %d\n",\
	__PRETTY_FUNCTION__,__FILE__,__LINE__);

 inline unsigned HRS(const double TinSEC)
                {return (unsigned) (TinSEC/3600.);}
 inline unsigned MINS(const double TinSEC) 
                {return unsigned ( ( (TinSEC/3600.) - (int)(TinSEC/3600.) )*60.);}

int CHECK_INDEX(const int i, const int iMax);
int CheckIndex(const int i, const int iMax);

// DOS or UNIX type line ends either with <lf><cr> or <lf>
enum EndLine { NOEnd = 0,  DOSEnd, UNIXEnd};

union Start {double fStart; int iStart;};
union End {double fEnd; int iEnd;};

  
// defined in stdfunc.cpp
#ifndef STDFUNC_CPP
extern char *szLorentzP[];
extern char *szLorentz12P[];
extern double RWLa21;
extern char *szCurieWeissP[];
extern char *szNonFermiResP[];
extern char *szStraightLineP[];
extern char *szPolynomP[];
#endif

#define FLOAT_NUM 0
#define INT_NUM 1        
#define UINT_NUM 2

#define TMPFILE_RETURN 1
#define TMPFILE_EXIT 2


inline double sign(double x) {return ( x < 0 ? -1 : ( x>0?1:0 )  ); }

  void Suspend(const double T);
  void Strupr(char *s);
  void RemoveCR_LF(char *szB);
   int DelSpace(char *szB);
   int ReplaceChar(char *szB, const int What, const int By);
   int CutRNum(char *szB); 
   int ReplaceLastChar(char *szB, const int What, const int By);
   int IsNumString(const char *s, const int NumFlag=FLOAT_NUM);
   int CountCols(char *szB, const char *szSeparator=" ");
   int GetCurrentTime(char *szB);
   int ShowFile(const char *szFile, FILE * fWhere);
  void Usage(const char *szPath, const char *szFile, FILE * fWhere);
  void Usage(const char *szPath, const char *szCmd, class String &Out);
  void PrintRCS(const char *szID, 
                const char *szPath, const char *szFile, FILE * fWhere);
unsigned AppendPath(const char *szP, const char *szF,
                    char * szR, const unsigned MaxL);
   int CopyAsciiFile(FILE *fpFrom, FILE * fpTo);
   int CreateTmpFile(char *szFile, char * szName, const int Flag=TMPFILE_EXIT);
   int GetFloatRange(char *szExp,double & lfStart, double & lfEnd, const char cDil=':');
   int GetIntRange(char *szExp,int & iStart, int & iEnd, const char cDil=':');
   int StrTok(const char *s, const char t, const int i, char *d);
   int GetFiles(const char *szArgs, char *szBuffer, const int iSize);
double fLinTrafo(double p_max,double p_min,double d_max,double d_min,double d);
double fBackTrafo(double p_max,double p_min, double b_max,double b_min,double b);
   int Scale(MDATA amin,MDATA amax, int max_res, MDATA &nmin, MDATA &nmax,
	     MDATA &ndelta);
int GetArrayFromLine(FILE *f, double *lfA, const int iLen,const char *szSep);
int GetArraysFromColumns(FILE *f, double **lfA, const int iRow, const int iCol,
                        const char *szSep);
int GetNumberFromString(String &S, double &lfA, const int iLen, const char cS);

double Horner(double *lfKoef, const int m, const double lfXw);

double Lorentz(const double x,const  int nP, double *p, const int nc, double *c);
double Lorentz12(const double x,const  int nP, double *p, const int nc, double *c);
double CurieWeiss(const double t, const int nP, double *p, const int nc, double *c);
double NonFermiRes(const double t, const int nP, double *p, const int nc, double *c);
double StraightLine(const double x,const  int nP, double *p, const int nc, double *c);
double Polynomal(const double x,const  int nP, double *p, const int nc, double *c);

inline double a2_2th(const double theta, const double Rwla21)
    {return 360./M_PI*asin(Rwla21*sin(theta*M_PI/180.));}


// *******************************************************
class Polynom
{
private:
    double * p;
      int    iDeg;
    double   Null;

public: 

  Polynom(void);
  Polynom(const int d, double *pi=NULL);
  Polynom(Polynom &N);
  ~Polynom();

        int  Degree(void) const {return iDeg;} 
   double &  operator [](const int i);
  Polynom &  operator=(Polynom &N);
     double  operator ()(const double x) const;
        int  Save(FILE *f, const double xS, const double xE, const int nPoints,
	          const char *szText, enum EndLine  LineT = NOEnd);
  
};
// **************************************************
#endif
