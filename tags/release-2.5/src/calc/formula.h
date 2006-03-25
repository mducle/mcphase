//File formula.h
// $Id: formula.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: formula.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#include <stdio.h>

#ifndef FORMULA_H
#define FORMULA_H

#include "stdinc.h"

// ****************
// Class Formula
// ****************
enum ArgType {NoARG=-1, ColARG, ValueARG, FuncARG};

union SrcType
 {    int iSrcCol;
      int iFunc;
   double lfVal;
 };

struct SrcArg
{union SrcType SType; 
  enum ArgType AType; 
  char cOp; 
};

#ifndef FORMULA_MODULE
extern MDATA (*FFuns[]) (double x);
#endif

#define MAX_SRC 2

class Formula
{
 protected:
  char szF[MAX_LINELENGTH+1];
 
  int iDestCol;
  struct SrcArg Src[MAX_SRC];

  char *pAct;

  static const struct SrcArg NoSrc;
  static const char *szFun[];
  static const int iMaxFun;

  int IsFNumber(double &lfV);
  int IsColumn(char *s);
  int SetSource(char *s, const int iArg);

 public:
  Formula(const char *szNF);
  //struct SrcArg GetSrc(const int i) const; 
  int GetDestCol(void) const {return iDestCol;}

  const char * GetFormula(void) const {return szF;}
  
  enum ArgType GetType(const int i) const;
           int GetSrc(const int i)  const; 
	  char GetOp(const int i)   const;
	   int GetFunc(const int i) const;
	double GetVal(const int i)  const;
	 
  int CheckArg(void) const;
  int IsFunc(void) const;

  int Analyze(void);
  void Print(FILE *f=stderr);
};

// ************
// Class Range
// ************
union IntFloat { int i; double lf;};
enum  SqType {NO_SQ=-1, INT_SQ, FLOAT_SQ};

struct RangeS
 {union IntFloat St;
  union IntFloat En;
 };
 
class Range
{ 
   char szExp[MAX_LINELENGTH+1];
 struct RangeS R;
    int nRanges;
    int * pR;

  enum SqType Type;

 public:
 Range(const char *, const enum SqType T);
 ~Range(){delete pR;}

 enum SqType GetType(void) const {return Type;}

         int GetiStart(void) const {return R.St.i;}  
      double GetfStart(void) const {return R.St.lf;}

         int GetiEnd(void) const {return R.En.i;}  
      double GetfEnd(void) const {return R.En.lf;}
     
         int GetRanges(void) const {return nRanges;}
         int GetRange(const int i) const;
         int IsValid(void) const {return nRanges>=0;}
        void Print(void) const;
};
#endif //FORMULA_H
