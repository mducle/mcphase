// File: simplex.h
// $Id: simplex.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: simplex.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef SIMPLEX_H
#define SIMPLEX_H 1

#include <math.h>

#include "cdata.h"
#include "stdfunc.h"

#define MAX_NAME_LEN 30
#define MAX_FN_NAME_LEN 80

#define STATUS_OK       0
#define STATUS_NOT_INIT 1
#define STATUS_ILL_DATA 2
#define STATUS_NO_MEM   3

#define STOP_MAXITER -1
#define STOP_LIMIT    1
#define ITERATE       0

#define MAX_FITFN 6

struct FitFn
   {   char FnName[MAX_FN_NAME_LEN+1];
        int nFixP,nSumP,nConstP;
     double (*FitFunc)(const double x, const int nP, double *p, const int nC=0, double *c=0);       
      char **szPName;
   };
//nSumP: number of pars to sum
//nFixP: number of fix pars  
//szPName: array of parameter names
// length: nSumP + nFixP
// nC: number of constant pars (not changed by fit)
// c: constant pars
// ****************************************************
class FitFunction
{protected:
  struct FitFn FF;

  double *lfSteps;

 public:
 FitFunction(const struct FitFn &ff);

 const char * GetFnName(void) const {return  FF.FnName;}
 const char * GetParName(const int iP) const;
 //(double (*f())) GetFn(void)     {return FF.FitFunc;} 
 
 struct FitFn GetFFStruct(void) const {return FF;}
 int GetNSumP(void) const {return FF.nSumP;} 
 int GetNFixP(void) const {return FF.nFixP;} 
 int GetNConstP(void) const {return FF.nConstP;} 

 virtual double GetStep(const int iP, const double lfP, const CRData *td=0) const = 0;

};
// ****************************************************
class LorentzFn : public FitFunction
{public:
 LorentzFn(void);

 //virtual 
  double GetStep(const int iP, const double lfP, const CRData *td=0) const;
};
// ****************************************************
class LorentzA12Fn : public FitFunction
{public:
 LorentzA12Fn(void);

 //virtual 
  double GetStep(const int iP, const double lfP, const CRData *td=0) const;
};
// ****************************************************
class CurieWeissFn : public FitFunction
{public:
 CurieWeissFn(void);

 //virtual 
  double GetStep(const int iP, const double lfP, const CRData *td=0) const;
};
// ****************************************************
class NonFermiResFn : public FitFunction
{public:
 NonFermiResFn(void);

 //virtual 
  double GetStep(const int iP, const double lfP, const CRData *td=0) const;
};
// ****************************************************
class StraightLineFn : public FitFunction
{public:
 StraightLineFn(void);

 //virtual 
  double GetStep(const int iP, const double lfP, const CRData *td=0) const;
};
// ****************************************************
class PolynomFn : public FitFunction
{public:
 PolynomFn(void);

 //virtual 
  double GetStep(const int iP, const double lfP, const CRData *td=0) const;
};
// ****************************************************
// ****************************************************
// CCCCCCCCCCCCCCCCCCCCCCCCCC
//    This procedure performs a simplex fitting algorithm (BYTE May 1984 p340)
//    for any number of (nonlinear) parameters
//   FitData:   pointer to data values
//   MaxIter:   nuber of max. iteration to be performed
//  *StSteps:  INPUT: 
//              pointer to array, holding the start values of parameter
// 	       OUTPUT:
//              StSteps[0]: residium (mean sum deviation) of fitted parameters
// 		step rate (size: nPars+1)
//    *Limit:  INPUT:
//              pointer to array holding max. error when iteration should
// 		stop (size: >= nPars)
//             OUTPUT:
//              estimated deviations of parameters
//     *Pars:   pointer to array of inital fitting parameters (size >= npar)
//     nPars:   number of Parameters
//   nPCount:   additional integer parameter of FitFunc
//       *FF:   pointer to structure which holds function to be fitted
//     nIter:   number of iterations really performed
//   CCCCCCCCCCCCCCCCCCCCCCCCC
// 
//  // *******************************************************

class Simplex
{
 double ALFA;    //=1.0    reflection cofficient >0
 double BETA;    //=0.5   contraction coefficient 0...1
 double GAMMA;   //=2.0   expansion coefficient >1

 CRData * FitData;
      int nData;
      int MaxIter;
  
      int nPars;
 double * StSteps;
 double * Limit;
 double * Pars;

      int ncPars;
 double * CPars;
 //double (*FitFunc)(const double x, const int nP, double *p); 
 const struct FitFn FF;

      int nSimp;
      int nPCount;
      int iStatus;

 double ** lfSimp;
  double * lfNext,*lfCenter,*lfError,*lfp,*lfq,*lfMean;      
     int * h,*l; // number of high,low parameters 
       int Done;
       int nIter;

//      char **PName;
        char *szParN;
      
     void Init(void);
     void Dealloc(void);
     void Assign(void);
          
   double SumResiduals(double *lfP);
     void Order(void);
     
public:
  Simplex(CRData *data, int max_iter, double *st_steps, double *limit,
	 int npar, double *pars, int p_count,
         const struct FitFn & fit_fn, int ncp=0, double *cpars=0);
  Simplex(CRData *data, int max_iter, double *st_steps, double *limit,
	 int npar, double *pars, int p_count,
         const FitFunction * Fit_fn_class,int ncp=0, double *cpars=0);
//	 double (*fit_fp)(const double x, const int pcount,double *p),
//	 const char *szN);
 ~Simplex() {Dealloc();}
 
          int Iterate(const int nMaxI=0);
         void NewIteration(const int nMaxI=0)
                          { nIter=0;  if(nMaxI>0)MaxIter=nMaxI; }
          int GetStatus(void)   const {return iStatus;}
          int GetNIter(void)    const {return nIter;}
          int GetMaxIter(void)  const {return MaxIter;}
 const char * GetPName(const int iP) const;	 
 const char * GetFnName(void) const {return FF.FnName;}	 
       double GetFitPar(const int iP);
       double GetErr(const int iP);
       double GetLimit(const int iP) const; 
         void GetMinMaxErr(double & MinErr, double & MaxErr,
                           int & iMi, int & iMa) const;
       double GetStStep(const int iP) const;
       double GetSigma(void) {return sqrt(GetErr(nSimp-1)/nData);}

         void SetLimit(const int iP, const double L);
         void SetStStep(const int iP, const double S);
//         void SetPName(const int iP, const char *szN);
};
//*************************************************************

#endif //SIMPLEX_H



