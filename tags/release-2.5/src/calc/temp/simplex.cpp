//File: simplex.cpp
//$Id: simplex.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
//$Log: simplex.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.2  1999/03/15 09:08:37  herbie
//*** empty log message ***
//

#include <math.h>
#include <string.h>

#ifndef SIMPLEX_H
#include "simplex.h"
#endif

#ifndef STDINC_H
#include "stdinc.h"
#endif

const char *NError={"Illegal Index"};

// To initialize derived classes of FitFunction
const struct FitFn FIT_FUNCS[MAX_FITFN]=
       {{     "Lorentzian",1,3,0, Lorentz,szLorentzP}, 
        {"Lorentz Alfa 12",2,2,1, Lorentz12,szLorentz12P},
	{    "Curie-Weiss",3,0,0, CurieWeiss,szCurieWeissP},
	{    "NonFermiRes",3,0,0, NonFermiRes,szNonFermiResP},
	{  "Straight line",2,0,0, StraightLine,szStraightLineP},
        {        "Polynom",0,1,0, Polynomal,szPolynomP}
       };

// *************************************************************
// Class FitFunction
// *************************************************************
// *************************************************************
FitFunction::FitFunction(const struct FitFn &ff)
{strncpy(FF.FnName,ff.FnName,MAX_FN_NAME_LEN);
 FF.FnName[MAX_FN_NAME_LEN]=0;
 FF.nSumP=FF.nFixP=0;
 if(ff.nSumP<0 || ff.nFixP<0)return;
 FF.nSumP=ff.nSumP;
 FF.nFixP=ff.nFixP;
 FF.nConstP=ff.nConstP;
 FF.FitFunc=ff.FitFunc; 
 int i,iLen=FF.nSumP+FF.nFixP;
 FF.szPName=NULL;

 FF.szPName=new char * [iLen+FF.nConstP];
 CHECK_POINTER_RET(FF.szPName)
  for(i=0;i<iLen+FF.nConstP;i++)FF.szPName[i]=0;
 

 for(i=0;i<iLen+FF.nConstP;i++)
   {FF.szPName[i]=new char [MAX_NAME_LEN+1];
    CHECK_POINTER_RET(FF.szPName[i]);
    if(ff.szPName)
      {strncpy(FF.szPName[i],ff.szPName[i],MAX_NAME_LEN); 
       FF.szPName[i][MAX_NAME_LEN]=0;
      }
    else
      {sprintf(FF.szPName[i],"Par[%d]",i);}
   }

 
 lfSteps = new  double [iLen];
 CHECK_POINTER_RET(lfSteps);

} 
// *************************************************************
 const char * FitFunction::GetParName(const int iP) const
{if(iP>=0 && iP<FF.nSumP+FF.nFixP+FF.nConstP)return FF.szPName[iP];
 else return NError;
}
// *************************************************************
// Class Lorentz Functions
// *************************************************************
LorentzFn::LorentzFn(void) : FitFunction(FIT_FUNCS[0])
{ }
// *************************************************************
double LorentzFn::GetStep(const int iP, const double lfP, const CRData *td) const
{
  switch(iP)
   {case 0: // amplitude
           return lfP*0.1;
    case 1: // position x0
           return (td ? 0.1*(td->GetMaxX()-td->GetMinX()): 0.1*lfP);
    case 2: // half band width
           return 0.1*lfP;
    case 3: // background
           return 0.1*lfP;
    default:CHECK_INDEX_RETURN(iP,FF.nSumP+FF.nFixP-1,0)
            return 0;
   }	   	   	   
}
// *************************************************************
// Class LorentzAlfa12 Functions
// *************************************************************
LorentzA12Fn::LorentzA12Fn(void) : FitFunction(FIT_FUNCS[1])
{ }
// *************************************************************
double LorentzA12Fn::GetStep(const int iP, const double lfP, const CRData *td) const
{if(!td){PRINT_DEBUG("Illegal (zero) Pointer td\n")
         return 0;
	}
 switch(iP)
   {case 0: // amplitude
           return lfP*0.1;
    case 1: // position x0
           return (td ? 0.1*(td->GetMaxX()-td->GetMinX()): 0.1*lfP);
    case 2: // half band width
           return 0.1*lfP;
    case 3: // background
           return 0.1*lfP;
    default:CHECK_INDEX_RETURN(iP,FF.nSumP+FF.nFixP-1,0)
            return 0;
   }	   	   	   
}
// *************************************************************
// Class CurieWeiss Functions
// *************************************************************
CurieWeissFn::CurieWeissFn(void) : FitFunction(FIT_FUNCS[2])
{ }
// *************************************************************
double CurieWeissFn::GetStep(const int iP, const double lfP, const CRData *td) const
{if(td);
  switch(iP)
   {case 0: // c
           return lfP*0.1;
    case 1: // chi 0
           return 0.1*lfP;
    case 2: // theta
           return (td ? 0.1*(td->GetMaxX()-td->GetMinX()): 0.1*lfP);
    default:CHECK_INDEX_RETURN(iP,FF.nSumP+FF.nFixP-1,0)
            return 0;
   }	   	   	   
}
// *************************************************************
// *************************************************************
// Class NonFermiRes Functions
// *************************************************************
NonFermiResFn::NonFermiResFn(void) : FitFunction(FIT_FUNCS[3])
{ }
double NonFermiResFn::GetStep(const int iP, const double lfP, const CRData *td) const
{if(td);
  switch(iP)
   {case 0: // rho 0
           return lfP*0.1;
    case 1: // N2(Ef)
           return 0.1*lfP;
    case 2: // exponent
           return 0.01*lfP;
    default:CHECK_INDEX_RETURN(iP,FF.nSumP+FF.nFixP-1,0)
            return 0;
   }	   	   	   
}
// *************************************************************
// *************************************************************
// Class StraightLine Functions
// *************************************************************
StraightLineFn::StraightLineFn(void) : FitFunction(FIT_FUNCS[4])
{ }
double StraightLineFn::GetStep(const int iP, const double lfP, const CRData *td) const
{if(td);
  switch(iP)
   {case 0: // k
           return lfP*0.1;
    case 1: // d
           return 0.1*lfP;
     default:CHECK_INDEX_RETURN(iP,FF.nSumP+FF.nFixP-1,0)
            return 0;
   }	   	   	   
}
// *************************************************************
// *************************************************************
// Class StraightLine Functions
// *************************************************************
PolynomFn::PolynomFn(void) : FitFunction(FIT_FUNCS[5])
{ }
double PolynomFn::GetStep(const int iP, const double lfP, const CRData *td) const
{if(td);
 CHECK_INDEX_RETURN(iP,FF.nSumP+FF.nFixP-1,0)
 return lfP*0.1;
}
// *************************************************************
// *************************************************************
// Class simplex Functions
// *************************************************************
// *************************************************************
void Simplex::Init(void)
{
    ALFA = 1.0;   /* reflection cofficient >0*/
    BETA = 0.5;   /* contraction coefficient 0...1*/
   GAMMA = 2.0;  /* expansion coefficient >1*/
 
   FitData=NULL;
   MaxIter=0;
   nPars=nSimp=nPCount=0;
   iStatus=STATUS_NOT_INIT;
   lfSimp=0;
   lfNext=lfCenter=lfError=lfp=lfq=lfMean=0;
   h=l=0;
}
// ************************************************************
void Simplex::Assign(void)
{
 int i;
     
// memory allocation for simplex 
 iStatus=STATUS_NO_MEM;
 szParN=new char [MAX_NAME_LEN+5];
 if(szParN == NULL) return;
 
 lfSimp = new double * [nSimp];
 if (lfSimp == NULL)return;

 for(i=0;i<nSimp;i++)
   {lfSimp[i] = new double[nSimp];
    if (lfSimp[i] == NULL)return;
   }

// number of high,low parameters 

 h = new int [nSimp];
 if (h == NULL){Dealloc();return;}

 l = new int [nSimp];
 if (h == NULL){Dealloc();return;}

//lfNext:    new vertex to be tested
 lfNext = new double [nSimp];
 if (lfNext == NULL){Dealloc();return;}

// lfCenter:  center of hyperplane described of all vertexes exept the worst
 lfCenter = new double [nSimp];
 if (lfCenter == NULL){Dealloc();return;}

// lfError: Calculating deviation
 lfError = new double  [nSimp];
 if (lfError == NULL){Dealloc();return;}

// lfp,lfq:     helps calculating first simplex
 lfp = new double [nSimp];
 if (lfp == NULL){Dealloc();return;}
 
 lfq = new double [nSimp];
 if (lfq == NULL){Dealloc();return;}

// lfMean:    meanvalue of fittet parameters
 lfMean = new double [nSimp];
 if (lfMean == NULL){Dealloc();return;}

 for(i=0;i<nPars;i++)lfSimp[0][i]=Pars[i];  // start values to simplex

 lfSimp[0][nSimp-1]=SumResiduals(Pars); // first vertex

//calculate offset of vertex of starting simplex
 for(i=0;i<nPars;i++)
    {lfp[i]=StSteps[i]*(sqrt(nSimp) + (double)nPars - 1.)/((double)nPars*M_SQRT2);
     lfq[i]=StSteps[i]*(sqrt(nSimp) - 1.)/((double)nPars*M_SQRT2);
   }
 int j;
 for(i=1;i<nSimp;i++)
    {for (j=0;j<nPars;j++)lfSimp[i][j]=lfSimp[0][j]+lfq[j];
     lfSimp[i][i-1]=lfSimp[0][i-1]+lfp[i-1];
     lfSimp[i][nSimp-1]=SumResiduals(lfSimp[i]);
   }

 for(i=0;i<nSimp;i++)l[i]=h[i]=0;

 Order();
 iStatus=STATUS_OK;

}
// *************************************************************
void Simplex::Dealloc(void)
{
if(lfSimp)
   {int i;
    for (i=0;i<nSimp; i++)delete lfSimp[i];
    delete lfSimp;
   }
 delete [] lfNext;
 delete [] lfCenter;
 delete [] lfError;
 delete [] lfp; 
 delete [] lfq;
 delete [] lfMean;
 
 delete [] h;
 delete [] l;
 delete szParN;
}
// *************************************************************
Simplex::Simplex(CRData *data, int max_iter, double *st_steps, double *limit,
	         int npar, double *pars, int p_count,
                 const struct FitFn & fit_fn, int ncp, double *cpars) : FF(fit_fn)
{Init();
 if(data==0 || data->GetSteps()==0 ||
    max_iter<=0 || npar<=0 || p_count <=0 ||
    st_steps==0 || limit==0 ||pars==0 )
              {iStatus=STATUS_ILL_DATA; return;}
 FitData=data;
 nData=FitData->GetSteps(); 
  
 MaxIter=max_iter;
 nPars=npar;
 
 StSteps=st_steps;
 Limit=limit;
 Pars=pars;
 
 nPCount=p_count;
 //FF=fit_fn;
 nIter=0;
 Done=0;
 
 nSimp = nPars+1;
 ncPars=ncp;
 CPars=cpars;   

     
 Assign();
}
// ***********************************************************
Simplex::Simplex(CRData *data, int max_iter, double *st_steps, double *limit,
	 int npar, double *pars, int p_count,
         const FitFunction * Fit_fn_class,int ncp, double *cpars) : FF(Fit_fn_class->GetFFStruct())
{Init();
 if(data==0 || data->GetSteps()==0 ||
    max_iter<=0 || npar<=0 || p_count <=0 ||
    st_steps==0 || limit==0 ||pars==0 )
              {iStatus=STATUS_ILL_DATA; return;}
 FitData=data;
 nData=FitData->GetSteps(); 
  
 MaxIter=max_iter;
 nPars=npar;
 
 StSteps=st_steps;
 Limit=limit;
 Pars=pars;
 
 nPCount=p_count;
 nIter=0;
 Done=0;
 
 nSimp = nPars+1;

 ncPars=ncp;
 CPars=cpars;   
 Assign();

} 
// *******************************************************
double Simplex::SumResiduals(double *lfP)
{int i;
 double s;

 for (i=0,s=0;i<nData;i++)
     s+=SQUARE( (*FF.FitFunc)(FitData->GetPointX(i),nPCount,lfP,ncPars,CPars ) - 
                            (double) FitData->GetPointY(i));
 return s;
}
// *******************************************************
void Simplex::Order(void)
{int i,j;

 for (j=0;j<nSimp;j++)
     {for(i=0;i<nSimp;i++)
	 {if(lfSimp[i][j] < lfSimp[l[j]][j])l[j]=i;
	  if(lfSimp[i][j] > lfSimp[h[j]][j])h[j]=i;
	 }/*for(i)*/
     }/*for(j)*/
}
// *******************************************************
int Simplex::Iterate(const int nMaxI)
{
//return: -1  nIter >= MaxIter  (STOP_MAXITER)
//         1  error less than limits (STOP_LIMIT)
//         0 if iteration should go on (ITERATE)
//           -> limits greater than start limits and nIter < MaxIter
//

if(nIter==0)
  {if(nMaxI > 0)MaxIter=nMaxI;}
  
int i,j;
// do /*iteration loop*/
   Done = TRUE;
   nIter++;

   for(i=0;i<nSimp;i++)lfCenter[i]=0.; /*compute centeroid*/
   for(i=0;i<nSimp;i++)
      {if(i != h[nSimp-1]){for(j=0;j<nPars;j++)lfCenter[j]+=lfSimp[i][j];}
      }/*for(i)*/

   for(i=0;i<nSimp;i++)
      {lfCenter[i]/=nPars;
       lfNext[i]=(1. + ALFA) * lfCenter[i] - ALFA * lfSimp[h[nSimp-1]][i];
       /* next vertex specular reflection of the worst */
      }/*for(i)*/

   lfNext[nSimp-1]=SumResiduals(lfNext);

   if(lfNext[nSimp-1] <= lfSimp[l[nSimp-1]][nSimp-1])
     {for(i=0;i<nSimp;i++)lfSimp[h[nSimp-1]][i]=lfNext[i]; /*new vertex*/

      for(i=0;i<nPars;i++)
         lfNext[i]=GAMMA*lfSimp[h[nSimp-1]][i]+(1.-GAMMA)*lfCenter[i];
	       /*expand*/
      lfNext[nSimp-1]=SumResiduals(lfNext);

      if(lfNext[nSimp-1]<=lfSimp[l[nSimp-1]][nSimp-1])
	      {for(i=0;i<nSimp;i++)lfSimp[h[nSimp-1]][i]=lfNext[i];} /*new vertex*/
     }/*if*/

    else  /*not better than best  */
      { /*else 1*/
       if(lfNext[nSimp-1] <= lfSimp[h[nSimp-1]][nSimp-1])
	 {for(i=0;i<nSimp;i++)lfSimp[h[nSimp-1]][i]=lfNext[i];} /*new vertex*/

       else /*worse than worst*/
	 {/* else 2*/
	  for(i=0;i<nPars;i++)
	     lfNext[i]=BETA*lfSimp[h[nSimp-1]][i]+(1.-BETA)*lfCenter[i];
	  lfNext[nSimp-1]=SumResiduals(lfNext);
	  if(lfNext[nSimp-1] <= lfSimp[h[nSimp-1]][nSimp-1])
	    {for(i=0;i<nSimp;i++)lfSimp[h[nSimp-1]][i]=lfNext[i];} /*new vertex*/
	   else /*if still bad */
	    {/*else 3*/
	     for(i=0;i<nSimp;i++)
                {for(j=0;j<nPars;j++)lfSimp[i][j]=(lfSimp[i][j]+lfSimp[l[nSimp-1]][j])*BETA;
		 lfSimp[i][nSimp-1]=SumResiduals(lfSimp[i]);
		}/*for(i)*/

	     }/*else 3*/
	  }/*else 2*/
        }/*else 1*/

   Order();

   /*check for convergence*/
   for(j=0;j<nSimp;j++)
      {lfError[j]=(lfSimp[h[j]][j] - lfSimp[l[j]][j]) / lfSimp[h[j]][j];
       if(Done){if(lfError[j] > Limit[j])
                   {Done=FALSE;
	            //fprintf(stderr,"False @ %d\n",j);
		   }
	       }
      }/*for (j)*/

   //if(iOCtrl)
   //  sprintf(szB,"ITERATION: %d of %d  ERROR: %6.2g @ %2d / %6.2g @ %2d",
   //		n_iter,*max_iter,min_err,imi+1,max_err,ima+1);

// } /*iteration loop*/
// while( !done && n_iter < *max_iter );
//
if(Done)return STOP_LIMIT;
if(nIter >= MaxIter)return STOP_MAXITER;

return ITERATE;

} 
// *****************************************************
// *****************************************************
// void Simplex::SetParsToFit(void)
// {
//  /*average parameter*/
//  int i,j;
//  for(i=0;i<nSimp;i++)
//     {lfMean[i]=0.;
//      for(j=0;j<nSimp;j++)lfMean[i]+=lfSimp[j][i];
//      lfMean[i]/=(double)nSimp;
//     }
//  for(i=0;i<nPars;i++)nPars[i]=lfMean[i];
//  lfStSteps[0]=lfMean[nSimp-1]; //square sum
// }
// *******************************************************
double Simplex::GetFitPar(const int iP)
{CHECK_INDEX_RETURN(iP,nSimp,0)
 // iP = nPars+1 = nSimp-1 holds the sum of squares of deviation  
 int j;
 lfMean[iP]=0.;
 for(j=0;j<nSimp;j++)lfMean[iP]+=lfSimp[j][iP];
 lfMean[iP]/=(double)nSimp;
 return lfMean[iP]; 
}
// *******************************************************
double Simplex::GetErr(const int iP)
{CHECK_INDEX_RETURN(iP,nSimp,0)
 return lfError[iP]; 
}
// *******************************************************
double Simplex::GetLimit(const int iP) const 
{CHECK_INDEX_RETURN(iP,nSimp,0)
 return Limit[iP]; 
}
// *******************************************************
double Simplex::GetStStep(const int iP) const 
{CHECK_INDEX_RETURN(iP,nPars,0)
 return StSteps[iP]; 
}
// *******************************************************
void Simplex::GetMinMaxErr(double & MinErr, double & MaxErr,
                           int & iMi, int & iMa) const
{int j;
 for(j=0;j<nSimp;j++)
    {if(j==0)
       {MaxErr=MinErr=lfError[j];iMi=iMa=0;}
     else
       {if(lfError[j]>MaxErr){MaxErr=lfError[j];iMa=j;}
	if(lfError[j]<MinErr){MinErr=lfError[j];iMi=j;}
       }
    }
}    
// *******************************************************
void Simplex::SetLimit(const int iP, const double L)
{CHECK_INDEX_RET(iP,nPars)
 Limit[iP]=L;
}
// ********************************************************
void Simplex::SetStStep(const int iP, const double S)
{CHECK_INDEX_RET(iP,nPars)
 StSteps[iP]=S;
}
// ********************************************************
//void Simplex::SetPName(const int iP, const char *szN)
// {if(iP>0 && iP<nPars)
//    {strncpy(PName[iP],szN,MAX_NAME_LEN-1);
//     PName[MAX_NAME_LEN-1]=0;
//    }
//  else fprintf(stderr,"Simplex::SetPName: Index %d out of range (0 ... %d)\n",
//               iP,nPars);
// }
// ******************************************************* 
const char * Simplex::GetPName(const int iP) const	 
{// iP=0 returns first parameter
 if(iP<0 || iP>=nPars)sprintf(szParN,"%s: %d",NError,iP);

 //const int iSP=SumPars(FF.nSum);
 //const int iFP=FixPars(FF.nSum);

 if(iP<nPCount*FF.nSumP)
   sprintf(szParN,"%s[%d]",FF.szPName[iP%FF.nSumP],iP/(nPCount+1)+1);
 else
   sprintf(szParN,"%s",FF.szPName[(FF.nSumP==0 ? iP: (iP%FF.nSumP)+FF.nSumP )]); 
 return szParN;  
 }
// *******************************************************
