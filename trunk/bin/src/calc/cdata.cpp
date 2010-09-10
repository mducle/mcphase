//File: cdata.cpp
// $Id: cdata.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: cdata.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.6  1999/04/29 09:45:06  herbie
// Remove-bug in MathOper corrected
//
// Revision 1.5  1999/03/15 09:08:37  herbie
// *** empty log message ***
//
// Revision 1.3  1999/02/22 13:39:44  herbie
// Export format changed to 15.8g
//
// Revision 1.2  1999/02/17 08:48:29  herbie
// Changed MAX_LINELENGTH to 1024.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "stdinc.h"
#include "stdfunc.h"
#include "cdata.h"
#include "formula.h"
#include "spline.h"

//#define TEST_CDATA
//#define DEBUG_CDATA

//extern int errno;
#include <errno.h>
// ******************************************************
// class CData Functions
// ******************************************************
// CData
// ******************************************************
// ******************************************************
 MDATA Null = 0;
 CData CNull;
CRData CRNull;
// ******************************************************
#if defined(DEBUG_CDATA)
void CData::PrintObj(const char *s, FILE *f)
{
  fprintf(f,"%s CData: %p, nSteps %d, pCol: %p\n",s,
           this,nSteps,pCol);  

}
#endif
// ******************************************************
CData::CData(void)
{
   pCol =  NULL;
 nSteps =  0;
   mMin =  0;
   mMax =  0;
   iMin = -1;
   iMax = -1;
#if defined (DEBUG_CDATA)
PrintObj("Def. Const");
#endif
}
// ******************************************************
CData::CData(const int iSteps, const MDATA v)
{// Error: nSteps =  = 0  Not enough memory
    pCol = NULL;
  nSteps = 0;
  
    pCol = new MDATA [iSteps];
  CHECK_POINTER_RET(pCol);

  nSteps = iSteps;

  for(int i = 0; i<nSteps; i++)pCol[i] = v;

    mMin = v;
    mMax = v;
    iMin = -1;
    iMax = -1;
#if defined (DEBUG_CDATA)
PrintObj("Const");
#endif
}
// ******************************************************
CData::CData(const CData &N)
{
  nSteps = 0;
    pCol =  new MDATA [N.nSteps];
  CHECK_POINTER_RET(pCol);

  nSteps = N.nSteps;
    mMin = N.mMin;
    mMax = N.mMax;
    iMin = N.iMin;
    iMax = N.iMax;

  for(int i=0; i<nSteps; i++) pCol[i]=N.pCol[i];
  NewMinMax();
#if defined (DEBUG_CDATA)
PrintObj("Copy Const");
fprintf(stderr,"From %p, nSteps %d, pCol %p\n",&N, N.nSteps, N.pCol);
#endif

}
// ******************************************************
CData::~CData(void)
{delete pCol;
#if defined (DEBUG_CDATA)
PrintObj("Destr");
#endif
 //pCol = NULL;
}
// ******************************************************
// Operators
// ******************************************************
MDATA & CData::operator[](const int i)
{if( !pCol || !nSteps)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return Null;
   }

 return pCol[CheckIndex(i,nSteps)];
}
// *********************************************************
CData & CData::operator=(const CData &N)
{ if(this != &N)
    {if(nSteps>0 && pCol)delete pCol;
       pCol = 0;
     nSteps = 0; 
     if(N.pCol && N.nSteps)pCol=new MDATA [N.nSteps];
     CHECK_POINTER_RETURN(pCol, *this)
     nSteps=N.nSteps;
     for(int i=0; i<nSteps; i++)pCol[i]=N.pCol[i];
     mMin = N.mMin;
     mMax = N.mMax;
     iMin = N.iMin;
     iMax = N.iMax;
#if defined (DEBUG_CDATA)
PrintObj("Op=");
#endif
    }
  return *this;
}
// ********************************************************
CData CData::operator+(const MDATA v)
{
 CData C(*this);  
 C.BasicCalc('+',v);
 return C;
}
// ********************************************************
CData CData::operator+(const CData & V)
{
 if(!pCol || nSteps==0 || !V.pCol || V.nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('+',V);
 return C;
}
// ********************************************************
CData CData::operator-(const CData & V)
{
 if(!pCol || nSteps==0 || !V.pCol || V.nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('-',V);
 return C;
}
// ********************************************************
CData CData::operator-(const MDATA v)
{
 if(!pCol || nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('-',v);
 return C;
}
// ********************************************************
CData CData::operator*(const CData & V)
{
 if(!pCol || nSteps==0 || !V.pCol || V.nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('*',V);
 return C;
}
// ********************************************************
CData CData::operator*(const MDATA v)
{
 if(!pCol || nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('*',v);
 return C;
}
// ********************************************************
CData CData::operator/(const CData & V)
{
 if(!pCol || nSteps==0 || !V.pCol || V.nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('/',V);
 return C;
}
// ********************************************************
CData CData::operator/(const MDATA v)
{
 if(!pCol || nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('/',v);
 return C;
}
// ********************************************************
CData CData::operator-()
{ if(!pCol || nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return *this;
   }
 CData C(*this);  
 C.BasicCalc('*',-1);
 
 return C;
}
// ********************************************************
// friend operators
// ********************************************************
CData operator+(const MDATA v,CData &V)
{
 CData C(V);
 C.BasicCalc(v,'+');
 return C;
}
// ********************************************************
CData operator-(const MDATA v, CData &V)
{
 CData C(V);
 C.BasicCalc(v,'-');
 return C;
}
// ********************************************************
CData operator*(const MDATA v, CData &V)
{
 CData C(V);
 C.BasicCalc(v,'*');
 return C;
}
// ********************************************************
CData operator/(const MDATA v, CData &V)
{
 CData C(V);
 C.BasicCalc(v,'/');
 return C;
}
// ********************************************************
MDATA  InProd(const CData &C1, const CData &C2)
{if( !C1.pCol || !C1.nSteps || !C2.pCol || !C2.nSteps)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   }
 if(C1.nSteps != C2.nSteps)
   {PRINT_DEBUG("Sizes are not equal %d != %d\n",C1.nSteps,C2.nSteps)
    return 0;
   }

 int i;
 MDATA s = 0;
 for(i=0; i<C1.nSteps; i++)s+=C1.pCol[i]*C2.pCol[i];
 return s;

}   
// *******************************************************
MDATA Abs(const CData & V)
{
 if(!V.pCol || V.nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   } 

 int i;
 MDATA s = 0;
 for(i=0; i<V.nSteps; i++)s+=V.pCol[i]*V.pCol[i];
 return sqrt(s);

}
// ***********************************************************
// ********************************************************
// Other member functions
// ********************************************************
void CData::Print(void)
{if(pCol==0 || nSteps==0){fprintf(stderr,"CData: Null\n"); return;}

 int i;
 fprintf(stderr,"Steps: %d\n",nSteps);
 for(i=0;i<nSteps;i++) fprintf(stderr,"%12.7g,",pCol[i]);
 fprintf(stderr,"\n");
} 
// ********************************************************
MDATA CData::At(const int iNo) const
{if( !pCol || !nSteps)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   }

 return pCol[CheckIndex(iNo,nSteps)];}
// ******************************************************
void CData::Set(const MDATA NewP, const int iNo, const int bNew)
{
  CHECK_INDEX_RET(iNo,nSteps)

  if( !pCol || !nSteps)
    {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
     return;
    }
  pCol[iNo]=NewP;
  if(bNew)NewMinMax();
}
// ******************************************************
int CData::Insert(const MDATA NewP, const int bNew)
{
 if(nSteps<=0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   } 

  MDATA *p = new MDATA[nSteps+1];
  CHECK_POINTER_RETURN(p,0)
  for(int i=0;i<nSteps;i++)p[i]=pCol[i];
  delete pCol;
  pCol = p;

  nSteps++;
  pCol[nSteps] = NewP;
  if(bNew)NewMinMax();
  return nSteps;
}
// ******************************************************
int CData::Remove(const int iNo, const int bNew)
{
  CHECK_INDEX_RETURN(iNo,nSteps,0)
  if( !pCol || nSteps<=1)
    {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
     return 0;
    }

  MDATA *p = new MDATA[nSteps-1];
  CHECK_POINTER_RETURN(p,0)

  int i,j=0;
  for(i=0;i<iNo;i++,j++)p[i]=pCol[i];
  for(i=iNo+1; i<nSteps; i++,j++)p[j]=pCol[i];

  delete pCol;
  pCol=p;
  nSteps--;

  if(bNew)NewMinMax();
  return nSteps;
}
// ******************************************************
int CData::NewRange(const int iStart, int iEnd, const int bNew)
{
 if(iEnd - iStart <= 0 || iStart < 0 )
   {PRINT_DEBUG("Illegal parameter. Cannot build NewRange\n")
    return 0;
   }

 CHECK_INDEX_RETURN(iEnd,nSteps,0);
 CHECK_INDEX_RETURN(iStart,nSteps,0);

  MDATA *p=new MDATA[iEnd-iStart+1];
  CHECK_POINTER_RETURN(p,0)

  int i,j;
    iEnd = (iEnd>nSteps?nSteps:iEnd);         //+1);
  for(i=iStart,j=0;i<=iEnd;i++,j++)p[j]=pCol[i];
  if(iEnd>nSteps)
    {for(i=nSteps;i<iEnd;i++)p[i]=0;}
    
  nSteps = iEnd-iStart+1;      //  +1;
  delete pCol;
    pCol = p;

  if(bNew)NewMinMax();
  return nSteps;
}
// ******************************************************
int CData::DelRange(const int iStart, int iEnd, const int bNew)
{
 if(iEnd - iStart <= 0 || iStart < 0 )
   {PRINT_DEBUG("Illegal parameter. Cannot build NewRange\n")
    return 0;
   }

 CHECK_INDEX_RETURN(iEnd,nSteps,0);
 CHECK_INDEX_RETURN(iStart,nSteps,0);

 int    iS1 = 0;
 int    iE1 = iStart+1;
 int    iS2 = iEnd+1;
 int    iE2 = nSteps;
 int iSteps = (iE1-iS1)+(iE2-iS2);

  MDATA *p = new MDATA[iSteps];
  CHECK_POINTER_RETURN(p,0)

  int i,j;
  for(i=iS1,j=0; i<iE1; i++,j++)p[j]=pCol[i];
  for(i=iS2; i<iE2; i++,j++)p[j]=pCol[i];

  nSteps = iSteps;
  delete pCol;
    pCol = p;

  if(bNew)NewMinMax();
  return nSteps;
}
// ******************************************************
int CData::Reduce(const int iR)
{
  if(iR<=1)
    {PRINT_DEBUG("Illegal parameter. Cannot reduce\n")
     return 0;
    } 
    int  i,j;
  MDATA *p = new MDATA[nSteps/iR + (nSteps%iR ? 1:0)];
  CHECK_POINTER_RETURN(p,0)

  for(i=0,j=0;i<=nSteps;i+=iR,j++)p[j]=pCol[i];

  nSteps = nSteps/iR + (nSteps%iR ? 1:0);
  delete pCol;
    pCol = p;

  NewMinMax();
  return nSteps;
}
// ******************************************************
int CData::Copy(CData *CSrc,int iStart, int iEnd)
{
 if(iEnd - iStart <= 0 || iStart < 0 || iEnd > CSrc->GetSteps())
   {PRINT_DEBUG("Illegal pars. Cannot copy\n")
    return 0;
   } 

 MDATA *p = new MDATA[iEnd-iStart+1];
 CHECK_POINTER_RETURN(p,0)

 int i,j;
 for(i=iStart,j=0;i<iEnd;i++,j++)p[j]=CSrc->At(i);

 nSteps = iEnd-iStart+1;
 delete pCol;
   pCol = p;


 NewMinMax();
 return nSteps;
}
// ******************************************************
int CData::NewMinMax(int iStart, int iEnd)
{/* Sucht Min und Max
	return  1: min/max gefunden
	return  0: iStart=iEnd
	return -1: xmin=xmax */

int i;
int iRet;

  if(pCol==0 || nSteps<=0)
    {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
     return 0;
    }
     
  if(iEnd==-1)iEnd=nSteps;
  mMin = pCol[iStart];
  mMax = pCol[iStart];
  iMin = iStart;
  iMax = iStart;
  if (iStart>iEnd)return 0;
  if (iStart<0 || iEnd<0)return 0;
  for(i=iStart;i<iEnd;i++)
     {if (pCol[i] > mMax){mMax=pCol[i];iMax=i;}
      if (pCol[i] < mMin){mMin=pCol[i];iMin=i;}
      }
  iRet=1;
  if(mMax==mMin)
    {mMax =  1.1*mMin;
     iRet = -1;
    }
  if(mMax==0 && mMin==0)
    {mMax =  0.1;
     mMin = -0.1;
     iRet = -1;}

 return iRet;
}
// *************************************************************
int CData::SortData(void)
{ int iXchCount;
  int i1,i2;

  if (pCol==0 || nSteps<=1)
     {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
      return 0;// empty
     } 
  do
   {iXchCount=0;
    for(i1=0,i2=1;i2<nSteps;i2++)
       {if(pCol[i1]>pCol[i2])
	  {SWAP(MDATA,pCol[i1],pCol[i2])
 	   iXchCount++;
	  }//if
    i1=i2;
       }//for
	 }//while
  while (iXchCount != 0);

 return 1;
}
// ***************************************************
int CData::BasicCalc(const char cOp,const MDATA Val)
{  int i;
 MDATA t;
 if(!pCol || nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   } 
 for(i=0;i<nSteps;i++)
    {switch (cOp)
	  {case '+':pCol[i]+=Val; break;
	   case '-':pCol[i]-=Val; break;
	   case '*':pCol[i]*=Val; break;
	   case '/':if(Val==0)
	              {PRINT_DEBUG("Division by zero (Index:%d)\n",i)
		       continue;
                      }
                    pCol[i]/=Val; break;
           case '^':
	   case '$':t=pow(pCol[i],Val);
                    if(errno)perror("CData::BasiCalc Illegal power value");
    	            else pCol[i]=t;
	   	    break;
    	   default: PRINT_DEBUG("Illegal operation char %c (%x)\n",cOp,cOp)
                    return 0;
	  }//switch 
     }//for

 NewMinMax();
 return nSteps;

}
// ********************************************************
int CData::BasicCalc(const MDATA Val,const char cOp)
{  int i;
 MDATA t;
 
 if(!pCol || nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   } 

 if(cOp=='+' || cOp=='*')return BasicCalc(cOp,Val);
	 
 for(i=0;i<nSteps;i++)
    {switch (cOp)
	  {case '-':pCol[i]=Val-pCol[i];
   		    break;
   	   case '/':if(pCol[i]==0)
	              {PRINT_DEBUG("Division by zero (Index:%d)\n",i)
		       continue;
                      }
		    pCol[i]=Val/pCol[i];
   		    break;
           case '^':
	   case '$':t=pow(Val,pCol[i]);
                    if(errno)perror("CData::BasiCalc Illegal power value");
    	            else pCol[i]=t;
		    break;
	   default: PRINT_DEBUG("Illegal operation char %c (%x)\n",cOp,cOp)
                    return 0;
	  }//switch
   }//for
 NewMinMax();
 return nSteps;
}
// **********************************************************
int CData::BasicCalc(const char cOp, const CData &V)
{  int i;
 MDATA t;
 
 if(!pCol || nSteps==0 || !V.pCol || V.nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   } 

 int iE = (nSteps<V.nSteps ? nSteps : V.nSteps);

 for(i=0;i<iE;i++)
    {switch (cOp)
	  {case '+':pCol[i]+=V.pCol[i]; break;
	   case '-':pCol[i]-=V.pCol[i]; break;
	   case '*':pCol[i]*=V.pCol[i]; break;
	   case '/':if(V.pCol[i]==0)
	              {PRINT_DEBUG("Division by zero (Index:%d)\n",i)
		       continue;
                      }
                    pCol[i]/=V.pCol[i]; break;
           case '^':
	   case '$':t=pow(pCol[i],V.pCol[i]);
                    if(errno)perror("CData::BasiCalc Illegal power value");
    	            else pCol[i]=t;
		    break;
       	   default: PRINT_DEBUG("Illegal operation char %c (%x)\n",cOp,cOp)
                    return 0;
	  }//switch 
     }//for

 NewMinMax();
 return iE;
}
// **********************************************************
int CData::BasicCalc(const char cOp, const CData &V1,const CData &V2)
{  int i;
 MDATA t;
 
 if(!pCol || nSteps==0 || !V1.pCol || V1.nSteps==0  
                       || !V2.pCol || V2.nSteps==0)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   } 

 int iE1 = (nSteps<V1.nSteps ? nSteps : V1.nSteps);
 int iE2 = (V1.nSteps<V2.nSteps ? V1.nSteps : V2.nSteps);
 int  iE = (iE1<iE2 ? iE1:iE2);

 for(i=0;i<iE;i++)
    {switch (cOp)
	  {case '+':pCol[i]=V1.pCol[i]+V2.pCol[i]; break;
	   case '-':pCol[i]=V1.pCol[i]-V2.pCol[i]; break;
	   case '*':pCol[i]=V1.pCol[i]*V2.pCol[i]; break;
	   case '/':if(V2.pCol[i]==0)
	              {PRINT_DEBUG("Division by zero (Index:%d)\n",i)
		       continue;
                      }
                    pCol[i]=V1.pCol[i]/V2.pCol[i]; break;
           case '^':
	   case '$':t=pow(V1.pCol[i],V2.pCol[i]);
                    if(errno)perror("CData::BasiCalc Illegal power value");
    	            else pCol[i]=t;
		    break;
       	   default: PRINT_DEBUG("Illegal operation char %c (%x)\n",cOp,cOp)
                    return 0;
	  }//switch 
     }//for

 NewMinMax();
 return iE;
}
// **********************************************************
int CData::ApplyFunc(MDATA (*fp)(double), double x)
{if( !pCol || !nSteps)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   }
 for(int i=0; i<nSteps; i++)pCol[i]=fp(x);

 NewMinMax();
 return nSteps;
}
// **********************************************************
int CData::Merge(CData *TD)
{
 int iE = TD->GetSteps();
 if(iE<=0 || TD==NULL)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   } 
 MDATA *p = new MDATA[nSteps+iE];
 CHECK_POINTER_RETURN(p,0)

    int i;
 for(i=0;i<nSteps;i++)p[i]=pCol[i];
 for(i=nSteps;i<nSteps+iE;i++)p[i]=TD->At(i-nSteps);

  nSteps+=iE;
  delete pCol;
  pCol = p;

  NewMinMax();
  return nSteps;
}
// ********************************************************
int CData::CalculateDl_l(const double lfL0)
{
 if(lfL0<=0)return 0;

 double fL = pCol[0];

 BasicCalc('-',fL);
 BasicCalc('/',lfL0);
 return 1;
}
// ******************************************************
int CData::CalculateGap(const double lfCDiam, const double lfCGap)
{ double fCKonst = 10.*0.0885*0.25*lfCDiam*lfCDiam*M_PI*-1e-3;
  // lfCDiam: Cell diameter in [cm]
  // lfCGap: Gap between Capacitor and ground in [cm] if -1 not used
  // Calculates gap in [mm]!!!
  if(lfCDiam<=0)return 0;
 if( !pCol || !nSteps)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   }
 if(lfCGap==-1.) BasicCalc(fCKonst,'/');
  else
  {double a,b,c,l1,l2,d;
   double r = 0.5*lfCDiam;
   for(int i=0;i<nSteps; i++)
      {a=pCol[i]/0.0885;
       b=0.22/0.0885*pCol[i]*lfCGap - r*r*M_PI - r*lfCGap*M_PI*(1.+lfCGap/lfCDiam);
       c=-0.22*r*r*M_PI*lfCGap;
       if( (d=b*b-4.*a*c)<0)return 0;
       l1=(-b+sqrt(d))/2/a;
       l2=(-b-sqrt(d))/2/a;
       if(l1<0 && l2<0)pCol[i]=0.;
       else {if(l1<0)pCol[i]=l2*10.;
       if(l2<0)pCol[i]=l1*10.;
      }
     }//for
  }//else
  return 1;

}
// **********************************************************************
int CData::SPrintStatistic(char * szB, const unsigned Flag)
{
 char *sp = szB;
 if(Flag & PR_ROW)sp+=sprintf(sp,"%d values : ",nSteps);

 if(Flag & PR_DATA)
   {NewMinMax();
    sp+=sprintf(sp,"%12.5g ... %12.5g\n",mMin,mMax);
   }
 return (int) (sp-&szB[0]);
}
// ******************************************************
void CData::SPrintStatistic(String &S, const unsigned Flag)
{
 if(Flag & PR_ROW)S.Addf("%d values : ",nSteps);

 if(Flag & PR_DATA)
   {NewMinMax();
    S.Addf("%12.5g ... %12.5g\n",mMin,mMax);
   }
}
// ******************************************************
// ******************************************************
// class CRData Functions
// ******************************************************
// CRData
// ******************************************************
// ******************************************************

// ******************************************************
#if defined(DEBUG_CDATA)
void CRData::PrintObj(const char *s, FILE *f)
{
 fprintf(f,"%s CRData: %p, nCols %d, pT: %p\n",s,
           this,nCols,pT);  
}
#endif
// ******************************************************
CRData::CRData(void)
{
    pT =  NULL;
 nCols =  0;
   iCx = -1;
   iCy = -1;
   iCz = -1;
#if defined (DEBUG_CDATA)
PrintObj("Def. Const");
#endif
}
// ******************************************************
CRData::CRData(const int iCols, const int iSteps, const MDATA v)
{// Error: nSteps==0  Not enough memory
      pT = NULL;
   nCols = 0;
  if(!Alloc(iCols,iSteps,v))return;
  nCols = iCols;
    iCx = -1;
    iCy = -1;
    iCz = -1;
#if defined (DEBUG_CDATA)
PrintObj("Const");
#endif
}
// ******************************************************
CRData::CRData(const CRData & N)
{ nCols = 0;
     pT = 0;
  if(!Alloc(N))return;
  nCols = CopyData(N);         
#if defined (DEBUG_CDATA)
PrintObj("Copy Const");
fprintf(stderr,"From %p, nCols %d, pT %p\n",&N, N.nCols, N.pT);
#endif
}
// ******************************************************
int CRData::CopyData(const CRData & N)
{int i;
 for (i=0;i<N.nCols;i++)*(pT[i])=*(N.pT[i]);
  iCx = N.iCx;
  iCy = N.iCy;
  iCz = N.iCz;
 return N.nCols;
}
// ******************************************************
int CRData::Alloc(const CRData & N)
{		
    pT = new CData * [N.nCols];
  CHECK_POINTER_RETURN(pT,0);

  for(int i=0; i<N.nCols; i++)
     {pT[i] = new CData(N.pT[i]->nSteps);
      CHECK_POINTER_RETURN(pT[i],0);
     }
 return 1;
}
// ********************************************************
int CRData::Alloc(const int iCols, const int iSteps, MDATA v)
{ if(iCols<=0 || iSteps<=0)
    {PRINT_DEBUG("Can not allocate %d columns or %d steps\n",iCols,iSteps);
     return 0;
    }
  pT = new CData * [iCols];
  CHECK_POINTER_RETURN(pT,0);

  for(int i=0; i<iCols; i++)
     {pT[i]=new CData(iSteps,v);
      CHECK_POINTER_RETURN(pT[i],0);
     }
  return 1;
}
// ******************************************************
void CRData::Dealloc(void)
{
  for(int i=0; i<nCols; i++)delete pT[i];

  delete pT;
  nCols = 0;

}
// *******************************************************
CData ** CRData::ReallocColData(const int nC, const int nR)
{CData **pD = 0;
 if(nC>0 && nR>1)
  {pD = new CData * [nC];
   CHECK_POINTER_RETURN(pD,0)
   int i;
   for(i=0;i<nC;i++)
      {pD[i] = new CData(nR);
       CHECK_POINTER_RETURN(pD[i],0)
      }//for
  }
 return pD;
 //nCols,nRows is NOT set to new values!!
}
// **************************************************
CRData::~CRData(void)
{Dealloc();
#if defined (DEBUG_CDATA)
PrintObj("Destr");
#endif
}
// ******************************************************
// Operators
// ******************************************************   
CData &  CRData::operator [](const int i)
{if( !pT || !nCols)
   {PRINT_DEBUG("Illegal (zero) data pointer pT\n")
    return CNull;
   }
 CHECK_INDEX_RETURN(i,nCols,CNull)
 return *(pT[i]);		 
} 
		
// ******************************************************   
CRData &  CRData::operator=(const CRData &N)
{ if(this != &N)
    {if(nCols>0 && pT)Dealloc();
     if(!Alloc(N)){pT=0; nCols=0; return *this;}
     nCols=CopyData(N);
#if defined (DEBUG_CDATA)
PrintObj("Op=");
#endif
    }
  return *this;
}
// ********************************************************
// Other Functions
// ******************************************************
int CRData::GetSteps(void) const 
{if(pT==0 || nCols<=0)return 0;
 int i1 = GetRows(0);
 int i;
 for(i=0; i<nCols; i++)
    {if(GetRows(i)!=i1)return -i1;}
 return i1;
}
// ******************************************************
void CRData::Print(void)
{if(nCols==0 || pT==0){fprintf(stderr,"CRData: Null\n"); return;}

 int i;
 fprintf(stderr,"nCols: %d\n",nCols);
 fprintf(stderr,"x:%d y:%d z:%d\n",iCx,iCy,iCz);
 for(i=0;i<nCols;i++) pT[i]->Print();
 fprintf(stderr,"\n\n");
} 
// ******************************************************
int CRData::GetRows(const int iCol) const
{CHECK_INDEX_RETURN(iCol,nCols,-1)
 return pT[iCol]->GetSteps();
}
// ******************************************************
void CRData::NewMinMax(int iStart, int iEnd)
{if(pT==0 || nCols==0)
   {PRINT_DEBUG("No data for Min/Max\n")
    return;
   }
 int i;
 for(i=0;i<nCols;i++)pT[i]->NewMinMax(iStart,iEnd);
}
// ******************************************************
MDATA CRData::GetColMax(int iCol) const
{CHECK_INDEX_RETURN(iCol,nCols,0.1)
 return pT[iCol]->mMax;
}
// ******************************************************
MDATA CRData::GetColMin(int iCol) const
{CHECK_INDEX_RETURN(iCol,nCols,-0.1)
 return pT[iCol]->mMin;
}
// *****************************************************
int CRData::GetColIMax(int iCol) const
{CHECK_INDEX_RETURN(iCol,nCols,-1)
 return pT[iCol]->iMax;
}
// ******************************************************
int CRData::GetColIMin(int iCol) const
{CHECK_INDEX_RETURN(iCol,nCols,-1)
 return pT[iCol]->iMin;
}
// ******************************************************
void CRData::GetMinMax(MPoint &PMin, MPoint & PMax) const
{CHECK_INDEX_RET(iCx,nCols)
 CHECK_INDEX_RET(iCy,nCols)

 PMin.x = pT[iCx]->mMin;
 PMin.y = pT[iCy]->mMin;
 PMax.x = pT[iCx]->mMax;
 PMax.y = pT[iCy]->mMax;
}
// ******************************************************			
int CRData::SetColX(const int iZ)
{CHECK_INDEX_RETURN(iZ,nCols,0)
 iCx = iZ;
 return 1;
}
// ******************************************************
int CRData::SetColY(const int iZ)
{CHECK_INDEX_RETURN(iZ,nCols,0)
 iCy = iZ;
 return 1;
}
// ******************************************************
int CRData::SetColZ(const int iZ)
{CHECK_INDEX_RETURN(iZ,nCols,0)
 iCz = iZ;
 return 1;
}
// ******************************************************
int CRData::AssignXYZ(const int ixyz, const int iCol)
{  switch(ixyz)
	{case 0:return SetColX(iCol);
	 case 1:return SetColY(iCol);
	 case 2:return SetColZ(iCol);
         default: return 0;
        }
}
// ******************************************************		
int CRData::ExchgCol(const int iFrom, const int iTo)
{
 CHECK_INDEX_RETURN(iFrom,nCols,0)
 CHECK_INDEX_RETURN(iTo,nCols,0)

  CData *t = pT[iFrom];
 pT[iFrom] = pT[iTo];
   pT[iTo] = t;

 return 1;
}
// ******************************************************
int CRData::Remove(const int iNo, const int bNew)
{if( !pT || !nCols)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   }

 int i;
 int ir = 1;
 for(i=0;i<nCols; i++)
    {if(!pT[i]->Remove(iNo,bNew))ir=0;}	 
 return ir;
}
// ******************************************************		
int CRData::InsertCol(const int iCol)
{
 if( !pT || !nCols)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   }

 if(GetRows(iCol)<=0)
   {PRINT_DEBUG("Cannot insert Column\n")
    return 0;
   }
 CHECK_INDEX_RETURN(iCol,nCols+1,0)
   
 CData ** p= new CData * [nCols+1];
 CHECK_POINTER_RETURN(p,0)

 int i,j;
 for(i=0,j=0;i<=nCols;i++)
    {if(i==iCol+1)
       {p[i]=new CData(GetRows(iCol));
        CHECK_POINTER_RETURN(p[i],0)
       }	
     else {p[i]=pT[j]; j++;}
    }

 delete pT;
 pT = p;
 nCols++;

 return 1;
}
// ******************************************************
CData * CRData::At(const int iCol)
{if( !pT || !nCols)
   {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
    return 0;
   }
 CHECK_INDEX_RETURN(iCol,nCols,0)
 return pT[iCol];		 
} 
// ******************************************************   
int CRData::SortData(const int iC)
{ int iXchCount;
  int i1,i2,j;

  if (nCols<=1 || pT==0)
     {PRINT_DEBUG("Illegal (zero) data pointer pT\n")
      return 0;// empty
     } 
  CHECK_INDEX_RETURN(iC,nCols,0)

  int iSteps = GetSteps();
  if(iSteps<=0)return 0;

  do
   {iXchCount=0;
    for(i1=0,i2=1;i2<iSteps;i2++)
       {if( pT[iC]->operator[](i1)> pT[iC]->operator[](i2))
	  {for(j=0;j<nCols;j++)SWAP(MDATA,pT[j]->operator[](i1),
                                          pT[j]->operator[](i2))
 	   iXchCount++;
	  }//if
        i1 = i2;
       }//for
   }//while
  while (iXchCount != 0);

 return 1;
}
// ***************************************************
int CRData::ApplyFunc(MDATA (*fp)(double), const int iCSrc, const int iCDest)
{ CHECK_INDEX_RETURN(iCSrc,nCols,0)
  CHECK_INDEX_RETURN(iCDest,nCols,0)

  if(pT[iCSrc]->GetSteps()!=pT[iCDest]->GetSteps())
    {PRINT_DEBUG("Steps of source and dest are not the same (%d != %d)\n",
                 pT[iCSrc]->GetSteps(),pT[iCDest]->GetSteps())
     return 0;
    }
     		 
  int i;
  for(i=0;i<pT[iCSrc]->GetSteps();i++)
      (*pT[iCDest])[i]=fp( (*pT[iCSrc])[i]);

  return pT[iCSrc]->GetSteps();
        
}
// ***************************************************
void CRData::SPrintStatistic(char * szB, const unsigned Flag)
{char *sp = szB;
 if(Flag & PR_COL)sp+=sprintf(sp,"%d data columns\n",nCols);
 int i;
 for (i=0;i<nCols;i++)
     {sp+=sprintf(sp,"Column %d",i+1);
      if(i==iCx)sp+=sprintf(sp,"X: ");
       else if(i==iCy)sp+=sprintf(sp,"Y: ");
          else if(i==iCz)sp+=sprintf(sp,"Z: ");
            else sp+=sprintf(sp," : ");
      sp+=pT[i]->SPrintStatistic(sp,Flag);
     }
}
// ******************************************************
void CRData::SPrintStatistic(String &S, const unsigned Flag)
{
 if(Flag & PR_COL)S.Addf("%d data columns\n",nCols);
 int i;
 for (i=0;i<nCols;i++)
     {S.Addf("Column %d",i+1);
      if(i==iCx)S.Addf("X: ");
       else if(i==iCy)S.Addf("Y: ");
          else if(i==iCz)S.Addf("Z: ");
            else S.Addf(" : ");
      pT[i]->SPrintStatistic(S,Flag);
     }
}
// ******************************************************
int CRData::ReadCols(FILE *InStream, const int iSteps, const int nC,
                    const int iX,const  int iY,const  int iZ)
{// Reads complete column data from file
 char szB[MAX_LINELENGTH+1],*s;

 if(pT && nCols>0)Dealloc();

 if(!Alloc(nC,iSteps))return 0;
 nCols = nC;

 int i,j;

 for(i=0; i<iSteps; i++)
   {if(fgets(szB,MAX_LINELENGTH,InStream)==NULL)break;
    szB[MAX_LINELENGTH]=0;
    while( (s=strchr(szB,'D'))!=NULL)*s='E';
    for(j=0;j<nC;j++)
       {if(j==0)s=strtok(szB," ");
	else s=strtok(NULL," ");
	if(s==NULL)break;
	(*this)[j][i]=atof(s);
       }//for j
    }//for i
 if(i<iSteps)return 0; //return i;

 for(i=0; i<nCols; i++)pT[i]->NewMinMax();

 iCx = (iX >= 0 && iX<nCols? iX:0);
 iCy = (iY >= 0 && iY<nCols? iY:0);
 iCz = (iZ >= 0 && iZ<nCols? iZ:0);

return iSteps;
}
// ***********************************************************
int CRData::ReadCols(LineString &LS, const int iSteps, const int nC,
                     const int iX,const  int iY,const  int iZ)
{// Reads complete column data from LineString
 // return 0: inittialisation error
 //        n>0: number of readed value lines
 //        n<0: error during read -n gives line

 if(LS.GetNLines()!=iSteps)return 0;
 if(pT && nCols>0)Dealloc();
 if(!Alloc(nC,iSteps))return 0;
 nCols = nC;

 int i,j;

 for(i=0; i<iSteps; i++)
   {String L(LS.Line(i));
    //fprintf(stderr,"Line: %d >%s\n",i,(const char *)L);
    for(j=0;j<nC;j++)
       {int iR=GetNumberFromString(L,(*this)[j][i],(j==0 ?nC:0),(j==0?' ':0));
        if(j==nC-1){if(iR!=0)return -i;}
        else if(iR<=0)return -i;        
       }//for j
    }//for i
 if(i<iSteps)return 0; //return i;

 for(i=0; i<nCols; i++)pT[i]->NewMinMax();

 iCx = (iX >= 0 && iX<nCols? iX:0);
 iCy = (iY >= 0 && iY<nCols? iY:0);
 iCz = (iZ >= 0 && iZ<nCols? iZ:0);

return iSteps;
}
// ***********************************************************
int CRData::ExportData(FILE * OutSt, const enum EndLine LineT, 
                       const int iCS, const int iCE,
		       const int iRS, const int iRE)
{
 int iCStart = iCS;
 int   iCEnd = (iCE==-1?nCols:iCE);
 int iRStart = iRS;
 int   iREnd = (iRE==-1?GetSteps():iRE);

 if(iCStart<0 || iCEnd<1 || iRStart<0 || iREnd<1 ||
    iCStart>=iCEnd || iRStart>=iREnd || iCEnd>nCols || iREnd>GetSteps())
   {PRINT_DEBUG("Illegal row/column selection\nStart(%d,%d) End:(%d,%d)\n",
                  iCStart, iRStart, iCEnd, iCStart)
   
    return 0;
   }
   
 if(pT==NULL || nCols==0)
   {PRINT_DEBUG("Illegal (zero) data pointer\n");
    return 0;
   }

 char szB[MAX_LINELENGTH+1];
 char szH[80];

 char szE[3] = {"\n"};
#if defined(LINUX)
  if(LineT==NOEnd || LineT==DOSEnd)strcpy(szE,"\r\n");
#endif 

 int i,j;

 szH[0] = 0;
 for(i=iRStart;i<iREnd;i++)
   {szB[0] = 0;
    for(j=iCStart;j<iCEnd;j++)
       {sprintf(szH,"%15.8g  ",(*this)[j][i]);
        strcat(szB,szH);
       }// for j
    fprintf(OutSt,"%s%s",szB,szE);
   }//for i
 return 1;
}
// ******************************************************
int CRData::Spline(const int iX, const int iY, const int iZ)
{ 

 int ix = (iX==-1?iCx:iX);
 int iy = (iY==-1?iCy:iY);
 int iz = (iZ==-1?iCz:iZ);

 CHECK_INDEX_RETURN(ix,nCols,0)
 CHECK_INDEX_RETURN(iy,nCols,0)
 CHECK_INDEX_RETURN(iz,nCols,0)

 int iSteps=GetSteps();
 if(iSteps<3)return 0;

  double p,qn,sig,un;
  double fD0 = 1.E30;
  double fDn = 1.E30;

 double *u = new double[iSteps];
 if(u==0)return 0;

 if(fD0>0.99E30)(*pT[iz])[0]=u[0]=0;
 else
 {(*pT[iz])[0]=-0.5;
  u[0]=(3.0/((*pT[ix])[1]-(*pT[ix])[0])) * (((*pT[iy])[1]-(*pT[iy])[0]) /
                                            ((*pT[ix])[1]-(*pT[ix])[0])-fD0);
 }

 for(int i=1;i<iSteps-1;i++)
    {  sig = ((*pT[ix])[i]-(*pT[ix])[i-1])/((*pT[ix])[i+1]-(*pT[ix])[i-1]);
         p = sig*(*pT[iz])[i-1]+2.0;
     (*pT[iz])[i] = (sig-1.0)/p;
      u[i] = ((*pT[iy])[i+1]-(*pT[iy])[i])/((*pT[ix])[i+1]-(*pT[ix])[i]) -
                   ((*pT[iy])[i]-(*pT[iy])[i-1])/((*pT[ix])[i]-(*pT[iy])[i-1]);
      u[i] = (6.0*u[i]/((*pT[ix])[i+1]-(*pT[ix])[i-1])-sig*u[i-1])/p;
    }
 if(fDn > 0.99e30) qn=un=0.0;
 else
 {qn = 0.5;
  un = (3.0/((*pT[ix])[iSteps-1]-(*pT[ix])[iSteps-2])) *
   (fDn-(*pT[iy])[iSteps-1]-(*pT[iy])[iSteps-2])/((*pT[ix])[iSteps-1]-(*pT[ix])[iSteps-2]);
 }
 (*pT[iz])[iSteps-1] = (un-qn*u[iSteps-2])/(qn*(*pT[iz])[iSteps-2]+1.0);

 for(int k=iSteps-2; k>0; k--)(*pT[iz])[k]=(*pT[iz])[k]*(*pT[iz])[k+1]+u[k];

 delete u;
 return 1;

}
// ******************************************************
int CRData::Calculate(Formula &F)
{//return number of calculated values (>0) if OK!!
 // negative or zero if an error occured

  int i = F.Analyze();
  if(i)return -i;
  F.Print();
  CHECK_INDEX_RETURN(F.GetDestCol()-1,nCols,0)

  //CHECK_INDEX_RETURN(F.GetSrc(0),nCols,0)
  //CHECK_INDEX_RETURN(F.GetSrc(1),nCols,0)

  if(!F.CheckArg())return -3;
  
  if(!F.IsFunc())
    {if(F.GetType(0)==NoARG)return 0;
     switch(F.GetType(0))
       {case ColARG:
           if(F.GetType(1)==ValueARG)
              return pT[F.GetDestCol()-1]->BasicCalc(F.GetOp(0),F.GetVal(1));
           if(F.GetType(1)==ColARG)
              return pT[F.GetDestCol()-1]->BasicCalc(F.GetOp(0),(*pT[F.GetSrc(0)-1]),
	                                                        (*pT[F.GetSrc(1)-1]));
           PRINT_DEBUG("Uuups, internal logical program error\n")		    	                  
           return 0;
	case ValueARG:
	   if(F.GetType(1)==ColARG)
              return pT[F.GetDestCol()-1]->BasicCalc(F.GetVal(0),F.GetOp(0));
           PRINT_DEBUG("Uuups, internal logical program error\n")		    	                  
           return 0;
        default:
	   PRINT_DEBUG("Uuups, internal logical program error\n")
           return 0;
      }// case first arg
     PRINT_DEBUG("Illegal Formula %s\n",F.GetFormula())		    	                  
     return 0;
    }// No function call  
  else
    {if(F.GetType(1)==FuncARG)
       {PRINT_DEBUG("Illegal formula %s\n",F.GetFormula())		    	                  
        return 0;
       }	
     if(F.GetType(1)==ValueARG)  
       return pT[F.GetDestCol()-1]->ApplyFunc(FFuns[F.GetFunc(0)],F.GetVal(1));  

     if(F.GetType(1)==ColARG)  
       return ApplyFunc(FFuns[F.GetFunc(0)],F.GetSrc(1)-1,F.GetDestCol()-1);  
     PRINT_DEBUG("Illegal Formula %s\n",F.GetFormula())		    	                  
     return 0;
    }//Function call
  PRINT_DEBUG("Illegal Formula %s\n",F.GetFormula())		    	                  
  return 0;
}
// ******************************************************
int CRData::LinIntpol(MDATA Xp, MDATA &Yp, const int iX, const int iY)
{
 int ix = (iX==-1?iCx:iX);
 int iy = (iY==-1?iCy:iY);

 CHECK_INDEX_RETURN(ix,nCols,0)
 CHECK_INDEX_RETURN(iy,nCols,0)

 int iSteps = GetSteps();
 if(iSteps<3)return 0;

 int    k;
 int    iLo = 0;
 int    iHi = iSteps-1;

 if(Xp > (*pT[ix])[iSteps-1] || Xp < (*pT[ix])[0])return 0;

 while(iHi-iLo>1)
	  {k=(iHi+iLo) >> 1;
		if( (*pT[ix])[k] > Xp)iHi=k;
		else iLo=k;
	  }
 if(iHi==iLo)return 0;

 MDATA Yk =( (*pT[iy])[iHi] - (*pT[iy])[iLo]) / 
	   ( (*pT[ix])[iHi] - (*pT[ix])[iLo])*(Xp - (*pT[ix])[iLo]);
 Yp = (*pT[iy])[iLo]+Yk;

 return 1;
}
// ******************************************************
int CRData::MathOper(CRData *TD,const  char cSy, const int iX, const int iY)
{
 // iCx, iCy of this must be set
 CHECK_INDEX_RETURN(iCx,nCols,0)
 CHECK_INDEX_RETURN(iCy,nCols,0)

 //int iSteps = GetSteps();
 if(GetSteps()<3)return -1;

 int iX2=(iX==-1 ? iCx: iX);
 int iY2=(iY==-1 ? iCy: iY);
 
 CHECK_INDEX_RETURN(iX2,TD->GetCols(),0);
 CHECK_INDEX_RETURN(iY2,TD->GetCols(),0);
 
  if(!SortData(iCx))return -1;
  if(!TD->SortData(iX2))return -1;

 NewMinMax();
 TD->NewMinMax();


 if(TD->GetColMin(iX2) > GetColMax(iCx) || GetColMin(iCx) > TD->GetColMax(iX2))
   return -1;
     		 
 MDATA Yn;
 int i,iSC;
 for(i=0,iSC=0;i<GetSteps();i++)
    {if(TD->LinIntpol((*pT[iCx])[i], Yn, iX2, iY2)==0)
       {iSC++;
	Remove(i,FALSE);
	i--;
       }//if
     else
       {switch (cSy)
	  {case '+': (*pT[iCy])[i]+=Yn;break;
	   case '-': (*pT[iCy])[i]-=Yn;break;
	   case '*': (*pT[iCy])[i]*=Yn;break;
	   case '/': if(Yn!=0)(*pT[iCy])[i]/=Yn;
	             else (*pT[iCy])[i]=HUGE_VAL;
                     break;
           default: return -1;
  	  }//switch
	}//else
     }//for
 SortData(iCx);
 NewMinMax();
 return iSC;
}
// ******************************************************
int CRData::Deriv(const int nP, const int iX, const int iY)
{
 int ix = (iX==-1?iCx:iX);
 int iy = (iY==-1?iCy:iY);

 CHECK_INDEX_RETURN(ix,nCols,-1)
 CHECK_INDEX_RETURN(iy,nCols,-1)

 int iSteps = GetSteps();
 if(iSteps<nP+1)return -1;
 if(nP<1)return -1;

 SortData(ix);
 int i,nErrs;
 MDATA Dx = 0;
 MDATA Dy = 0;

  for(i=0,nErrs=0;i<iSteps-nP;i++)
     {Dx=Dy=0.;
      for (int j=i+nP;j>=i+1;j--)
	  {Dx+=(*pT[ix])[j] - (*pT[ix])[i];
	   Dy+=(*pT[iy])[j] - (*pT[iy])[i];
	  }
	  Dx/=nP;
	  Dy/=nP;
	  if(Dx==0)nErrs++;
	  else (*pT[iy])[i] = Dy/Dx;
	 }
 for(i=0;i<nCols; i++)pT[i]->NewRange(0,iSteps-nP-1);
 return nErrs;
}
// **************************************************************
int CRData::Integrate(const int iX, const int iY)
{
 int ix = (iX==-1?iCx:iX);
 int iy = (iY==-1?iCy:iY);

 CHECK_INDEX_RETURN(ix,nCols,0)
 CHECK_INDEX_RETURN(iy,nCols,0)

 int iSteps = GetSteps();
 if(iSteps<3)return 0;

 SortData(ix);
 MDATA SumI = 0;
 MDATA   Dx,y; 

 int i;
 for(i=1;i<iSteps;i++)
   { Dx = (*pT[ix])[i]-(*pT[ix])[i-1];
      y = 0.5*((*pT[iy])[i-1]+(*pT[iy])[i])*Dx;
    SumI+=y;
    (*pT[iy])[i-1] = SumI;
   }

 for(i=0;i<nCols; i++)pT[i]->NewRange(0,iSteps-1);
 return 1;
}
// **************************************************************
int CRData::NewColRange(const int iCStart, const int iCEnd)
{
 if(iCEnd - iCStart <= 0 || iCStart < 0 )
   {PRINT_DEBUG("Illegal parameter. Cannot build NewRange\n")
    return 0;
   }

 CHECK_INDEX_RETURN(iCEnd,nCols,0)
 CHECK_INDEX_RETURN(iCStart,nCols,0)

     int iCols = iCEnd-iCStart+1;
  CData ** npT = new CData * [iCols];
  CHECK_POINTER_RETURN(npT,0);

  int i,j;

  for(i=iCStart,j=0;i<=iCEnd;i++,j++)npT[j]=pT[i];

  for(i=0;i<nCols; i++)
     {if(i<iCStart ||i>iCEnd)delete pT[i];}

  nCols = iCols;
  delete pT;
     pT = npT;

  return nCols;
}
// ******************************************************
int CRData::NewRowRange(const int iRStart, const int iREnd)
{int ir = 1;
 int  i;
 for(i=0;i<nCols; i++){if(!pT[i]->NewRange(iRStart,iREnd))ir=0;}
 return ir;
}
// **************************************************************
int CRData::SelVal(const int iC, const double lfS, const double lfE)
{// Select those rows which have the values lfS-lfE in column iC

 CHECK_INDEX_RETURN(iC,nCols,-1)
 
 int i,j,icc;

 if(GetSteps()<=0)
   {PRINT_DEBUG("Illegal number of rows\n")
    return -1;
   }

 for(i=0,icc=0; i<GetSteps(); i++)
    {double t = (*pT[iC])[i];
     if(t>=lfS && t<=lfE)icc++;
    }

 if(icc==0)
   {PRINT_DEBUG("Values %f:%f not found\n",lfS,lfE)
    return -1;
   }
   int  iCols = nCols;       
 CData   **pm = ReallocColData(iCols,icc);
 CHECK_POINTER_RETURN(pm,0)
 int  ii;
 for(i=0,ii=0; i<GetSteps(); i++)
    {if((*pT[iC])[i]>=lfS && (*pT[iC])[i]<=lfE)
       {for(j=0;j<nCols;j++)(*pm[j])[ii]=(*pT[j])[i];
	ii++;   
       }//if
    }//for i

 Dealloc();
    pT = pm;
 nCols = iCols;
 NewMinMax();
 return icc;   
}    
// **************************************************************		
int CRData::DelRange(const int iCs, const int iCe, const int iRs, const int iRe)
{
 
  if(GetSteps()<=0 || iCe<iCs)
   {PRINT_DEBUG("Illegal range selection rows: %d, %d cols: %d, %d\n",
                iRs,iRe,iCs,iCe)

    return -1;
   }       		  		 


 int ice1 = iCs;
 int ics1 = 0;
 if(iCe==-1 && iCs==-1)
   {ics1 = 0;
    ice1 = nCols;}

 int ics2 = iCe+1;
 int ice2 = nCols;
 if(iCe==-1 && iCs==-1)
    {ics2 = 0;
     ice2 = 0;}

 CHECK_INDEX_RETURN(ics1,nCols,-1)
 CHECK_INDEX_RETURN(ice1-1,nCols,-1)
 CHECK_INDEX_RETURN(ics2-1,nCols,-1)
 CHECK_INDEX_RETURN(ice2-1,nCols,-1)

 int irs = iRs;
 int ire = iRe;
 if(iRs==-1 && iRe==-1){irs=0; ire=0;}

 CHECK_INDEX_RETURN(irs,GetSteps(),-1)
 CHECK_INDEX_RETURN(ire,GetSteps(),-1)

 int iCols = (ice1-ics1)+(ice2-ics2);       
 int iRows = ( irs==0 && ire==0 ? GetSteps() : GetSteps()-(ire-irs+1) );
 int     i;

 if(! (irs==0 && ire==0 ) )
   {for(i=0; i<nCols; i++)pT[i]->DelRange(irs,ire,FALSE);}


 
 if(! (iCs==-1 && iCe==-1))
   { CData **pm = new CData * [iCols];
     CHECK_POINTER_RETURN(pm,0);
     int  ii;
     for(i=ics1,ii=0; i<ice1; i++,ii++)pm[ii]=pT[i];

     for(i=ics2; i<ice2; i++,ii++)pm[ii]=pT[i];

     for (i=iCs; i<=iCe; i++)delete pT[i];
 
     nCols = iCols;
        pT = pm;
   }

 NewMinMax();
 return iRows*iCols;   
}    
// **************************************************************		
int CRData::DelVal(const int iC, const double lfS, const double lfE)
{// Delete those rows which have the values lfS-lfE in column iC

 CHECK_INDEX_RETURN(iC,nCols,-1)
 
 int i,icc;

 if(GetSteps()<=0)
   {PRINT_DEBUG("Illegal number of rows\n")
    return -1;
   }

 for(i=0,icc=0;  i<GetSteps(); i++)
    {if((*pT[iC])[i]>=lfS && (*pT[iC])[i]<=lfE)
       {Remove(i,FALSE); 
        icc++;
         i--;
       }
    }//for i

 if(icc==0)
   {PRINT_DEBUG("Warning: values %f:%f not found\n",lfS,lfE)
    return 1;
   }

 NewMinMax();
 return GetSteps()-icc;   
}    
// **************************************************************		
int CRData::LookTable(CRData *Table, const int iInsert, 
                                     const int iX, const int iY)
{
 // this->iCx: 1st x-column of DataFile (data values)
 // this->iCy: 1st y-column of DataFile (interpolated values)
 // Table->iCx: x-column of  TableFile (data values)
 // Table->iCy: y-column of  TableFile (table values)
 // iInsert: TRUE to be inserted after column iCy; else replaced

 if(GetSteps()<=0)
   {PRINT_DEBUG("Invalid data table\n")
    return 0;
   } 

 if(Table->GetSteps()==0)
   {PRINT_DEBUG("Invalid ppline table\n")
    return 0;
   } 

 SplineTab *SpT=new SplineTab(Table);
 if(!SpT || SpT->GetLength()<=0)
   {PRINT_DEBUG("Illegal spline table");
    return -1;
   }

 int ix = (iX==-1?iCx:iX);
 int iy = (iY==-1?iCy:iY);

 CHECK_INDEX_RETURN(ix,nCols,0)
 CHECK_INDEX_RETURN(iy,nCols,0)
 CHECK_INDEX_RETURN(Table->GetColX(),Table->GetCols(),0)
 CHECK_INDEX_RETURN(Table->GetColY(),Table->GetCols(),0)

 if(iInsert)
   {if(!InsertCol(iy))
      {PRINT_DEBUG("Out of memory\n");
       return 0;
      }
    if(ix>=iy)ix++; 
    }

 double lfT;
 for(int i=0; i<pT[iCy]->GetSteps();i++)
    {if(!SpT->Intpol((*pT[ix])[i],lfT))
       {PRINT_DEBUG("Spline value (%d: %f) out of table range\n",i,(*pT[ix])[i])
       }
     (*pT[(iInsert ? iy+1:iy)])[i]=lfT;
    }//for

  NewMinMax();
  delete SpT;
  return 1;
}
// ********************************************************************
int CRData::ZeroShift(MDATA xz, const int iX, const int iY)
{
 if(GetSteps()<=0)
   {PRINT_DEBUG("Invalid data\n")
    return 0;
   } 

 int ix = (iX==-1?iCx:iX);
 int iy = (iY==-1?iCy:iY);

 CHECK_INDEX_RETURN(ix,nCols,0)
 CHECK_INDEX_RETURN(iy,nCols,0)

 SortData(ix);

 MDATA y0;
 if(!LinIntpol(xz,y0,ix,iy))
   {PRINT_DEBUG("x-value %f out of range\n",xz)
    return 0;
   }
 pT[iy]->BasicCalc('-',y0);
 NewMinMax();
 return 1;
}
// ******************************************************
int CRData::Reduce(const int iR)
{int ir=1;
 int  i;
 if(GetSteps()<=0 || !nCols)
   {PRINT_DEBUG("Reduce not possible\n")
    return 0;
   }
	  
 for(i=0;i<nCols; i++){if( !(ir=pT[i]->Reduce(iR)) ) return 0;}
 return ir;
}
// **************************************************************
int CRData::MergeRows( CRData & C)
{int i,ii;
 if(GetSteps()<=0 || C.GetSteps()<=0 ||
    !nCols || !C.GetCols() || nCols != C.GetCols())
   {PRINT_DEBUG("Merge not possible\n")
    return 0;
   }

  for(i=0,ii=0; i<nCols; i++)
    {if(!pT[i]->Merge(&C[i]) )ii++;}//for i

 NewMinMax();
 fprintf(stderr,"Steps: %d Error:%d\n",GetSteps(),ii);
 return (ii ? 0:GetSteps());   
}
// **************************************************************
int CRData::PolyFit(Polynom &P, const int iOrder,const int iX, const int iY)
{ double  lfA[MAX_ORDER+3][MAX_ORDER+3];
  double lfxx[MAX_ORDER+3];
  double  lfXh,lfYh,lfFx,lfam,lfF;
  //double lfDm,lfDiff;
  int i,j,k,m,n;
  int m1,m2;
  int kk;

  int nPoints=GetSteps();
 if(nPoints<=3)
   {PRINT_DEBUG("Invalid data\n")
    return 0;
   } 

 int ix = (iX==-1?iCx:iX);
 int iy = (iY==-1?iCy:iY);

 CHECK_INDEX_RETURN(ix,nCols,0)
 CHECK_INDEX_RETURN(iy,nCols,0)

 if(iOrder <= 0 || iOrder > MAX_ORDER)return 0;
 SortData(ix);
 m = iOrder;

 double *lfKoff = new double [iOrder+1];
 CHECK_POINTER_RETURN(lfKoff,0)

 for(i=0;i<m;i++)lfKoff[i] = lfxx[i] = 0;
 lfxx[m] = 0;

 for(i=1,n=nPoints,m1=m+1,m2=m+2;i<=m1;i++)
        lfA[1][i] = lfA[i][m1] = lfA[i][m2] = 0;

  for(i=0;i<nPoints;i++)
     { for(j=1,lfFx=1,lfXh=(*pT[ix])[i],lfYh=(*pT[iy])[i];j<=m1;j++)
	  {lfA[1][j]+=lfFx; lfFx*=lfXh; lfA[j][m2]+=lfYh; lfYh*=lfXh;}
  	   for(j=2;j<=m1;j++){lfA[j][m1]+=lfFx; lfFx*=lfXh;}
	  }

  for(i=2;i<=m1;i++)
     {for(j=1;j<=m;j++) lfA[i][j]=lfA[i-1][j+1];}

  /*GAUSS*/
  for(k=1,n=m1;k<=n;k++)
     { lfam=fabs(lfA[k][k]);
       for(i=k,kk=k;i<=n;i++)
	  {if(lfam<fabs(lfA[i][k])){lfam=fabs(lfA[i][k]);kk=i;}
	  }
       if(lfam<1e-12)return -1;//failed

       if(kk!=k)
	 { for(j=k;j<=n+1;j++){SWAP (double,lfA[k][j],lfA[kk][j])}
	 }
       for(i=1;i<=n;i++)
	  { if(i!=k)
	      {if(lfA[k][k]==0)return -1; //failed
	       lfF=lfA[i][k]/lfA[k][k];
	       for(j=k+1;j<=n+1;j++)lfA[i][j]=lfA[i][j]-lfF*lfA[k][j];
	      }
	   }

      }// for(k)
  for(k=1;k<=n;k++)
     {if (lfA[k][k]==0)return -1; //failed
      lfxx[k]=lfA[k][n+1]/lfA[k][k];
     }
  /*GAUSS end*/
 for(i=0;i<=m;i++)lfKoff[i]=lfxx[i+1];
 Polynom p(iOrder,lfKoff);
 P=p;
 return 1;
}
//**********************************************************************


#if defined (TEST_CDATA)

main(void)
{CRData A(3,5,3.14159);
 CRData B(3,5,2);
 Formula F("C2=POW(C1)");
 //Formula G("C3=C1/C2");
 A.Print();
 A.Calculate(F);
 //A.Print();
 //A.Calculate(G);
 A.Print();
  
}

#endif
