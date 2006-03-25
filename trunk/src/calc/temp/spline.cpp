// File: spline.cpp
//$Id: spline.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
//$Log: spline.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.2  1999/03/15 09:08:37  herbie
//*** empty log message ***
//

#include <stdio.h>
#include <string.h>

#ifndef SPLINE_H
#include "spline.h"
#endif

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif 

#define SWAP(TYP,X,Y) \
	{TYP D=X;X=Y;Y=D;}

//#define SPLINE_DEBUG
#define SPLINE_DEBUGFILE "spline.bug"
#if defined(SPLINE_DEBUG)
FILE* fdd;
#endif

// ******************************************************
// class SplineTab Functions
// ******************************************************
// Splinetab
// ******************************************************

// ******************************************************
SplineTab::SplineTab(double *pX, double *pY, const int nP)
{ plfX=plfY=plfD2=NULL;

  if(nP>3)
    {if(pX!=NULL)
       {if( (plfX=new double[nP]) == NULL){nPoints=0;return;} 
       }

     if(pY!=NULL)
      {if( (plfY=new double[nP]) == NULL)
	 {nPoints=0; delete plfX; return;}
      }

    if( (plfD2=new double[nP]) == NULL)
	{nPoints=0; delete plfX; delete plfY; plfX=plfY=NULL; return;}
    for (int i=0;i<nP;i++){plfX[i]=pX[i];plfY[i]=pY[i];}

    nPoints=nP;
#if defined(SPLINE_DEBUG)
 fdd=fopen(SPLINE_DEBUGFILE,"a");
 if(fdd!=NULL)
  {fprintf(fdd,"New Object SPLINE (Spline): %p: X:%p Y:%p D:%p Size:%d\n",
           this,plfX,plfY,plfD2,nPoints);
   fclose(fdd);
  }
#endif

	SortData();
	Spline();
   }
}
//*********************************************************
SplineTab::SplineTab(CRData *TD)
{ plfX=plfY=plfD2=NULL;
  nPoints=0;
  int nP=TD->GetSteps();
  if(TD==0 || nP<=0)return;

  if( (plfX=new double[nP]) == NULL)return;

  if( (plfY=new double[nP]) == NULL) { delete plfX; return;}

  if( (plfD2=new double[nP]) == NULL)
		 { delete plfX; delete plfY; plfX=plfY=NULL; return;}

  int ix=TD->GetColX();
  int iy=TD->GetColY();
  CHECK_INDEX_RET(ix,TD->GetCols())
  CHECK_INDEX_RET(iy,TD->GetCols())

  nPoints=nP;
  
  for (int i=0;i<nPoints;i++)
      {plfX[i]=(*TD)[ix][i];
       plfY[i]=(*TD)[iy][i];
      }

#if defined(SPLINE_DEBUG)
 fdd=fopen(SPLINE_DEBUGFILE,"a");
 if(fdd!=NULL)
  {fprintf(fdd,"New Object SPLINE (Spline): %p: X:%p Y:%p D:%p Size:%d\n",
			  this,plfX,plfY,plfD2,nPoints);
	fclose(fdd);
  }
#endif

	SortData();
	Spline();
}
//*********************************************************
int SplineTab::ReadTable(const char *szFileName)
{// Reads complete data from file
 // return -2 : Error opening file
 //        -1 : Illegal Header
 //         0 : Out of memory
 //       > 3 : Ok # of data points

 FILE *f;

 if ((f=fopen(szFileName,"r")) == NULL)
	{nPoints=0;return(-2);}

 fgets(szHead,80,f);
 RemoveCR_LF(szHead);
 if(strstr(szHead,"SPLINETABLE")==NULL)return -1;
 fscanf(f,"%d",&nPoints);
 if(nPoints<3)return -1;

 if( (plfX=new double[nPoints]) == NULL)
     {nPoints=0;return 0;}

 if( (plfY=new double[nPoints]) == NULL)
	 {nPoints=0; delete plfX; return 0;}

 if( (plfD2=new double[nPoints]) == NULL)
	 {nPoints=0; delete plfX; delete plfY; return 0;}

#if defined(SPLINE_DEBUG)
 fdd=fopen(SPLINE_DEBUGFILE,"a");
 if(fdd!=NULL)
  {fprintf(fdd,"New Object SPLINE (ReadTable): %p: X:%p Y:%p D:%p Size:%d\n",
           this,plfX,plfY,plfD2,nPoints);
   fclose(fdd);
  }
#endif

 for(int i=0; i<nPoints; i++)
     fscanf(f,"%lf %lf", &(plfX[i]),&(plfY[i]));
 SortData();
 fclose(f);
 return nPoints;
}
// ******************************************************
void SplineTab::GetHeader(char *szH)
{
  strncpy(szH,szHead,80);
}
// ******************************************************
int SplineTab::SortData(void)
{ int iXchCount;
  int i1,i2;

  if (nPoints<3)return 0;// empty

  do
	 {iXchCount=0;

	  for(i1=0,i2=1;i2<nPoints;i2++)
	  {if(plfX[i1]>plfX[i2])
		     {SWAP(double,plfX[i1],plfX[i2])
			  SWAP(double,plfY[i1],plfY[i2])
		      iXchCount++;
		     }
	   i1=i2;
	  }
     }
  while (iXchCount != 0);

 return(1);
}
//**********************************************************************
void SplineTab::PrintTable(int bPrD2)
{
 printf("%s",szHead);

 for(int i=0; i<nPoints; i++)
	{if(!bPrD2)printf("%f  %f\n", plfX[i],plfY[i]);
	   else printf("%f  %f  %f\n", plfX[i],plfY[i],plfD2[i]);
	 if(i!=0 && i%18==0)
	   {printf("Type any key and <CR> to continue ...");
	    getchar();
	   }
    }
 printf("Table printed.\n");

}
// ******************************************************
int SplineTab::SaveTable(const char * szFileName)
{
 FILE *f;

 if(nPoints<3)return -1;

 if ((f=fopen(szFileName,"w")) == NULL)return -2;

 fprintf(f,"%s",szHead);
 fprintf(f,"%d\n",nPoints);

 for(int i=0; i<nPoints; i++)
	fprintf(f,"%f  %f  %f\n", plfX[i],plfY[i],plfD2[i]);
 fclose(f);
 return nPoints;
}
// ******************************************************
int SplineTab::Spline(float fD0, float fDn)
{
 double p,qn,sig,un;

 if(nPoints<3)return 0;

 double *u= new double[nPoints];
 if(u==0)return 0;

 if(fD0>0.99E30)plfD2[0]=u[0]=0;
 else
 {plfD2[0]=-0.5;
  u[0]=(3.0/(plfX[1]-plfX[0])) * ((plfY[1]-plfY[0])/(plfX[1]-plfX[0])-fD0);
 }

 for(int i=1;i<nPoints-1;i++)
    {sig=(plfX[i]-plfX[i-1])/(plfX[i+1]-plfX[i-1]);
     p=sig*plfD2[i-1]+2.0;
     plfD2[i]=(sig-1.0)/p;
     u[i]=(plfY[i+1]-plfY[i])/(plfX[i+1]-plfX[i]) -
                         (plfY[i]-plfY[i-1])/(plfX[i]-plfX[i-1]);
     u[i]=(6.0*u[i]/(plfX[i+1]-plfX[i-1])-sig*u[i-1])/p;
    }
 if(fDn > 0.99e30) qn=un=0.0;
 else
 {qn=0.5;
  un=(3.0/(plfX[nPoints-1]-plfX[nPoints-2])) *
   (fDn-plfY[nPoints-1]-plfY[nPoints-2])/(plfX[nPoints-1]-plfX[nPoints-2]);
 }
 plfD2[nPoints-1]=(un-qn*u[nPoints-2])/(qn*plfD2[nPoints-2]+1.0);

 for(int k=nPoints-2; k>0; k--)plfD2[k]=plfD2[k]*plfD2[k+1]+u[k];

#if defined(SPLINE_DEBUG)
 fdd=fopen(SPLINE_DEBUGFILE,"a");
 if(fdd!=NULL)
  {fprintf(fdd,"Splining Object SPLINE (Spline): %p: X:%p Y:%p D:%p Size:%d\n",
           this,plfX,plfY,plfD2,nPoints);
   fclose(fdd);
  }
#endif

 delete u;
 return 1;
}
// ******************************************************
SplineTab::~SplineTab()
{
#if defined(SPLINE_DEBUG)
 fdd=fopen(SPLINE_DEBUGFILE,"a");
 if(fdd!=NULL)
  {fprintf(fdd,"Deleting Object SPLINE (~Spline): %p: X:%p Y:%p D:%p Size:%d\n",
           this,plfX,plfY,plfD2,nPoints);
   fclose(fdd);
  }
#endif

   if(plfX!=NULL)delete plfX;
   if(plfY!=NULL)delete plfY;
   if(plfD2!=NULL)delete plfD2;
}
// ******************************************************
int SplineTab::Intpol(const double lfX,double &lfY)
{
 int iLo=0,k;
 int iHi=nPoints-1;

 if(nPoints<3)return 0;
 if(lfX>plfX[nPoints-1]){lfY=plfY[nPoints-1];return 0;}
 if(lfX<plfX[0]){lfY=plfY[0];return 0;}

 while(iHi-iLo>1)
	  {k=(iHi+iLo) >> 1;
	   if(plfX[k] > lfX)iHi=k;
	   else iLo=k;
	  }
 double h=plfX[iHi]-plfX[iLo];
 if(h==0.0)return 0;
 double a=(plfX[iHi]-lfX)/h;
 double b=(lfX-plfX[iLo])/h;
 lfY=a*plfY[iLo]+b*plfY[iHi]+
	 ((a*a*a-a)*plfD2[iLo]+(b*b*b-b)*plfD2[iHi])*(h*h)/6.0;
 return 1;
}
// ******************************************************
