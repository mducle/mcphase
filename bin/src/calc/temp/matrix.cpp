//File: matrix.cpp
// $Id: matrix.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: matrix.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <bits/nan.h>

#include "stdinc.h"
#include "stdfunc.h"
#include "matrix.h"

#define TEST_MATRIX

extern int errno;
 MDATA Null = 0;
 Matrix  MNull;
// ********************************************************
// ******************************************************
// class Matrix Functions
// ******************************************************
// Matrix
// ******************************************************
// ******************************************************

// ******************************************************
Matrix::Matrix(void)
{
    pT = NULL;
 nRows = 0;
 nCols = 0;
 indx  = 0;
}
// ******************************************************
Matrix::Matrix(const unsigned iRows, const unsigned iCols, const MDATA v)
{// Error: nSteps==0  Not enough memory
      pT = NULL;
    indx = 0;
   nRows = 0;
   nCols = 0;
  if(!Alloc(iRows,iCols,v))return;
  nRows = iRows;
  nCols = iCols;
}
// ******************************************************
// Matrix::Matrix(MDATA **p, const unsigned iRows, const unsigned iCols)
// {     pT = NULL;
//    nRows = 0;
//    nCols = 0;
//   if(!Alloc(iRows,iCols))return;
//   nRows = iRows;
//   nCols = iCols;
//   CopyData(p);
// }
// ******************************************************
Matrix::Matrix(double *p, const unsigned iRows, const unsigned iCols)
{     pT = NULL;
    indx = 0;
   nRows = 0;
   nCols = 0;
  if(!Alloc(iRows,iCols))return;
  nRows = iRows;
  nCols = iCols;
  unsigned i,j;
  for(i=0;i<iRows;i++)
     for(j=0;j<iCols;j++)
        pT[i][j]=*(p++);
}
// ******************************************************
Matrix::Matrix(const Matrix &N)
{ nRows = 0;
  nCols = 0;
     pT = 0;
   indx = 0;
  if(!Alloc(N))return;
  nRows=N.nRows;
  nCols=N.nCols;
  CopyData(N);         
}
// ******************************************************
int Matrix::CopyData(const Matrix &N)
{ return CopyData(N.pT);}
// ******************************************************
int Matrix::CopyData(double **N)
{unsigned i,j;
 for (i=0;i<nRows;i++)
       for(j=0;j<nCols;j++)pT[i][j]=N[i][j];
 indx=0;
 return (int)nRows;
}
// ******************************************************
int Matrix::Alloc(const Matrix & N)
{		
    pT = new MDATA * [N.nRows];
  CHECK_POINTER_RETURN(pT,0);

  for(unsigned i=0; i<N.nRows; i++)
     {pT[i] = new MDATA[N.nCols];
      CHECK_POINTER_RETURN(pT[i],0);
     }
 return 1;
}
// ********************************************************
int Matrix::Alloc(const unsigned iRows, const unsigned iCols, MDATA v)
{ if(!iCols || !iRows)
    {PRINT_DEBUG("Can not allocate %d columns or %d rows\n",iCols,iRows);
     return 0;
    }
  pT = new MDATA * [iRows];
  CHECK_POINTER_RETURN(pT,0);

  for(unsigned i=0; i<iRows; i++)
     {pT[i]=new MDATA [iCols];
      CHECK_POINTER_RETURN(pT[i],0);
      unsigned j;
      for(j=0;j<iCols;j++)pT[i][j]=v;
     }

  return 1;
}
// ******************************************************
void Matrix::Dealloc(void)
{
  if(pT)for(unsigned i=0; i<nRows; i++)delete pT[i];

  delete pT;
  nCols = 0;
  nRows = 0;
  delete indx;
}
// *******************************************************
MDATA ** Matrix::ReallocColData(const unsigned nR, const unsigned nC)
{MDATA **pD = 0;
 if(nC && nR)
  {pD = new MDATA * [nC];
   CHECK_POINTER_RETURN(pD,0)
   unsigned i;
   for(i=0;i<nR;i++)
      {pD[i] = new MDATA[nC];
       CHECK_POINTER_RETURN(pD[i],0)
      }//for
  }
 return pD;
 //nCols,nRows is NOT set to new values!!
}
// **************************************************
Matrix::~Matrix(void)
{Dealloc();}
// ******************************************************
// Operators
// ******************************************************   
const MDATA * Matrix::operator [](const unsigned i)
{if( !pT || !nCols || !nRows)
   {PRINT_DEBUG("Illegal (zero) data pointer\n")
    return NULL;
   }
 CHECK_INDEX_RETURN(i,nCols,NULL)
 return pT[i];		 
} 
// ******************************************************   
MDATA & Matrix::operator ()(const unsigned iRow, const unsigned iCol)
{if( !pT || !nCols || !nRows)
   {PRINT_DEBUG("Illegal (zero) data pointer\n")
    return Null;
   }
 CHECK_INDEX_RETURN(iCol,nCols,Null)
 CHECK_INDEX_RETURN(iRow,nRows,Null)
 return pT[iRow][iCol];		 
}
// ******************************************************   
Matrix &  Matrix::operator=(const Matrix &N)
{ if(this != &N)
    {if(nCols && pT)Dealloc();
     if(!Alloc(N)){pT=0; nCols=0; nRows=0; return *this;}
     nRows=N.nRows;
     nCols=N.nCols;
     CopyData(N);
    }
  return *this;
}
// ******************************************************   
// ********************************************************
// Other Functions
// ******************************************************
void Matrix::Print(String &S)
{
 S.Setf("nRows: %d\n",nRows);
 S.Addf("nCols: %d\n",nCols);

 if(!nCols || !nRows || !pT){S.Addf("Matrix: Null\n"); return;}

 unsigned i,j;
 for(i=0;i<nRows;i++)
    {for(j=0;j<nCols;j++)S.Addf("(%u,%u):%f ", i,j,pT[i][j]);
     S.Add("\n");
    }
 S.Add("\n\n");
} 
// ******************************************************
// CData * CRData::At(const int iCol)
// {if( !pT || !nCols)
//    {PRINT_DEBUG("Illegal (zero) data pointer pCol\n")
//     return 0;
//    }
//  CHECK_INDEX_RETURN(iCol,nCols,0)
//  return pT[iCol];		 
// } 
// // ******************************************************   
// int CRData::ReadCols(FILE *InStream, const int iSteps, const int nC,
//                     const int iX,const  int iY,const  int iZ)
// {// Reads complete column data from file
//  char szB[MAX_LINELENGTH+1],*s;
// 
//  if(pT && nCols>0)Dealloc();
// 
//  if(!Alloc(nC,iSteps))return 0;
//  nCols = nC;
// 
//  int i,j;
// 
//  for(i=0; i<iSteps; i++)
//    {if(fgets(szB,MAX_LINELENGTH,InStream)==NULL)break;
//     szB[MAX_LINELENGTH]=0;
//     while( (s=strchr(szB,'D'))!=NULL)*s='E';
//     for(j=0;j<nC;j++)
//        {if(j==0)s=strtok(szB," ");
// 	else s=strtok(NULL," ");
// 	if(s==NULL)break;
// 	(*this)[j][i]=atof(s);
//        }//for j
//     }//for i
//  if(i<iSteps)return 0; //return i;
// 
//  for(i=0; i<nCols; i++)pT[i]->NewMinMax();
// 
//  iCx = (iX >= 0 && iX<nCols? iX:0);
//  iCy = (iY >= 0 && iY<nCols? iY:0);
//  iCz = (iZ >= 0 && iZ<nCols? iZ:0);
// 
// return iSteps;
// }
// // ***********************************************************
// int CRData::ReadCols(LineString &LS, const int iSteps, const int nC,
//                      const int iX,const  int iY,const  int iZ)
// {// Reads complete column data from LineString
//  // return 0: inittialisation error
//  //        n>0: number of readed value lines
//  //        n<0: error during read -n gives line
// 
//  if(LS.GetNLines()!=iSteps)return 0;
//  if(pT && nCols>0)Dealloc();
//  if(!Alloc(nC,iSteps))return 0;
//  nCols = nC;
// 
//  int i,j;
// 
//  for(i=0; i<iSteps; i++)
//    {String L(LS.Line(i));
//     //fprintf(stderr,"Line: %d >%s\n",i,(const char *)L);
//     for(j=0;j<nC;j++)
//        {int iR=GetNumberFromString(L,(*this)[j][i],(j==0 ?nC:0),(j==0?' ':0));
//         if(j==nC-1){if(iR!=0)return -i;}
//         else if(iR<=0)return -i;        
//        }//for j
//     }//for i
//  if(i<iSteps)return 0; //return i;
// 
//  for(i=0; i<nCols; i++)pT[i]->NewMinMax();
// 
//  iCx = (iX >= 0 && iX<nCols? iX:0);
//  iCy = (iY >= 0 && iY<nCols? iY:0);
//  iCz = (iZ >= 0 && iZ<nCols? iZ:0);
// 
// return iSteps;
// }

//**********************************************************************
Matrix operator+(Matrix &N1,Matrix &N2)
{if(N1.nCols!=N2.nCols || N1.nRows!=N2.nRows || !N1.pT || !N2.pT)
   {PRINT_DEBUG("Illegal (zero) data pointer or matrices dont match\n")
    return N1;
   }
 Matrix M(N1.nCols,N1.nRows);
 unsigned i,j;

 for(i=0; i<N1.nRows; i++)
    for(j=0;j<N1.nCols;j++)M(i,j)=N1(i,j)+N2(i,j);
     
 return M;
}
// ******************************************************   
Matrix operator-(Matrix &N1,Matrix &N2)
{if(N1.nCols!=N2.nCols || N1.nRows!=N2.nRows || !N1.pT || !N2.pT)
   {PRINT_DEBUG("Illegal (zero) data pointer or matrices dont match\n")
    return N1;
   }
 Matrix M(N1.nCols,N1.nRows);
 unsigned i,j;

 for(i=0; i<N1.nRows; i++)
    for(j=0;j<N1.nCols;j++)M(i,j)=N1(i,j)-N2(i,j);
     
 return M;

}
// ******************************************************   
Matrix operator*(Matrix &N1,Matrix &N2)
{
if(N1.nCols!=N2.nRows || !N1.pT || !N2.pT)
   {PRINT_DEBUG("Illegal (zero) data pointer or matrices dont match\n")
    return N1;
   }
 Matrix M(N1.nRows,N2.nCols);
 unsigned i,j,k;

 for(i=0; i<N1.nRows; i++)
    {for(k=0;k<N2.nCols;k++)
       {for(j=0;j<N1.nCols;j++)
	  {M(i,k)+=N1(i,j)*N2(j,k);}
       }
    } 
 
 return M;

}
// ******************************************************   
Matrix operator*(const MDATA N, Matrix &N2)
{if(!N2.nCols || !N2.nRows || !N2.pT)
   {PRINT_DEBUG("Illegal (zero) data pointer\n")
    return MNull;
   }
 Matrix M(N2.nCols,N2.nRows);
 unsigned i,j;

 for(i=0; i<N2.nRows; i++)
    for(j=0;j<N2.nCols;j++)M(i,j)=N2(i,j)*N;
     
 return M;
}
// ******************************************************
int Matrix::LUDcmp(double *d)
{
 // LU decomposition ALGORITHM FROM 
 // W.H. Press, et al "Numerical Recipies in C" Cambridge University Press 1992
if(nCols==0 || nRows==0 || pT==0)
   {PRINT_DEBUG("Illegal (zero) data pointer\n")
    return 0;
   }
if(nCols!=nRows)
   {PRINT_DEBUG("Cannot perform LU decomosition of non square matrix\n")
    return 0;
   }

    int i,imax,j,k;
 double big,dum,sum,temp;
 const double tiny=1.0e-20;

 double *vv=new double [nRows];
 CHECK_POINTER_RETURN(vv,0)

 delete indx;
 indx=new int[nRows];
 CHECK_POINTER_RETURN(indx,0)

*d=1.0;
for(i=0;i<(int)nRows;i++)
   {big=0.0;
    for(j=0;j<(int)nRows;j++)
       if((temp=fabs(pT[i][j]))>big)big=temp;
    if(big==0.0){PRINT_DEBUG("Singular matrix\n")
                 return 0;
                }
    vv[i]=1.0/big;
   }
for(j=0;j<(int)nRows;j++)
   {for(i=0;i<j;i++)
       {sum=pT[i][j];
        for(k=0;k<i;k++)
             sum-=pT[i][k]*pT[k][j];
        pT[i][j]=sum;
       }
    big=0.0;
    for(i=j;i<(int)nRows;i++)
       {sum=pT[i][j];
        for(k=0;k<j;k++)
                sum-=pT[i][k]*pT[k][j];
        pT[i][j]=sum;
        if((dum=vv[i]*fabs(sum))>=big)
          {big=dum; imax=i;}
       
       }
    if(j!=imax)
      {for(k=0;k<(int)nRows;k++)
          {dum=pT[imax][k];
           pT[imax][k]=pT[j][k];
           pT[j][k]=dum;
          }
        *d=-(*d);
        vv[imax]=vv[j];
      }

    indx[j]=imax;
    if(pT[j][j]==0.0)pT[j][j]=tiny;
    if(j!=(int)nRows-1)
      {dum=1.0/(pT[j][j]);
       for(i=j+1;i<(int)nRows;i++)pT[i][j]*=dum;
      }
   }
 delete vv;
 return 1;
}
// ******************************************************
int  Matrix::LUBakSub(double *b)
{// *b: Input right-hand side vector B
 // LU back substitution ALGORITHM FROM 
 // W.H. Press, et al "Numerical Recipies in C" Cambridge University Press 1992
 // indx must be set from LUDcmp
if(nCols==0 || nRows==0 || pT==0)
   {PRINT_DEBUG("Illegal (zero) data pointer\n")
    return 0;
   }
if(nCols!=nRows)
   {PRINT_DEBUG("Cannot perform LU decomosition of non square matrix\n")
    return 0;
   }
if(!indx){PRINT_DEBUG("LUDcmp not performed\n")
          return 0;
         }
int i,ii=-1,ip,j;
double sum;

for(i=0;i<(int)nRows;i++)
   {ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii!=-1)
      {for(j=ii;j<=i-1;j++)
            sum-=pT[i][j]*b[j];
      }
    else if(sum)ii=i;
    b[i]=sum;
   }
for(i=nRows-1;i>=0;i--)
   {sum=b[i];
    for(j=i+1;j<(int)nRows;j++)sum-=pT[i][j]*b[j];
    b[i]=sum/pT[i][i];
   }
return 1;
}
// ******************************************************
MDATA Matrix::Det(void)
{double d;
 if(LUDcmp(&d)==0)return NAN;
 int i;
 for(i=0;i<(int)nRows;i++)d*=pT[i][i];
 return d;
}
// ******************************************************
int Matrix::Inv(void)
{double d;
 Matrix H(nRows,nRows);

 if(LUDcmp(&d)==0)return 0;
 unsigned j,i;

 double *col=new double[nRows];
 CHECK_POINTER_RETURN(col,0)
 
 for(j=0;j<nRows;j++)
    {for(i=0;i<nRows;i++)col[i]=0.0;
     col[j]=1;
     LUBakSub(col);
     for(i=0;i<nRows;i++)H(i,j)=col[i];
    }
 *this=H;
 delete col;
 return 1;
}
// ******************************************************
int Matrix::Transp(void)
{ 
 if(nCols==0 || nRows==0 || pT==0)
   {PRINT_DEBUG("Illegal (zero) data pointer\n")
    return 0;
   }
 Matrix H(nCols,nRows);
 unsigned i,j;
 for(i=0;i<nRows;i++)
    {for(j=0;j<nCols;j++)H(j,i)=pT[i][j];
    }
 *this=H;
 return 1;
}
// ******************************************************
#if defined (TEST_MATRIX)

main(void)
{double c[3][3]={{3,2,1},{1,0,2},{4,1,3}};
 double a[3][3]={{1,0,0},{0.75,1,0},{0.25,-0.2,1}};
 double b[3][3]={{4,1,3},{0,1.25,-1.25,},{0,0,1}};
 //double c[4][4]={{2,9,9,4},{2,-3,12,8},{4,8,3,-5},{1,2,6,4}};
 Matrix A((double *)&(a[0]),3,3);
 Matrix B((double *)&(b[0]),3,3);
 Matrix C((double *)&(c[0]),3,3);

 String S;
 A.Print(S);
 fprintf(stdout,"A: %s",(const char *)S);
 B.Print(S);
 fprintf(stdout,"B: %s",(const char *)S);
 C.Print(S);
 fprintf(stdout,"C: %s",(const char *)S);

//  Matrix D(A);
//  D.Print(S);
//  fprintf(stdout,"D: %s",(const char *)S);
 Matrix D;
 D=A*B;
 D.Print(S);
 fprintf(stdout,"D=AB: %s",(const char *)S);

 //fprintf(stdout,"Det=%f\n",B.Det()); 
 C.Inv();
 C.Print(S);
 fprintf(stdout,"inv(C) %s",(const char *)S);

 C.Transp();  
 C.Print(S);
 fprintf(stdout,"transp(C) %s",(const char *)S);
}
#endif
