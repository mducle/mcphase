// File: formula.cpp
// $Id: formula.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: formula.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// *** empty log message ***
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#ifndef FORMULA_H
#include "formula.h"
#endif

#ifndef STDINC_H
#include "stdinc.h"
#endif

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

#define FORMULA_MODULE

//#define FORMULA_TEST
//#define RANGE_TEST

// **************************************************
// class formula functions
// **************************************************

const struct SrcArg Formula::NoSrc={{-1},NoARG,0};
const char *Formula::szFun[]=
                  {"ABS","COS","EXP","LOG","SIN","TAN","SEQ",
                   "ACOS","ASIN","ATAN","SQRT"};
MDATA (*FFuns[]) (double x)=
                  {fabs, cos, exp, log10, sin, tan, 0 ,
		   acos, asin, atan, sqrt};

const int Formula::iMaxFun=11;
//
// Allowed formulas
// #: Columnumber : 1 - n
// c# = 1.23 <+-*/^> c#     1
// c# = c# <+-*/^> 1.234    2
// c# = c# <+-*/^> c#       3
// c# = FUN(c#)             4
// c# = FUN(1.234)          5
//
// FUN:  ABS, COS, EXP, LOG, SIN,  TAN, SEQ, ACOS, ASIN, ATAN, SQRT,
// iFun:   0    1   2    3    4    5    6    7     8     9     10  

Formula::Formula(const char *szNF)
{
  strncpy(szF,szNF,MAX_LINELENGTH);
  szF[MAX_LINELENGTH]=0;
  Strupr(szF);
  iDestCol=-1;

  int i;
  for(i=0;i<MAX_SRC;i++)Src[i]=NoSrc;
     
  pAct=&szF[0];

}
// **************************************************
// struct SrcArg Formula::GetSrc(const int i) const
// {if(i<0 || i>=MAX_SRC)return NoSrc;
//  else return Src[i];
// }
// **************************************************
enum ArgType Formula::GetType(const int i) const
{if(i<0 || i>=MAX_SRC)return NoARG;
 else return Src[i].AType;
}
// **************************************************
int Formula::GetSrc(const int i)  const
{if(i<0 || i>=MAX_SRC)return -1;
 else return Src[i].SType.iSrcCol;
}
// **************************************************
char Formula::GetOp(const int i)   const
{if(i<0 || i>=MAX_SRC)return 0;
 else return Src[i].cOp;
}
// **************************************************
int Formula::GetFunc(const int i) const
{if(i<0 || i>=MAX_SRC)return -1;
 else return Src[i].SType.iFunc;
}
// **************************************************
double Formula::GetVal(const int i) const
{if(i<0 || i>=MAX_SRC)return 0;
 else return Src[i].SType.lfVal;
}
// **************************************************
int Formula::IsFNumber(double &lfV)

{
 char *cpStart, szNumStr[80];
 int bDecimal, iLen, iNumLen;

 while (*pAct == ' ')pAct++;

 if (*pAct == 0)return 0;

 if (strchr("+-0123456789.", *pAct))
	 {cpStart = pAct;
	  iLen = 0;
	 bDecimal = 0;
	  while ((isdigit(*pAct)) ||
		((*pAct == '.') && (!bDecimal)) || (*pAct=='+') || (*pAct=='-'))
		{ if (*pAct == '.')bDecimal = 1;
		  pAct++;
		  iLen++;
		 }
	  if ((iLen == 1) && (cpStart[0] == '.'))return 0;
     if (*pAct == 'E')
	{pAct++;
	 iLen++;
	 if (strchr("+-", *pAct) != NULL){pAct++;iLen++;}
	 iNumLen = 0;
	 while ((isdigit(*pAct)) && (++iNumLen <= 3)){pAct++;iLen++;}
	}
	  strncpy(szNumStr, cpStart, iLen);
      szNumStr[iLen] = 0;
      lfV = atof(szNumStr);
      if (errno == ERANGE){errno=0;return 0;}
      return 1;
    }
return 0;
}
// **********************************************
int Formula::IsColumn(char *s)
{//returns Column number; negative if invalid 
 if(*s!='C')return -1; // no valid column
 s++;
 if(!IsNumString(s,UINT_NUM))return 0;
 return atoi(s);
}
// ***********************************************
int Formula::SetSource(char *s, const int iArg)
{if(*s=='C') //column follows
   {Src[iArg].SType.iSrcCol=IsColumn(s);
    if(Src[iArg].SType.iSrcCol<=0)return 0; // Illegal Col#
    Src[iArg].AType=ColARG;    
    //if(*s)return 0; // Illegal formula;
   }// column
   else  // number arg in function
    {if(!IsFNumber(Src[iArg].SType.lfVal))return 0; // illegal Function Arg
     Src[iArg].AType=ValueARG;
    }
 return 1;
}       
// *************************************************
 int Formula::Analyze(void)
{
 // return 0 if syntax OK
  
  char sc; //Save character at running pointer
  char *psc; // running pointer
 //char *pAct Actual string to analyze

  int i;

  pAct=szF;  
  if( (psc=strchr(pAct,'='))==NULL)
      { PRINT_DEBUG("Formula syntax error near [1]\n")
        return 1; // NO = char
      }	
  if( strchr( (psc+1),'='))
      { PRINT_DEBUG("Formula syntax error near [%d]\n",(int)(psc-szF))
       return psc-szF; // two = char
      }
  sc=*psc;
  *psc=0;
  
  iDestCol=IsColumn(pAct);

  if(iDestCol<=0)
        { PRINT_DEBUG("Formula syntax error near [%d]\n",(int)(psc-szF))
          return psc-szF;     // Illegal Col#
        }
 
  *psc=sc;
  
  pAct=psc+1;
  
  if( (psc=strchr(pAct,'(') ) )   //Analyze 1st after = for function
    {if(!strchr(pAct,')'))
        { PRINT_DEBUG("Formula syntax error near [%d]\n",(int)(psc-szF))
          return psc-szF;     // Unbalanced ()
        }
    
     //Its a Function: FUN(c1) or FUN(1.2)
     sc=*psc;
     *psc=0;
     for(i=0;i<iMaxFun;i++)
        {if(strstr(szFun[i],pAct))
           {Src[0].SType.iFunc=i;
            Src[0].AType=FuncARG;
            break;
           }// if
        }//for
     if(i==iMaxFun)
        { PRINT_DEBUG("Formula syntax error near [%d]\n",(int)(psc-szF))
          return psc-szF;     // Illegal function
        }

     *psc=sc;
     pAct=psc+1;

     // a number or c# has to follow
     psc=strchr(pAct,')');
     sc=*psc;
     *psc=0;

     i=SetSource(pAct,1);  
     *psc=sc;
     pAct=psc+1;
     if(!i)
       { PRINT_DEBUG("Formula syntax error near [%d]\n",(int)(psc-szF))
          return psc-szF;     //Invalid source arg
       }

     fprintf(stdout,"remain:%s\n",pAct);
     return 0;
    }//if() analyze for function

  psc=pAct+1;

  while(*psc)
    {switch(*psc)
           {case '+':
            case '-':
            case '*':
            case '/':
            case '^':
            case '$':Src[0].cOp=*psc; break;
           }
     if(Src[0].cOp!=0)break;
     psc++;
     }

  if(Src[0].cOp==0)
    { PRINT_DEBUG("Formula syntax error near [%d]\n",(int)(psc-szF))
      return psc-szF;     // Illegal Op char
    }
 
  sc=*psc;
  *psc=0;
  // a number or c# has to follow
  //Ananalyze 1st part after = for Column   
  i=SetSource(pAct,0);  
  if(!i)
    { PRINT_DEBUG("Formula syntax error near [%d]\n",(int)(psc-szF))
      return psc-szF;     // Illegal Col#
    }

  *psc=sc;
  pAct=psc+1;
  
  i=SetSource(pAct,1);  
	         
  return 0; 
}
// ***************************************************
int Formula::CheckArg(void) const
{ 
 if(iDestCol<=0)return 0;
 for(int i=0; i<MAX_SRC; i++)
    {if(Src[i].AType==NoARG)return 0;
     if(Src[i].AType!=ValueARG)
       {if(Src[i].SType.iSrcCol<0)return 0;}
    }
 return 1;
}
// ***************************************************
int Formula::IsFunc(void) const
{ 
 for(int i=0; i<MAX_SRC; i++)
    {if(Src[i].AType==FuncARG)return 1;
    }
 return 0;
}
// ***************************************************
void Formula::Print(FILE *f)
{fprintf(f,"Formula: %s\n",szF);
 fprintf(f,"Destination Column: %d\n",iDestCol);
 int i;
 for(i=0;i<MAX_SRC;i++)
    fprintf(f,"Source Argument %d: Func/SrcCol: %d, Value: %f, Type: %d, Op: %c(%x)\n",
                   i+1,Src[i].SType.iSrcCol, Src[i].SType.lfVal, (int)Src[i].AType,
                   Src[i].cOp,Src[i].cOp);
 fprintf(f,"\n");
}
// ***************************************************
// class Range functions
// ***************************************************
Range::Range(const char *S, const enum SqType T)
{R.St.i=R.En.i=-1;
 R.St.lf=R.En.lf=0;
 Type=NO_SQ;
 nRanges=-1;
 pR=0;
 szExp[0]=0;
 if(!S || S[0]==0)return;

 strncpy(szExp,S,MAX_LINELENGTH);
 szExp[MAX_LINELENGTH]=0;
 
 if(T==FLOAT_SQ) // szExp: "1.23 : 345"
   {
    if(GetFloatRange(szExp,R.En.lf,R.En.lf))return;
    nRanges=0;
    Type=FLOAT_SQ;    
   }

 if(T==INT_SQ)
   {char *pc=strchr(szExp,':');
    char szB[MAX_LINELENGTH+1];
    if(pc){if(strchr(++pc,':'))return;} // more than one ":"
  
    int i=0,id=0;
      while(1)
      { if(!StrTok(szExp,',',i,szB))break;
        fprintf(stderr,"Item %d: %s (%d)\n",i,szB,id);
        if(strchr(szB,':'))id++;
        i++;
      } 
    fprintf(stderr,"%d (%d) Items\n",i-id,i);
     if(i>0)
       {pR=new int[i-id];
        CHECK_POINTER_RET(pR)
       }     
 
     strncpy(szB,szExp,MAX_LINELENGTH);
     szB[MAX_LINELENGTH]=0;           
     pc=szB;
 
     int j;
     double a,b;
   
     for(j=0,id=0; j<i; j++)
        {if(!StrTok(szExp,',',j,szB))return;
         if(strchr(szB,':'))
           {if(GetFloatRange(szB,a,b))return;
            R.St.i=(int)a;
            R.En.i=(int)b;
	    id++;
           }
         else
           {if(!IsNumString(szB))return;
	    pR[j-id]=strtol(szB,NULL,10);}
          fprintf(stderr,"Analyse Item: %s\n",szB);
        }
     Type=INT_SQ;
     nRanges=i-id;
     return;
    }
 	
}	
// ***************************************************
int Range::GetRange(const int i) const
{CHECK_INDEX_RETURN(i,nRanges,0);
 return nRanges;
}
// ***************************************************
void Range::Print(void) const
{ int i; 
 fprintf(stderr,"Range expression: %s\n",szExp);
 switch(Type)
       {case INT_SQ:
            fprintf(stderr,"Type: INT (%d)\n",(int)Type);
            fprintf(stderr,"Start: %i\n",R.St.i);
            fprintf(stderr,"  End: %i\n",R.En.i);
            fprintf(stderr,"nRanges: %d\n",nRanges);
            for(i=0; i<nRanges; i++)
                fprintf(stderr,"Range %d: %d\n",i,pR[i]);    
            break;
         case FLOAT_SQ:
            fprintf(stderr,"Type: FLOAT (%d)\n",(int)Type);
            fprintf(stderr,"Start: %f\n",R.St.lf);
            fprintf(stderr,"  End: %f\n",R.En.lf);
            fprintf(stderr,"nRanges: %d\n",nRanges);
            break;
         case NO_SQ:
            fprintf(stderr,"Type: ERORR (%d)\n",(int)Type);
            fprintf(stderr,"Start: %i\n",R.St.i);
            fprintf(stderr,"  End: %i\n",R.En.i);
            fprintf(stderr,"nRanges: %d\n",nRanges);
            break;
        }
}
// ***************************************************
// ***************************************************

#if defined (FORMULA_TEST)
main(int iArgC, char ** szArgV)
{

 if(iArgC!=2)
   {fprintf(stderr,"Usage: %s formula\n",szArgV[0]);
    exit(EXIT_FAILURE);
   }
  DelSpace(szArgV[1]);
    
 fprintf(stdout,"Testing Formula: %s\n",szArgV[1]);

 Formula F(szArgV[1]);
 int i=F.Analyze();
 fprintf(stderr," Analyze result :%d\n",i); 
 F.Print();
 exit(EXIT_SUCCESS);
}
#endif 
// **********************************************
#if defined (RANGE_TEST)
main(int iArgC, char ** szArgV)
{

 if(iArgC!=2)
   {fprintf(stderr,"Usage: %s range-exp\n",szArgV[0]);
    exit(EXIT_FAILURE);
   }
   
 fprintf(stdout,"Testing range: %s\n",szArgV[1]);

 Range R(szArgV[1],INT_SQ);
 int i=R.IsValid();
 fprintf(stderr," Analyze result :%d\n",i); 
 R.Print();
 exit(EXIT_SUCCESS);
}
#endif 
