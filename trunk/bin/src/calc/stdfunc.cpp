//File: stdfunc.cpp
// $Id: stdfunc.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: stdfunc.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.3  1999/02/19 09:32:07  herbie
// Added GetArraysFrom...
//
// Revision 1.2  1999/02/17 08:48:03  herbie
// Changed MAX_LINELENGTH to 1024.
//

#define STDFUNC_CPP

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>

#include "stdfunc.h"
#include "stdinc.h"
#include "simplex.h"

// ************************************************
int CHECK_INDEX(const int i, const int iMax)
{if(i<0 || i>=iMax)
   {fprintf(stderr,"Index %d out of range (0 ...%d)\n",i,iMax-1);
    return 0;
   }
 return 1;
} 
// ************************************************
int CheckIndex(const int i, const int iMax)
{if(i<0 || i>=iMax)
   {fprintf(stderr,"Index %d out of range (0 ...%d)\n",i,iMax-1);
    return (i<0?0:iMax-1);
   }
 return i;
} 
// **************************************************
void Suspend(const double T)
{clock_t start,end;
 
 start=clock();
 while(1)
      {end=clock();
       if(end==-1)break;
       if( (double) (end-start)/CLOCKS_PER_SEC >= T)break;
      }
}
// ************************************************
void Strupr(char *s)
{ if(!s)return;
  while (*s){*s=toupper(*s); s++;}
}
// ************************************************
void RemoveCR_LF(char *szB)
{
  //   Replaces \r\n combinations at first occurence by \0
  //
  if (!szB)return;    
  char *s;
  if( (s=strchr(szB,'\n')) )*s=0;
  if( (s=strchr(szB,'\r')) )*s=0;
}
// ************************************************
int DelSpace(char *szB)
{// deletes space in 0-terminated string
  int iC = 0; 
 char *p = strchr(szB,' ');
 
 while(p)	
  {char *src=p+1;
   while(*p)
    {*p=*src; p++; src++;}
   p=strchr(szB,' ');
   iC++;
  }
 return iC;       	
}		
// ************************************************
int  ReplaceChar(char *szB, const int What, const int By)
{
  //   Replaces every What by By
 if(!szB)return 0;
  int  i = 0;
 char *s;
 while( (s=strchr(szB,What))!=NULL){*s=By;i++;}
 return i;
}
// ************************************************
int CutRNum(char *szB)
{if(!szB)return 0;
 char *p=szB;
 while(*p==' ')p++;
 while(*p)
      {if(strchr("0123456789+-.Ee",*p)==0)
         {*p=0; break;}
       p++;
      }
 return p-szB;
}
// ************************************************
int ReplaceLastChar(char *szB, const int What, const int By)
{if(!szB)return 0;
 char *p=strrchr(szB,What);
 if(!p)return 0;
 *p=By;
 return p-szB;
}
// ************************************************
int  IsNumString(const char *s, const int NumFlag)
{
  //   returns 0 if string contains non-numerical chars
  //   numerical are +-0123456789Ee.
  char szT[]="0123456789+-.Ee";

  if(NumFlag==UINT_NUM)szT[10]=0;
  if(NumFlag==INT_NUM)szT[12]=0; 

  if(!s || *s==0)return 0;

  while(*s)
    {if(*s==' '){s++;continue;}
    if(strchr(szT,*s)==NULL)return 0;
    s++;
   }
  return 1;
}
// ***********************************************
int  CountCols(char *szB, const char *szSeparator)
{
  // Counts the number of columns
 if(!szB || !szSeparator || szB[0]==0)return 0;
 char    *s = strtok(szB,szSeparator);
  int iColC = 1;
 if(s==NULL)return 0;
 while(strtok(NULL,szSeparator)!=NULL)iColC++;
 return iColC;
}
// ********************************************************************
int GetCurrentTime(char *szB)
{     time_t ct = time(NULL);
  struct tm *lt = localtime(&ct);
  
 strcpy(szB,asctime(lt));
 return strlen(szB);
}
// ********************************************************************
int ShowFile(const char *szFile, FILE * fWhere)
{if(!(fWhere==stderr || fWhere==stdout) )
   {PRINT_DEBUG("Illegal file pointer\n")
    return 0;
   }
 FILE *fIn=fopen(szFile,"r");
 CHECK_FILE_POINTER_RETURN(fIn,szFile,0)
 
 char szB[MAX_LINELENGTH+1];
  int iC = 0;
  
 while(1)
      {if(!fgets(szB,MAX_LINELENGTH,fIn))break;
       szB[MAX_LINELENGTH]=0;
       fprintf(fWhere,"%s",szB);
       iC++;
      }    
 fclose(fIn);
 return iC;
}
// ********************************************************************
void Usage(const char *szPath, const char *szFile, FILE * fWhere)
{String B;
 char *p=strrchr(szFile,'/');
 const char *s=(p ? p+1 : szFile);
 if(szPath)
  {B.Set(szPath);
   if(szPath[strlen(szPath)-1]!='/')B.Add("/");
   B.Add(s);
  }
 else B.Set(s);	 
 B.Add(".man");
fprintf(fWhere,szPath);

// ShowFile(B.GetBuf(),fWhere);
 //fprintf(stderr,"\a");
 exit(EXIT_FAILURE);
}
// ********************************************************************  
void Usage(const char *szPath, const char *szCmd, String &Out)
{String B;
 char *p=strrchr(szCmd,'/');
 const char *s=(p ? p+1 : szCmd);
 if(szPath)
  {B.Set(szPath);
   if(szPath[strlen(szPath)-1]!='/')B.Add("/");
   B.Add(s);
  }
 else B.Set(s);	 
 B.Add(".man");

 FILE *fIn=fopen(B.GetBuf(),"r");
 CHECK_FILE_POINTER_RET(fIn,B.GetBuf())
 
 char szB[MAX_LINELENGTH+1];
 while(1)
      {if(!fgets(szB,MAX_LINELENGTH,fIn))break;
       Out.Add(szB);
      }
 fclose(fIn);
}
// ********************************************************************
void PrintRCS(const char *szID, const char *szPath, const char *szFile, FILE * fWhere)
{char szB[MAX_LINELENGTH+1];
 if(!AppendPath(szPath,szFile,szB,MAXPATH))return;
 strncat(szB,".man",MAX_LINELENGTH);
 szB[MAX_LINELENGTH]=0;
 FILE *fIn=fopen(szB,"r");
 CHECK_FILE_POINTER_RET(fIn,szB)
 
 int bFound=0;
 while(!bFound)
      {if(!fgets(szB,MAX_LINELENGTH,fIn))break;
       if(strstr(szB,"$Id:"))bFound=1;
      }
 fprintf(fWhere,"Man-ID: %s",(bFound ? szB : "not found\n"));
 fprintf(fWhere,"Cmd-ID: %s\n",szID);
 fclose(fIn); 
}
// ********************************************************************
unsigned AppendPath(const char *szP, const char *szF, char * szR, const unsigned MaxL)
{//Apppends szP (Path), szF (File) to szR -> complete Path/Filename
 if(MaxL==0)return 0;
 String B;
 char *p=strrchr(szF,'/');
 const char *s=(p ? p+1 : szF);
 if(szP)
  {B.Set(szP);
   if(szP[strlen(szP)-1]!='/')B.Add("/");
   B.Add(s);
  }
 else B.Set(s);	 
 strncpy(szR,B.GetBuf(),MaxL);
 szR[MaxL-1]=0;
 return strlen(szR);
}
// ********************************************************************
int CopyAsciiFile(FILE *fpFrom, FILE * fpTo)
{char szB[MAX_LINELENGTH+1];
 while(fgets(szB,MAX_LINELENGTH,fpFrom))
      {if(fputs(szB,fpTo)==EOF)return 0;}
 return 1;
} 	   		
// ********************************************************************
int CreateTmpFile(char *szFile, char * szName, const int Flag)
{// creates a temporary file with uniqe name
 if(isatty(STDIN_FILENO))
   {fprintf(stderr,"ERROR %s: No input file\n",szName);
    if(Flag==TMPFILE_RETURN)return 0;
    else exit(EXIT_FAILURE);
   }  
 sprintf(szFile,"%sXXXXXX",szName);
     if(!mktemp(szFile))
       {perror("ERROR Calc: Cannot create temporary file");
        if(Flag==TMPFILE_RETURN)return 0;
        else exit(EXIT_FAILURE);
       }
 strncat(szFile,".tmp",MAXPATH);
 FILE *ft=fopen(szFile,"w");
 if(!ft)
   {fprintf(stderr,"ERROR Calc: Error opening temporary file %s\n",szFile);
    if(Flag==TMPFILE_RETURN)return 0;
    else exit(EXIT_FAILURE);
   }
 CopyAsciiFile(stdin,ft);
 fclose(ft);
 return 1;
} 
// ******************************************************
int GetFloatRange(char *szExp,double & lfStart, double & lfEnd, const char cDil)
{//Get start end from 1.234:456.7
 // return 0 if OK
 char szB[MAX_LINELENGTH];
 char *pc = strchr(szExp,cDil);
 
 DelSpace(szExp);	
 if(!pc)return 1;  // no ":" in string
 if(strchr(++pc,cDil))return 2; // more than one ":"

 if(!StrTok(szExp,cDil,0,szB))return 3;
 if(!IsNumString(szB))return 4;

 errno=0;
 lfStart=strtod(szB,NULL);    
 if(errno){perror("GetFloatRange: Illegal input"); return 5;}
 if(!StrTok(szExp,cDil,1,szB))return 6;
 
 if(!IsNumString(szB))return 7;
 lfEnd=strtod(szB,NULL);    
 if(errno){perror("GetFloatRange: Illegal input"); return 8;}
 return 0;    
}
// *****************************************************
int GetIntRange(char *szExp,int & iStart, int & iEnd, const char cDil)
{//Get start end from 1.234:456.7
 // return 0 if OK
 char szB[MAX_LINELENGTH];
 char *pc = strchr(szExp,cDil);
 
 DelSpace(szExp);	
 if(!pc)return 1;  // no ":" in string
 if(strchr(++pc,cDil))return 2; // more than one ":"

 if(!StrTok(szExp,cDil,0,szB))return 3;
 if(!IsNumString(szB))return 4;

 errno=0;
 iStart=strtol(szB,(char**)NULL,10);    
 if(errno){perror("GetFloatRange: Illegal input"); return 5;}
 if(!StrTok(szExp,cDil,1,szB))return 6;
 
 if(!IsNumString(szB))return 7;
 iEnd=strtol(szB,(char**)NULL,10);    
 if(errno){perror("GetFloatRange: Illegal input"); return 8;}
 return 0;    
}
// *****************************************************
int StrTok(const char *s, const char t,const int i, char *d)
{if(i<0 || !s)return 0;
        int   j = 0;
 const char *pc = s;

 if(i==0)
   {char *p;
     strncpy(d,s,strlen(s));
     d[strlen(s)]=0;
     p=strchr(d,t);
     if(p)*p=0;
     return 1;           
   }

 for(j=0; j<i; j++)
    {pc=strchr(pc,t);
     if(!pc)return 0;
     pc++;
     strcpy(d,pc);
     char *ec = strchr(d,t);
     if(ec)*ec=0;
     /* if(!*d)return 0;*/
     if(!ec && !*d)return 0;     	   
    }
 return j;    	
}
// *****************************************************
#if defined(LINUX)
int GetFiles(const char *szArg, char *szBuffer, const int iSize)
// gets all files in szBuffer that start with szArg
{ 
  if(szBuffer==NULL || iSize<=0)return-1;
  if(szArg[0]==0 || szArg==NULL){szBuffer[0]=0;return 0;}

            DIR *Entry,*SubEntry;
  struct dirent *Found;

      const int  LINE_LENGTH = 256;
     const char *s = szArg;
           int   j = strlen(s);

  char szH[LINE_LENGTH+1], szPath[LINE_LENGTH+1], szFile[LINE_LENGTH+1];

  szBuffer[0] = 0;

  strcpy(szH,s);

  while(j>0 && s[j]!=' ')j--;
  s=&s[j+1];

  char *n = strchr(szH,' ');
  char *p = strrchr(s,'/');

  if(n)*n=0;
  //fprintf(stderr,"Arg:%s, s:%s, p:%s szH:%s *n:%s*\n",szArg,s,p,szH,n);

  if(!p){sprintf(szPath,"./");
         strncpy(szFile,s,LINE_LENGTH);
        }
  else {if(p==s){sprintf(szPath,"/");
                 strncpy(szFile,p+1,LINE_LENGTH);
                }
         else{strncpy(szPath,s,LINE_LENGTH);
              szPath[p-s+1]=0;
              p++;
              strncpy(szFile,p,LINE_LENGTH);
	     }     
	}
 
  //fprintf(stderr,"Path:%s File:%s\n",szPath,szFile); 

  if(!(Entry=opendir(szPath)))
    {snprintf(szBuffer,iSize,"Error opening %s\a\n",szPath);
     PRINT_DEBUG("%s\n",szBuffer)
     return 0;
    }
  int  iCount = 0;
  int    iLen = 0;
  int iRemain = iSize;

  strcpy(szBuffer,szPath);

  while( (Found=readdir(Entry)) )
    {//fprintf(stderr,"%s\n",Found->d_name);
     if(strstr(Found->d_name,szFile)==Found->d_name)
       {//fprintf(stderr,"%s\n",Found->d_name);
        iRemain-=strlen(Found->d_name)+1;
        if(iRemain){strcat(szBuffer,Found->d_name);
                    strcat(szBuffer,"\n");
                    iLen+=strlen(Found->d_name)+1;
	           }
        iCount++;
       }
    }
  
  if(iCount==1)
    {char szAnswer[LINE_LENGTH];
     szBuffer[strlen(szBuffer)-1]=0;
     if( (SubEntry=opendir(szBuffer)))
          {strcat(szBuffer,"/");
           closedir(SubEntry);
          }
     sprintf(szAnswer,"%s %s",szH,szBuffer);
     strcpy(szBuffer,szAnswer);
     iLen=1;
    } 
  else strcat(szBuffer,"\a");

  //fprintf(stderr,"stdfunc: %s\n",szBuffer);  

  closedir(Entry); 

  return iLen;
}
#endif
// ************************************************************************
double fLinTrafo(double p_max,double p_min,double d_max,double d_min,double d)

{double diffz, diffn;
 diffz= p_max-p_min;
 diffn=d_max-d_min;
 return  (diffz/diffn*d + (p_max-diffz/diffn*d_max));
}
//********************************************************************
double fBackTrafo(double p_max,double p_min, double b_max,double b_min,double b)
{double diffz, diffn;
 double h = p_max;
 diffz= p_max-p_min;
 diffn=b_max-b_min;
 if(diffn==0 || diffz==0)return 0;
 return  ((b-(h-(diffz/diffn)*b_max))*diffn/diffz);
}
//********************************************************************

int Scale(MDATA amin,MDATA amax, int max_res, MDATA &nmin, MDATA &nmax,
	  MDATA &ndelta)

{
  float  sf[3]={1.,2.,5.};
  float  delta,range,min_t,mult;
    int  i,j,k;
 double  h;
 double  n1,n2;
 double  tick;


  if(amin>=amax)return 0;
  if(max_res<=0)return 0;

  delta=amax-amin;
  min_t=delta/(float)max_res;
  mult=(min_t<sf[0] ? 0.1 : 10.);
  j=-1;
  if(min_t<sf[0])
    {range=sf[2];
     sf[2]=sf[0];
     sf[0]=range;
     for(i=0;;i++)
        {if(sf[i%3]<=delta && j==-1)j=i;
	 if(sf[i%3]<min_t)break;
	 sf[i%3]*=mult;
	}
     sf[0]=5.,sf[1]=2.,sf[2]=1.;
    }
  else
    {for(i=0;;i++)
        {if(sf[i%3]>=min_t && j==-1)j=i;
	 if(sf[i%3]>=delta)break;
	 sf[i%3]*=mult;
	}
      sf[0]=1.,sf[1]=2.,sf[2]=5.;
     }

  float *s = new float [abs(j-i)];
  CHECK_POINTER_RETURN(s,0)

  for(k=0;k<i;k++)
     {if(k>=j && k<=i)s[k-j]=sf[k%3];
      sf[k%3]*=mult;
     }

  tick=s[(mult<1?abs(i-j)-1:0)];

  h=amin/tick;
  if(!(h==0 || floor(h)==0))
    {if(fabs(ceil(h)/h-1)<=8e-7 || fabs(h/floor(h)-1)<=8e-7)h=floor(h+.5);}
  n1=floor(h);

  h=amax/tick;
  if(!(h==0 || floor(h)==0))
    {if(fabs(ceil(h)/h-1)<=8e-7 || fabs(h/floor(h)-1)<=8e-7)h=floor(h+.5);}
  n2=ceil(h);

  nmin=n1*tick;
  nmax=n2*tick;
  ndelta=tick;
  delete s;
  if(nmin==nmax)return 0;
  return 1;
}
// *****************************************************************
int GetArrayFromLine(FILE *f, double *lfA, const int iLen, const char *szSep)
{//return nuber of readed values OK; return 0 Error
 char *pc;
 char szLBuffer[MAX_LINELENGTH+1];
 if(fgets(szLBuffer,MAX_LINELENGTH,f)==0)return 0;
 szLBuffer[MAX_LINELENGTH]=0;
 pc=strtok(szLBuffer,szSep);
 lfA[0]=atof(pc);
 int i;
 for(i=1;i<iLen;i++)
	 {if((pc=strtok(NULL," ,"))==NULL)break;
	 lfA[i]=atof(pc);
     }
 return i-1;
}
// **********************************************************************
int GetArraysFromColumns(FILE *f, double **lfA, const int iRow, const int iCol,
                        const char *szSep)
{char *pc;
 char szLBuffer[MAX_LINELENGTH+1];
 int i,iCnt;
 for(iCnt=0;iCnt<iRow;iCnt++)
    {if(fgets(szLBuffer,MAX_LINELENGTH,f)==NULL){iCnt=0;break;}
     szLBuffer[MAX_LINELENGTH]=0;
     pc=strtok(szLBuffer,szSep);
     lfA[0][iCnt]=atof(pc);
	 for(i=1;i<iCol;i++)
	    {if((pc=strtok(NULL,szSep))==NULL)break;
	     lfA[i][iCnt]=atof(pc);
	     }//for i
     }//for iCnt
 return iCnt-1;
}
// **********************************************************************
int GetNumberFromString(String &S, double &lfA, 
		        const int iLen, const char cS)
{// First call: iLen>0 : Max. no of numbers (tokens) in line
 //              cS: seperator
 // Next calls  iLen=0
 //            szSep=0
 // lfA: returns double of converted number
 // return: >0: number of scan (token)
 //         =0: if no more token is found
 //         <0: max number of tokens (specified with 1st call) is exceed
 static char *ppp=0;
 static char *pp;
 static int iL=0,iC=0;
 static char cSep=0;
 static String SO;
 //static char DeLim;

 if(iLen>0 && cS)
   {if(ppp)
      {//fprintf(stderr,"deleting %p\n",ppp);
       delete ppp;
       ppp=0;
      }
    if(ppp==0){ppp=new char [S.GetSize()+1];
               //fprintf(stderr,"alocating %p\n",ppp);
              }
    pp=ppp;
    strncpy(pp,(const char *)S,S.GetSize());
    pp[S.GetSize()]=0;
    iL=iLen;
    iC=0;
    SO=S;
    cSep=cS;
    //DeLim.Set("");
    //SetDelim(S,DeLim);
   }
 else {if(pp==0)return -iC;}

 if((const char *)S==0 || (const char *)S[0]==0)return -1;

 lfA=strtod(pp,&pp);
 if(errno==ERANGE)PRINT_DEBUG("Not a number\n");
 while (*pp==' ')pp++;
 //if(*pp==0 || !strcmp(pp,(const char *)DeLim))return 0;
 if(*pp==0)return 0;
 iC++;
 if(cSep!=' ')
   {if(*pp==cSep)pp++;
    else return 0;
   }             
 if(!strcmp(SO,pp))return -iC;
 SO.Set(pp);
 //if(! (*pp && strcmp(pp,(const char *)DeLim)) )return 0;
 if(!*pp)return 0;
  //if((*pp==0 || !strcmp(pp,(const char *)DeLim)) && iC){i=-(iC+1);break;}
 return (iC>iL ? -iC : iC);
  
}
// **********************************************************************
double Horner(double *lfKoef,const int m,const double lfXw)

{double lfSum;
    int i;
 for(i=m-1,lfSum=lfKoef[m];i>=0;i--)lfSum=lfSum*lfXw+lfKoef[i];
 return lfSum;
}
// ************************************************************
// *********************************************************************
char *szLorentzP[]={"Amp","X0","Hbw","Backg"};

double Lorentz(const double x, const int nP, double *p, const int nC, double *c)
{// Sum of Lorentz functions: BG + A / { 1 + [ 2* (x - x0) / HBW]**2 } 
 //   A: Amplitude
 //  x0: Position of peak (maximum)
 // HBW: Half band width
 //  BG: Background 
 // nP: number of curves (NOT the number of pars)
 // p[0]: Amplitude[0] -> [n*3]
 // p[1]: x0[0]
 // p[2]: HBW[0]
 // p[3]: Amlidude[1]
 // p[4]: x0[1]
 // p[5]: HBW[1]
 // ... 
 // p[3*nP] Background
 
    int i,j;
 double s;
 for(i=0,s=0; i<nP; i++) 
    {j=3*i;
     s+=( p[j] / (1 + SQUARE(2. * ( x - p[j+1])/p[j+2]) ) );
    }
 return s+p[3*nP];    
}
// ********************************************************************
// Global variables to pass wavelength and peak ratio of alfa1 and alfa2
// can be set by caller
// extern double Ra1a2;
// extern double RWLa21;
double  Ra1a2 = 0.5;
double RWLa21 = 1.79278/1.78892;
//
char *szLorentz12P[]={"Amp","Theta0","Hbw","Backg","RWl"};
double Lorentz12(const double x, const int nP, double *p, const int nC, double *c)
{//A sum of lorentzians can be generated
 // nP holds the number of curves to be summed up
 // sequence: a1,deta1,a2,deta2,...,hbw,bgr
 //  hbw and bgr the same for all curves
 //  nP*2+1 parameters
 //  nP*2+1 : bgr
 //  
    int i;
 double s;

 for(i=0,s=0; i<nP*2; i+=2)
    s+=( p[i] / (1 + SQUARE(2.*( x - p[i+1]) / p[nP*2])) +
    Ra1a2*p[i] / (1 +SQUARE(2.*( x - (a2_2th(0.5*p[i+1],c[0])))/p[nP*2]) ) );

 return(s+p[nP*2+1]);
}
// ********************************************************************
char *szCurieWeissP[]={"C","Chi0","Theta"};
double CurieWeiss(const double t, const int nP, double *p, const int nC, double *c)
{// p[0]: c  p[1]:chi0  p[3]:theta
 if(nP!=3){PRINT_DEBUG("Illegal # of pars %d (3)\n",nP)
           return 0;
	  }
 return p[1] +	(p[0]/ (t - p[2]));
}  
// ********************************************************************	   
char *szNonFermiResP[]={"Rho0","N2(Ef)","Exponent"};
double NonFermiRes(const double t, const int nP, double *p, const int nC, double *c)
{// p[0]: rho0  p[1]:N2(Ef)  p[3]:exponent
 if(nP!=3){PRINT_DEBUG("Illegal # of pars %d (3)\n",nP)
           return 0;
	  }
 return p[0] +	p[1]*pow(t,p[2]);
}  
// ********************************************************************	   
char *szStraightLineP[]={"Slope","d"};
double StraightLine(const double t, const int nP, double *p, const int nC, double *c)
{// p[0]: k (slope)  p[1]: d 
 if(nP!=2){PRINT_DEBUG("Illegal # of pars %d (2)\n",nP)
           return 0;
	  }
 return p[0]*t + p[1];
}  
// **********************************************************	   
char *szPolynomP[]= {"C"};
double Polynomal(const double t, const int nP, double *p, const int nC, double *c)
{// p[0]: x^0, p1:x^1, p[2]:x^2
 return Horner(p,nP-1,t);
}
// **********************************************************
// class Polynom functions
// **********************************************************
Polynom::Polynom(void):Null(0)
{   p = NULL;
 iDeg = 0;
}
// **********************************************************
Polynom::Polynom(const int d, double *pi): Null(0)
{   p = NULL;
 iDeg = -1;

 if(d<0)return;

 p=new double [d+1];
 CHECK_POINTER_RET(p);
 iDeg=d;
 for(int i=0; i<=iDeg; i++)p[i]=(pi ? pi[i]: 0.0);
}
// ********************************************************* 
Polynom::Polynom(Polynom &N):Null(0)
{   p = NULL;
 iDeg = -1;

 if(N.iDeg<0)return;

 p=new double [N.iDeg+1];
 CHECK_POINTER_RET(p);

 iDeg=N.iDeg;
 for(int i=0; i<=iDeg; i++)p[i]=N.p[i];
}
// *********************************************************
Polynom::~Polynom()
{if(iDeg>=0)delete [] p;

}
// **********************************************************
double & Polynom::operator[](const int i)
{
 CHECK_INDEX_RETURN(i,iDeg+1,Null)
 return p[i];
}
// *********************************************************
Polynom & Polynom::operator=(Polynom &N)
{ if(this != &N)
    {if(iDeg>=0 && p)delete p;
     p=new double [N.iDeg+1];
     iDeg=N.iDeg;
     for(int i=0; i<=iDeg; i++)p[i]=N.p[i];
    }
  return *this;
}
// ********************************************************
double Polynom::operator()(const double x) const
{if(iDeg<0 || p==0)
   {PRINT_DEBUG("Undefined Polynom (deg: %d)\n",iDeg)
    return 0; 
   }
 double w = p[iDeg];
 for(int i=iDeg-1; i>=0; i--)w=w*x+p[i];
 return w;
}
// ***********************************************************
int Polynom::Save(FILE *f, const double xS, const double xE, const int nPoints,
	          const char * szText, enum EndLine  LineT)
{if(xS>=xE)
   {PRINT_DEBUG("xStart (%14.7g) >= xEnd (%14.7g)\n",xS,xE)
    return 0;
   } 

 if(nPoints<=2)
   {PRINT_DEBUG("Illegal number of points (%d)\n",nPoints)
    return 0;
   } 

 double lfD = (xE-xS)/nPoints;  

 char szB[MAX_LINELENGTH],szE[3];
 GetCurrentTime(szB);

 strcpy(szE,"\n");
#if defined(LINUX)
  if(LineT==NOEnd || LineT==DOSEnd)strcpy(szE,"\r\n");
#endif 

 fprintf(f,"File created from XAusw (C)HM %s",szB);
 if(szText)fprintf(f,"%s%s",szText,szE);
 int i;
 for(i=0; i<=Degree(); i++)fprintf(f,"x^%d: %-14.7g%s",i,p[i],szE);

 double lfX = xS;
 while (lfX<xE)
       {fprintf(f,"%14.7g  %14.7g%s",lfX,this->operator()(lfX),szE);
        lfX+=lfD;
       }
 return 1; 
}
// *****************************************************************
