//File: filepar.cpp
//$Log: filepar.new.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//$Id: filepar.new.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>

#ifndef SFILEPAR_H
#include "filepar.h"
#endif

#ifndef STDINC_H
#include "stdinc.h"
#endif

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

//#define TEST_FILEPAR

char * FilePar::szMsg[]={"FilePar:Ok",
                        "FilePar:File not found", "FilePar:Topic not found",
                        "FilePar:Name not found", "FilePar:Out of memory",
                        "FilePar:Bad input",      "FilePar:Invalid file"
                       };

//**********************************************************************
FilePar::FilePar(const char *szFileName,const int iClose)
{FILE *fp;
 nLines=0;
 iStatus=PFILE_NFOUND;
 if ((fp=fopen(szFileName,"r")) == NULL)return;
 ReadFile(fp,iClose);
} 
//**********************************************************************
FilePar::FilePar(FILE *f, const int iClose)
{
 nLines=0;
 iStatus=FOUND;
 ReadFile(f,iClose);
}
//**********************************************************************
void FilePar::ReadFile(FILE *f,const int iClose)
{char szB[MAX_LINELENGTH+2],*n;
 int iPCnt=0;
 Pars.Set("");
 while((n=fgets(szB,MAX_LINELENGTH,f))!=NULL)
 {szB[MAX_LINELENGTH+1]=0;
  if(iPCnt==0)
    {DeLim.Set(strchr(szB,'\n'));
     if(!DeLim.GetBuf()){iStatus=INVALID_FILE;return;}
    }
  if(!strcmp(szB,"@EOH"))break;
  if(szB[0]==';' || szB[0]=='\n'|| szB[0]=='#')continue;
  if(szB[0]=='[' || strchr(szB,'=')!=0)
    {if(!Pars.Add(szB)){iStatus=NO_MEM;fclose(f);return;}
     iPCnt++;
    }
 }
 if(iPCnt<2){iStatus=BAD_INPUT; fclose(f); return;}

 nLines=Pars.Lines(DeLim);
 iStatus=FOUND;

 if(iClose==F_CLOSE)fclose(f);
}
// **********************************************************************
FilePar::~FilePar(void)
{ }
// *********************************************************************
const char* FilePar::FindPar(const char *szTopic, const char *szName)
{

String NewPar(Pars(szTopic));
if(!NewPar.GetSize()){iStatus=TOPIC_NFOUND; return 0;}

char szB[MAX_LINELENGTH+1];
sprintf(szB,"%s[",DeLim.GetBuf());
const char *pNext=NewPar(szB);
const char *pPar=NewPar.Find(szName,0,DeLim);

if(pPar==0){iStatus=TOPIC_NFOUND; return 0;}
if(pNext){if(pPar>pNext){iStatus=TOPIC_NFOUND; return 0;}}

iStatus=FOUND;
return pPar;

}
//**********************************************************************
int FilePar::GetPar(const char *szTopic, const char *szName, char *szFormat, ...)
{
 va_list arg_ptr;
 char szB[MAX_LINELENGTH+1];
 const char *Act=FindPar(szTopic,szName);
 if(Act==0 || iStatus!=FOUND)return iStatus;

 va_start(arg_ptr, szFormat);
 if( strchr(szFormat,'S')!=0)
	{char *s=strchr(szFormat,'%'),*e=strchr(szFormat,'S');
	 if(!s || !e)return BAD_INPUT;
	 int i;
         if(s+1 == e)i=0;
	 else i=atoi(s+1);
	 if(i<0){iStatus=BAD_INPUT; return iStatus;}
	 s=strchr(Act,'=');
         if(i==0)i=(int)strlen(s);
	 if(i<(int)strlen(s))*(s+i)=0;
	 strncpy(va_arg(arg_ptr,char *),s+1,i);
	 va_end(arg_ptr);
	 return iStatus;
	}
 sprintf(szB,"%s=%s",szName,szFormat);
 vsscanf(Act, szB, arg_ptr);
 va_end(arg_ptr);

 return iStatus;

}
//**********************************************************************
void FilePar::Print(FILE *f) const
{
 fprintf(f,"%s",Pars.GetBuf());
}
//**********************************************************************
#ifdef TEST_FILEPAR
main(int iArgC, char ** szArgV)
{
 if(iArgC<1)exit(EXIT_FAILURE);

 if(iArgC>2)
   {fprintf(stderr,"Usage %s filename\n",szArgV[0]);
    exit(EXIT_FAILURE);
   }

 FilePar FP(szArgV[1]);
 if(FP.GetStatus()!=0)
   {fprintf(stderr,"Error %s\n",FP.GetMsg());
    exit(EXIT_FAILURE);
   }
 FP.Print(stdout);

 char szB[MAX_LINELENGTH+1]={"?"};

 int i=-1;
 i=FP.GetPar("Cell_M","CalFile","%S",szB);
 fprintf(stdout,"%d: %s\n",i,szB);

 unsigned u=0;
 i=FP.GetPar("User","ernst","%x",&u);
 fprintf(stdout,"%d: %x\n",i,u);

 double d=0;
 i=FP.GetPar("SxsPars","Alpha","%lf",&d);
 fprintf(stdout,"%d: %f\n",i,d);

 exit(EXIT_SUCCESS); 
}

#endif








