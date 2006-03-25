//File: filepar.cpp
//$Log: filepar.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.3  1999/03/15 09:08:37  herbie
//*** empty log message ***
//
//$Id: filepar.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
#define _GNU_SOURCE 1
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>

#include "filepar.h"
#include "stdinc.h"
#include "stdfunc.h"

//#define FPAR_DEBUG

#if defined (FPAR_DEBUG)
FILE *fd=NULL;
#endif

//**********************************************************************
FilePar::FilePar(const char *szFileName,const int iClose)
{FILE *fp;
 iLen=0;
 iStatus=PFILE_NFOUND;
 Lines=NULL;

#if defined (FPAR_DEBUG)
if(!fd){fd=fopen("fpar.bug","a");fprintf(fd,"File Opend(Name)\n");}
else fflush(fd);
#endif

 if ((fp=fopen(szFileName,"r")) == NULL)return;
 ReadFile(fp,iClose);
}
//**********************************************************************
FilePar::FilePar(FILE *f, const int iClose)
{
 Lines=NULL;
 iLen=0;

#if defined (FPAR_DEBUG)
if(!fd){fd=fopen("fpar.bug","a");fprintf(fd,"File Opend (FILE)\n");}
else fflush(fd);
#endif

ReadFile(f,iClose);
}
//**********************************************************************
void FilePar::ReadFile(FILE *f,const int iClose)
{char szB[MAX_PAR_LEN],*n;
 int iPCnt=0;
 while((n=fgets(szB,MAX_PAR_LEN-2,f))!=NULL)
 {RemoveCR_LF(szB);
  if(!strcmp(szB,"@EOH"))break;
  if(szB[0]==';' || szB[0]=='\n'|| szB[0]=='#')continue;
  if(szB[0]=='[' || strchr(szB,'=')!=0)
	 { n=strchr(szB,'\n');
	  if(n!=NULL)n[0]=0;
	  if(!Insert(szB)){iStatus=NO_MEM;fclose(f);return;}
	  iPCnt++;
	 }
 }
 if(iPCnt<2){iStatus=BAD_INPUT; fclose(f); return;}

  iLen=iPCnt;
  iStatus=FOUND;

 if(iClose==F_CLOSE)fclose(f);
}
// **********************************************************************
FilePar::~FilePar(void)
{DeleteLines();}
// *********************************************************************

void FilePar::DeleteLines(void)
{FLINE *Act,*Next;
 Act=Lines;
 while(Act)
		{if(Act->szLine)
#if defined (FPAR_DEBUG)
if(fd)fprintf(fd,"%p: : Del Sz %p: %s Next -> %p\n",Act,Act->szLine,Act->szLine,Act->Next);
#endif
			delete Act->szLine;Act=Act->Next;
		 }
 Act=Next=Lines;
 while(Next){Next=Act->Next;
#if defined (FPAR_DEBUG)
if(fd)fprintf(fd,"Del %p: Next -> %p\n",Act,Act->Next);
#endif
				 delete Act; Act=Next;}
#if defined (FPAR_DEBUG)
if(fd){fprintf(fd,"File Closed\n");fclose(fd);fd=NULL;}
#endif
Lines=NULL;
iLen=0;
}
//**********************************************************************
FLINE *FilePar::FindPar(const char *szTopic, const char *szName)
{int bFound=FALSE;
 FLINE *Act=Lines;

//  if(iStatus!=1 || Act->Next==NULL){iStatus=TOPIC_NFOUND; return NULL;}
  if(Act->Next==NULL){iStatus=TOPIC_NFOUND; return NULL;}
  while(Act->Next)
		 {if(Act->szLine[0]=='[' && strstr(Act->szLine,szTopic)!=0)
			 {bFound=TRUE;break;}
		  Act=Act->Next;
		 }
 if(!bFound){iStatus=TOPIC_NFOUND; return NULL;}

 bFound=FALSE;
 Act=Act->Next;
 while(Act->szLine[0]!='[')
 {if(strstr(Act->szLine,szName)!=0)
	{bFound=TRUE;break;}
  if(Act->Next==NULL)break;
  Act=Act->Next;
 }

 if(!bFound){iStatus=NAME_NFOUND; return NULL;}
 else iStatus=FOUND;

 return Act;
}
//**********************************************************************
int FilePar::GetPar(const char *szTopic, const char *szName, char *szFormat, ...)
{
 va_list arg_ptr;
 char szB[2*MAX_PAR_LEN];
 FLINE *Act;

 Act=FindPar(szTopic,szName);
 if(Act==NULL || iStatus!=FOUND)return iStatus;

 va_start(arg_ptr, szFormat);
 if( strchr(szFormat,'S')!=0)
	{char *s=strchr(szFormat,'%'),*e=strchr(szFormat,'S');
	 if(!s || !e)return BAD_INPUT;
	 int i;
         if(s+1 == e)i=0;
	 else i=atoi(s+1);
	 if(i<0)return BAD_INPUT;
	 s=strchr(Act->szLine,'=');
         if(i==0)i=(int)strlen(s);
	 if(i<(int)strlen(s))*(s+i)=0;
	 strncpy(va_arg(arg_ptr,char *),s+1,i);
	 va_end(arg_ptr);
	 return iStatus;
	}
 sprintf(szB,"%s=%s",szName,szFormat);
 vsscanf(Act->szLine, szB, arg_ptr);
 va_end(arg_ptr);

 return iStatus;

}
//**********************************************************************
int FilePar::Insert(const char *NewE)
{FLINE *Act,*New;
 int iLen=strlen(NewE)+1;

 if( !(New=new FLINE) )return 0;
#if defined (FPAR_DEBUG)
 if(fd)fprintf(fd,"Alloc New %p: ",New);
#endif
 if(!(New->szLine=new char [iLen]) )return 0;
 strcpy(New->szLine,NewE);
 New->Next=NULL;
#if defined (FPAR_DEBUG)
 if(fd)fprintf(fd,"New  Sz:(%p.%d) %s Next -> %p\n",New->szLine,iLen,New->szLine,New->Next);
#endif

 if(Lines==NULL){Lines=New; return 1;}
 Act=Lines;
 while(Act->Next)Act=Act->Next;
 Act->Next=New;
#if defined (FPAR_DEBUG)
 if(fd)fprintf(fd,"Insert%p: Sz:(%p) %s Next -> %p\n",Act,Act->szLine,Act->szLine,Act->Next);
 fflush(fd);
#endif
 return 1;
}
//**********************************************************************
void FilePar::Print(void)
{FLINE *Act=Lines;
 while(Act->Next){printf("%s\n\r",Act->szLine);Act=Act->Next;}
}
//**********************************************************************
#if defined (FPAR_DEBUG)
void FilePar::PrintTree(void)
{FLINE *Act=Lines;
 if(!fd)return;
 fprintf(fd,"Top: %p\n",Lines);
 while(Act->Next){fprintf(fd,"%p: %s Next -> %p\n",Act,Act->szLine,Act->Next);Act=Act->Next;}
}
#endif









