//File: strings.cpp
// $Id: strings.cpp,v 1.2 2001/10/16 16:44:04 xausw Exp xausw $
// $Log: strings.cpp,v $
// Revision 1.2  2001/10/16 16:44:04  xausw
// Calib and Gap ignore empty lines and stop reading after the next header
//
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.1  1999/06/29 16:30:34  herbie
// Initial revision
//
#define _GNU_SOURCE 1

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
#include <regex.h>

#ifndef STRINGS_H
#include "strings.h"
#endif

#ifndef STDFUNC__H
#include "stdfunc.h"
#endif

//#define STRINGS_TEST
// ***************************************************
int SetDelim(const char *s, String &Delim)
{ if(strstr(s,"\r\n")){Delim.Set("\r\n"); return 1;}
  if(strstr(s,"\n\r")){Delim.Set("\n\r"); return 1;} 
  if(strstr(s,"\n")){Delim.Set("\n"); return 1;} 
  if(strstr(s,"\r")){Delim.Set("\r"); return 1;} 
  return 0;
}
// **************************************************

// ***************************************************
// ***************************************************
// class String
// ***************************************************
char String::Null=0;
// ***************************************************
String::~String()
{if(iSize>0 && b)delete b;} 
// ***************************************************
int String::Copy(const String & cb)
{
 if(cb.b)
      {b=new char [cb.iSize+1];
       CHECK_POINTER_RETURN(b,0)
       strcpy(b,cb.b);
       iSize=strlen(b);
      }
    else {b=0; iSize=0;}
 return 1;
}
// **************************************************
int String::Copy(const char *cb)
{
 if(cb)
      {b=new char [strlen(cb)+1];
       CHECK_POINTER_RETURN(b,0)
       strcpy(b,cb);
       iSize=strlen(b);
      }
    else {b=0; iSize=0;}
 return 1;
}
// **************************************************
int String::Add(const char *c, const int iLen )
{if(!c)return 0;
 if(c[0]==0 || iLen<0)return 0;

 char *h;
 int iL=(iLen==0 ? strlen(c) : iLen);

 h=new char [iSize+iL+1];
 CHECK_POINTER_RETURN(h,0)

 h[0]=0;
 if(b)strcpy(h,b);
 if(b)delete b;
 b=h;
 strncat(b,c,iL);
 iSize+=iL;
 b[iSize]=0;
 return iSize;
} 
// *************************************************** 
int String::Set(const char *c, const int iLen )
{if(iSize>0)delete b;
 b=NULL;
 iSize=0;
 int iL=(iLen==0 ? strlen(c) : iLen);
 b=new char [iL+1];
 CHECK_POINTER_RETURN(b,0)
 if(iL>0)strncpy(b,c,iL);
 b[iL]=0;
 iSize=strlen(b);
 return iSize;
}

// **********************************************************
int String::Setf(const char *szFormat, ...)
{
 va_list arg_ptr;
 char szB[MAX_LINELENGTH+1];

 va_start(arg_ptr, szFormat);
 vsnprintf(szB,MAX_LINELENGTH,szFormat,arg_ptr);
 szB[MAX_LINELENGTH]=0;
 va_end(arg_ptr);
 return String::Set(szB);
 
}
// *********************************************************
int String::Cut(const unsigned iBehind)
{if(!b || !iSize)return 0;
 if((int)iBehind>=iSize)return 0;

 char *nb=new char [iBehind+1];
 CHECK_POINTER_RETURN(nb,0);
 b[iBehind]=0;
 strcpy(nb,b);
 delete b;
 b=nb;
 iSize=iBehind+1;
 return 1;
}
// *********************************************************
int String::Select(const char * szD, const unsigned iStart,int iEnd)
{//szD: delimiter
 //iStart: first occurence of delimiter
 //iEnd: last occurence of delimiter (-1: end)
 //  ----------------------------------
 // |         |           |            |
 // -----------------------------------
 // 0         1           2            3 (iStart, iEnd)
 if(!b || !iSize)return 0;
 if(iEnd==-1)iEnd=Count(szD);
 if((int)iStart>=iEnd)return 0;
 int i;
 char *ps=b;
 for(i=0;i<(int)iStart;i++)
    {ps=strstr(ps,szD);
     if(!ps)return 0;
     ps+=strlen(szD);
    }
 //if(iStart>0)ps+=strlen(szD);

 char *pe=ps;
 for(i=(int)iStart;i<iEnd;i++)
    {pe=strstr(pe,szD);
     if(!pe)return 0;
     pe+=strlen(szD);
    }
 //pe+=strlen(szD);
 
 int iL=pe-ps+1;
 char *nb=new char [iL];
 CHECK_POINTER_RETURN(nb,0);
 *pe=0;
 strcpy(nb,ps);
 delete b;
 b=nb;
 iSize=iL;
 return pe-ps;
}
// *********************************************************
int String::Count(const char *szD)
{if(!szD || !*szD)return 0;
 int i=0;
 const char *p=b;
 while( (p=strstr(p,szD)) ){i++; p+=strlen(szD);;}
 return i;
}
// *********************************************************
int String::Addf(const char *szFormat, ...)
{
 va_list arg_ptr;
 char szB[MAX_LINELENGTH+1];

 va_start(arg_ptr, szFormat);
 vsnprintf(szB,MAX_LINELENGTH,szFormat,arg_ptr);
 szB[MAX_LINELENGTH]=0;
 va_end(arg_ptr);
 return String::Add(szB);
}
// ***************************************************
int String::Strupr(void)
{ if(!b || !iSize)return 0;
  int i;
  for(i=0;i<iSize;i++)b[i]=toupper(b[i]); 
  return i-1;
}
// ************************************************
void String::RemoveCR_LF(void)
{
  //   Replaces \r\n combinations by \0
  //
  if (!b || !iSize)return;    
  char *p=strchr(b,'\r');
  char *q=strchr(b,'\n');

  if( !(p-b>iSize || p==0) )*p=0;
  if( !(q-b>iSize || q==0) )*q=0;
  if((int)strlen(b)==iSize)return;

  char *n=new char [strlen(b)+1];
  strcpy(n,b);
  delete b;
  b=n;
  iSize=strlen(b); 
}
// ************************************************
int String::DelSpace(void)
{// deletes space in buffer
  int iC = 0; 
  int i,j;
  if (!b || !iSize)return 0;    
  for(i=0;i<iSize;i++)if(b[i]==' ')iC++;
  if(iC==0)return 0;

  char *n=new char [iSize-iC+1];

  for(i=0,j=0;j<iSize;i++,j++)
     {if(b[j]!=' ')n[i]=b[j];
      else {i--;}
     }
  delete b;
  b=n;
  iSize-=iC; 
  return iC;       	
}		
// ************************************************
char String::operator[](const int i)
{if(!CHECK_INDEX(0,iSize))return Null;
 return  b[i];
}
// ***************************************************
int String::Read(const char *szFile, const int iAppend)
{// iAppend==1 : Contens of file is appended on string
 // iAppend==0 : String is set to contens of file
 FILE *f=fopen(szFile,"r");
 if(!f)
   {PRINT_DEBUG("Can not open %s\n",szFile);
    perror("Reason");
    return 0;
   }
 int i=Read(f,iAppend);
 fclose(f);
 return i;
}
// ***************************************************
int String::Read(FILE *f, const int iAppend)
{if(!f)return 0;
 char szB[MAX_LINELENGTH+1];
 int i=0;
 while(1)
      {if(!fgets(szB,MAX_LINELENGTH,f))return 1;
       if(i==0)
         {if(iAppend)Add(szB);
          else Set(szB);
         }//i==0
       else Add(szB);
       i++;
      }
}
// ***************************************************
int String::Write(const char *szFile) const
{FILE *f=fopen(szFile,"w");
 if(!f)
   {PRINT_DEBUG("Can not open %s\n",szFile);
    perror("Reason");
    return 0;
   }
 int i=Write(f);
 fclose(f);
 return i;
}
// ***************************************************
int String::Write(FILE *f) const
{if(fputs(b,f)==EOF)return 0;
 return 1;
}
// ***************************************************
int String::Contains(const char *s) const
{if(!b)return 0;
 const char *p=strstr(b,s);
 return (p ? p-b+1 : 0);
}
// ***************************************************
const char * String::operator()(const int i, const enum String::Dir d) const
{if(i>iSize || i<=0)return 0;
 switch(d)
       {case FROM_LEFT: return &(b[i-1]);
        case FROM_RIGHT: return &(b[iSize-i]);
       }
 return 0;
}
// ***************************************************
const char * String::operator()(const char *s) const
{if(!*s || !s || iSize<=0 || b==0 )return 0;
 return strstr(b,s);
}
// ***************************************************
const char * String::Find(const char *s, const int iNext) const
{
 //iNext=0 :   pointer 1st occurence of s is returned
 //iNext=1 :   pointer to next occurence of previous search is returned

 if(!*s || !s || iSize<=0 || b==0)return 0;
 char *p;
 static char *pFound;
 if(iNext==0)pFound=0; 
 if(iNext && pFound>=b && pFound<b+iSize)p=strstr(++pFound,s);
 else p=strstr(b,s);
 pFound=p;
 return pFound;
}
// ***************************************************
// Class LineString
// **************************************************
String LineString::Result("");
// **************************************************
void LineString::Init(const char *s)
{
 //SetDelim((const char *)s,DeLim);
 nLines=CountLines(s);
 if(!Alloc(nLines)){nLines=0;iPlceFor=0;}
 if(nLines<=0)return;

 int iC=0;
 const char *p=s;
 const char *e;
 while(iC<nLines)
      {e=strstr(p,(const char *)DeLim);
       if(!e)break;
       Lines[iC]->Set(p,e-p);
       p=e+strlen(DeLim);
       iC++;       
      }
}
// **************************************************
int LineString::Alloc(const int iL)
{if(iL<=0)return 0;
 Lines=new String * [iL];
 CHECK_POINTER_RETURN(Lines,0)
 for (int i=0;i<iL;i++)Lines[i]=new String;
 nLines=iPlceFor=iL;
 return nLines;
}
// **************************************************
void LineString::Dealloc(void)
{if(nLines<=0)return;
 int i;
 for(i=0;i<nLines;i++)delete Lines[i];
 delete Lines;
 Lines=0;
 nLines=iPlceFor=0;
}
// **************************************************
int LineString::Realloc(const int iNewS)
{if(iNewS<=0)return 0;
 String **NL=new String * [iNewS];
 CHECK_POINTER_RETURN(NL,0)
 int i;
 for (i=0; i<nLines;i++)NL[i]=Lines[i];
 if(iNewS>nLines)
   {for(i=nLines; i<iNewS; i++)NL[i]=new String;
    iPlceFor=iNewS;
   }
 else {iPlceFor=nLines=iNewS;}
 delete Lines;
 Lines=NL;  
 return 1;
}
// **************************************************
int LineString::CountLines(const char *s)
{if(*s==0 || s==0)return -1;
 if(!SetDelim(s,DeLim))return -1;
 int iC=0;
 const char *p=s;
 while(*p)
      {if( (p=strstr(p,(const char *)DeLim)) )
          {iC++;p+=strlen((const char *)DeLim);}
      }	  
 return iC;
}
// **************************************************
const char * LineString::Line(const int i, const int WithDel)
{// Index of first line is 0
 if(Lines==0)return 0;
 if(DeLim.IsEmpty() || DeLim.IsZero())return 0;

 if(!WithDel)return Lines[i]->GetBuf();
 
 Result.Set(Lines[i]->GetBuf());
 Result.Add((const char *)DeLim);
 return (const char *)Result;
}
// **************************************************
int LineString::InsertLine(const char * szL, int iBefore )
{// Index of first line is 0
 // iBefore==0: insert at begin (before first)
 // iBefore==-1 or = nLines insert at end
if(DeLim.IsEmpty() || DeLim.IsZero())return 0;

if(iBefore==-1)iBefore=nLines;
if(iBefore<0 || iBefore>nLines)return 0;

if(iPlceFor<nLines+1)
  {if(!Realloc(nLines+100))return 0;}

int i;
 for(i=nLines;i>iBefore;i--)Lines[i]=Lines[i-1];
 if(iBefore<nLines)Lines[i]=new String(szL);
 else Lines[i]->Set(szL);
 nLines++;
 return nLines;
}
// **************************************************
int LineString::ReplaceLine(const char * szL, const int iL)
{
 if(iL<0 || iL>=nLines || szL==0)return 0;
 if(DeLim.IsEmpty() || DeLim.IsZero())return 0;
 String *NS=new String(szL);
 CHECK_POINTER_RETURN(NS,0)
 delete Lines[iL];
 Lines[iL]=NS;
 return 1;
}
// **************************************************
int LineString::RemoveLine(const int iL)
{
 if(iL<0 || iL>nLines)return 0;
 if(DeLim.IsEmpty() || DeLim.IsZero())return 0;
 int i;

 delete Lines[iL];

 for(i=iL;i<nLines-1;i++)Lines[i]=Lines[i+1];

 nLines--;
 return 1;
}

// **************************************************
int LineString::Read(const char *szFile, const char *szEnd=0,
                     const enum Search es)
{// szEnd==0 : Read unil eof
 // szEnd!=0 : read until szEnd is found
 // se: controls what is "found"
 FILE *f=fopen(szFile,"r");
 if(!f)
   {PRINT_DEBUG("Can not open %s\n",szFile);
    perror("Reason");
    return 0;
   }
 int i=Read(f,szEnd,es);
 fclose(f);
 return i;
}
// **************************************************
int LineString::Read(FILE *f, const char *szEnd=0, 
                     const enum Search es=S_CONTAINS)
{
 // Reads file until a condition controlled by es occurs
 //   es==S_EQUAL: until line is equal szEnd
 //                return: number of lines 
 //                        >0 string was not found
 //                        <0 string was found
 //es==S_CONTAINS: until line contains szEnd
 //                return: number of lines 
 //                        >0 string was not found
 //                        <0 string was found
 //  es==S_NUMBER: until line is a correct series of numbers
 //                with szEnd[0] as separator
 //                return: No. of columns of first line detected
 //   es==S_ISCOL: while Line is a correct series of numbers
 //                with szEnd[0] as separator
 //                return: No. of data rows
 // szEnd==0: Complete file is readed
 //           return: No. of readed line
 if(!f)return 0;
 char szB[MAX_LINELENGTH+1];
 int i=0;
 int iC;
 char *ppp=0;
 int iRet=-1;

 if(es==S_NUMBER || es==S_ISCOL || es==S_ISCOL_BREAK)
                        ppp=new char [MAX_LINELENGTH+1];
 char *pp=0;

 Dealloc();
 while(1)
      {if(!fgets(szB,MAX_LINELENGTH,f)){iRet=i; break;}
       szB[MAX_LINELENGTH]=0;
       if(i==0)SetDelim(szB,DeLim);
       if(szEnd && *szEnd)
         {if(es==S_CONTAINS)
            {if(strstr(szB,szEnd))break;}
          if(es==S_EQUAL)
            {if(!strcmp(szB,szEnd))break;}
          if(es==S_NUMBER)
            { iC=0;
              pp=ppp;
              char szO[MAX_LINELENGTH+1]; 
              strncpy(szO,szB,MAX_LINELENGTH);
              szO[MAX_LINELENGTH]=0;
              strncpy(pp,szB,MAX_LINELENGTH);
              pp[MAX_LINELENGTH]=0;
              do{double d=strtod(pp,&pp);
                 if(errno==ERANGE)PRINT_DEBUG("Not a number\n");
                 while (*pp==' ')pp++;
                 if(!strcmp(pp,(const char *)DeLim) || *szB==0)break;
                 iC++;
                 if(*szEnd!=' ')
                   {if(*pp==*szEnd)pp++;
                    else break;
                   }             
                 if(!strcmp(szO,pp))break;
                 strcpy(szO,pp);
                }while(*pp && strcmp(pp,(const char *)DeLim));
             if((*pp==0 || !strcmp(pp,(const char *)DeLim)) && iC){i=-(iC+1);break;}
            }// if(es==S_NUMBER)
           if(es==S_ISCOL || es==S_ISCOL_BREAK)
            { iC=0;
              pp=ppp;
              char szO[MAX_LINELENGTH+1]; 
              strncpy(szO,szB,MAX_LINELENGTH);
              szO[MAX_LINELENGTH]=0;
              strncpy(pp,szB,MAX_LINELENGTH);
              pp[MAX_LINELENGTH]=0;
              do{double d=strtod(pp,&pp);
                 if(errno==ERANGE)PRINT_DEBUG("Not a number\n");
                 while (*pp==' ')pp++;
                 if(!strcmp(pp,(const char *)DeLim) || *szB==0)break;
                 iC++;
                 if(*szEnd!=' ')
                   {if(*pp==*szEnd)pp++;
                    else break;
                   }             
                 if(!strcmp(szO,pp))
                   {if(es==S_ISCOL_BREAK)
                      { delete ppp;
                        return i;
                      }
                    else break;
                   }
                 strcpy(szO,pp);
                }while(*pp && strcmp(pp,(const char *)DeLim));
             if(! ((*pp==0 || !strcmp(pp,(const char *)DeLim)) && iC))
                {continue;}
            }//if(es==S_ISCOL)
         }//if(szEnd && *szEnd)
       szB[strlen(szB)-strlen((const char *)DeLim)]=0;
       InsertLine(szB,-1);
       i++;
      }//while(1)
 delete ppp;
 return (iRet<0 ? -i : iRet);
}
// **************************************************
const char * LineString::Find(int &iL, const char *s,
                              const unsigned iStart, const unsigned iNext)
{
 //returns Line where *s is found
 //iStart: First line to search 
 // iNext: number of matches after that should be returned
 //
 iL=0;
 if(!*s || !s )return 0;
 if(Lines==0)return 0;
 if(DeLim.IsEmpty() || DeLim.IsZero())return 0;
 if(!iNext)return 0;
 if((int)iStart>nLines)return 0;
 int i;
 unsigned j=0;
 //const char *p;
 for(i=iStart;i<nLines;i++)
    {if(strstr(Lines[i]->GetBuf(),s))j++;
     if(j==iNext)break;
    }
 iL=(i<nLines?i:0);
 return (i<nLines?Lines[i]->GetBuf():0);
}
// **************************************************
int LineString::GetPar(const char *szTopic, const char *szName, char *szFormat, ...)
{
 const char *pp=FindPar(szTopic,szName);
 if(!pp)return 0;
 va_list arg_ptr;
 char szB[MAX_LINELENGTH+1];

 va_start(arg_ptr, szFormat);
 if( strchr(szFormat,'S')!=0)
	{char *s=strchr(szFormat,'%'),*e=strchr(szFormat,'S');
	 if(!s || !e)return 0;
	 int i;
         if(s+1 == e)i=0;
	 else i=atoi(s+1);
	 if(i<0)return 0;
	 s=strchr(pp,'=');
         if(!s)return 0;
         if(i==0)i=(int)strlen(s);
	 if(i<(int)strlen(s))*(s+i)=0;
	 strncpy(va_arg(arg_ptr,char *),s+1,i);
	 va_end(arg_ptr);
	 return 1;
	}
 sprintf(szB,"%s=%s",szName,szFormat);
 vsscanf(pp, szB, arg_ptr);
 va_end(arg_ptr);

 return 1;
}
// **************************************************
const char *LineString::FindPar(const char *szTopic, const char *szName,
		                const int iWithDel)
{int i,j,k=0;
//const char *p;
regex_t re;
regcomp( &re, "^$", REG_EXTENDED );

if(szTopic)
  {for(i=0;i<nLines;i++)
      {if(Lines[i]->GetBuf()[0]=='#' || Lines[i]->GetBuf()[0]==';'
          || (regexec( &re, Lines[i]->GetBuf(), (size_t)0, (regmatch_t *)NULL, 0)==0)  )continue;
     int l1=strlen(Lines[i]->GetBuf());
     int l2=strlen(szTopic);
     if(!strncmp((Lines[i]->GetBuf()),szTopic,(l1<l2 ? l1 :l2))){k=i;break; }
//     if(!strncmp((Lines[i]->GetBuf()+1),szTopic,(l1<l2 ? l1 :l2)))break;
//       if(strstr(Lines[i]->GetBuf(),szTopic))break; 
      }
    
   if(i==nLines)return 0;
  }
else i=0;

 for(j=i;j<nLines;j++)
    {if(Lines[j]->GetBuf()[0]=='#' || Lines[j]->GetBuf()[0]==';'
         || (regexec( &re, Lines[j]->GetBuf(), (size_t)0, (regmatch_t *)NULL, 0)==0)  )continue;
    if(Lines[j]->GetBuf()[0]=='[' && j>k )return 0; // stop at next header
//     if((p=strstr(Lines[j]->GetBuf(),szName)))break;
     int l1=strlen(Lines[j]->GetBuf());
     int l2=strlen(szName);
     if(!strncmp(Lines[j]->GetBuf(),szName,(l1<l2 ? l1 :l2)))break;
    }

 if(j==nLines)return 0;

 if(!iWithDel)return Lines[j]->GetBuf();
 
 Result.Set(Lines[j]->GetBuf());
 Result.Add((const char *)DeLim);
 return (const char *)Result;
}
//***************************************************
int LineString::Write(const char *szFile) const
{ FILE *f=fopen(szFile,"w");
 if(!f)
   {PRINT_DEBUG("Can not open %s\n",szFile);
    perror("Reason");
    return 0;
   }
 int i=Write(f);
 fclose(f);
 return i;
}
// **************************************************
int LineString::Write(FILE *f) const
{if(nLines<=0)return 0;
 int i;
 for(i=0; i<nLines; i++)
 fprintf(f,"%s%s",Lines[i]->GetBuf(),DeLim.GetBuf());
 return i-1;
}
// **************************************************
int LineString::Write(String &S, const int iAdd) const
{if(nLines<=0)
   {if(iAdd)S.Set("No text available\n");
    else S.Set("No text available\n");
   }

 int i;
 for(i=0; i<nLines; i++)
    if(i==0 && !iAdd)S.Setf("%s%s",Lines[i]->GetBuf(),DeLim.GetBuf());
    else S.Addf("%s%s",Lines[i]->GetBuf(),DeLim.GetBuf());
 return i-1;
}
// **************************************************
#ifdef STRINGS_TEST
main(int iArgC, char ** szArgV)
{LineString LS;//("ABC\ndef\n123\n");
  LS.Write(stdout); 
  fprintf(stdout,"---------------\n");
  fprintf(stdout,"%d Lines\n",LS.GetNLines());
  int i;
  for(i=0;i<LS.GetNLines();i++)fprintf(stdout,"Line %d:%s\n",i,LS.Line(i));
 
  fprintf(stdout,"Replace 2---------------\n");
  const char *s={"NewLine"};
  LS.ReplaceLine(s,2);
  for(i=0;i<LS.GetNLines();i++)fprintf(stdout,"Line %d:%s\n",i,LS.Line(i));

  FILE *f=fopen(szArgV[1],"r");
  if(!f){fprintf(stderr,"Error opening file\n"); exit(EXIT_FAILURE);}
  
//   i=LS.Read(f," ",LineString::S_NUMBER);
//    if(!i)
//      {fprintf(stderr,"Erorr reading file\n");
//       exit(EXIT_FAILURE);
//      }
//    fprintf(stdout,"%d columns\n",i);
//    LS.Write(stdout);


// while(1)
//    {i=LS.Read(f,"@EOH",LineString::S_CONTAINS);
//     fprintf(stdout,"Until %d\n",i);
//     if(!i)break;
// 
//     i=LS.Read(f," ",LineString::S_ISCOL_BREAK);
//     if(!i)
//       {fprintf(stderr,"Error reading file\n");
//        exit(EXIT_FAILURE);
//       }
//     fprintf(stdout,"Result %d\n",i);
//     if(!i)break;  
//     LS.Write(stdout);
//    }
    int j,iset=4;
      for(j=1;j<iset;j++)
         {i=LS.Read(f,"@EOH",LineString::S_CONTAINS);
          fprintf(stdout,"Until %d\n",i);
         }
      if(!i)exit(EXIT_FAILURE);
      if(iset>1)i=LS.Read(f," ",LineString::S_ISCOL_BREAK);
      fprintf(stdout,"i=%d\n",i);
      LS.Write(stdout);
      i=LS.Read(f," ",LineString::S_NUMBER);
      fprintf(stdout,"i=%d\n",i);
      LS.Write(stdout);
fclose(f);
//   const char *p=LS.Find(i,"SampleSensor",1);
//   fprintf(stdout,"found %d:%s\n",i,p);
//  
//   fprintf(stdout,"Par: %s\n",LS.FindPar("[FileHeader]","NoofCols"));
//   double d=0;
//   LS.GetPar("[FileHeader]","DataPoints","%lf",&d);
//   fprintf(stdout,"Par: %f\n",d);
 
}
#endif

// ***************************************************
