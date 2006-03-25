//File:dfile.cpp
//$Log: dfile.orig.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.3  1999/03/15 09:08:37  herbie
//*** empty log message ***
//
//$Id: dfile.orig.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>

#ifndef STDINC_H
#include "stdinc.h"
#endif

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

#ifndef THEFUNC_H
#include "thefunc.h"
#endif

#ifndef CDATA_H
#include "cdata.h"
#endif

#ifndef FILEPAR_H
#include "filepar.h"
#endif

#define DFILE_MODULE
#include "dfile.h"


const char * TheTYPES[] =
 { "THE: None", "THE: RCP", "THE: CapV1", "THE: Membrane", "THE: CapOld"};

const char * TransTYPES[]=
 { "TRANSPORT: None", "TRANSPORT: Datap", "TRANSPORT: Measure"};

const char * XDifTYPES[]=
 { "XDIF: None", "XDIF: raw", "XDIF: dif (peak)"};

const char * XSxsTYPES[]=
 { "SXS: None", "SXS: old", "SXS: new"};

const char szNoText[]={"No more textlines available"};

const char *szLineType[]={"None","DOS","UNIX"};

const SXS_PAR SxSFile::SPar=SxSFile::ReadPars(PARAMETER_FILE);

 // ******************************************************
 // General Functions
 // ******************************************************
int CheckFileType(const char *szFile, enum EndLine &LineT )
{FILE *f;
 f=fopen(szFile,"r");
 CHECK_FILE_POINTER_RETURN(f,szFile,FtNONE)

 LineT=NOEnd;

   //int iSupType=FtILLG;
 char szB[MAX_LINELENGTH+1];

 // Get 1. line
 fgets(szB,MAX_LINELENGTH-1,f);

 szB[MAX_LINELENGTH]=0;
 int iL=strlen(szB);
 if(iL==MAX_LINELENGTH ||iL<=0)return FtERROR;
 
 char *s=&szB[iL-2];
 if(strcmp(s,"\r\n")==0)LineT=DOSEnd;
 else
   {if(*(s+1)=='\n')LineT=UNIXEnd;}
   
 RemoveCR_LF(szB);
 szB[MAX_LINELENGTH]=0;
 if(strncmp(szB,"RAW",3)==0) {fclose (f);return FtXDIF_RAW;}
 if(strncmp(szB,"PEAK",4)==0){ fclose (f);return FtXDIF_DIF;}
 if(strncmp(szB,"SPLINETABLE",11)==0){ fclose (f);return FtSPLINE;}
 if(strncmp(szB,"[SXSFileHeader]",15)==0){ fclose (f);return FtSXS_NEW;}
 
 if(szB[0]=='{')
   { // Looks like Martins RCP-Files
     do
       {
        //Get 2. Line
	 if( !fgets(szB,MAX_LINELENGTH-1,f) ){fclose (f);return FtERROR;}
       }
     while(szB[0]!='}');
     
     // Get next lines
     while(fgets(szB,MAX_LINELENGTH-1,f))
       {RemoveCR_LF(szB);
       if(!IsNumString(szB)){fclose (f);return FtILLG;}
       }  
      fclose (f);
      return FtTHE_RCP;
   }

 // 1.Line
 if(strcmp(szB,"[FileHeader]")==0)
   { //Newest file header
     fclose (f);
     return FtTHE_CAPV1;
   }
 // 1.Line
 if(strncmp(szB,"File created",12)==0 || strncmp(szB,"File from",9)==0 ||
    strncmp(szB,"\"File created",13)==0 )
   { // DATAP Files or older THE-Files
     // Read 2. line
     if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
     RemoveCR_LF(szB);

     // Old Thermal Expansion Membran file
     if(strncmp(szB,"l0=",3)==0 || strncmp(szB,"L0=",3)==0 ||
	  strncmp(szB,"\"L0=",4)==0){fclose (f);return FtTHE_MEMBRAN;}

     // May be an older Version of THE Capacitor files
     if(szB[0]=='@')
       {//Get 3. line
	 if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
        RemoveCR_LF(szB);
        if(szB[0]=='@'){fclose (f);return FtTHE_CAPOLD;}
       fclose (f);
       return FtILLG;
       }
     
     // Now it may be a Transport DATAP-File
     // 2. Line
     if(IsNumString(szB)){fclose (f);return FtILLG;}
     if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
        RemoveCR_LF(szB);
	//3. line 
        if(IsNumString(szB)){fclose (f);return FtTRANS_DATAP;}
        return FtILLG;
   }
 // 1.Line
 if(IsNumString(szB)){fclose (f);return FtILLG;}
 //Check for Transport Measurement or He3 measurement or different

 // 2.Line
 if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
 RemoveCR_LF(szB);
 if(IsNumString(szB))
   {// Check for SXS-Data file:
    // 1.Line Text
    // 2.Line (uint) Startposition (SXS_POS_MIN ... SXS_POS_MAX)
    unsigned uSt,uEn,u,uDe;
    sscanf(szB,"%u",&uSt);
    if( ! (uSt>SXS_POS_MIN && uSt<=SXS_POS_MAX) ){fclose (f);return FtILLG;}

    // 3.Line (uint) Endposition (SXS_POS_MIN ... SXS_POS_MAX)

    if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
    RemoveCR_LF(szB);
    if(!IsNumString(szB)){fclose (f);return FtILLG;}
    sscanf(szB,"%u",&uEn);
    if( ! (uEn>=SXS_POS_MIN && uEn<=SXS_POS_MAX) ){fclose (f);return FtILLG;}
    if(uSt>=uEn){fclose (f);return FtILLG;}

    // 4.Line (uint) No of Steps (multiple of 5)
    if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
    RemoveCR_LF(szB);
    if(!IsNumString(szB)){fclose (f);return FtILLG;}
    sscanf(szB,"%u",&u);
    if( ! (u>SXS_POS_MIN && u<=SXS_POS_MAX)){fclose (f);return FtILLG;}
    if( u%5 || (uEn-uSt)%u ){fclose (f); return FtILLG;}

    // 5.Line (uint) Delta
    if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
    RemoveCR_LF(szB);
    if(!IsNumString(szB)){fclose (f);return FtILLG;}
    sscanf(szB,"%u",&uDe);
    if( (uEn-uSt)/u != uDe){fclose (f); return FtILLG;}

    int i;
    // 6.Line (uint) No of Spectra
    if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
    RemoveCR_LF(szB);
    if(!IsNumString(szB)){fclose (f);return FtILLG;}
    sscanf(szB,"%d",&i);
    if (i<=0){fclose (f);return FtILLG;}

    float fT;
    // 7.Line (float) Timebase (sec)
    if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
    RemoveCR_LF(szB);
    if(!IsNumString(szB)){fclose (f);return FtILLG;}
    sscanf(szB,"%f",&fT);
    if(fT<=0){fclose (f);return FtILLG;}

    // 8.Line (int) Max Counts
    if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
    RemoveCR_LF(szB);

    // next lines: No of seps intensity values (float) in blocks of 5 vals

   unsigned iu;
   for(iu=0; iu<u; iu+=5)
      {if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
       if(CountCols(szB," ")!=5){fclose (f);return FtILLG;}
      }
   fclose (f);
   return FtSXS_OLD;
  }//Check for SXS data file


 // 3.Line
 if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
 RemoveCR_LF(szB);
 if(IsNumString(szB)){fclose (f);return FtTRANS_MEAS;}

 // 4.Line
 if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f); return FtERROR;}
 RemoveCR_LF(szB);
 if(IsNumString(szB)){fclose (f); return FtILLG;}

 // 5.Line
 if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f); return FtERROR;}
 RemoveCR_LF(szB);
 if(!IsNumString(szB)){fclose (f); return FtILLG;}
 if(CountCols(szB)==12) {fclose (f); return FtHE3_MEAS;}
 fclose (f);
 return FtILLG;
}
// ******************************************************
int GetTextLines(FILE *f, const char * szDelim, char* szA)
{ 
 // return 1: a textline is read
 // return 0: when EOF is found
 // return -#: # columns no text in line
 // iStart=1: initializes Line count
 // szLine: Text Line  Buffer
 // szLine=NULL: no Text is returned

  char szB[MAX_LINELENGTH+1];
  char szH[MAX_LINELENGTH+1];

 if(!fgets(szB,MAX_LINELENGTH-1,f))return 0;
 RemoveCR_LF(szB);
 int j=0;
 while(szB[j])
   {if(szB[j]==' '){j++; continue;}
       if(strchr("+-0123456789.Ee",szB[j])==NULL)
	 {
          if(szA)strcpy(szA,szB);
          return 1;
         }
     j++;
    }//while
  if(j==0 && szB[j]==0)return 1; //Leerzeile
  
  strcpy(szH,szB);
  char *s=strtok(szH,szDelim);
   if(s==NULL)
     {PRINT_DEBUG("Unexpected Error analyzing file\n");
      return 0;
     }
   int iColC=1;
   while(strtok(NULL,szDelim)!=NULL)iColC++;
   if(iColC<=1)
      {PRINT_DEBUG("GetTextLines: # of columns (%d) < 2\n",iColC);
       return 0;
      }
  if(szA)strcpy(szA,szB);
  return -iColC;
}
// ******************************************************
int CheckAsciiFile(const char *szFile, int &nTextL, int &nDataL,
                   const int CheckMode)
{
 nTextL=0; nDataL=0;
 FILE *f=fopen(szFile,"r");
 CHECK_FILE_POINTER_RETURN(f,szFile,0)

 int iTxt=0,iData=1;
 int i=0;

 while ( (i=GetTextLines(f," ",NULL))>0)iTxt++;

 if(i==0){fclose(f); return -1;}

 int j;
 while( (j=GetTextLines(f," ",NULL)) )
   {if(j!=i )
      {if(CheckMode==STRICT){iData*=-1;break;}
       else {PRINT_DEBUG("Warning: Illegal data line %d in %s\n",iData,szFile); 
	    } 
      }
    iData++;
   }

 //fprintf(stderr,"File has %d columns\n",-i);
 nTextL=iTxt;
 nDataL=iData;
 fclose(f);
 return -i;
}
// ******************************************************
// DataFile Functions
// ******************************************************
DataFile::DataFile(void)
{
 TD=NULL;
 nSets=0;

 szCID=NULL;
 nCols=0;

 iFType=FtILLG;
 iFGroup=NoGROUP;
 iFLEnd=NOEnd;
 strcpy(szID,"None");
 szName[0]=0;
 nTxtLines=0;

}
// ******************************************************
int DataFile::Init(const char *szFile,
                   const int iG, const int iT, const char *szI)
{
 strncpy(szID,szI,MAX_IDLENGTH);
 szID[MAX_IDLENGTH]=0;

 strncpy(szName,szFile,MAXPATH);
 szName[MAXPATH]=0;

 iFGroup=(enum FileGroup)iG;
 iFType= (enum FileType) (iG+iT);

 if(!CheckType())
   {PRINT_DEBUG("Illegal File Group (%d) or File Type (%d)\n",iG,iT)
    return 0;
   }

 return 1;

}
// ***********************************************************
int DataFile::CheckType(void)
{switch (iFGroup)
   {case TheGROUP:if(iFType>TheGROUP && iFType<=FtTHE_CAPOLD)return 1;
                  break;
    case SxSGROUP:if(iFType>SxSGROUP && iFType<=FtSXS_NEW)return 1;
                  break;
    case XDifGROUP:if(iFType>XDifGROUP && iFType<=FtXDIF_DIF)return 1;
                  break;
   case TransGROUP:if(iFType>TransGROUP && iFType<=FtTRANS_MEAS)return 1;
                  break;
   case He3GROUP:if(iFType>He3GROUP &&iFType<=FtHE3_MEAS)return 1;
                  break;
   case SplineGROUP:if(iFType==FtSPLINE)return 1;
                 break;
   case NoGROUP:if(iFType==FtASCII)return 1;
                break;
   }
 iFGroup=NoGROUP;
 iFType=FtERROR;
 return 0;
}
// ******************************************************
int DataFile::AllocColID(const int iC)
{
 DeleteColID();
 if(iC>0)
   {szCID = new char* [iC];
    CHECK_POINTER_RETURN(szCID,0)

#if defined (DFILE_DEBUG)
 if(fd){fprintf(fd,"Alloc szCID %p:%d\n",szCID,iC);
        fflush(fd);}
#endif
    int i;
    for(i=0;i<iC;i++)
       {szCID[i]=new char[MAX_IDLENGTH+1];
        CHECK_POINTER_RETURN(szCID[i],0)

#if defined (DFILE_DEBUG)
 if(fd){fprintf(fd,"Alloc szCID[%d]: %p, %d\n\n",i,szCID[i],MAX_IDLENGTH);
        fflush(fd);}
#endif
        }//for
 nCols=iC;
  }//if
return 1;
}
// ******************************************************
void DataFile::DeleteColID(void)
{if(nCols>0 && szCID!=NULL)
   {int i;
    for(i=0;i<nCols;i++)
       {
#if defined (DFILE_DEBUG)
        if(fd!=NULL){fprintf(fd,"Delete szCID[%d]: %p\n",i,szCID[i]);
                     fflush(fd);}
#endif
         if(szCID[i])delete szCID[i];
        }
#if defined (DFILE_DEBUG)
 if(fd!=NULL){fprintf(fd,"\n");
              fprintf(fd,"Delete szCID: %p\n\n",szCID);
              fflush(fd);
             }
#endif
      delete szCID;

    }
 szCID=NULL;
 nCols=0;
}
// ********************************************************
int DataFile::InsertColID(const char *szN, const int iCol)
{
    
 if( !szCID || !nCols)
   {PRINT_DEBUG("Cannot insert Column Id\n\n")
    return 0;
   }

 CHECK_INDEX_RETURN(iCol,nCols+1,0)
   
 char **p= new char * [nCols+1];
 CHECK_POINTER_RETURN(p,0)

 int i,j;
 for(i=0,j=0;i<=nCols;i++)
    {if(i==iCol+1)
       {p[i]=new char[strlen(szN)+1];
        CHECK_POINTER_RETURN(p[i],0)
        strcpy(p[i],szN);  
       }	
     else {p[i]=szCID[j]; j++;}
    }

 delete szCID;
 szCID=p;
 nCols++;

 return 1;

}
// ********************************************************
const char * DataFile::GetColID(int iC) const
{
 CHECK_INDEX_RETURN(iC,nCols,NULL)
 return szCID[iC];
}
// ********************************************************
int DataFile::SetColID(const char *szI, const int iC)
{
 CHECK_INDEX_RETURN(iC,nCols,0)
 strncpy(szCID[iC],szI,MAX_IDLENGTH);
 szCID[MAX_IDLENGTH]=0;
 return 1;
}
// ******************************************************
int DataFile::SetColRange(const int iCS, const int iCE)
{ CHECK_INDEX_RETURN(iCS,nCols,0)
  CHECK_INDEX_RETURN(iCE,nCols,0)
  if(iCE<=iCS)return 0;

  int iC=iCE-iCS+1;
  char **szNCID = new char* [iC];
    CHECK_POINTER_RETURN(szNCID,0)

  int i,j;
    for(i=iCS,j=0;i<=iCE;i++,j++)szNCID[j]=szCID[i];
 
  for(i=0;i<nCols; i++)
     {if(i<iCS ||i>iCE)delete szCID[i];}

  szCID=szNCID;
  nCols=iC;
  return nCols;
}
// ******************************************************
char *DataFile::SPrintInfo(char * szB) const
{
 char *sp=szB;

 sp += sprintf(sp,"File: %s (%s)\n",szName, szLineType[(int)iFLEnd]);

 if(TD)
   {sp+=sprintf(sp,"%d data points, %d columns\n",
                TD->GetSteps(),TD->GetCols());
   }
 else {sp+=sprintf(sp,"NO data columns\n");}
 *sp=0;
 return sp;
 // return number of printed chars
}
// *****************************************************
void DataFile::SPrintInfo(String &B) const
{
 B.Addf("File: %s (%s)\n",szName, szLineType[(int)iFLEnd]);

 if(TD)
   B.Addf("%d data points, %d columns\n",TD->GetSteps(),TD->GetCols());
 else B.Add("NO data columns\n");
}
// *****************************************************
char *DataFile::SPrintColInfo(char * szB) const
{ char *sp=szB;

 if(iFType==FtERROR || iFType==FtILLG)
   {sp+=sprintf(sp,"Illegal file type %d; cannot print data\n",iFType);
    return sp;
   }

 if(!TD)
   {sp+=sprintf(sp,"No data to print\n");
    return sp;
   }

  TD->NewMinMax();

  for(int j=0;j<TD->GetCols();j++)
     {if(j==TD->GetColX())sp+=sprintf(sp,"X:");
	   else {if(j==TD->GetColY())sp+=sprintf(sp,"Y:");
	           else {if(j==TD->GetColZ())sp+=sprintf(sp,"Z:");
	                  else sp+=sprintf(sp,"  ");
                         }
                 }
      sp+=sprintf(sp,"Column %d:%-14.14s",j+1,(szCID?szCID[j]:" "));
      sp+=sprintf(sp,"%12.5g ... %12.5g\n",TD->GetColMin(j),TD->GetColMax(j));
     }

 // return number of printed chars
 *sp=0;
 return sp;
}
// *****************************************************
void DataFile::SPrintColInfo(String &B) const
{
 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d; cannot print data\n",iFType);
    return;
   }

 if(!TD)
   {B.Add("No data to print\n");
    return;
   }

  TD->NewMinMax();

  for(int j=0;j<TD->GetCols();j++)
     {if(j==TD->GetColX())B.Add("X:");
	   else {if(j==TD->GetColY())B.Add("Y:");
	           else {if(j==TD->GetColZ())B.Add("Z:");
	                  else B.Add("  ");
                         }
                 }
      B.Addf("Column %d:%-14.14s",j+1,(szCID?szCID[j]:" "));
      B.Addf("%12.5g ... %12.5g\n",TD->GetColMin(j),TD->GetColMax(j));
     }
}
// *****************************************************
int DataFile::SaveData(const char *szText, const char *szFile)
{
  if(!TD || iFType==FtERROR || iFType==FtILLG)return 0;


  FILE *f=fopen(szFile,"w");
  CHECK_FILE_POINTER_RETURN(f,szFile,0)
  int i=SaveData(szText, f);
  fclose(f);
  strncpy(szName,szFile,MAXPATH);
  return i;
}
// ******************************************************
int DataFile::SaveData(const char *szText, FILE *f)
{ 
  SetEndChr();
  fprintf(f,"%s%s",szText,szE);
  TD->ExportData(f,iFLEnd);
  return 1;
}
// *****************************************************
int DataFile::SaveSplineTable(FILE *f, const int iX, const int iY)
{
  if(!TD || iFType==FtERROR || iFType==FtILLG || TD->GetSteps()<=0)return 0;

  CHECK_INDEX_RETURN(iX,TD->GetCols(),0);
  CHECK_INDEX_RETURN(iY,TD->GetCols(),0);

  SetEndChr();
  if(!f)return 0;

  fprintf(f,"SPLINETABLE from %s (x:%d, y:%d)%s",szName,iX+1,iY+1,szE);
  fprintf(f,"%d%s",TD->GetSteps(),szE);

  int i,j;
  double lfx,lfy;
  for(i=0;i<TD->GetSteps();i++)
     {for(j=0;j<TD->GetCols();j++)
         {if(j==iX)lfx=(*TD)[j][i];
          if(j==iY)lfy=(*TD)[j][i];
         }// for j
      fprintf(f,"%13.7g  %13.7g%s",lfx,lfy,szE);
     }//for i

  return 1;

}		
// *****************************************************		
void DataFile::SetEndChr(void)
{ strcpy(szE,"\n");
#if defined(LINUX)
  if(iFLEnd==NOEnd || iFLEnd==DOSEnd)strcpy(szE,"\r\n");
#endif 
}
// *****************************************************
DataFile::~DataFile()
{
 if(TD)delete TD;
 DeleteColID();
}

// ********************************************************
// TheFile
// ******************************************************
 // Header of new Data File:
 // [FileHeader]
 //Idn=1.0:THE (orig.meas.data)
 //Cell=N
 //p_Cell=10
 //p_Heater=20
 //H_Field=0.0
 //KMon=0.3097
 //DMon=-0.2664
 //HeliPot=4.2
 //CalConst=926.77
 //CalFile=NULL
 //SampleSensor=GaAs9257
 //HeaterSensor=9243
 //L0=1.23
 //SampleID=test 3 col
 //K0=0
 //;
 //StartTime=Fri Mar 01 15:42:19 1996
 //EndTime=Fri Mar 01 15:52:22 1996
 //NoofCols=#
 //DataPoints=6
 //ColID1=Tsample[K]
 //ColID2=C[pF]
 //ColId3=Tsample[V]
 //...
 //ColId#=LastV[??]
 //@EOH
 // N Data Points C1 C2 C3 ... C# separeted by blancs
// **********************************************************************
// **********************************************************************
TheFile::TheFile(const char *szFile) : DataFile()
{
  //Inital values will be overwritten if read from file
  Head.szPCell[0]=Head.szPHeat[0]=0;
  Head.fMField=Head.fL0=Head.fKMon=Head.fDMon=-1;
  Head.fHeliPot=Head.fCalConst=Head.lfK0=-1.;
  Head.nCols=0;
  Head.cInsert='?';
  Head.cCell=0;
  Head.szCalFile[0]=0;
  Head.szSSens[0]=0;
  Head.szHSens[0]=0;
  Head.szSTime[0]=0;
  Head.szETime[0]=0;
  Head.nSteps=0;

  int iT=CheckFileType(szFile,iFLEnd);

  if(iT==FtERROR || iT==FtILLG)
    {iFType=(enum FileType)iT;
     iFGroup=NoGROUP;
     PRINT_DEBUG("Error analyzing file %s\n",szFile)
     return;
    }

  DataFile::Init(szFile,TheGROUP,iT-TheGROUP,TheTYPES[iT-TheGROUP]);
  nTxtLines=1;
  nSets=1;
}
// ********************************************************
TheFile::~TheFile() { }
// ********************************************************
int TheFile::ReadHeader(FILE *InStream)
{// Read Header from InStream
 // return 0: Illegal Header
 //        1: Ok
// ******************************************************
// Header of new Data File see above
// SOME header example
// File from THECAP Version 1.0 on Sat Dec 10 11:07:01 1994
// @TSlope=0.002000, Cell=A, p(Cell)=5, p(Heater)=5
// @H=2.5, CF=FILNAME ,  L0=2.970000
// Text Specifing the sample
// N Data Points X Y  separeted by blancs

 if(iFType==FtERROR || iFType==FtILLG || iFGroup != TheGROUP)return 0;

 char szB[MAX_LINELENGTH+1];
 int i;
 char *c;

 switch (iFType)
   {//case   FtERROR: return 0;
    //case    FtILLG: return 0;
    //case   FtASCII: return 0;
   case   FtTHE_RCP:
   case FtTHE_CAPV1:
       {fgets(Head.szHead,MAX_LINELENGTH-1,InStream);
        RemoveCR_LF(Head.szHead);
        if(!(strstr(Head.szHead,"File")!=NULL || Head.szHead[0]=='{'))
	       return 0;
        else
         {if(Head.szHead[0]=='{')
            {fgets(Head.szHead,MAX_LINELENGTH-1,InStream);
             RemoveCR_LF(Head.szHead);
            }
         }

        if(strcmp(Head.szHead,"[FileHeader]")) return 0;
        int iErr=0;
        rewind(InStream);
        FilePar *FP=new FilePar(InStream,F_NCLOSE);
        if(FP==0)return 0;
        if(FP->GetStatus()!=FOUND){delete FP;return 0;}
        FP->GetPar("FileHeader","Insert","%c",&(Head.cInsert));
        if(FP->GetPar("FileHeader","Idn","%79S",Head.szHead)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","Cell","%c",&(Head.cCell))!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","p_Cell","%19S",Head.szPCell)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","p_Heater","%19S",Head.szPHeat)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","T_Slope","%f",&(Head.fTSlope))!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","H_Field","%f",&(Head.fMField))!=FOUND)iErr=1;

        FP->GetPar("FileHeader","K0","%lf",&(Head.lfK0));
        FP->GetPar("FileHeader","KMon","%f",&(Head.fKMon));
        FP->GetPar("FileHeader","DMon","%f",&(Head.fDMon));
        FP->GetPar("FileHeader","HeliPot","%f",&(Head.fHeliPot));
        FP->GetPar("FileHeader","CalConst","%f",&(Head.fCalConst));

        if(FP->GetPar("FileHeader","CalFile","%s",Head.szCalFile)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","SampleSensor","%s",Head.szSSens)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","HeaterSensor","%s",Head.szHSens)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","L0","%f",&(Head.fL0))!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","SampleID","%79S",Head.szText)!=FOUND)iErr=1;

        if(FP->GetPar("FileHeader","StartTime","%25S",Head.szSTime)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","EndTime","%25S",Head.szETime)!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","NoofCols","%i",&(Head.nCols))!=FOUND)iErr=1;
        if(FP->GetPar("FileHeader","DataPoints","%i",&Head.nSteps)!=FOUND)iErr=1;

        if(!AllocColID(Head.nCols))return 0;
        for(i=0; i<Head.nCols; i++)
           {sprintf(szB,"ColID%d",i+1);
	        if(FP->GetPar("FileHeader",szB,"%15S",szCID[i])!=FOUND)iErr=1;
           }
         delete FP;
         if(iErr)return 0;

        if(iFType==FtTHE_RCP)fgets(szB,MAX_LINELENGTH-1,InStream);
        return Head.nSteps;
       }

  case FtTHE_MEMBRAN:
     fgets(Head.szHead,MAX_LINELENGTH-1,InStream);
     RemoveCR_LF(Head.szHead);
     fgets(szB,MAX_LINELENGTH-1,InStream);
     RemoveCR_LF(szB);
     if( (c=strstr(szB,"l0=")) !=0 )sscanf(c,"l0=%f",&(Head.fL0));
     if( (c=strstr(szB,"L0=")) !=0 )sscanf(c,"L0=%f",&(Head.fL0));
	 Head.fMField=0.0;
     c=strchr(szB,'>');
     if(c!=0)sprintf(Head.szText,"%s",c+1);
     Head.nCols=2;
     if(!AllocColID(Head.nCols))return 0;
     for(i=0; i<Head.nCols; i++)sprintf(szCID[i],"Col%d",i+1);
     return 1;
  case FtTHE_CAPOLD:
     fgets(Head.szHead,MAX_LINELENGTH-1,InStream);
     RemoveCR_LF(Head.szHead);
     fgets(szB,MAX_LINELENGTH-1,InStream);
     RemoveCR_LF(szB);
     sscanf(szB,"@TSlope=%f, Cell=%c, p(Cell)=%3s, p(Heater)=%3s",
		           &(Head.fTSlope),&(Head.cCell),Head.szPCell,Head.szPHeat);
     fgets(szB,81,InStream);
     RemoveCR_LF(szB);
	 switch (szB[0])
		{case '@':sscanf(szB,"@H=%f, CF=%s , L0=%f",
				 &(Head.fMField),Head.szCalFile,&(Head.fL0));
			  fgets(Head.szText,MAX_LINELENGTH-1,InStream);
			  //return 1;
			  break;
		 default:sscanf(szB,"L0=%f",&(Head.fL0));
			 Head.fMField=0.0;
			 break;
                }
      Head.nCols=2;
     if(!AllocColID(Head.nCols))return 0;
     for(i=0; i<Head.nCols; i++)sprintf(szCID[i],"Col%d",i+1);
      return 1;

    default: return 0;

   }
 return 0;
}
// ******************************************************
int TheFile::ReadData(const int iSet)
{// Reads complete data from file
 // return -1 : Error opening file
 //         0 : Illegal Header
 //         1 : Ok

 if(iFType==FtERROR || iFType==FtILLG || iFGroup != TheGROUP)return 0;
 // Multiple sets not possible
 if(iSet!=1)
   {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
    return 0;}

 FILE *f;
 char szB[MAX_LINELENGTH+1];
 int i,j=0;

  if ((f=fopen(szName,"r")) == NULL)
    {PRINT_DEBUG("Error opening file %s\n",szName)
     iFType=FtERROR;
     iFGroup=NoGROUP;
     return -1;
    }

  i=ReadHeader(f);
  if(!i)
    {fclose(f);
     PRINT_DEBUG("Error reading file header %s\n",szName)
     iFType=FtERROR;
     iFType=FtERROR;
     iFGroup=NoGROUP;
     return 0;
    }
  if(i<=1)
    {while (fgets(szB,MAX_LINELENGTH-1,f))j++;
    }
  else j=i;

 TD=new CRData(Head.nCols,i<=1?j:Head.nSteps);
 CHECK_POINTER_RETURN(TD,-1)

 rewind(f);
 ReadHeader(f);

 if(TD->ReadCols(f,j,Head.nCols)==0)
   {PRINT_DEBUG("Error reading columns from file %s\n",szName)
    iFType=FtERROR;
    iFGroup=NoGROUP;
    return 0;
   }
 fclose(f);
 return 1;
}
// *******************************************************
void TheFile::SetFileTime(const int iStartEnd)
{ long lSec_Now;
  char *cpDate;
  time(&lSec_Now);
  cpDate = ctime(&lSec_Now);
  if(iStartEnd & SET_START)sprintf(Head.szSTime,"%s",cpDate);
  if(iStartEnd & SET_END)sprintf(Head.szETime,"%s",cpDate);
}
// *******************************************************
char *TheFile::SPrintInfo(char * szB, const unsigned Flag)
{
 char *sp=szB;

 if(iFType==FtERROR || iFType==FtILLG || iFGroup != TheGROUP)
   {sp+=sprintf(sp,"Illegal file type %d ;cannot print data\n",iFType);
    return sp;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(sp);

 if(Flag & PR_TEXT)
   sp+=sprintf(sp,"Text: %s\n",Head.szText);
 if(Flag & PR_TYPE)
   sp+=sprintf(sp,"       File Type: %s (%1d)\n",szID,iFType);

 if(Flag & PR_HEAD)
   {sp+=sprintf(sp,"              L0: %5.3f [mm]\n",Head.fL0);

   if(iFType==FtTHE_RCP || iFType==FtTHE_CAPV1)
     {sp+=sprintf(sp,"         T Slope: %6.5f [K/s]\n",Head.fTSlope);
      sp+=sprintf(sp,"         p(Cell): %s\n",Head.szPCell);
      sp+=sprintf(sp,"       p(Heater): %s\n",Head.szPHeat);
      sp+=sprintf(sp,"         Cell ID: %c       \n",Head.cCell);
      sp+=sprintf(sp,"  Mag. Field: %6.2f [Tesla]\n",Head.fMField);
      //sp+=sprintf(sp,"     Monitor: %6.4f [mV/A]\n",Head.fKMon);
      //sp+=sprintf(sp," Heli Pot Setting: %4.2f\n",Head.fHeliPot);
      //sp+=sprintf(sp,"Field Calibration: %8.3f [Gauss/A]\n",Head.fCalConst);
      sp+=sprintf(sp," Calibration File:%s\n",Head.szCalFile);
      sp+=sprintf(sp,"    Sample Sensor:%s\n",Head.szSSens);
      sp+=sprintf(sp,"    Heater Sensor:%s\n",Head.szHSens);
     }
   }
  if(Flag & PR_TIME && iFType==FtTHE_CAPV1)
    {
	 sp+=sprintf(sp,"Start Time: %s\n",Head.szSTime);
	 sp+=sprintf(sp,"  End Time: %s\n",Head.szETime);
    }//if

 if(Flag & PR_COLS)sp=DataFile::SPrintColInfo(sp);
 // return number of printed chars
 *sp=0;
 return sp;
}
// ******************************************************
void TheFile::SPrintInfo(String &B, const unsigned Flag)
{
 
 if(iFType==FtERROR || iFType==FtILLG || iFGroup != TheGROUP)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(Flag & PR_BASIC)
   DataFile::SPrintInfo(B);

 if(Flag & PR_TEXT)
   B.Addf("Text: %s\n",Head.szText);
   
 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",szID,iFType);

 if(Flag & PR_HEAD)
 {B.Addf("              L0: %5.3f [mm]\n",Head.fL0);

  if(iFType==FtTHE_RCP || iFType==FtTHE_CAPV1)
    {B.Addf("         T Slope: %6.5f [K/s]\n",Head.fTSlope);
     B.Addf("         p(Cell): %s\n",Head.szPCell);
     B.Addf("       p(Heater): %s\n",Head.szPHeat);
     B.Addf("         Cell ID: %c       \n",Head.cCell);
     B.Addf("  Mag. Field: %6.2f [Tesla]\n",Head.fMField);
     //B.Addf("     Monitor: %6.4f [mV/A]\n",Head.fKMon);
     //B.Addf(" Heli Pot Setting: %4.2f\n",Head.fHeliPot);
     //B.Addf("Field Calibration: %8.3f [Gauss/A]\n",Head.fCalConst);
     B.Addf(" Calibration File:%s\n",Head.szCalFile);
     B.Addf("	 Sample Sensor:%s\n",Head.szSSens);
     B.Addf("	 Heater Sensor:%s\n",Head.szHSens);
     B.Addf("Start Time: %s\n",Head.szSTime);
     B.Addf("  End Time: %s\n",Head.szETime);
    }//if
  }  
 if(Flag & PR_TIME && iFType==FtTHE_CAPV1)
     {
      B.Addf("Start Time: %s\n",Head.szSTime);
      B.Addf("  End Time: %s\n",Head.szETime);
     }//if
 if(Flag & PR_COLS)DataFile::SPrintColInfo(B);
}
// ******************************************************
int TheFile::SaveHeader(FILE *OutStream, const int nCS, const int nCE,
                                         const int nR)
{
 SetEndChr();

 switch(iFType)
   {case FtTHE_RCP:
    case FtTHE_CAPV1:
      {fprintf(OutStream,"[FileHeader]%s",szE);
       fprintf(OutStream,"Idn=%s%s",Head.szHead,szE);
       fprintf(OutStream,"Insert=%c%s",Head.cInsert,szE);
       fprintf(OutStream,"Cell=%c%s",Head.cCell,szE);
       fprintf(OutStream,"p_Cell=%-19.19s%s",Head.szPCell,szE);
       fprintf(OutStream,"p_Heater=%-19.19s%s",Head.szPHeat,szE);
       fprintf(OutStream,"T_Slope=%f%s",Head.fTSlope,szE);
       fprintf(OutStream,"H_Field=%f%s",Head.fMField,szE);
       fprintf(OutStream,"K0=%f%s",Head.lfK0,szE);

       fprintf(OutStream,"KMon=%f%s",Head.fKMon,szE);
       fprintf(OutStream,"DMon=%f%s",Head.fDMon,szE);
       fprintf(OutStream,"HeliPot=%f%s",Head.fHeliPot,szE);
       fprintf(OutStream,"CalConst=%f%s",Head.fCalConst,szE);

       fprintf(OutStream,"CalFile=%s%s",Head.szCalFile,szE);
       fprintf(OutStream,"SampleSensor=%s%s",Head.szSSens,szE);
       fprintf(OutStream,"HeaterSensor=%s%s",Head.szHSens,szE);
       fprintf(OutStream,"L0=%f%s",Head.fL0,szE);
       fprintf(OutStream,"SampleID=%s%s",Head.szText,szE);
       fprintf(OutStream,";%s",szE);


       if(strchr(Head.szSTime,'\n')==0)
	{fprintf(OutStream,"StartTime=%s%s",Head.szSTime,szE);
	 fprintf(OutStream,"EndTime=%s%s",Head.szETime,szE);
	}
       else
	{fprintf(OutStream,"StartTime=%s",Head.szSTime);
	 fprintf(OutStream,"EndTime=%s",Head.szETime);
	}

        int nCStart=0;
        if(nCS>0 && nCS>nCE)nCStart=nCS;

        int nCEnd=TD->GetCols();
        if(nCE<TD->GetCols() && nCE>nCStart)nCEnd=nCE;
        int nRR=(nR==-1?TD->GetSteps():nRR);  
       fprintf(OutStream,"NoofCols=%i%s",nCEnd-nCStart,szE);
       fprintf(OutStream,"DataPoints=%i%s",nRR,szE);
       int i;
       for(i=nCStart; i<nCEnd;i++)
           fprintf(OutStream,"ColID%d=%s%s",i-nCStart+1,szCID[i],szE);

       char szB[MAX_LINELENGTH];
       GetCurrentTime(szB);
       fprintf(OutStream,";XAusw: CreationTime=%s%s",szB,szE);
       fprintf(OutStream,"@EOH%s",szE);
      }
       return 1;
   default: return 0;
   }
}
// ********************************************************
int TheFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG || iFGroup != TheGROUP)return 0;
 
 FILE *f;

 char szB[MAXPATH+1];

 if(szFile)strncpy(szB,szFile,MAXPATH);
 else strcpy(szB,szName);
 f=fopen(szB,"w");
 CHECK_FILE_POINTER_RETURN(f,szB,0)

 strncpy(szName,szB,MAXPATH);
 int iRet=TheFile::SaveData(f);
 fclose(f);
 return iRet; 
}
// ********************************************************
int TheFile::SaveData(FILE *fOut)
{
 int i=SaveHeader(fOut);
 i|=TD->ExportData(fOut,iFLEnd);
 return i;
}
// ********************************************************
const char * TheFile::GetInfoText(const int iL)
{
 if(iL==-1)return Head.szText;
 else return szNoText;
}
// *******************************************************
int TheFile::SetInfoText(const char *szT, const int iL)
{if(iL==-1){strncpy(Head.szText,szT,MAX_LINELENGTH);
            return 1;
           }
 else return 0;
}
// *******************************************************
// AsciiFile
// *******************************************************
AsciiFile::AsciiFile(const char *szFile, const int nL,
                     const int CheckMode) : DataFile()
{
  szTextLine=0;
  iLStart=iLEnd=nLines=0;

   // nL==0 automatic analysis
  // nL>0 Reads max. nL textlines

  int nTxtL=0,nDataL=0;

  int nCs=CheckAsciiFile(szFile,nTxtL,nDataL,CheckMode);

  if(nCs==-1 || nDataL <= 0)
    {iFType=FtERROR;
     iFGroup=NoGROUP;
     PRINT_DEBUG("Error analyzing file %s in line %d\n",szFile,-nDataL)
     return;
    }
  nCols=nCs;
  nDataLines=nDataL;
  nLines=nTxtL;
  iLStart=1;
  iLEnd=nL;
  if(nL==0 || nL>nTxtL)iLEnd=nTxtL;

  DataFile::Init(szFile,NoGROUP,FtASCII,"ASCII File");

  iFType=FtASCII;
  iFGroup=NoGROUP;
  nSets=1;
}
// ********************************************************
// AsciiFile(const char *szFile, CRData *td, const int nL, char **szLs)
//{
//}
// ********************************************************
int  AsciiFile::ReadData(const int iSet)
{// Reads complete data from file
 // return -1 : Error opening file
 //         0 : Illegal Header
 //         1 : Ok

 if(iFType==FtERROR || iFType==FtILLG)return 0;
 // Multiple sets not possible
 if(iSet!=1)
  {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
   return 0;}
 FILE *f;
 int i;

  if ((f=fopen(szName,"r")) == NULL)
    {PRINT_DEBUG("Error opening file %s\n",szName)
     iFType=FtERROR;
     iFGroup=NoGROUP;
     return -1;
    }

 if(szTextLine)
   {for(i=0;i<iLEnd-iLStart+1;i++)delete szTextLine[i];
    delete szTextLine;
   }


 szTextLine= new char * [iLEnd-iLStart+1];
 CHECK_POINTER_RETURN(szTextLine,-1)

 for(i=0;i<iLEnd-iLStart+1;i++)
   {szTextLine[i]=new char [MAX_LINELENGTH+1];
    CHECK_POINTER_RETURN(szTextLine[i],-1)
   }

 for(i=0; i<nLines; i++)
    {if(i+1>=iLStart && i+1<=iLEnd)
       {if(!fgets(szTextLine[i-iLStart+1],MAX_LINELENGTH-1,f))return -1;
        RemoveCR_LF(szTextLine[i-iLStart+1]);
       }
    }

 if(TD)delete TD;

 TD=new CRData(nCols,nDataLines);
 CHECK_POINTER_RETURN(TD,-1)

 if(TD->ReadCols(f,nDataLines,nCols)==0)
   {PRINT_DEBUG("Error reading columns from file %s\n",szName)
    iFType=FtERROR;
    iFGroup=NoGROUP;
    return 0;
   }

 nLines=iLEnd-iLStart+1;
 nTxtLines=nLines;

 //DeleteColID();
 if(!AllocColID(nCols))return 0;
 for(i=0; i<nCols; i++)sprintf(szCID[i],"Col%d",i+1);


 fclose(f);
 return 1;
}
// *******************************************************
int  AsciiFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG)return 0;

 FILE *f;

 char szB[MAXPATH+1];

 if(szFile)strncpy(szB,szFile,MAXPATH);
 else strcpy(szB,szName);
 f=fopen(szB,"w");
 CHECK_FILE_POINTER_RETURN(f,szB,0)

 strncpy(szName,szB,MAXPATH);
 int iRet=AsciiFile::SaveData(f);
 fclose(f);
 return iRet;
 }
// *******************************************************
int  AsciiFile::SaveData(FILE *fOut)
{int j;
 SetEndChr();
 for(j=iLStart; j<=iLEnd; j++)
    {if(!fprintf(fOut,"%s%s",szTextLine[j-1],szE))return 0;}
 int i=TD->ExportData(fOut,iFLEnd);
 if(i==0)return 0;
 return i;
}
// *******************************************************
char *AsciiFile::SPrintInfo(char * szB, const unsigned Flag)
{char *sp=szB;

 if(iFType==FtERROR || iFType==FtILLG)
   {sp+=sprintf(sp,"Illegal file type %d ;cannot print data\n",iFType);
    return sp;
   }

 if(iFType==FtERROR || iFType==FtILLG)
   {sp+=sprintf(sp,"Illegal file type %d ;cannot print data\n",iFType);
    return sp;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(sp);

 int i;

 if(Flag & PR_TYPE)
   sp+=sprintf(sp,"       File Type: %s (%1d)\n",szID,iFType);
 if(Flag & PR_TEXT && szTextLine)
   {for(i=0;i<nLines;i++)
     sp+=sprintf(sp,"%s\n",szTextLine[i]);
    sp+=sprintf(sp,"... %d textlines...\n",nLines);

   }
 if(Flag & PR_COLS)sp=SPrintColInfo(sp);
 *sp=0;
 // returns number of printed chars
 return sp;
}
// ******************************************************
void AsciiFile::SPrintInfo(String &B, const unsigned Flag)
{
 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(B);

 int i;

 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",szID,iFType);
 if(Flag & PR_TEXT && szTextLine)
   {for(i=0;i<nLines;i++)
     B.Addf("%s\n",szTextLine[i]);
     B.Addf("... %d textlines...\n",nLines);
   }
 if(Flag & PR_COLS)DataFile::SPrintColInfo(B);
}
// ******************************************************
const char * AsciiFile::GetLine(const int iL) const
{CHECK_INDEX_RETURN(iL,nLines,szNoText)
 return szTextLine[iL];
}
// ******************************************************
int AsciiFile::SetLine(const char * szT, const int iL)
{CHECK_INDEX_RETURN(iL,nLines,0)
 strncpy(szTextLine[iL],szT,MAX_LINELENGTH);
 return 1;
}
// ******************************************************
AsciiFile::~AsciiFile()
{if(szTextLine)
   {int i;
    for(i=0;i<nLines;i++)
     {if(szTextLine[i])delete szTextLine[i];}
    delete szTextLine;
   }
}
// ******************************************************
// *******************************************************
// SxSFile
// *******************************************************
// ******************************************************
// Old Header SPEC_HEAD
//------------------------------------------------------
//        Text: A1 rhombo kV*mA 1, 4V*3,4A
//   Start Pos:27000
//     End Pos:32000
// No of steps:250
//       Delta:20
//No of spectra:2
// Timebase [s]:5.000000
//   Max Counts:0
// followed by intensity 5 values / line from Start Pos to End Pos in steps of
// Delta
//	 104            112             97             94             92
//...
//-------------------------------------------------------
// New Header SXS_HEAD
//[SXSFileHeader]
//Idn=File from SXSC (orig.meas.data) 2.0
//StartTime=Mon Jun 8 14:34:24 1998
//EndTime=Non Jun 8 18:35:33 1998
//NoofCols=3
//DataPoints=250
//ColID1=Energy [eV]
//ColID2=Counts
//ColId3=Position
//
//Greating=600
//SampleVacuum=1e-7
//IXRay=1.5
//UXRay=3.5
//UFilament=1.8
//IFilament=3.8
//UWehnelt=0.5
//UCEM=1.5kV
//CEMId=Galileo 4839
//UPhoto=125
//PhotoId=KaliumIodid
//Gain=1.5
//Operator=NN
//SampleID=text describing measurement
//SampleID1=....
//SampleId2=....
//
//StartPos=27000
//EndPos=32000
//Delta=20
//NoofSpectra=2
//Timebase=5.000000
//
//NoofMeas=1
//TotalMeas=5
//@EOH
//[######] [######] [######]
//Energy   Counts   Position
// *****************************************************
SxSFile::SxSFile(const char *szFile) : DataFile()
{
  Head.SpHead.szText[0]=0;
  Head.SpHead.uStartP=0;
  Head.SpHead.uEndP=0;
  Head.SpHead.uSteps=0;
  Head.SpHead.uDeltaP=0;
  Head.SpHead.unSpecs=0;
  Head.SpHead.fTimeBase=0;
  Head.SpHead.iMaxCounts=0;

 Head.szHead[0]=0;
 Head.szSTime[0]=Head.szETime[0]=0;
 Head.nCols=0;
 Head.szCID=0;
 Head.iFType=-1;

 Head.iGreating=0;
 Head.fSampVac=0;
 Head.fXI=0;
 Head.fXV=0;
 Head.fUFil=0;
 Head.fIFil=0;
 Head.fUWehn=0;
 Head.fUCEM=0;
 Head.szCEMId[0]=0;
 Head.fUPhoto=0;
 Head.fAGain=0;
 Head.szPhoto[0]=0;
 Head.szOperator[0]=0;
 Head.iMeas=0;
 Head.iTotMeas=0;
 Head.szFileName[0]=0;

  int iT=CheckFileType(szFile,iFLEnd);
  
  if(iT==FtERROR || iT==FtILLG )
    {iFType=(enum FileType)iT;
     iFGroup=NoGROUP;
     PRINT_DEBUG("Error analyzing file %s\n",szFile)
     return;
    }


  DataFile::Init(szFile,SxSGROUP,iT-SxSGROUP,"SxS Data File");
  if(!FindSets())
     {iFType=FtILLG;
      iFGroup=NoGROUP;
      PRINT_DEBUG("Error finding data set %s\n",szFile)
      return;
     } 
 
  if(!SPar.lfPhi)
    {PRINT_DEBUG("Cannot calculate energies; SxS pars missing\n")
    } 
  nTxtLines=1;
}
// ********************************************************
int SxSFile::FindSets(void)
{ nSets=0;
  switch(iFType)
    {case FtSXS_OLD:nSets=1;
	            return 1;
     case FtSXS_NEW:
	  {FilePar *FP=new FilePar(szName);
	   if(FP==0)return 0;
	   if(FP->GetStatus()!=FOUND){delete FP;return 0;}
           int i;
	   if(FP->GetPar("SXSFileHeader","TotalMeas","%i",&i)!=FOUND)return 0;
           nSets=i;
	   delete FP;
	   return nSets;
	  }
  default: return 0;
 }

}
// ********************************************************
int SxSFile::ReadHeader(FILE *InStream,FilePar *FP)
{// Read Header from InStream
 // return 0: Illegal Header
 //        1: Ok
 switch(iFType)
   {case FtSXS_OLD:
      fgets(Head.SpHead.szText,81,InStream);
      fscanf(InStream,"%u",&Head.SpHead.uStartP);
      fscanf(InStream,"%u",&Head.SpHead.uEndP);
      fscanf(InStream,"%u",&Head.SpHead.uSteps);
      fscanf(InStream,"%u",&Head.SpHead.uDeltaP);
      fscanf(InStream,"%u",&Head.SpHead.unSpecs);
      fscanf(InStream,"%f",&Head.SpHead.fTimeBase);
      fscanf(InStream,"%i",&Head.SpHead.iMaxCounts);
      return 1;
    case FtSXS_NEW:
      if(FP==0)return 0;
      if(FP->GetStatus()!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","Idn","%79S",Head.szHead)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","StartTime","%25S",Head.szSTime)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","EndTime","%25S",Head.szETime)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","NoofCols","%i",&(Head.nCols))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","DataPoints","%u",&(Head.SpHead.uSteps))!=FOUND)return 0;
      if(!AllocColID(Head.nCols))return 0;
      char szB[80];
      int i;
      for(i=0; i<Head.nCols; i++)
      	{sprintf(szB,"ColID%d",i+1);
	 if(FP->GetPar("SXSFileHeader",szB,"%25S",szCID[i])!=FOUND)return 0;
	}
      if(FP->GetPar("SXSFileHeader","Greating","%i",&(Head.iGreating))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","SampleVacuum","%f",&(Head.fSampVac))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","IXRay","%f",&(Head.fXI))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","UXRay","%f",&(Head.fXV))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","UFilament","%f",&(Head.fUFil))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","IFilament","%f",&(Head.fIFil))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","UWehnelt","%f",&(Head.fUWehn))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","UCEM","%f",&(Head.fUCEM))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","CEMId","%14S",Head.szCEMId)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","UPhoto","%f",&(Head.fUPhoto))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","Gain","%f",&(Head.fAGain))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","PhotoId","%14S",Head.szPhoto)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","Operator","%14S",Head.szOperator)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","SampleID","%80S",Head.SpHead.szText)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","StartPos","%u",&Head.SpHead.uStartP)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","EndPos","%u",&Head.SpHead.uEndP)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","Delta","%u",&Head.SpHead.uDeltaP)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","NoofSpectra","%u",&Head.SpHead.unSpecs)!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","Timebase","%f",&Head.SpHead.fTimeBase)!=FOUND)return 0;
      Head.SpHead.iMaxCounts=0;
      if(FP->GetPar("SXSFileHeader","NoofMeas","%i",&(Head.iMeas))!=FOUND)return 0;
      if(FP->GetPar("SXSFileHeader","TotalMeas","%i",&(Head.iTotMeas))!=FOUND)return 0;
      return 1;
  default: return 0;
 }
}
// ********************************************************
int  SxSFile::ReadData(const int iSet)
{// Reads complete data from file
 // return -1 : Error opening file
 //         0 : Illegal Header
 //         1 : Ok file complete
 //        less than -1: Number of data sets and data points is decoded
 //                      in the negative number 
 //                      -( (SetNo)*SXS_SET_OFFS + NoOfSteps)
 //                      SXS_SET_OFFS see dfile.h
 // iSet==0 Reads incomplete Datafile to show progress of measurement
 // or last data set if file is complete
 
 if(iFType==FtERROR || iFType==FtILLG || iFGroup != SxSGROUP)return 0;
 if(iFType==FtSXS_OLD && iSet != 1)
   {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
    return 0;}

 FILE *f;

  if ((f=fopen(szName,"r")) == NULL)
    {PRINT_DEBUG("Error opening file %s\n",szName)
     iFType=FtERROR;
     iFGroup=NoGROUP;
     return -1;
    }
 if(iFType==FtSXS_OLD)
   {if(!ReadHeader(f))return 0;;
    TD=new CRData(3,Head.SpHead.uSteps);
    CHECK_POINTER_RETURN(TD,-1)

    double pCh[5];
    unsigned i,j;
    for(i=0; i<Head.SpHead.uSteps; i+=5)
      {if(fscanf(f,"%lf  %lf  %lf  %lf  %lf",
	   &(pCh[0]),&(pCh[1]),&(pCh[2]),&(pCh[3]),&(pCh[4])) == EOF)
          {PRINT_DEBUG("Error reading columns from file %s\n",szName)
           iFType=FtERROR;
           iFGroup=NoGROUP;
           return 0;
          }
       for(j=0; j<5; j++)(*TD)[1][i+j]=pCh[j];
      }
    for(i=0; i<Head.SpHead.uSteps; i++)
      {j=Head.SpHead.uStartP+i*Head.SpHead.uDeltaP;
       (*TD)[2][i]=j;
       (*TD)[0][i]=(j>ZERO_ORDER && SPar.lfPhi ? SpecFunc((double)j) : j) ;
      }
    nCols=3;
    if(!AllocColID(nCols))return 0;

    sprintf(szCID[0],"Energy[eV]");
    sprintf(szCID[1],"Counts [%4.1f * %2d sec]",
                     Head.SpHead.fTimeBase,Head.SpHead.unSpecs);
    sprintf(szCID[2],"POS [Steps]");
    TD->SetColX(0);
    TD->SetColY(1);
    TD->SetColZ(2);
    TD->SortData(0);
    fclose(f);
    return 1;
   }// FtSXS_OLD

 if(iFType==FtSXS_NEW)
   {int iiSet=(iSet==0 ? nSets : iSet);
    if(iiSet<1 || iiSet>nSets)
     {PRINT_DEBUG("Illegal set # %d / %d\n",iSet,nSets)
      return 0;
     }
    FilePar *FP=new FilePar(f,F_NCLOSE);
    int i,iRet=1;
    unsigned j;
    float e,ii,p;
    for(i=0;i<iiSet-1;i++)
      {if(!ReadHeader(f,FP)){fclose(f);delete FP; return 0;}
       if(iSet==0)
         {TD=new CRData(3,Head.SpHead.uSteps);
          CHECK_POINTER_RETURN(TD,-1)
	 }
       for(j=0; j<Head.SpHead.uSteps; j++)
          {if(fscanf(f,"%f  %f  %f",&e,&ii,&p)==EOF)
	      {iRet=-((i+1)*SXS_SET_OFFS+j);break;}
           if(iSet==0)//Data file may be incomplete
	     {(*TD)[0][j]=e;
              (*TD)[1][j]=ii;
	      (*TD)[2][j]=p;
	     } 
          }
        if(iRet<0)
	  {nCols=Head.nCols;
           sprintf(szCID[1],"Counts [%4.1f * %2d sec]",
                             Head.SpHead.fTimeBase,Head.SpHead.unSpecs);
	   fclose(f);
	   delete FP;
	   return iRet;
	  }
        if(iSet==0)delete TD;           
        FP->Update(f);
        if(FP->GetStatus()!=FOUND){fclose(f);delete FP;return 0;}
       }


    if(!ReadHeader(f,FP)){delete FP; return 0;}
    TD=new CRData(3,Head.SpHead.uSteps);
    CHECK_POINTER_RETURN(TD,0)

    for(j=0; j<Head.SpHead.uSteps; j++)
       {if(fscanf(f,"%f  %f  %f",&e,&ii,&p)==EOF)
               {iRet=-((i+1)*SXS_SET_OFFS+j);break;}
        (*TD)[0][j]=e;
        (*TD)[1][j]=ii;
        (*TD)[2][j]=p;
       }
    if(iSet==0 && iRet<0)
      {for(unsigned jj=j; jj<Head.SpHead.uSteps; jj++)
          {unsigned s=Head.SpHead.uStartP+jj*Head.SpHead.uDeltaP;
           (*TD)[2][jj]=s;
           (*TD)[0][jj]=(s>ZERO_ORDER && SPar.lfPhi ? SpecFunc((double)s) : s);
           (*TD)[1][jj]=0;		
          }//for jj
      }//if

    nCols=Head.nCols;
    sprintf(szCID[1],"Counts [%4.1f * %2d sec]",
                      Head.SpHead.fTimeBase,Head.SpHead.unSpecs);

    TD->SetColX(0);
    TD->SetColY(1);
    TD->SetColZ(2);
    TD->SortData(0);
    
    fclose(f);
    delete FP;
    return iRet;
   }// FtSXS_NEW
 return 0;
}
// *******************************************************
double SxSFile::SpecFunc(double lfPos)
{
 double Alpha=M_PI/180.*SPar.lfAlpha;
 double Phi=M_PI/180.*SPar.lfPhi;
 double z=232.22857-lfPos*9.07143E-3;
 double d2q2=SPar.lfD*SPar.lfD-SPar.lfRQ*SPar.lfRQ;
 double Temp=SPar.lfRQ*sqrt(z*z+d2q2)/d2q2;
 Temp=.5*atan(SPar.lfD*z/d2q2+Temp);
 double Beta=Phi-Temp;
 double Scosa=SPar.lfSigma*cos(Alpha);
 double Lambda=1.e7*(Scosa-SPar.lfSigma*cos(Beta));
 return 12400./Lambda;
}
// ******************************************************
int  SxSFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG || iFGroup != SxSGROUP)return 0;

 FILE *OutStream;

 char szB[MAXPATH+1];

 if(szFile)strncpy(szB,szFile,MAXPATH);
 else strcpy(szB,szName);
 OutStream=fopen(szB,"w");
 CHECK_FILE_POINTER_RETURN(OutStream,szB,0)

 strncpy(szName,szB,MAXPATH);

 int iRet=SxSFile::SaveData(OutStream);
 fclose(OutStream);
 return iRet;
}
// ******************************************************
int  SxSFile::SaveData(FILE *OutStream)
{
 SetEndChr();
 Head.iFType=SXS_FTYPE_NEW;
 
 Head.SpHead.uSteps=(unsigned)TD->GetRows(0);
 Head.iMeas=1;
 Head.iTotMeas=1;
 Head.nCols=TD->GetCols();

 fprintf(OutStream,"[SXSFileHeader]%s",szE);
 fprintf(OutStream,"Idn=File from SXSC (orig.meas.data) 2.0%s",szE);
 fprintf(OutStream,"StartTime=%s%s",Head.szSTime,szE);
 fprintf(OutStream,"EndTime=%s%s",Head.szETime,szE);
 fprintf(OutStream,"NoofCols=%i%s",Head.nCols,szE);
 fprintf(OutStream,"DataPoints=%i%s",Head.SpHead.uSteps,szE);
 int i;
 for(i=0; i<Head.nCols; i++)
   {fprintf(OutStream,"ColID%1d=%s%s",i+1,szCID[i],szE);
   }
 fprintf(OutStream,"Greating=%i%s",Head.iGreating,szE);
 fprintf(OutStream,"SampleVacuum=%7.1e%s",Head.fSampVac,szE);
 fprintf(OutStream,"IXRay=%f%s",Head.fXI,szE);
 fprintf(OutStream,"UXRay=%f%s",Head.fXV,szE);
 fprintf(OutStream,"UFilament=%f%s",Head.fUFil,szE);
 fprintf(OutStream,"IFilament=%f%s",Head.fIFil,szE);
 fprintf(OutStream,"UWehnelt=%f%s",Head.fUWehn,szE);
 fprintf(OutStream,"UCEM=%f%s",Head.fUCEM,szE);
 fprintf(OutStream,"CEMId=%s%s",Head.szCEMId,szE);
 fprintf(OutStream,"UPhoto=%f%s",Head.fUPhoto,szE);
 fprintf(OutStream,"Gain=%f%s",Head.fAGain,szE);
 fprintf(OutStream,"PhotoId=%s%s",Head.szPhoto,szE);
 fprintf(OutStream,"Operator=%s%s",Head.szOperator,szE);
 fprintf(OutStream,"SampleID=%s%s",Head.SpHead.szText,szE);
 fprintf(OutStream,"StartPos=%u%s",Head.SpHead.uStartP,szE);
 fprintf(OutStream,"EndPos=%u%s",Head.SpHead.uEndP,szE);
 fprintf(OutStream,"Delta=%u%s",Head.SpHead.uDeltaP,szE);
 fprintf(OutStream,"NoofSpectra=%u%s",Head.SpHead.unSpecs,szE);
 fprintf(OutStream,"Timebase=%f%s",Head.SpHead.fTimeBase,szE);
 fprintf(OutStream,"NoofMeas=%i%s",Head.iMeas,szE);
 fprintf(OutStream,"TotalMeas=%i%s",Head.iTotMeas,szE);
 fprintf(OutStream,"@EOH%s",szE);

 i=TD->ExportData(OutStream,iFLEnd);
 if(i==0)return 0;

 return i;

}
// *******************************************************
char *SxSFile::SPrintInfo(char * szB, const unsigned Flag)
{char *sp=szB;

 if(iFType==FtERROR || iFType==FtILLG || iFGroup != SxSGROUP)
   {sp+=sprintf(sp,"Illegal file type %d ;cannot print data\n",iFType);
    return sp;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(sp);

 if(Flag & PR_TYPE)
   sp+=sprintf(sp,"       File Type: %s (%1d)\n",szID,iFType);
 if(Flag & PR_TEXT && Head.SpHead.szText)
     sp+=sprintf(sp,"%s\n",Head.SpHead.szText);
 if(Flag & PR_HEAD)
   {if(iFType==FtSXS_OLD)
      {sp+=sprintf(sp,"Start Position: %5d \n",Head.SpHead.uStartP); 
       sp+=sprintf(sp,"  End Position: %5d \n",Head.SpHead.uEndP); 
       sp+=sprintf(sp,"         Steps: %5d \n",Head.SpHead.uSteps);
       sp+=sprintf(sp,"         Delta: %5d \n",Head.SpHead.uDeltaP); 
       sp+=sprintf(sp,"  # of Spectra: %5d \n",Head.SpHead.unSpecs); 
       sp+=sprintf(sp," Time base [s]: %f\n",  Head.SpHead.fTimeBase);
       sp+=sprintf(sp,"   Max. Counts: %5d \n",Head.SpHead.iMaxCounts);
      } 
//       sp+=sprintf(sp,"   Columns: %i\n",Head.nCols);
//       sp+=sprintf(sp,"DataPoints:%i\n",Head.SpHead.uSteps);
//        int i;  
//        for(i=0; i<Head.nCols; i++)
//           {sp+=sprintf(sp,"ColID%1d: %s\n",i+1,szCID[i]);
//           }
    if(iFType==FtSXS_NEW)
      {sp+=sprintf(sp,"    Greating: %i\n",Head.iGreating);
       sp+=sprintf(sp,"SampleVacuum: %7.1e [torr]\n",Head.fSampVac);
       sp+=sprintf(sp,"       IXRay: %4.2f [mA]\n",Head.fXI);
       sp+=sprintf(sp,"       UXRay: %4.2f [kV]\n",Head.fXV);
       sp+=sprintf(sp,"   UFilament: %4.1f [V]\n",Head.fUFil);
       sp+=sprintf(sp,"   IFilament: %4.1f [A]\n",Head.fIFil);
       sp+=sprintf(sp,"    UWehnelt: %4.1f [V]\n",Head.fUWehn);
       sp+=sprintf(sp,"        UCEM: %4.1f [kV]\n",Head.fUCEM);
       sp+=sprintf(sp,"       CEMId: %s\n",Head.szCEMId);
       sp+=sprintf(sp,"      UPhoto: %3.0f [V]\n",Head.fUPhoto);
       sp+=sprintf(sp,"        Gain: %3.1f\n",Head.fAGain);
       sp+=sprintf(sp,"     PhotoId: %s\n",Head.szPhoto);
       sp+=sprintf(sp,"    Operator: %s\n",Head.szOperator);
       sp+=sprintf(sp,"    SampleID: %s\n",Head.SpHead.szText);
       sp+=sprintf(sp,"    StartPos: %u\n",Head.SpHead.uStartP);
       sp+=sprintf(sp,"      EndPos: %u\n",Head.SpHead.uEndP);
       sp+=sprintf(sp,"       Delta: %u\n",Head.SpHead.uDeltaP);
       sp+=sprintf(sp,"# of Spectra: %u\n",Head.SpHead.unSpecs);
       sp+=sprintf(sp,"    Timebase: %5.1f [s]\n",Head.SpHead.fTimeBase);
       sp+=sprintf(sp,"Run %d of %d\n",Head.iMeas,Head.iTotMeas);
      } 

 
   }
 if(Flag & PR_TIME && iFType==FtSXS_NEW)
   {sp+=sprintf(sp," StartTime: %s\n",Head.szSTime);
    sp+=sprintf(sp,"   EndTime: %s\n",Head.szETime);
   }

 if(Flag & PR_COLS)sp=SPrintColInfo(sp);
 *sp=0;
 // returns number of printed chars
 return sp;
}
// ******************************************************
void SxSFile::SPrintInfo(String &B, const unsigned Flag)
{

 if(iFType==FtERROR || iFType==FtILLG || iFGroup != SxSGROUP)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 DataFile::SPrintInfo(B);

 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",szID,iFType);
 if(Flag & PR_TEXT && Head.SpHead.szText)
     B.Addf("%s\n",Head.SpHead.szText);
 if(Flag & PR_HEAD)
   {if(iFType==FtSXS_OLD)
      {B.Addf("Start Position: %5d \n",Head.SpHead.uStartP); 
       B.Addf("  End Position: %5d \n",Head.SpHead.uEndP); 
       B.Addf("         Steps: %5d \n",Head.SpHead.uSteps);
       B.Addf("         Delta: %5d \n",Head.SpHead.uDeltaP); 
       B.Addf("  # of Spectra: %5d \n",Head.SpHead.unSpecs); 
       B.Addf(" Time base [s]: %f\n",  Head.SpHead.fTimeBase);
       B.Addf("   Max. Counts: %5d \n",Head.SpHead.iMaxCounts);
      } 
    if(iFType==FtSXS_NEW)
      {
//        B.Addf("   Columns: %i\n",Head.nCols);
//        B.Addf("DataPoints:%i\n",Head.SpHead.uSteps);
//        int i;  
//        for(i=0; i<Head.nCols; i++)
//           {B.Addf("ColID%1d: %s\n",i+1,szCID[i]);
//           }
       B.Addf("    Greating: %i\n",Head.iGreating);
       B.Addf("SampleVacuum: %7.1e [torr]\n",Head.fSampVac);
       B.Addf("       IXRay: %4.2f [mA]\n",Head.fXI);
       B.Addf("       UXRay: %4.2f [kV]\n",Head.fXV);
       B.Addf("   UFilament: %4.1f [V]\n",Head.fUFil);
       B.Addf("   IFilament: %4.1f [A]\n",Head.fIFil);
       B.Addf("    UWehnelt: %4.1f [V]\n",Head.fUWehn);
       B.Addf("        UCEM: %4.1f [kV]\n",Head.fUCEM);
       B.Addf("       CEMId: %s\n",Head.szCEMId);
       B.Addf("      UPhoto: %3.0f [V]\n",Head.fUPhoto);
       B.Addf("        Gain: %3.1f\n",Head.fAGain);
       B.Addf("     PhotoId: %s\n",Head.szPhoto);
       B.Addf("    Operator: %s\n",Head.szOperator);
       B.Addf("    SampleID: %s\n",Head.SpHead.szText);
       B.Addf("    StartPos: %u\n",Head.SpHead.uStartP);
       B.Addf("      EndPos: %u\n",Head.SpHead.uEndP);
       B.Addf("       Delta: %u\n",Head.SpHead.uDeltaP);
       B.Addf("# of Spectra: %u\n",Head.SpHead.unSpecs);
       B.Addf("    Timebase: %5.1f [s]\n",Head.SpHead.fTimeBase);
       B.Addf("Run %d of %d\n",Head.iMeas,Head.iTotMeas);
      } 
   }
if(iFType==FtSXS_NEW && (Flag & PR_COLS))
  {B.Addf(" StartTime: %s\n",Head.szSTime);
   B.Addf("   EndTime: %s\n",Head.szETime);
  }

if(Flag & PR_COLS)SPrintColInfo(B);
}
// ******************************************************
SXS_PAR SxSFile::ReadPars(const char * szFile)
{ 
  SXS_PAR S={0,0,0,0,0};

  const char *szConfPath=getenv("AUSW_CONF_PATH");
  char szB[MAXPATH+1];
  AppendPath(szConfPath,szFile,szB,MAXPATH);
  FilePar *FP=new FilePar(szB);
  if(FP==0)return S;
  if(FP->GetStatus()!=FOUND){delete FP;return S;}

 //[deg] Spektrometerkonstante
  FP->GetPar("SxsPars","Phi","%lf",&(S.lfPhi));
 // Rowlandkreisradius [mm]
  FP->GetPar("SxsPars","RowlandR","%lf",&(S.lfD));
  //[mm]_?Spindelsteigung_RQ
  FP->GetPar("SxsPars","ScrewSlope","%lf",&(S.lfRQ));
  // 1/600 1/d:Gitterparameter  Linien/mm
  FP->GetPar("SxsPars","InvLinesPmm","%lf",&(S.lfSigma));
  // ? Einfallswinkel Alpha ?
  FP->GetPar("SxsPars","Alpha","%lf",&(S.lfAlpha));
  delete FP;
  return S;
}
// ********************************************************
const char * SxSFile::GetInfoText(const int iL)
{
 if(iL==-1)return Head.SpHead.szText;
 else return szNoText;
}
// *******************************************************
int SxSFile::SetInfoText(const char *szT, const int iL)
{if(iL==-1){strncpy(Head.SpHead.szText,szT,MAX_SXSHEAD_LENGTH);
            Head.SpHead.szText[MAX_SXSHEAD_LENGTH]=0;
            return 1;
           }
 else return 0;
}
// *******************************************************
int SxSFile::AddSets(void)
{

 if(iFType==FtSXS_OLD)return -1;
 SxSFile *SF=new SxSFile(szName);
 
 int i,iAddC;
 float fTBase=0;
 for(i=0,iAddC=0;i<SF->GetNSets();i++)
    {
     if(i==0)
       {if(!ReadData(i+1))
          {PRINT_DEBUG("Error reading data\n")
           return 0;
          }
        if(TD==NULL || TD->GetSteps()<=0)
          {PRINT_DEBUG("Data Block empty\n")
           return 0;
          }
 
	fTBase=Head.SpHead.unSpecs*Head.SpHead.fTimeBase;
	iAddC++;  
       }	  
     else
       {if(!SF->ReadData(i+1))
          {PRINT_DEBUG("Error reading data\n")
           return 0;
          }

        if(Head.SpHead.uStartP != SF->GetSxSHead()->SpHead.uStartP ||
             Head.SpHead.uEndP != SF->GetSxSHead()->SpHead.uEndP ||
	    Head.SpHead.uSteps != SF->GetSxSHead()->SpHead.uSteps)continue;
 	(*TD)[1].BasicCalc('+',(*SF->GetCRData())[1]);
	fTBase+=SF->GetSxSHead()->SpHead.unSpecs*SF->GetSxSHead()->SpHead.fTimeBase;
	iAddC++;
       } 
    }    

 //TD->ExchgCol(0,2);   
 TD->SetColX(0);
 TD->SetColY(1);
 TD->SetColZ(2);

 Head.SpHead.fTimeBase=fTBase;
 Head.SpHead.unSpecs=1;
 Head.iTotMeas=1;
 Head.iMeas=1;

 char szC[MAX_IDLENGTH+1];
 sprintf(szC,"Counts [%*.0f sec]",(int)log10(Head.SpHead.fTimeBase)+1, 
              Head.SpHead.fTimeBase);
 SetColID(szC,1); 
 sprintf(szC," %d of %d spectra added",iAddC,SF->GetNSets() );
 strncat(Head.SpHead.szText,szC,MAX_SXSHEAD_LENGTH);
 Head.SpHead.szText[MAX_SXSHEAD_LENGTH]=0;
 delete SF;
 return iAddC;
}
// ******************************************************
SxSFile::~SxSFile()
{ }
// ******************************************************
// *******************************************************
// TransFile
// *******************************************************
TransFile::TransFile(const char *szFile) : AsciiFile(szFile)
{

 if(iFType==FtILLG || iFType==FtERROR)return; 

 int iT=CheckFileType(szFile,iFLEnd);

 if( (int) (iT/10)*10 != TransGROUP)
   {iFType=FtILLG;
    iFGroup=NoGROUP;
    return;
   }

  DataFile::Init(szFile,TransGROUP,iT-TransGROUP,TransTYPES[iT-TransGROUP]);
 nSets=1;
}
// ********************************************************
int  TransFile::ReadData(const int iSet)
{// Reads complete data from file
 // return -1 : Error opening file
 //         0 : Illegal Header
 //         1 : Ok
 if(iSet!=1)
  {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
   return 0;}

 if(iFGroup==TransGROUP)return AsciiFile::ReadData(iSet);
 // nTxtLines=2 updated from AsciiFile::ReadData
 return 0;
}
// *******************************************************
int  TransFile::SaveData(const char *szFile)
{
 if(iFGroup==TransGROUP)return AsciiFile::SaveData(szFile);
 return 0;
}
// *******************************************************
int  TransFile::SaveData(FILE *fOut)
{
 if(iFGroup==TransGROUP)return AsciiFile::SaveData(fOut);
 return 0;
}
// *******************************************************
char *TransFile::SPrintInfo(char * szB, const unsigned Flag)
{ if(iFGroup==TransGROUP) return AsciiFile::SPrintInfo(szB,Flag);
  return 0;
}
// ******************************************************
void TransFile::SPrintInfo(String &B, const unsigned Flag)
{ if(iFGroup==TransGROUP) AsciiFile::SPrintInfo(B,Flag);
}
// ******************************************************
const char * TransFile::GetInfoText(const int iL)
{
 if(iFType==FtTRANS_DATAP)
  {if(iL==-1)return AsciiFile::GetInfoText(1);}
 return  AsciiFile::GetInfoText(iL==-1 ? 0 : iL);
}
// ****************************************************** 
int TransFile::SetInfoText(const char *szT, const int iL)
{if(iFType==FtTRANS_DATAP)
  {if(iL==-1) return AsciiFile::SetInfoText(szT,1);}
 return AsciiFile::SetInfoText(szT,iL);
}
// *****************************************************
TransFile::~TransFile()
{ }
// ******************************************************
// ********************************************************
// XDifFile
// ******************************************************
// Header of RAW/DIF File:
// File is a binary file
//    char Key[4] -> "RAW" or "PEAK"
//    INT4 nSteps -> Number of steps
//  FLOAT4 fTBase -> Time base [s]
//  FLOAT4 fStepSize -> 2 theta step size
//    INT4 DACOflag
//    INT4 SamplePos
//  FLOAT4 fTheta2Start,fThetaStart,fChiStart,fPhiStart -> Start angles
//    char SampleID[32]
//  FLOAT4 fWLKa1,fWLKa2  -> Wavelength
//    char Unused[68]
//    INT4 nPeaks -> number of peaks
//    INT4 nRanges -> number of ranges
// ********************************************************
XDifFile::XDifFile(const char *szFile) : DataFile()
{
  //Inital values will be overwritten if read from file
  Head.Key[0]=Head.SampleID[0]=0;
  Head.DACOflag=Head.SamplePos=Head.nPeaks=Head.nRanges=0;
  Head.nSteps=0;
  Head.fTBase=Head.fStepSize=0;
  Head.fTheta2Start = Head.fThetaStart=0;
  Head.fChiStart=Head.fPhiStart=0;
  Head.fWLKa1=Head.fWLKa2=0;

  int iT=CheckFileType(szFile,iFLEnd);

  if(iT==FtERROR || iT==FtILLG)
    {iFType=(enum FileType)iT;
     iFGroup=NoGROUP;
     PRINT_DEBUG("Error analyzing file %s\n",szFile)
     return;
    }

  DataFile::Init(szFile,XDifGROUP,iT-XDifGROUP,XDifTYPES[iT-XDifGROUP]);
  nTxtLines=1;
  
  if(iFType==FtXDIF_RAW)nSets=1;
  if(iFType==FtXDIF_DIF)nSets=Head.nRanges+1;
}
// ********************************************************
XDifFile::~XDifFile() { }
// ********************************************************
void XDifFile::SetRange(const int iR)
{if(iFType==FtXDIF_DIF && iR>0 && iR<Head.nRanges)ReadRange=iR;
 else 
    {PRINT_DEBUG("Range %d out of range (1 ... %d)\n",iR,Head.nRanges)
    }
}
// ********************************************************
int XDifFile::ReadHeader(FILE *InStream)
{// Read Header from InStream
 // return 0: Illegal Header
 //        1: Ok
// ******************************************************
// Header of XDif File see above

 if(iFGroup != XDifGROUP)return 0;
 int i;
 switch (iFType)
  {case    FtERROR: return 0;
   case     FtILLG: return 0;
   case    FtASCII: return 0;
   case FtXDIF_RAW: if(!AllocColID(2))return 0;
                    fseek(InStream,0,SEEK_SET);
                    i=fread(&Head,sizeof(Head),1,InStream);
                    sprintf(szCID[0],"2 Theta");
                    sprintf(szCID[1],"Counts/%5.2f [s]",Head.fTBase);
		    nSets=1;
                    if(!i)return 0;
                    return Head.nSteps;
   case FtXDIF_DIF: if(!AllocColID(2))return 0;
                    fseek(InStream,0,SEEK_SET);
                    i=fread(&Head,sizeof(Head),1,InStream);
                    sprintf(szCID[0],"d");
                    sprintf(szCID[1],"Intensity");
                    if(!i)return 0;
		    nSets=Head.nRanges;
                    return Head.nPeaks;
   default: return 0;
  }

 //return 0;
}
// ******************************************************
int XDifFile::ReadData(const int iSet)
{// Reads complete data from file
 // return -1 : Error opening file
 //         0 : Illegal Header
 //         1 : Ok

 if(iFType==FtERROR || iFType==FtILLG || iFGroup != XDifGROUP)return 0;

 FILE *f;
 if(iSet!=1 && iFType==FtXDIF_RAW)
  {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
   return 0;}

  if ((f=fopen(szName,"r")) == NULL)
    {PRINT_DEBUG("Error opening file %s\n",szName)
     iFType=FtERROR;
     iFGroup=NoGROUP;
     return -1;
    }

  int i=ReadHeader(f);
  if(!i)
    {fclose(f);
     PRINT_DEBUG("Error reading file header %s\n",szName)
     iFType=FtERROR;
     iFType=FtERROR;
     iFGroup=NoGROUP;
     fclose(f);
     return 0;
    }

 if(iFType==FtXDIF_RAW)
   { TD=new CRData(2,Head.nSteps);
     CHECK_POINTER_RETURN(TD,-1)

     int kstart=sizeof(Head);
     int k,iC=0;
     FLOAT4 c;
     for(k=kstart;k<kstart+Head.nSteps*4;k+=4)
        {fseek(f,k,SEEK_SET);
         fread(&c,sizeof (c),1,f);
         (*TD)[0][iC]=Head.fTheta2Start+Head.fStepSize*iC;
	 (*TD)[1][iC]=c;
         //lambda = 2.d.sin(theta)
         //TD->SetPointV(iC)=c
	     iC++;
       }/*for(k)*/
     TD->SetColX(0);
     TD->SetColY(1);
     TD->NewMinMax();
     TD->SortData(0);
     fclose(f);
     return 1;
    }
 if(iFType==FtXDIF_DIF)
   {if(iSet <=0 || iSet>Head.nRanges)
      {PRINT_DEBUG("Illegal data set %d / %d\n",iSet, Head.nRanges)
       return 0;}
    int ii,k,steps;
    for(ii=0,steps=0;ii<iSet;ii++)
       {k=sizeof(Head)*ii+2*sizeof(FLOAT4)*steps;
        fseek(f,k,SEEK_SET);
        fread(&Head,sizeof(Head),1,f);
        steps+=Head.nPeaks;
       }/*for(ii)*/
      
       
     TD=new CRData(2,Head.nPeaks);
     CHECK_POINTER_RETURN(TD,-1)

     int kstart=k+sizeof(Head);
     int iC=0;
     FLOAT4 c,d;
     for(k=kstart;k<kstart+Head.nPeaks+2*(int)sizeof(FLOAT4);k+=2*sizeof(FLOAT4))
        {fseek(f,k,SEEK_SET);
         fread(&d,sizeof (d),1,f);
         fread(&c,sizeof (c),1,f);
         (*TD)[0][iC]=d;
	 (*TD)[1][iC]=c;
	 iC++;
        }/*for(k)*/
     TD->SetColX(0);
     TD->SetColY(1);
     TD->NewMinMax();
     TD->SortData(0);
     fclose(f);
     return 1;
   }
  


 fclose(f);
 return 0;
}
// *******************************************************
char *XDifFile::SPrintInfo(char * szB, const unsigned Flag)
{
 char *sp=szB;

 if(iFType==FtERROR || iFType==FtILLG || iFGroup != XDifGROUP)
   {sp+=sprintf(sp,"Illegal file type %d ;cannot print data\n",iFType);
    return sp;
   }

 if(Flag & PR_BASIC)
   DataFile::SPrintInfo(sp);

 if(Flag & PR_TEXT)
   sp+=sprintf(sp,"Sample ID: %s\n",Head.SampleID);
 if(Flag & PR_TYPE)
   sp+=sprintf(sp,"       File Type: %s (%1d)\n",szID,iFType);

 if(Flag & PR_HEAD)
   {
    sp+=sprintf(sp,"        Key: %-4.4s\n",Head.Key);
    sp+=sprintf(sp," # of steps: %-ld\n",(long)Head.nSteps);
    sp+=sprintf(sp,"  Time base: %-6.2f\n",Head.fTBase);
    sp+=sprintf(sp,"  Step size: %-6.3f\n",Head.fStepSize);
    sp+=sprintf(sp,"          DACO flag: %0lx\n",(unsigned)Head.DACOflag);
    sp+=sprintf(sp,"        Sample Pos.: %-lu\n",Head.SamplePos);
    sp+=sprintf(sp," 2 THETA Start: %-7.3f\n",Head.fTheta2Start);
    sp+=sprintf(sp,"   THETA Start: %-7.3f\n",Head.fThetaStart);
    sp+=sprintf(sp,"     CHI Start: %-7.3f\n",Head.fChiStart);
    sp+=sprintf(sp,"     PHI Start: %-7.3f\n",Head.fPhiStart);
    sp+=sprintf(sp,"     Sample ID: %-32.32s\n",Head.SampleID);
    sp+=sprintf(sp,"   Wavelength Ka1: %-7.5f\n",Head.fWLKa1);
    sp+=sprintf(sp,"   Wavelength Ka2: %-7.5f\n",Head.fWLKa2);
    sp+=sprintf(sp,"  # of peaks: %-ld\n",(long)Head.nPeaks);
    sp+=sprintf(sp," # of ranges: %-ld\n",(long)Head.nRanges);
 }
 if(Flag & PR_COLS)sp=DataFile::SPrintColInfo(sp);
 // return number of printed chars
 *sp=0;
 return sp;
}
// ******************************************************
void XDifFile::SPrintInfo(String &B, const unsigned Flag)
{
 if(iFType==FtERROR || iFType==FtILLG || iFGroup != XDifGROUP)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(B);

 if(Flag & PR_TEXT)
   B.Addf("Sample ID: %s\n",Head.SampleID);
 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",szID,iFType);

 if(Flag & PR_HEAD)
   {
    B.Addf("        Key: %-4.4s\n",Head.Key);
    B.Addf(" # of steps: %-ld\n",(long)Head.nSteps);
    B.Addf("  Time base: %-6.2f\n",Head.fTBase);
    B.Addf("  Step size: %-6.3f\n",Head.fStepSize);
    B.Addf("          DACO flag: %0lx\n",(unsigned)Head.DACOflag);
    B.Addf("        Sample Pos.: %-lu\n",Head.SamplePos);
    B.Addf(" 2 THETA Start: %-7.3f\n",Head.fTheta2Start);
    B.Addf("   THETA Start: %-7.3f\n",Head.fThetaStart);
    B.Addf("     CHI Start: %-7.3f\n",Head.fChiStart);
    B.Addf("     PHI Start: %-7.3f\n",Head.fPhiStart);
    B.Addf("     Sample ID: %-32.32s\n",Head.SampleID);
    B.Addf("   Wavelength Ka1: %-7.5f\n",Head.fWLKa1);
    B.Addf("   Wavelength Ka2: %-7.5f\n",Head.fWLKa2);
    B.Addf("  # of peaks: %-ld\n",(long)Head.nPeaks);
    B.Addf(" # of ranges: %-ld\n",(long)Head.nRanges);
 }
 if(Flag & PR_COLS)DataFile::SPrintColInfo(B);
}
// ******************************************************
void XDifFile::SaveHeader(FILE *OutStream)
{ char szB[80];
 SetEndChr();

 GetCurrentTime(szB);
 fprintf(OutStream,FILE_ID,VERSION,szB);
 fprintf(OutStream,"%s%s",Head.SampleID,szE);
}
// ********************************************************
int XDifFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG || iFGroup != XDifGROUP)return 0;

 FILE *f;

 char szB[MAXPATH+1];

 if(szFile)strncpy(szB,szFile,MAXPATH);
 else strcpy(szB,szName);
 f=fopen(szB,"w");
 CHECK_FILE_POINTER_RETURN(f,szB,0)

 strncpy(szName,szB,MAXPATH);
 PRINT_DEBUG("File %s saved as ASCII; not in original data format\n",szName)
 int iRet=XDifFile::SaveData(f);

 fclose(f);
 return iRet;
}
// ********************************************************
int XDifFile::SaveData(FILE *fOut)
{ SaveHeader(fOut);
  return TD->ExportData(fOut,iFLEnd);
}
// ********************************************************
const char * XDifFile::GetInfoText(const int iL)
{
 if(iL==-1)return Head.SampleID;
 else return szNoText;
}
// *******************************************************
int XDifFile::SetInfoText(const char *szT, const int iL)
{if(iL==-1){strncpy(Head.SampleID,szT,31);
            return 1;
           }
 else return 0;
}
// *****************************************************
// *******************************************************
// SplineFile
// *******************************************************
// SPLINETABLE Sensor:GaAs SerialNo.:9243
// 85 -> #of data points
// 1.2413254 45.657465
// .....
// ******************************************************
SplineFile::SplineFile(const char *szFile) : DataFile()
{
 nDataLines=0;

 iFType=FtERROR;
 iFGroup=NoGROUP;

 int iT=CheckFileType(szFile,iFLEnd);

  if(iT==FtERROR || iT==FtILLG)
    {iFType=(enum FileType)iT;
     iFGroup=NoGROUP;
     PRINT_DEBUG("Error analyzing file %s\n",szFile)
     return;
    }

  DataFile::Init(szFile,SplineGROUP,iT-SplineGROUP,"Spline File");
  nTxtLines=1;

   iFType=FtSPLINE;
  iFGroup=NoGROUP;
  nSets=1;
}
// ********************************************************
int  SplineFile::ReadData(const int iSet)
{// Reads complete data from file
 // return -1 : Error opening file
 //         0 : Illegal Header
 //         1 : Ok

 if(iFType==FtERROR || iFType==FtILLG)return 0;
 int iErr=0;

 if(iSet!=1)
  {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
   return 0;}


 FILE *f=fopen(szName,"r");
 CHECK_FILE_POINTER_RETURN(f,szName,0)
 char szB[MAX_LINELENGTH];

 if(!fgets(szTextLine,MAX_LINELENGTH-1,f))iErr=1;
 RemoveCR_LF(szTextLine);
 if(strncmp(szTextLine,"SPLINETABLE",11))iErr=1;

 if(!fgets(szB,MAX_LINELENGTH-1,f))iErr=1;
 RemoveCR_LF(szB);
 int j=0;
 while(szB[j])
   {if(szB[j]==' '){j++; continue;}
    if(strchr("0123456789",szB[j])==NULL){iErr=1; break;}
    j++;
    }//while

 int nPoints=atoi(szB);;
 if(nPoints<3)iErr=1;
 
 if(iErr)
   {PRINT_DEBUG("Error reading from file %s\n",szName)
    return 0;
   }
 
 if(TD)delete TD;
 TD = new CRData(3,nPoints);
 // j=TD->GetSteps();
 if(! TD || TD->GetSteps()<=3)
   {PRINT_DEBUG("Error allocating memory\n")
    return 0;
   }
 nDataLines=nPoints;

 MDATA x,y;
 for(j=0; j<nDataLines; j++)
    { fscanf(f,"%lf %lf", &x,&y);
      (*TD)[0][j]=x;
      (*TD)[1][j]=y;
//      *(TD)[2][j]=y;
    }
 TD->SetColX(0);
 TD->SetColY(1);
 TD->SetColZ(2);
 
 TD->Spline();
 TD->SortData(0);
 TD->NewMinMax();


 DeleteColID();
 if(!AllocColID(3))return 0;

 sprintf(szCID[0],"X");
 sprintf(szCID[1],"Y");
 sprintf(szCID[2],"Spline");

 fclose(f);
 return 1;
}
// *******************************************************
int  SplineFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG)return 0;

 FILE *f;

 char szB[MAXPATH+1];

 if(szFile)strncpy(szB,szFile,MAXPATH);
 else strcpy(szB,szName);
 f=fopen(szB,"w");
 CHECK_FILE_POINTER_RETURN(f,szB,0)

 strncpy(szName,szB,MAXPATH);
 int iRet=SaveData(f);
 if(iRet!=0)fclose(f);
 return iRet;
}
// *******************************************************
int  SplineFile::SaveData(FILE *fOut)
{
 SetEndChr();

 fprintf(fOut,"%s%s",szTextLine,szE);
 fprintf(fOut,"%d%s",nDataLines,szE);

 return TD->ExportData(fOut,iFLEnd,0,1);
}
// *******************************************************
char *SplineFile::SPrintInfo(char * szB, const unsigned Flag)
{char *sp=szB;

 if(iFType==FtERROR || iFType==FtILLG)
   {sp+=sprintf(sp,"Illegal file type %d ;cannot print data\n",iFType);
    return sp;
   }

 if(Flag & PR_BASIC)
   DataFile::SPrintInfo(sp);


 if(Flag & PR_TYPE)
   sp+=sprintf(sp,"       File Type: %s (%1d)\n",szID,iFType);
 if(Flag & PR_TEXT && szTextLine)
   {sp+=sprintf(sp,"%s\n",szTextLine);
   }
 if(Flag & PR_COLS)sp=SPrintColInfo(sp);
 *sp=0;
 // returns number of printed chars
 return sp;
}
// ******************************************************
void SplineFile::SPrintInfo(String &B, const unsigned Flag)
{
 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(B);

 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",szID,iFType);
 if(Flag & PR_TEXT && szTextLine)
   {B.Addf("%s\n",szTextLine);
   }
 if(Flag & PR_COLS)SPrintColInfo(B);
}
// ******************************************************
const char * SplineFile::GetInfoText(const int iL=0)
{if(iL==0)return szTextLine;
 else return szNoText; 
}
// ******************************************************
int SplineFile::SetInfoText(const char *szT, const int iL=0)
{if(iL==0){strncpy(szTextLine,szT,MAX_LINELENGTH);
           return 1;
          }
 return 0;
}
// ******************************************************
SplineFile::~SplineFile()
{ }
// ******************************************************

