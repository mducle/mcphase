//File:dfile.cpp
//$Log: dfile.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//$Id: dfile.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>

#include "stdinc.h"
#include "stdfunc.h"
#include "thefunc.h"
#include "cdata.h"

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

#define XAUSW_ID "XAusw for LINUX V 1.1"

 // ******************************************************
 // General Functions
 // ******************************************************
int CheckFileType(const char *szFile, enum EndLine &LineT )
{FILE *f;
 f=fopen(szFile,"r");
 CHECK_FILE_POINTER_RETURN(f,szFile,FtNONE)

 LineT=NOEnd;

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
//      if(IsNumString(szB)){fclose (f);return FtILLG;}
//      if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
//         RemoveCR_LF(szB);
// 	//3. line 
//         if(IsNumString(szB)){fclose (f);return FtTRANS_DATAP;}
//         return FtILLG;
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


//  // 3.Line
//  if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f);return FtERROR;}
//  RemoveCR_LF(szB);
// // if(IsNumString(szB)){fclose (f);return FtTRANS_MEAS;}
//  if(IsNumString(szB)){fclose (f);return FtILLG;}
// 
//  // 4.Line
//  if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f); return FtERROR;}
//  RemoveCR_LF(szB);
//  if(IsNumString(szB)){fclose (f); return FtILLG;}
// 
//  // 5.Line
//  if(!fgets(szB,MAX_LINELENGTH-1,f)){fclose (f); return FtERROR;}
//  RemoveCR_LF(szB);
//  if(!IsNumString(szB)){fclose (f); return FtILLG;}
//  if(CountCols(szB)==12) {fclose (f); return FtHE3_MEAS;}
 fclose (f);
 return FtILLG;
}
// ******************************************************
// ******************************************************
// DataFile Functions
// ******************************************************
void DataFile::Init(void)
{
 TD=NULL;
 nSets=0;

 ColID=0;
 nCols=0;

 iFType=FtILLG;
 iFLEnd=NOEnd;
 FId.Set("None");

 //FPars=0;
}
// ******************************************************
int DataFile::Init(const char *szFile,const int iT, const char *szI)
{
 FId.Set(szI);
 FName.Set(szFile);

 iFType= (enum FileType) (iT);

 if(!CheckType())
   {PRINT_DEBUG("Illegal File Type (%d)\n",iT)
    return 0;
   }

 return 1;

}
// ***********************************************************
int DataFile::CheckType(void)
{
 if(iFType>(FileType) TheGROUP && iFType<=(FileType) FtTHE_CAPOLD)return 1;
 if(iFType>(FileType) SxSGROUP && iFType<=(FileType) FtSXS_NEW)return 1;
 if(iFType>(FileType) XDifGROUP && iFType<=(FileType) FtXDIF_DIF)return 1;
// if(iFType>(FileType) TransGROUP && iFType<=(FileType) FtTRANS_MEAS)return 1;
// if(iFType>(FileType) He3GROUP &&iFType<=(FileType) FtHE3_MEAS)return 1;
 if(iFType==(FileType) FtSPLINE)return 1;
 if(iFType==(FileType) FtASCII)return 1;
 iFType=FtERROR;
 return 0;
}
// ******************************************************
int DataFile::InsertColID(const char *szN, const int iCol)
{
    
 if( !ColID || !nCols)
   {PRINT_DEBUG("Cannot insert Column Id\n\n")
    return 0;
   }

CHECK_INDEX_RETURN(iCol,nCols+1,0)
nCols++;
return ColID->InsertLine(szN,iCol);

}
// ********************************************************
const char * DataFile::GetColID(int iC) const
{
 CHECK_INDEX_RETURN(iC,nCols,NULL)
 return ColID->Line(iC);
}
// ********************************************************
int DataFile::SetColID(const char *szI, const int iC)
{
 CHECK_INDEX_RETURN(iC,nCols,0)
 ColID->ReplaceLine(szI,iC);
 return 1;
}
// ******************************************************
int DataFile::SetColRange(const int iCS, const int iCE)
{ CHECK_INDEX_RETURN(iCS,nCols,0)
  CHECK_INDEX_RETURN(iCE,nCols,0)
  if(iCE<=iCS)return 0;

  String NS;

  int i,j;
  for (i=iCS,j=0;i<=iCE;i++,j++)
      {if(j==0)NS.Setf("%s%s",ColID->Line(i),(const char*)DeLim);
       else NS.Addf("%s%s",ColID->Line(i),(const char*)DeLim);;
      }

  LineString *NCID = new LineString(NS);
  CHECK_POINTER_RETURN(NCID,0)

  delete ColID;
  ColID=NCID;
  
  nCols=iCE-iCS+1;
  return nCols;
}
// ******************************************************
void DataFile::SPrintInfo(String &B) const
{
 B.Addf("File: %s (%s)\n",(const char *)FName, szLineType[(int)iFLEnd]);

 if(TD)
   B.Addf("%d data points, %d columns\n",TD->GetSteps(),TD->GetCols());
 else B.Add("NO data columns\n");
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
      B.Addf("Column %d:%-14.14s",j+1,(ColID?ColID->Line(j):" "));
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
  FName.Set(szFile);
  return i;
}
// ******************************************************
int DataFile::SaveData(const char *szText, FILE *f)
{ 
  SetEndChr();
  fprintf(f,"%s",szText);
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

  fprintf(f,"SPLINETABLE from %s (x:%d, y:%d)%s",(const char*)FName,
                                                iX+1,iY+1,(const char*)DeLim);
  fprintf(f,"%d%s",TD->GetSteps(),(const char*)DeLim);

  int i,j;
  double lfx,lfy;
  for(i=0;i<TD->GetSteps();i++)
     {for(j=0;j<TD->GetCols();j++)
         {if(j==iX)lfx=(*TD)[j][i];
          if(j==iY)lfy=(*TD)[j][i];
         }// for j
      fprintf(f,"%13.7g  %13.7g%s",lfx,lfy,(const char*)DeLim);
     }//for i

  return 1;

}		
// *****************************************************		
void DataFile::SetEndChr(void)
{ DeLim.Set("\n");
#if defined(LINUX)
  if(iFLEnd==NOEnd || iFLEnd==DOSEnd)DeLim.Set("\r\n");
#endif 
}
// *****************************************************		
int DataFile::ReplaceLine(const char *SearchT, const char * WithT)
{int i;
 if(!FText.Find(i,SearchT))
   {PRINT_DEBUG("Can not find line (Text:%s)\n",SearchT);
    return 0;
   }
 return FText.ReplaceLine(WithT,i);
}
// *****************************************************
DataFile::~DataFile()
{
 if(TD)delete TD;
}
// *******************************************************
// *******************************************************
// AsciiFile
// *******************************************************
// *******************************************************
AsciiFile::AsciiFile(const char *szFile) : DataFile()
{
  DataFile::Init(szFile,FtASCII,"ASCII File");

  iFType=FtASCII;
  nSets=1;
}
// ********************************************************
int  AsciiFile::ReadData(const int iSet)
{// Reads complete data from file
 // return  0 : Ok
 //         1 : Error opening file
 //         2 : Illegal Header
 //         3 : No memory
 if(iFType==FtERROR || iFType==FtILLG)return 2;
 // Multiple sets not possible
 if(iSet!=1)
  {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
   return 2;}
 FILE *f;
 int i;

  if ((f=fopen((const char *)FName,"r")) == NULL)
    {PRINT_DEBUG("Error opening file %s\n",(const char *)FName)
     iFType=FtERROR;
     return 1;
    }

 i=FText.Read(f," ",LineString::S_NUMBER);
 if(i<=0){fclose(f); return 1;}

 DeLim.Set(FText.GetDelim());

 if(TD)delete TD;

 LineString Columns;
 
 rewind(f);

 if(!Columns.Read(f," ",LineString::S_ISCOL))
    {PRINT_DEBUG("Error reading from file %s\n",(const char *)FName)
     iFType=FtERROR;
     return 1;
    }
 nDataLines=Columns.GetNLines();
 nCols=i;
 TD=new CRData(nCols,nDataLines);
 CHECK_POINTER_RETURN(TD,3)

 i=TD->ReadCols(Columns,nDataLines,nCols);
 if(i<=0)
   {PRINT_DEBUG("Error reading columns from file %s line %d\n",
                (const char *)FName,i)
    iFType=FtERROR;
    return 2;
   }

 String CS;
 for(i=0; i<nCols; i++)CS.Addf("Col%d\n",i+1);
 
 delete ColID;
 ColID=new LineString(CS);

 fclose(f);
 return 0;
}
// *******************************************************
int  AsciiFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG)return 0;

 FILE *f;
 String FN;

 if(szFile)FN.Set(szFile);
 else FN=FName;
 f=fopen((const char *)FN,"w");
 CHECK_FILE_POINTER_RETURN(f,(const char *)FN,0)

 FName=FN;
 int iRet=AsciiFile::SaveData(f);
 fclose(f);
 return iRet;
 }
// *******************************************************
int  AsciiFile::SaveData(FILE *fOut)
{int j;
 SetEndChr();
 for(j=0; j<FText.GetNLines(); j++)
    {if(!fprintf(fOut,"%s%s",FText.Line(j),(const char *)DeLim))return 0;}
 int i=TD->ExportData(fOut,iFLEnd);
 if(i==0)return 0;
 return i;
}
// *******************************************************
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
   B.Addf("       File Type: %s (%1d)\n",(const char *)FId,iFType);
 if(Flag & PR_TEXT)
   {for(i=0;i<FText.GetNLines();i++)
     B.Addf("%s\n",FText.Line(i));
     B.Addf("... %d textlines...\n",FText.GetNLines());
   }
 if(Flag & PR_COLS)DataFile::SPrintColInfo(B);
}
// ******************************************************
const char * AsciiFile::GetLine(const int iL, const int iWithE)
{if(iL<0 || iL>=FText.GetNLines())return szNoText;
 return FText.Line(iL,iWithE);
}
// ******************************************************
int AsciiFile::SetLine(const char * szT, const int iL)
{return FText.ReplaceLine(szT,iL);
}
// ******************************************************
const char * AsciiFile::GetInfoText(const int iL, const int iWithE)
{return GetLine((iL==-1?0:iL),iWithE);
}
// ******************************************************
int AsciiFile::SetInfoText(const char *szT, const int iL=0)
{return SetLine(szT,iL);}
// *******************************************************
AsciiFile::~AsciiFile()
{ }

#ifdef MORE_FILETYPES
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
TheFile::TheFile(const char *szFile) 
{
 int iT=CheckFileType(szFile,iFLEnd);

 if(iT==FtERROR || iT==FtILLG)
   {iFType=(enum FileType)iT;
    PRINT_DEBUG("Error analyzing file %s\n",szFile)
    return;
  }
 
 iFType=(FileType)iT;
 nSets=1;
 DataFile::Init(szFile,iT,TheTYPES[iT-TheGROUP]);
}
// ********************************************************
TheFile::~TheFile() { }
// ********************************************************
int TheFile::ReadData(const int iSet)
{// Reads complete data from file

 int i=AsciiFile::ReadData(iSet);
 
 char szB[MAX_LINELENGTH+1];
 String S;
 
 int j;
 for(j=0;j<nCols;j++)
    {S.Setf("ColID%1d",j+1);
     if(!FText.GetPar(0,(const char *)S,"%S",szB))
       {/*PRINT_DEBUG("Data file integrity has changed (ColID)\n")*/
        break;
       }
     if(ColID)ColID->ReplaceLine(szB,j);
    }
    
 return i;
}
// *******************************************************
void TheFile::SPrintInfo(String &B, const unsigned Flag)
{
 
 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(Flag & PR_BASIC)
   DataFile::SPrintInfo(B);

 if(Flag & PR_TEXT)
   B.Addf("%s\n",FindPar("[FileHeader]","Id"));
   
   
 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",(const char *)FId,iFType);

 if(Flag & PR_HEAD)
   FText.Write(B,1);

 if(Flag & PR_TIME && iFType==FtTHE_CAPV1)
     {
      B.Addf("%s\n",FindPar("[FileHeader]","StartTime"));
      B.Addf("%s\n",FindPar("[FileHeader]","EndTime"));
     }//if
 if(Flag & PR_COLS)DataFile::SPrintColInfo(B);
}
// ******************************************************
int TheFile::SaveHeader(void)
{
 int iRet=0;
 SetEndChr();

 if(iFType==FtTHE_RCP || iFType==FtTHE_CAPV1)
      {int i;

       String S;
       S.Setf("NoofCols=%i",TD->GetCols());
       if(!ReplaceLine("NoofCols=",(const char *)S))
         {PRINT_DEBUG("Data file integrity is corrupt\n");
	  return iRet;
	 }

       S.Setf("DataPoints=%i",TD->GetRows(0));
       if(!ReplaceLine("DataPoints=",(const char *)S))
         {PRINT_DEBUG("Data file integrity is corrupt\n");
	  return iRet;
	 }


       int j=1;
       while(1)
         {S.Setf("ColID%1d",j);
          if(!FText.Find(i,(const char *)S))break;
          FText.RemoveLine(i);
          j++;
         }

       if(!FText.Find(i,"NoofCols="))
         {PRINT_DEBUG("Data file integrity is corrrupt\n");
	  return iRet;
	 }

       for(j=0;j<nCols;j++,i++)
          {S.Setf("ColID%1d=%s",j+1,ColID->Line(j));
           FText.InsertLine((const char *)S,i);
          }

       char szB[MAX_LINELENGTH];
       GetCurrentTime(szB);
       char *p=strchr(szB,'\n');
       if(p)*p=0;
       S.Setf("Idn=File from %s: CreationTime=%s",XAUSW_ID,szB);
       if(!ReplaceLine("Idn=",(const char *)S))
         {PRINT_DEBUG("Data file integrity is corrupt\n");
	  return iRet;
	 }

       iRet=1;
    }
return iRet;
}
// ********************************************************
int TheFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG)return 0;
 
 FILE *f;

 String FN;

 if(szFile)FN.Set(szFile);
 else FN=FName;
 f=fopen((const char *)FN,"w");
 CHECK_FILE_POINTER_RETURN(f,(const char *)FN,0)

 FName=FN;

 int iRet=TheFile::SaveData(f);
 fclose(f);
 return iRet; 
}
// ********************************************************
int TheFile::SaveData(FILE *fOut)
{
  int i=SaveHeader();
  FText.Write(fOut);
  i|=TD->ExportData(fOut,iFLEnd);
  return i;
}
// ********************************************************
const char * TheFile::GetInfoText(const int iL, const int iWithE)
{
 if(iL==-1)return FText.FindPar("[FileHeader]","SampleID",iWithE);
 if(iL<0 || iL>FText.GetNLines())return szNoText;
 return FText.Line(iL,iWithE);
}
// *******************************************************
int TheFile::SetInfoText(const char *szT, const int iL)
{if(iL==-1)
   { String S;
     S.Setf("SampleID=%s",szT);
     if(!ReplaceLine("SampleID=",(const char *)S))
       {PRINT_DEBUG("Data file integrity is corrupt\n");
        return 0;
       }
     return 1;
   }  
 return 0;
}
// ******************************************************
// *******************************************************
// SxSFile
// *******************************************************
// ******************************************************
// Old Header SPEC_HEAD
//------------------------------------------------------
//       Text: A1 rhombo kV*mA 1, 4V*3,4A
//  Start Pos:27000
//    End Pos:32000
//No of steps:250
//      Delta:20
//No of spectra:2
//Timebase [s]:5.000000
//  Max Counts:0
//followed by intensity 5 values / line from Start Pos to End Pos in steps of
//Delta
//	 104            112             97             94             92
//...
//-------------------------------------------------------
//New Header SXS_HEAD
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
SxSFile::SxSFile(const char *szFile)
{DataFile::Init();
 
 int iT=CheckFileType(szFile,iFLEnd);

 if(iT==FtERROR || iT==FtILLG)
   {iFType=(enum FileType)iT;
    PRINT_DEBUG("Error analyzing file %s\n",szFile)
    return;
  }
 
 iFType=(FileType)iT;
 nSets=1;
 DataFile::Init(szFile,iT,"SxS Data File");

 if(!FindSets())
   {iFType=FtILLG;
    PRINT_DEBUG("Error finding data set %s\n",szFile)
    return;
   } 
 
  if(!SPar.lfPhi)
    {PRINT_DEBUG("Cannot calculate energies; SxS pars missing\n")
    } 
}
// ********************************************************
int SxSFile::FindSets(void)
{ nSets=0;
  switch(iFType)
    {case FtSXS_OLD:nSets=1;
	            return 1;
     case FtSXS_NEW:
	  {FILE *f=fopen((const char *)FName,"r");
	   if (!f)
             {PRINT_DEBUG("Error opening file %s\n",(const char *)FName)
              iFType=FtERROR;
              return 0;
             }

           int i;
           i=FText.Read(f," ",LineString::S_NUMBER);
            fclose(f);
            if(i<=0)return 0;
	   
            DeLim.Set(FText.GetDelim());
	  
 	   if(!FText.GetPar("[SXSFileHeader]","TotalMeas","%i",&i))
             {PRINT_DEBUG("Error opening file %s\n",(const char *)FName)
              iFType=FtERROR;
              return 0;
             }
           nSets=i;
	   return nSets;
	  }
  default: return 0;
 }

}
// ********************************************************
int SxSFile::ReadHeader(FILE *InStream, const int iSet)
{// Read Header from InStream
 // return 0: Illegal Header
 //        >0: Ok (no of columns)
 if(iSet<0 || iSet>nSets)return 0;
 char szB[MAX_LINELENGTH+1];
 String S;
 int i;

 switch(iFType)
   {case FtSXS_OLD:
//1.Line: Text
//2.Line: uStartPosition
//3.Line: uEndPosition
//4.Line: uSteps
//5.Line: uDeltaPos
//6.Line: unSpecs
//7.Line: fTimeBase
//8.Line: iMaxCounts

      for(i=0;i<OLD_SXS_HEADL;i++)
         {if(!fgets(szB,MAX_LINELENGTH,InStream))return 0;
          szB[MAX_LINELENGTH]=0;
          if(i==0)S.Set(szB);
          else S.Add(szB);
         }
      FText.Set(S);
      DeLim.Set(FText.GetDelim());
      return 3;
    case FtSXS_NEW:
     for(i=1;i<iSet;i++)
         {if(!FText.Read(InStream,"@EOH",LineString::S_CONTAINS))return 0;
         }
      if(iSet>1)i=FText.Read(InStream," ",LineString::S_ISCOL_BREAK);
      if(i<=0)return 0;
      i=FText.Read(InStream," ",LineString::S_NUMBER);
      if(iSet>1)FText.InsertLine("[SXSFileHeader]",0);
      DeLim.Set(FText.GetDelim());
      return i;
   default: return 0;
 }
 return 0;
}
// ********************************************************
int  SxSFile::ReadData(const int iSet)
{// Reads complete data from file
// Reads complete data from file
 // return  0 : Ok
 //         1 : Error opening file
 //         2 : Illegal Header or data
 //         3 : No memory 
 //        less than -1: Number of data sets and data points is decoded
 //                      in the negative number 
 //                      -( (SetNo)*SXS_SET_OFFS + NoOfSteps)
 //                      SXS_SET_OFFS see dfile.h
 // iSet==0 Reads incomplete Datafile to show progress of measurement
 // or last data set if file is complete
 
 if(iFType==FtERROR || iFType==FtILLG)return 0;
 if(iFType==FtSXS_OLD && iSet != 1)
   {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
    return 0;}

 FILE *f;
 int iCols;

  if ((f=fopen((const char *)FName,"r")) == NULL)
    {PRINT_DEBUG("Error opening file %s\n",(const char *)FName)
     iFType=FtERROR;
     return 1;
    }

  if(iFType==FtSXS_OLD)
   {if(!ReadHeader(f,1)){fclose(f);return 0;}
    unsigned uSteps=atoi(FText.Line(3));
    if(uSteps<5)return 2;
    iCols=3;
    TD=new CRData(iCols,uSteps);
    CHECK_POINTER_RETURN(TD,3)

    double pCh[5];
    unsigned i,j;
    for(i=0; i<uSteps; i+=5)
      {if(fscanf(f,"%lf  %lf  %lf  %lf  %lf",
	   &(pCh[0]),&(pCh[1]),&(pCh[2]),&(pCh[3]),&(pCh[4])) == EOF)
          {PRINT_DEBUG("Error reading columns from file %s\n",(const char *)FName)
           iFType=FtERROR;
           fclose(f);
	   return 2;
          }
       for(j=0; j<5; j++)(*TD)[1][i+j]=pCh[j];
      }
    fclose(f);
    unsigned uStartP=atoi(FText.Line(1));
    unsigned uDeltaP=atoi(FText.Line(4));
    if(uStartP>SXS_POS_MAX || uDeltaP<1 || uDeltaP>SXS_POS_MAX/2)return 2;

    for(i=0; i<uSteps; i++)
      {j=uStartP+i*uDeltaP;
       (*TD)[2][i]=j;
       (*TD)[0][i]=(j>ZERO_ORDER && SPar.lfPhi ? SpecFunc((double)j) : j) ;
      }
    nCols=iCols;

    String CS;
    float fTimeBase=atof(FText.Line(6));
    unsigned unSpecs=atoi(FText.Line(5));
    CS.Addf("Energy[eV]\n");
    CS.Addf("Counts [%4.1f * %2d sec]\n",fTimeBase,unSpecs);
    CS.Addf("POS [Steps]\n");

    delete ColID;
    ColID=new LineString(CS);

    TD->SetColX(0);
    TD->SetColY(1);
    TD->SetColZ(2);
    TD->SortData(0);
    return 0;
   }// FtSXS_OLD

 if(iFType==FtSXS_NEW)
   {int iiSet=(iSet==0 ? nSets : iSet);
    if(iiSet<1 || iiSet>nSets)
     {PRINT_DEBUG("Illegal set # %d / %d\n",iSet,nSets)
      fclose(f);
      return 0;
     }
    int i,iRet=0;
    unsigned j;
    nCols=0;
    iCols=ReadHeader(f,iiSet);
    if(!iCols){fclose(f); return 0;}
    rewind(f);
    LineString Columns;
    for(i=0;i<iiSet;i++)
         {if(!Columns.Read(f,"@EOH",LineString::S_CONTAINS))break;
         }
 
    if(i<=iiSet && iSet==0)
      {//Last set should be readed
       int j;
       for(j=1;j<i;i++)
          {if(!Columns.Read(f,"@EOH",LineString::S_CONTAINS))
              {fclose(f);return 0;}
           iRet=-(i+1)*SXS_SET_OFFS;
          }
       }

    if(!Columns.Read(f," ",LineString::S_ISCOL_BREAK))
      {PRINT_DEBUG("Error reading from file %s\n",(const char *)FName)
       iFType=FtERROR;
       fclose(f);
       return 1;
      }
    fclose(f);
    nDataLines=Columns.GetNLines();

    nCols=iCols;
    TD=new CRData(nCols,nDataLines);
    CHECK_POINTER_RETURN(TD,3)

    i=TD->ReadCols(Columns,nDataLines,nCols);
    if(i<=0)
      {PRINT_DEBUG("Error reading columns from file %s line %d\n",
                (const char *)FName,i)
       iFType=FtERROR;
       return 0;
      }

    unsigned uSteps=(unsigned) GetIntPar("[SXSFileHeader]","DataPoints",i);
    if(uSteps==0 || !i)return 2;

    unsigned uStartP=(unsigned) GetIntPar("[SXSFileHeader]","StartPos",i);
    if(uStartP>=SXS_POS_MAX || !i) return 2;

    unsigned uDeltaP=(unsigned) GetIntPar("[SXSFileHeader]","Delta",i);;
    if(uDeltaP>=SXS_POS_MAX << !i)return 2;

    if(iSet==0 && (int)uSteps!=nDataLines)
      {for(unsigned jj=j; (int)jj<nDataLines; jj++)
          {unsigned s=uStartP+jj*uDeltaP;
           (*TD)[2][jj]=s;
           (*TD)[0][jj]=(s>ZERO_ORDER && SPar.lfPhi ? SpecFunc((double)s) : s);
           (*TD)[1][jj]=0;
          }//for jj
       iRet-=nDataLines;
      }//if

    String S,SC;
    char szB[MAX_LINELENGTH+1];

      for(i=0;i<nCols;i++)
         {S.Setf("ColID%1d",i+1);
          if(!FText.GetPar(0,(const char *)S,"%S",szB))
            {PRINT_DEBUG("Data file integrity is corrupt (ColID)\n")
             break;
            }
          SC.Addf("%s\n",szB);
         }
 
     delete ColID;
     ColID=new LineString(SC);

    float fTimeBase=(float) GetDblPar("[SXSFileHeader]","Timebase",i);
    if(fTimeBase<0 || !i)return 2;

    unsigned unSpecs=(unsigned) GetIntPar("[SXSFileHeader]","NoofSpectra",i);
    if(unSpecs<=0 || !i)return 2;

    sprintf(szB,"Counts [%4.1f * %2d sec]",fTimeBase,unSpecs);
    ColID->ReplaceLine(szB,1);

    TD->SetColX(0);
    TD->SetColY(1);
    TD->SetColZ(2);
    TD->SortData(0);

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
 if(iFType==FtERROR || iFType==FtILLG)return 0;
 return AsciiFile::SaveData(szFile);
}
// ******************************************************
int  SxSFile::SaveData(FILE *OutStream)
{
 String S;
 if(iFType==FtSXS_NEW)
   {S.Setf("Idn=File from %s",XAUSW_ID);
    if(!ReplaceLine("Idn=",(const char *)S))PRINT_DEBUG("Data file integrity is corrupt\n")    
   }
 else 
   {S.Setf("File from %s",XAUSW_ID);
    FText.InsertLine((const char *)S,0);
   }

 return AsciiFile::SaveData(OutStream);
}
// *******************************************************
void SxSFile::SPrintInfo(String &B, const unsigned Flag)
{
 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 DataFile::SPrintInfo(B);

 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",(const char *)FId,iFType);
 if(Flag & PR_TEXT)
    {if(iFType==FtSXS_OLD)B.Addf("%s\n",FText.Line(0));
     else  B.Addf("%s\n",FindPar("[SXSFileHeader]","SampleIDId"));
    }
 if(Flag & PR_HEAD)
   {if(iFType==FtSXS_OLD)
      {B.Addf("Start Position: %s\n",FText.Line(1)); 
       B.Addf("  End Position: %s\n",FText.Line(2)); 
       B.Addf("         Steps: %s\n",FText.Line(3));
       B.Addf("         Delta: %s\n",FText.Line(4)); 
       B.Addf("  # of Spectra: %s\n",FText.Line(5)); 
       B.Addf(" Time base [s]: %s\n",FText.Line(6));
       B.Addf("   Max. Counts: %s\n",FText.Line(7));
      } 
    if(iFType==FtSXS_NEW)
      {if(Flag & PR_HEAD)
       FText.Write(B,1);

       if(Flag & PR_TIME && iFType==FtSXS_NEW)
         {B.Addf("%s\n",FindPar("[SXSFileHeader]","StartTime"));
          B.Addf("%s\n",FindPar("[SXSFileHeader]","EndTime"));
         }//if
      }
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
//   FilePar *FP=new FilePar(szB);
//   if(FP==0)return S;
//   if(FP->GetStatus()!=FOUND){delete FP;return S;}
  LineString *FP=new LineString();
  if(!FP->Read(szB)){delete FP;  return S;}
 //[deg] Spektrometerkonstante
  FP->GetPar("[SxsPars]","Phi","%lf",&(S.lfPhi));
 // Rowlandkreisradius [mm]
  FP->GetPar("[SxsPars]","RowlandR","%lf",&(S.lfD));
  //[mm]_?Spindelsteigung_RQ
  FP->GetPar("[SxsPars]","ScrewSlope","%lf",&(S.lfRQ));
  // 1/600 1/d:Gitterparameter  Linien/mm
  FP->GetPar("[SxsPars]","InvLinesPmm","%lf",&(S.lfSigma));
  // ? Einfallswinkel Alpha ?
  FP->GetPar("[SxsPars]","Alpha","%lf",&(S.lfAlpha));
  delete FP;
  return S;
}
// ********************************************************
const char * SxSFile::GetInfoText(const int iL, const int iWithE)
{
 if(iL==-1)
   {if(iFType==FtSXS_OLD)return FText.Line(0,iWithE);
     else  return FindPar("[SXSFileHeader]","SampleID",iWithE);
    }
 if(iL<0 || iL>FText.GetNLines())return szNoText;
 return FText.Line(iL,iWithE);
}
// *******************************************************
int SxSFile::SetInfoText(const char *szT, const int iL)
{if(iL==-1)
   {if(iFType==FtSXS_OLD)FText.ReplaceLine(szT,0);
    else
     {String S;
      S.Setf("SampleID=%s",szT);
      if(!ReplaceLine("SampleID=",(const char *)S))
        {PRINT_DEBUG("Data file integrity is corrupt\n");
         return 0;
        }
      }
    return 1;
   }  
 return 0;
}
// *******************************************************
int SxSFile::AddSets(void)
{
 if(iFType==FtSXS_OLD)return -1;
 SxSFile *SF=new SxSFile((const char *)FName);
 
 int i,iAddC,iErr;
 float fTBase=0;

 float fTimeBase=(float) GetDblPar("[SXSFileHeader]","Timebase",iErr);
 if(fTimeBase<0 || !iErr)return 0;

 unsigned unSpecs=(unsigned) GetIntPar("[SXSFileHeader]","NoofSpectra",iErr);
 if(!iErr)return 0;

 unsigned uStartP=(unsigned) GetIntPar("[SXSFileHeader]","StartPos",iErr);
 if(uStartP>SXS_POS_MAX || !iErr)return 0;

 unsigned uEndP=(unsigned) GetIntPar("[SXSFileHeader]","EndPos",iErr);
 if(uStartP>=uEndP || uEndP>SXS_POS_MAX || !iErr)return 0;

 unsigned uSteps=(unsigned) GetIntPar("[SXSFileHeader]","DataPoints",iErr);
 if(uSteps<1 || uSteps>SXS_POS_MAX || !iErr)return 0;

 for(i=0,iAddC=0;i<SF->GetNSets();i++)
    {
     if(i==0)
       {if(ReadData(i+1))
          {PRINT_DEBUG("Error reading data\n")
           return 0;
          }
        if(TD==NULL || TD->GetSteps()<=0)
          {PRINT_DEBUG("Data Block empty\n")
           return 0;
          }
 

     fTBase=unSpecs*fTimeBase;
      iAddC++;  
       }	  
     else
       {if(SF->ReadData(i+1))
          {PRINT_DEBUG("Error reading data\n")
           return 0;
          }

        if((int) uStartP != SF->GetIntPar("[SXSFileHeader]","StartPos",iErr) ||
           (int)uEndP != SF->GetIntPar("[SXSFileHeader]","EndPos",iErr) ||
	   (int)uSteps != SF->GetIntPar("[SXSFileHeader]","DataPoints",iErr))continue;
 	(*TD)[1].BasicCalc('+',(*SF->GetCRData())[1]);
	fTBase+=SF->GetIntPar("[SXSFileHeader]","NoofSpectra",iErr) * 
	        SF->GetDblPar("[SXSFileHeader]","Timebase",iErr);
	iAddC++;
       } 
    }    

 //TD->ExchgCol(0,2);   
 TD->SetColX(0);
 TD->SetColY(1);
 TD->SetColZ(2);


 String S;
 S.Set("NoofMeas=1");
 if(!ReplaceLine("NoofMeas=",(const char *)S))
   {PRINT_DEBUG("Data file integrity is corrupt\n");
    return 0;
   }

 S.Set("TotalMeas=1");
 if(!ReplaceLine("TotalMeas=",(const char *)S))
   {PRINT_DEBUG("Data file integrity is corrupt\n");
    return 0;
   }

 S.Set("NoofSpectra=1");
 if(!ReplaceLine("NoofSpectra=",(const char *)S))
   {PRINT_DEBUG("Data file integrity is corrupt\n");
    return 0;
   }
	
 S.Setf("Counts [%*.0f sec]",(int)log10(fTBase)+1,fTBase);
 ColID->ReplaceLine((const char *)S,1); 

 S.Setf("ColID2=%s",ColID->Line(1));
 if(!ReplaceLine("ColID2=",(const char *)S))
   {PRINT_DEBUG("Data file integrity is corrupt\n");
    return 0;
   }

 S.Setf("Timebase=%f",fTBase);
 if(!ReplaceLine("Timebase=",(const char *)S))
   {PRINT_DEBUG("Data file integrity is corrupt\n");
    return 0;
   }
  
 S.Setf("%s; %d of %d spectra added",
        FText.FindPar("[SXSFileHeader]","SampleID"),iAddC,SF->GetNSets() );
 if(!ReplaceLine("SampleID=",(const char *)S))
   {PRINT_DEBUG("Data file integrity is corrupt\n");
    return 0;
   }

 delete SF;
 return iAddC;
}
// ******************************************************
SxSFile::~SxSFile()
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
     PRINT_DEBUG("Error analyzing file %s\n",szFile)
     return;
    }

  DataFile::Init(szFile,iT,XDifTYPES[iT-XDifGROUP]);
  
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
 // return -1: Illegal Header
 //       >0: Ok # of data points
// ******************************************************
// Header of XDif File see above

 if(!(iFType==FtXDIF_DIF || iFType==FtXDIF_RAW))return 1;
 int i, iRet=-1;
 String S;
 switch (iFType)
  {case    FtERROR: return 0;
   case     FtILLG: return 0;
   case    FtASCII: return 0;
   case FtXDIF_RAW: fseek(InStream,0,SEEK_SET);
                    i=fread(&Head,sizeof(Head),1,InStream);
		    nSets=1;
                    if(!i)return -1;
                    iRet=Head.nSteps;
                    S.Add("2 Theta\n");
                    S.Addf("Counts/%5.2f [s]\n",Head.fTBase);
                    break;
   case FtXDIF_DIF: fseek(InStream,0,SEEK_SET);
                    i=fread(&Head,sizeof(Head),1,InStream);
                   if(!i)return -1;
		    nSets=Head.nRanges;
                    iRet=Head.nPeaks;
                    S.Add("d\n");
                    S.Add("Intensity\n");
                    break;
   default: return -1;
  }

 delete ColID;
 ColID=new LineString(S);

  String B;
  char szB[80];
  GetCurrentTime(szB);
 
  B.Setf("File from %s %s",XAUSW_ID,szB);
  B.Addf("Sample ID: %s\n",Head.SampleID);
  B.Addf("Key: %-4.4s\n",Head.Key);
  B.Addf("No of steps: %-ld\n",(long)Head.nSteps);
  B.Addf("Time base: %-6.2f\n",Head.fTBase);
  B.Addf("Step size: %-6.3f\n",Head.fStepSize);
  B.Addf("DACO flag: %0lx\n",(unsigned)Head.DACOflag);
  B.Addf("Sample Pos.: %-lu\n",Head.SamplePos);
  B.Addf("2 THETA Start: %-7.3f\n",Head.fTheta2Start);
  B.Addf("THETA Start: %-7.3f\n",Head.fThetaStart);
  B.Addf("CHI Start: %-7.3f\n",Head.fChiStart);
  B.Addf("PHI Start: %-7.3f\n",Head.fPhiStart);
  B.Addf("Sample ID: %-32.32s\n",Head.SampleID);
  B.Addf("Wavelength Ka1: %-7.5f\n",Head.fWLKa1);
  B.Addf("Wavelength Ka2: %-7.5f\n",Head.fWLKa2);
  B.Addf("No of peaks: %-ld\n",(long)Head.nPeaks);
  B.Addf("No of ranges: %-ld\n",(long)Head.nRanges);

  FText.Set(B);

 return iRet;
}
// ******************************************************
int XDifFile::ReadData(const int iSet)
{// Reads complete data from file
 // return 2 : Error opening file
 //        1 : Illegal Header
 //        0 : Ok

 if(iFType==FtERROR || iFType==FtILLG )return 1;

 FILE *f;
 if(iSet!=1 && iFType==FtXDIF_RAW)
  {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
   return 1;}

  if ((f=fopen((const char *)FName,"r")) == NULL)
    {PRINT_DEBUG("Error opening file %s\n",(const char *)FName)
     iFType=FtERROR;
     return 2;
    }

  int i=ReadHeader(f);
  if(i<=0)
    {fclose(f);
     PRINT_DEBUG("Error reading file header %s\n",(const char *)FName)
     iFType=FtERROR;
     fclose(f);
     return 1;
    }

 switch(iFType)
  {case FtXDIF_RAW:
   { TD=new CRData(2,Head.nSteps);
     CHECK_POINTER_RETURN(TD,2)

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
     break;
    }
   case FtXDIF_DIF:
   {if(iSet <=0 || iSet>Head.nRanges)
      {PRINT_DEBUG("Illegal data set %d / %d\n",iSet, Head.nRanges)
       return 1;}
    int ii,k,steps;
    for(ii=0,steps=0;ii<iSet;ii++)
       {k=sizeof(Head)*ii+2*sizeof(FLOAT4)*steps;
        fseek(f,k,SEEK_SET);
        fread(&Head,sizeof(Head),1,f);
        steps+=Head.nPeaks;
       }/*for(ii)*/
      
       
     TD=new CRData(2,Head.nPeaks);
     CHECK_POINTER_RETURN(TD,2)

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
     break;
   }
   default: return 1;
  }
 
 fclose(f);
 return 0;
}
// *******************************************************
void XDifFile::SPrintInfo(String &B, const unsigned Flag)
{
 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(B);

 if(Flag & PR_TEXT)
   B.Addf("Sample ID: %s\n",Head.SampleID);
 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",(const char *)FId,iFType);

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
{
 SetEndChr();
 FText.Write(OutStream);
}
// ********************************************************
int XDifFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG)return 0;
 
 FILE *f;
 String FN;

 if(szFile)FN.Set(szFile);
 else FN=FName;
 f=fopen((const char *)FN,"w");
 CHECK_FILE_POINTER_RETURN(f,(const char *)FN,0)

 FName=FN;
 PRINT_DEBUG("File %s saved as ASCII; not in original data format\n",
             (const char *)FName)
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
const char * XDifFile::GetInfoText(const int iL, const int iWithE)
{
 if(iL==-1)return Head.SampleID;
 if(iL<0 || iL>FText.GetNLines())return szNoText;
 return FText.Line(iL,iWithE);
}
// *******************************************************
int XDifFile::SetInfoText(const char *szT, const int iL)
{if(iL==-1){strncpy(Head.SampleID,szT,31);
            return 1;
           }
 else return FText.ReplaceLine(szT,iL);
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
SplineFile::SplineFile(const char *szFile) : AsciiFile()
{


 DataFile::Init();
 
 int iT=CheckFileType(szFile,iFLEnd);

 if(iT==FtERROR || iT==FtILLG)
   {iFType=(enum FileType)iT;
    PRINT_DEBUG("Error analyzing file %s\n",szFile)
    return;
  }
 
 iFType=(FileType)iT;
 nSets=1;
 DataFile::Init(szFile,iT,"Spline File");

}
// ********************************************************
int  SplineFile::ReadData(const int iSet)
{// Reads complete data from file
 // return  1 : Error opening file
 //         2 : Illegal Header
 //         3 : No memory
 //         0 : Ok

 if(iFType==FtERROR || iFType==FtILLG)return 2;
 //int iErr=0;

 if(iSet!=1)
  {PRINT_DEBUG("Unsupported # of data sets %d\n",iSet)
   return 2;}


 FILE *f=fopen((const char *)FName,"r");
 CHECK_FILE_POINTER_RETURN(f,(const char *)FName,1)


 char szB[MAX_LINELENGTH];
 String S;
 if(!fgets(szB,MAX_LINELENGTH,f))return 2;
 szB[MAX_LINELENGTH]=0;
 if(strstr(szB,"SPLINETABLE")!=szB)return 2;
 S.Set(szB);

 if(!fgets(szB,MAX_LINELENGTH,f))return 2;
 szB[MAX_LINELENGTH]=0;

 int nPoints=atoi(szB);
 if(nPoints<3)return 2;

 S.Add(szB);
 
 FText.Set(S);
 DeLim.Set(FText.GetDelim());

 if(TD)delete TD;
 TD = new CRData(3,nPoints);
 if(! TD || TD->GetSteps()<=3)
   {PRINT_DEBUG("Error allocating memory\n")
    return 0;
   }

 MDATA x,y;
 int j;
 for(j=0; j<nPoints; j++)
    { fscanf(f,"%lf %lf", &x,&y);
      (*TD)[0][j]=x;
      (*TD)[1][j]=y;
    }
 TD->SetColX(0);
 TD->SetColY(1);
 TD->SetColZ(2);
 
 TD->Spline();
 TD->SortData(0);
 TD->NewMinMax();


 String CS;
 CS.Addf("X");
 CS.Addf("Y");
 CS.Addf("Spline");
 
 delete ColID;
 ColID=new LineString(CS);

 fclose(f);
 return 1;
}
// *******************************************************
int  SplineFile::SaveData(const char *szFile)
{
 if(iFType==FtERROR || iFType==FtILLG)return 0;
 return AsciiFile::SaveData(szFile);
}

// *******************************************************
int  SplineFile::SaveData(FILE *fOut)
{
 return AsciiFile::SaveData(fOut);
}
// *******************************************************
void SplineFile::SPrintInfo(String &B, const unsigned Flag)
{
 if(iFType==FtERROR || iFType==FtILLG)
   {B.Addf("Illegal file type %d ;cannot print data\n",iFType);
    return;
   }

 if(Flag & PR_BASIC)
    DataFile::SPrintInfo(B);

 if(Flag & PR_TYPE)
   B.Addf("       File Type: %s (%1d)\n",(const char *)FId,iFType);
 if(Flag & PR_TEXT)
   {B.Addf("%s\n",FText.Line(0));
   }
 if(Flag & PR_COLS)SPrintColInfo(B);
}
// ******************************************************
const char * SplineFile::GetInfoText(const int iL, const int iWithE)
{return AsciiFile::GetInfoText(iL,iWithE);
}
// ******************************************************
int SplineFile::SetInfoText(const char *szT, const int iL)
{return AsciiFile::SetInfoText(szT,iL);
}
// ******************************************************
SplineFile::~SplineFile()
{ }
// ******************************************************
#endif // MORE_FILETYPES
