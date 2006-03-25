#include <string.h>
#include <stdio.h>

#ifndef DFILE_H
#define DFILE_H

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

#ifndef STRINGS_H
#include "strings.h"
#endif

#ifndef STDINC_H
#include "stdinc.h"
#endif

#ifndef CDATA_H
#include "cdata.h"
#endif

#ifndef AUSW_H
#include "ausw.h"
#endif

#ifndef FILEPAR_H
#include "filepar.h"
#endif

#ifndef THEFUNC_H
#include "thefunc.h"
#endif

// $Id: dfile.orig.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: dfile.orig.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

//#define FILEDATA_DEBUG

#define   MAX_IDLENGTH 30
#define       MAXTYPES 10
//#define MAX_LINELENGTH 256

// ***File Header definitions***

#define SXS_FTYPE_NONE -1
#define SXS_FTYPE_OLD   0
#define SXS_FTYPE_NEW   1

// **********************************************************
// Old SXS File Header
typedef struct spec_file_header
   {	char szText[82];
    unsigned uStartP; //Start position
    unsigned uEndP;   //End position
    unsigned uSteps;  //no of steps
    unsigned uDeltaP; // number of spectra
    unsigned unSpecs; //time base in sec
       float fTimeBase;
	 int iMaxCounts;
   } SPEC_HEAD;

// ***********************************************************
#define MAX_SXSHEAD_LENGTH 81
typedef struct sxs_file_header
    {    char szHead[MAX_SXSHEAD_LENGTH+1];
         char szSTime[30];
	 char szETime[30];
          int nCols;
      char ** szCID;
         int iFType;

       int iGreating; //Lines / mm to identify greating
     float fSampVac;  //pressure in sample chamber in [torr]
     float fXI;     //x-ray current in mA
     float fXV;     //x-ray voltage in kV
     float fUFil;   //Filament Voltage in V
     float fIFil;   //Filament current in A
     float fUWehn;  //wehnelt voltage in V
     float fUCEM;  //Voltage on channeltron in kV
      char szCEMId[16]; //Channeltron ID
     float fUPhoto; //Voltage on photo cathode in V
     float fAGain;  //Gain of amplifier for channeltron; helipot position 1-10
      char szPhoto[16]; //Type of Photo cathode
      char  szOperator[16]; // Name of experimentator

      int iMeas;	 //no of measurement in file
      int iTotMeas;	 //total no of measurments
      char szFileName[MAXPATH+1];

     SPEC_HEAD SpHead;
     } SXS_HEAD;

// ***********************************************************
// SXS Pars to calculate energy from steps
typedef struct spect_pars
      {double lfPhi;  //[deg] Spektrometerkonstante
       double lfD; // Rowlandkreisradius
       double lfRQ; // ??
       double lfSigma; // 1/600 1/d:Gitterparameter  Linien/mm
       double lfAlpha; // ?

       //spect_pars() {lfPhi=0; lfD=0; lfRQ=0; lfSigma=0; lfAlpha=0;}
      } SXS_PAR;
// ********************************************************

#if defined (LINUX)
#define INT4 int
#define UNSIGNED4 unsinged
#define FLOAT4 float
#else
#define INT4 long int
#define UNSIGNED4 long unsinged
#define FLOAT4 float
#endif

typedef struct dif_file_header
   {  char Key[4];
      INT4 nSteps;
    FLOAT4 fTBase;
    FLOAT4 fStepSize;
      INT4 DACOflag;
      INT4 SamplePos;
    FLOAT4 fTheta2Start,fThetaStart,fChiStart,fPhiStart;
      char SampleID[32];
    FLOAT4 fWLKa1,fWLKa2;
      char Unused[68];
      INT4 nPeaks,nRanges;
      } RAW_DIF_HEAD;

// *******************************************************

#define SET_START 0x01
#define   SET_END 0x02

// #define PR_TEXT  0x0001U
// #define PR_TYPE  0x0002U
// #define PR_HEAD  0x0004U
// #define PR_COLS  0x0008U

#define PR_BASIC 0x0001U
#define PR_TEXT  0x0002U
#define PR_TYPE  0x0004U
#define PR_HEAD  0x0008U
#define PR_TIME  0x0010U
#define PR_COLS  0x0020U
#define PR_ALL  UINT_MAX

#define STRICT  1
#define MILD    2

#define SXS_POS_MIN 1U
#define SXS_POS_MAX 55000U
#define  ZERO_ORDER 7900U

#define SXS_SET_OFFS 1000000
#define GET_SXS_SET(I) ((int) ( -(I)/SXS_SET_OFFS) )
#define GET_SXS_POINTS(I) (-(I%SXS_SET_OFFS))

// Group IDNumbers for identifying files
enum FileGroup {    NoGROUP = 0,    TheGROUP = 10, SxSGROUP = 20,
                  XDifGROUP = 30, TransGROUP = 40, He3GROUP = 50,
                SplineGROUP = 60};

enum FileType {FtNONE=-2, FtERROR, FtILLG, FtASCII,

                   FtTHE_RCP = TheGROUP + 1,
                 FtTHE_CAPV1 = TheGROUP + 2,
               FtTHE_MEMBRAN = TheGROUP + 3,
                FtTHE_CAPOLD = TheGROUP + 4,

                  FtSXS_OLD = SxSGROUP + 1,
                  FtSXS_NEW = SxSGROUP + 2,

                  FtXDIF_RAW = XDifGROUP + 1,
                  FtXDIF_DIF = XDifGROUP + 2,

               FtTRANS_DATAP = TransGROUP + 1,
                FtTRANS_MEAS = TransGROUP + 2,

	          FtHE3_MEAS = He3GROUP + 1,

                    FtSPLINE = SplineGROUP + 1};

// ********************************************************
int CheckFileType(const char *szFile, enum EndLine &LineT);

int GetTextLines(FILE *f, const int iStart,
                 const char * szDelim, char* szA=NULL);

int CheckAsciiFile(const char *szFile, int &nTextL, int &nDataL, 
                   const int CheckMode=STRICT);
// ********************************************************
class DataFile
{
 protected:
   char szName[MAXPATH+1];
   char szID[MAX_IDLENGTH+1];

   enum FileGroup iFGroup;
   enum FileType  iFType;
   enum EndLine   iFLEnd;
   char szE[3];
  
  char ** szCID;
      int nCols;

  CRData *TD;
    int nSets;

  int nTxtLines;

   int  Init(const char *szFile,const int iG, const int iT, const char *szI);
   int  AllocColID(const int nC);
  void  DeleteColID(void);
   int  CheckType(void);
 
   void SetEndChr(void);

  DataFile(void);

 public:

    const char * GetID(void)       const {return &szID[0];}
    const char * GetFileName(void) const {return &szName[0];}
  enum FileGroup GetFGroup(void)   const {return iFGroup;}
   enum FileType GetFType(void)    const {return iFType;}
    enum EndLine GetLineType(void) const {return iFLEnd;}

             int GetError()        const {return iFType==FtERROR;}
        CRData * GetCRData(void)   const {return TD;} 
	     int GetNSets(void)    const {return nSets;}
	     int GetNoTxtL(void)   const {return nTxtLines;}
       
    const char * GetColID(int iC)  const ;
             int GetCols(void)     const {return nCols;}

             int SetColID(const char *szI, const int iC);
             int SetColRange(const int iCS, const int iCE);
             int InsertColID(const char *szN, const int iCol);
            void SetLineType(const enum EndLine L){iFLEnd=L;}  
        
          char * SPrintInfo(char * szB)    const; 
            void SPrintInfo(String& B)    const; 
          char * SPrintColInfo(char *szB)  const;
            void SPrintColInfo(String &B) const;
             int SaveData(const char *szText, FILE *f);
             int SaveData(const char *szText, const char *szFile);
             int SaveSplineTable(FILE *f, const int iX, const int iY);

 virtual              int HasManySets(void)=0;
 virtual              int ReadData(const int iSet=1)=0;
 virtual              int SaveData(FILE *fOut)=0;
 virtual              int SaveData(const char *szFile=NULL)=0;
 virtual           char * SPrintInfo(char *szB, const unsigned Flag)=0;
 virtual            void  SPrintInfo(String &B, const unsigned Flag)=0;
 virtual     const char * GetInfoText(const int iL=-1)=0;
 virtual              int SetInfoText(const char *szT, const int iL=-1)=0;
                     // If iL==-1 default infotext line is handled
                     // Else other textlines if available
 virtual              int AddSets(void)=0;

 virtual     const THE_HEAD * GetTheHead(void)=0;
 virtual     const SXS_HEAD * GetSxSHead(void)=0;
 virtual const RAW_DIF_HEAD * GetDifHead(void)=0;

virtual ~DataFile();

};
// *******************************************************
class TheFile : public DataFile
{
 private:
  THE_HEAD Head;

 int ReadHeader(FILE *InStream);
 int SaveHeader(FILE *OutStream, const int nCS=0, const int nCE=-1,const int nR=-1);

 public:
 TheFile(const char *szFile);

 inline  void   SetFileTime(const int iStartEnd);
 inline float   GetLength(void) const {return Head.fL0;}
 inline float   GetField(void)  const {return Head.fMField;}

 //int CalculateGap(const double lfCDiam, const double lfCGap=-1.0);
 //int CalculateDl_l(const double lfL0);
//VIRTUALS
              int HasManySets(void)   {return 0;}
              int ReadData(const int iSet=0);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
           char * SPrintInfo(char * szB, const unsigned Flag);
            void SPrintInfo(String &B, const unsigned Flag);
     const char * GetInfoText(const int iL=-1);
              int SetInfoText(const char *szT, const int iL=-1);
              int AddSets(void){return -1;}

     const THE_HEAD * GetTheHead(void){return &Head;}
     const SXS_HEAD * GetSxSHead(void){return NULL;}
 const RAW_DIF_HEAD * GetDifHead(void){return NULL;}

                  ~TheFile();

};
// *******************************************************
class AsciiFile : public DataFile
{
 private:
  char ** szTextLine;
      int nLines;

  int iLStart,iLEnd;
  int nDataLines;

 public:
  AsciiFile(const char *szFile, const int nL=0, const int CheckMode=STRICT );
  AsciiFile(const char *szFile, CRData *td, const int nL, char **szLs);

 const char * GetLine(const int iL) const;
          int SetLine(const char *szT, const int iL);

//VIRTUALS
              int HasManySets(void)   {return 0;}
              int ReadData(const int iSet=1);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
           char * SPrintInfo(char * szB, const unsigned Flag);
             void SPrintInfo(String &B, const unsigned Flag);
     const THE_HEAD * GetTheHead(void){return NULL;}
     const SXS_HEAD * GetSxSHead(void){return NULL;}
 const RAW_DIF_HEAD * GetDifHead(void){return NULL;}
     const char * GetInfoText(const int iL=-1) {return GetLine(iL==-1?0:iL);}
              int SetInfoText(const char *szT, const int iL=0)
                                     {return SetLine(szT,iL);}
              int AddSets(void){return -1;}
                  ~AsciiFile();

};
// *******************************************************
class SxSFile : public DataFile
{
 private:

  SXS_HEAD Head;

 protected:
 double SpecFunc(double lfPos);
    int FindSets(void);
    int ReadHeader(FILE *InStream, FilePar *FP=NULL);

 public:
  SxSFile(const char *szFile);
  static const SXS_PAR SPar;
        static SXS_PAR ReadPars(const char * szFile);
//VIRTUALS
              int HasManySets(void)   {return 1;}
              int ReadData(const int iSet=0);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
           char * SPrintInfo(char * szB, const unsigned Flag);
             void SPrintInfo(String &B, const unsigned Flag);
    const THE_HEAD * GetTheHead(void){return NULL;}
     const SXS_HEAD * GetSxSHead(void){return &Head;}
 const RAW_DIF_HEAD * GetDifHead(void){return NULL;}
     const char * GetInfoText(const int iL=-1);
              int SetInfoText(const char *szT, const int iL=-1);
              int AddSets(void);
                  ~SxSFile();

};
// *********************************************************
class TransFile : public AsciiFile
{

 public:
  TransFile(const char *szFile);

//VIRTUALS
              int HasManySets(void)   {return 0;}
              int ReadData(const int iSet=0);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
           char * SPrintInfo(char * szB, const unsigned Flag);
             void SPrintInfo(String &B, const unsigned Flag);
     const THE_HEAD * GetTheHead(void){return NULL;}
     const SXS_HEAD * GetSxSHead(void){return NULL;}
 const RAW_DIF_HEAD * GetDifHead(void){return NULL;}
     const char * GetInfoText(const int iL=-1);
              int SetInfoText(const char *szT, const int iL=-1);
              int AddSets(void){return -1;}
                  ~TransFile();

};
// *******************************************************
class XDifFile : public DataFile
{
 private:
  RAW_DIF_HEAD Head;
  int ReadRange;

 void SaveHeader(FILE *OutStream);

 public:
 XDifFile(const char *szFile);
 void SetRange(const int iR);
  int GetRange(void) const {return ReadRange;}
  int GetMaxRange(void) const {return Head.nRanges;}
  int ReadHeader(FILE *InStream);

//VIRTUALS
              int HasManySets(void) {return (iFType==FtXDIF_RAW ? 0:1);}
              int ReadData(const int iSet=0);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
           char * SPrintInfo(char * szB, const unsigned Flag);
             void SPrintInfo(String &B, const unsigned Flag);
    const THE_HEAD * GetTheHead(void){return NULL;}
     const SXS_HEAD * GetSxSHead(void){return NULL;}
 const RAW_DIF_HEAD * GetDifHead(void){return &Head;}
     const char * GetInfoText(const int iL=-1);
              int SetInfoText(const char *szT, const int iL=-1);
              int AddSets(void){return -1;}
                  ~XDifFile();
};
// *******************************************************
class SplineFile : public DataFile
{
 private:
  char szTextLine[MAX_LINELENGTH+1];
  int nDataLines;

 public:
  SplineFile(const char *szFile);
  int GetnPoints(){return nDataLines;}
//VIRTUALS
              int HasManySets(void)  {return 0;}
              int ReadData(const int iSet=0);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
           char * SPrintInfo(char * szB, const unsigned Flag);
             void SPrintInfo(String &B, const unsigned Flag);
     const THE_HEAD * GetTheHead(void){return NULL;}
     const SXS_HEAD * GetSxSHead(void){return NULL;}
 const RAW_DIF_HEAD * GetDifHead(void){return NULL;}
     const char * GetInfoText(const int iL=0);
              int SetInfoText(const char *szT, const int iL=0);
              int AddSets(void){return -1;}
                  ~SplineFile();

};
// *******************************************************

#endif



