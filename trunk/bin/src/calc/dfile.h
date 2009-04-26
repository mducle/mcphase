//File: dfile.h
// $Id: dfile.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: dfile.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef DFILE_H
#define DFILE_H 1

#include <string.h>
#include <stdio.h>

#include "stdfunc.h"
#include "strings.h"
#include "stdinc.h"
#include "cdata.h"
#include "thefunc.h"

#define MORE_FILETYPES

#define PR_BASIC 0x0001U
#define PR_TEXT  0x0002U
#define PR_TYPE  0x0004U
#define PR_HEAD  0x0008U
#define PR_TIME  0x0010U
#define PR_COLS  0x0020U
#define PR_ALL  UINT_MAX


// Group IDNumbers for identifying files
enum FileGroup {    NoGROUP = 0,    TheGROUP = 10, SxSGROUP = 20,
                  XDifGROUP = 30,   He3GROUP = 40, SplineGROUP = 50};

enum FileType {FtNONE=-2, FtERROR, FtILLG, FtASCII,

                   FtTHE_RCP = TheGROUP + 1,
                 FtTHE_CAPV1 = TheGROUP + 2,
               FtTHE_MEMBRAN = TheGROUP + 3,
                FtTHE_CAPOLD = TheGROUP + 4,

                  FtSXS_OLD = SxSGROUP + 1,
                  FtSXS_NEW = SxSGROUP + 2,

                  FtXDIF_RAW = XDifGROUP + 1,
                  FtXDIF_DIF = XDifGROUP + 2,

//	          FtHE3_MEAS = He3GROUP + 1,

                    FtSPLINE = SplineGROUP + 1};
// *******************************************************
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


// ********************************************************
int CheckFileType(const char *szFile, enum EndLine &LineT);
// ********************************************************
class DataFile
{
 protected:
   String FName;
   String FId;
   String DeLim;

   enum FileType  iFType;
   enum EndLine   iFLEnd;
  
   LineString *ColID;
   int nCols;

  CRData *TD;
    int nSets;

  LineString FText;
//  LineString *FPars;

   void Init(void);
   int  Init(const char *szFile,const int iT, const char *szI);
  int  CheckType(void);
 
  void SetEndChr(void);
   int ReplaceLine(const char *SearchT, const char *WithT);

  DataFile(void){Init();}

 public:

    const char * GetID(void)       const {return (const char *)FId;}
    const char * GetFileName(void) const {return (const char *)FName;}
             int GetFGroup(void)   const {return (int)iFType/10*10;}
   enum FileType GetFType(void)    const {return iFType;}
    enum EndLine GetLineType(void) const {return iFLEnd;}

             int GetError()        const {return iFType==FtERROR;}
        CRData * GetCRData(void)   const {return TD;} 
	     int GetNSets(void)    const {return nSets;}
	     int GetNoTxtL(void)   const {return FText.GetNLines();}
       
    const char * GetColID(int iC)  const ;
             int GetCols(void)     const {return nCols;}
    LineString & GetHead(void) {return FText;}

             int SetColID(const char *szI, const int iC);
             int SetColRange(const int iCS, const int iCE);
             int InsertColID(const char *szN, const int iCol);
            void SetLineType(const enum EndLine L){iFLEnd=L;}  
        
            void SPrintInfo(String& B)    const; 
            void SPrintColInfo(String &B) const;
             int SaveData(const char *szText, FILE *f);
             int SaveData(const char *szText, const char *szFile);
             int SaveSplineTable(FILE *f, const int iX, const int iY);

    int GetIntPar(const char *szTopic, const char *szName, int &iR)
                   {int i;
		    iR= FText.GetPar(szTopic,szName,"%i",&i);
		    return i;
		   }
 double GetDblPar(const char *szTopic, const char *szName, int &iR)
                   {double d;
		    iR= FText.GetPar(szTopic,szName,"%lf",&d);
		    return d;
		   }

   const char *FindPar(const char *szTopic, const char *szName, const int iWithE=0)
             {return FText.FindPar(szTopic,szName,iWithE);}


 virtual              int HasManySets(void)=0;
 virtual              int ReadData(const int iSet=1)=0;
 virtual              int SaveData(FILE *fOut)=0;
 virtual              int SaveData(const char *szFile=NULL)=0;
 virtual             void SPrintInfo(String &B, const unsigned Flag)=0;
 virtual     const char * GetInfoText(const int iL=-1, const int iWithE=0)=0;
 virtual              int SetInfoText(const char *szT, const int iL=-1)=0;
                     // If iL==-1 default infotext line is handled
                     // Else other textlines if available
 virtual              int AddSets(void)=0;
 virtual              int AddLogLine(const char * szText)=0;
 virtual    RAW_DIF_HEAD * GetDifHead()=0;
 virtual ~DataFile();

};
// *******************************************************
class AsciiFile : public DataFile
{
 protected:
  int nDataLines;
  AsciiFile(void){DataFile::Init();}
 public:
  AsciiFile(const char *szFile);

 const char * GetLine(const int iL, const int iWithE=0);
          int SetLine(const char *szT, const int iL);

//VIRTUALS
              int HasManySets(void)   {return 0;}
              int ReadData(const int iSet=1);
              int SaveData(const char *szFile=0);
              int SaveData(FILE *fOut);
             void SPrintInfo(String &B, const unsigned Flag);
     const char * GetInfoText(const int iL=-1, const int iWithE=0);
              int SetInfoText(const char *szT, const int iL=0);
              int AddSets(void){return -1;}
              int AddLogLine(const char * szText) {return (szText?0:0);}
   RAW_DIF_HEAD * GetDifHead(){return 0;}
                  ~AsciiFile();

};
// *******************************************************
// *******************************************************
// ***File Header definitions***
#ifdef MORE_FILETYPES

#define SXS_POS_MIN 1U
#define SXS_POS_MAX 55000U
#define  ZERO_ORDER 7900U

#define SXS_SET_OFFS 1000000
#define GET_SXS_SET(I) ((int) ( -(I)/SXS_SET_OFFS) )
#define GET_SXS_POINTS(I) (-(I%SXS_SET_OFFS))

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
class TheFile : public AsciiFile
{
 protected:
 int SaveHeader(void);
  
 public:
 TheFile(const char *szFile);

//VIRTUALS
              int HasManySets(void)   {return 0;}
              int ReadData(const int iSet=0);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
           char * SPrintInfo(char * szB, const unsigned Flag);
             void SPrintInfo(String &B, const unsigned Flag);
     const char * GetInfoText(const int iL=-1, const int iWithE=0);
              int SetInfoText(const char *szT, const int iL=-1);
              int AddSets(void){return -1;}
              int AddLogLine(const char * szText) {return (szText?0:0);}
   RAW_DIF_HEAD * GetDifHead(){return 0;}
        	  ~TheFile();

};
// *******************************************************
#define OLD_SXS_HEADL 8
#define PARAMETER_FILE "ausw.conf"
 class SxSFile : public AsciiFile
{
 private:

//  SXS_HEAD Head;

 protected:
 double SpecFunc(double lfPos);
    int FindSets(void);
    int ReadHeader(FILE *f, const int iSet);

 public:
  SxSFile(const char *szFile);
  static const SXS_PAR SPar;
        static SXS_PAR ReadPars(const char * szFile);
//VIRTUALS
              int HasManySets(void)   {return 1;}
              int ReadData(const int iSet=0);
              int SaveData(const char *szFile=NULL);
              int SaveData(FILE *fOut);
             void SPrintInfo(String &B, const unsigned Flag);
     const char * GetInfoText(const int iL=-1, const int iWithE=0);
              int SetInfoText(const char *szT, const int iL=-1);
              int AddSets(void);
              int AddLogLine(const char * szText) {return (szText?0:0);}
   RAW_DIF_HEAD * GetDifHead(){return 0;}
        	  ~SxSFile();

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
             void SPrintInfo(String &B, const unsigned Flag);
     const char * GetInfoText(const int iL=-1, const int iWithE=0);
              int SetInfoText(const char *szT, const int iL=-1);
              int AddSets(void){return -1;}
              int AddLogLine(const char * szText) {return (szText?0:0);}
   RAW_DIF_HEAD * GetDifHead(){return &Head;}
                  ~XDifFile();
};
// *******************************************************
class SplineFile : public AsciiFile
{
 public:
  SplineFile(const char *szFile);

//VIRTUALS
	      int HasManySets(void)  {return 0;}
	      int ReadData(const int iSet=0);
	      int SaveData(const char *szFile=NULL);
	      int SaveData(FILE *fOut);
	     void SPrintInfo(String &B, const unsigned Flag);
     const char * GetInfoText(const int iL=-1, const int iWithE=0);
	      int SetInfoText(const char *szT, const int iL=0);
	      int AddSets(void){return -1;}
              int AddLogLine(const char * szText) {return (szText?0:0);}
   RAW_DIF_HEAD * GetDifHead(){return 0;}
		  ~SplineFile();
};
// // *******************************************************

#endif //MORE_FILETYPES

#endif //DFILE_H

   
