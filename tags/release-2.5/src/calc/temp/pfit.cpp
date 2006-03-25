// $Id: pfit.cpp,v 1.3 2001/10/04 12:17:34 blau Exp blau $
// $Log: pfit.cpp,v $
// Revision 1.3  2001/10/04 12:17:34  blau
// *** empty log message ***
//
// Revision 1.1  2001/09/11 19:34:21  root
// Initial revision
//

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "dfile.h"
#include "stdfunc.h"
#include "simplex.h"

#define MAN_PATH "AUSW_MANPATH"
#define MAX_PARLENGTH 1024

const FitFunction  * FF[MAX_FITFN]=
      { new LorentzFn(),     new LorentzA12Fn(),  new CurieWeissFn(),
        new NonFermiResFn(), new StraightLineFn(), new PolynomFn()
      };

const char *szRCSID={"$Id"};

 extern char *optarg;
 extern int optind;

 char szTFile[MAXPATH+1]="";
 char szAFile[MAXPATH+1]="ausw.fit";
 char szEFile[MAXPATH+1]="errors";
 char szLFile[MAXPATH+1]="lastfit";
 void RmTmpFile(void);
 char szAnswer [2048];
 DataFile *DF;     //=new DataFile(szArgV[optind]);

 int ReadPars(LineString &FP);
 void PIterate(int nMaxI, double *plfStp, double *plfL, int pnPars, double *plfP, int pnCurves, int fn, int ncPars, double *lfC);

int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iPrint=0, plot=0;
 int ix=1, iy=2;
 const char *szManPath=getenv(MAN_PATH);
 char szPFile[MAXPATH+1]="";
   
 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "F:A:L:E:x:y:vpVh");
  if(iO=='F'){strncpy(szPFile,optarg,MAXPATH); szPFile[MAXPATH]=0;}
  if(iO=='A'){strncpy(szAFile,optarg,MAXPATH); szAFile[MAXPATH]=0;}
  if(iO=='L'){strncpy(szLFile,optarg,MAXPATH); szLFile[MAXPATH]=0;}
  if(iO=='E'){strncpy(szEFile,optarg,MAXPATH); szEFile[MAXPATH]=0;}  
  if(iO=='x'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               ix=strtol(optarg,(char **)NULL, 10);
             } 
  if(iO=='y'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               iy=strtol(optarg,(char **)NULL, 10);
             }
  if(iO=='v')iPrint=1;
  if(iO=='p')plot=1;  
  if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
    }	
 while (iO!=EOF);

 char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Fit", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }

  LineString *FP=0;

  // if no file parameter -> try filename.conf
  if(!szPFile[0]) {
  char *b;
   if ( (b = strchr(szArgV[optind],'.')) != 0 ) {
    strcpy(b,".conf");
    strncpy(szPFile,szArgV[optind],MAXPATH); szPFile[MAXPATH]=0;
   }
   else {
    String B;
    B.Set(szArgV[optind]);
    B.Addf( ".conf");   
    strncpy(szPFile,B,MAXPATH); szPFile[MAXPATH]=0;	   
   }	   
  }
  FP=new LineString();
  if(FP==0 || !FP->Read(szPFile) || FP->GetNLines()==0)
  {FP=0;
   PRINT_DEBUG("Error Pfit: Cannot open parameter file %s\a\n",szPFile)
   exit(EXIT_FAILURE);
  }

  int iType;
      iType=CheckFileType(szFile,LineT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Fit: ERROR-exit\n");
     exit(EXIT_FAILURE);
    }

       switch ((FileGroup)(int)(iType/10)*10)
        {case TheGROUP:
              DF=new TheFile(szFile);
              break;
         case SxSGROUP:
              DF=new SxSFile(szFile);
              break;
        case XDifGROUP:
             DF=new XDifFile(szFile);
             break;
        case SplineGROUP:
             DF=new SplineFile(szFile);
             break;
        case NoGROUP:
             DF=new AsciiFile(szFile);
             break;
       default:fprintf(stderr,"ERROR Fit: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Fit: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    
  if( DF->ReadData()>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Fit: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
  if(!DF->GetCRData()->SetColX(ix-1))
    {fprintf(stderr,"ERROR Fit: Illegal x column %d\n",ix);
     exit(EXIT_FAILURE);
    }
   if(!DF->GetCRData()->SetColY(iy-1))
    {fprintf(stderr,"ERROR Fit : Illegal y column %d\n",iy);
     exit(EXIT_FAILURE);
    }

 ReadPars(*FP);
 if (plot==1) {
  String S;
  S.Addf("echo -e  plot \\'lastfit\\' using 1:2:\\(sqrt\\($ 2\\)\\) with yerrorbars pointsize 0,\\'%s\\' with lines | /usr/bin/gnuplot -persist",szFile);
  system((const char *)S);
 }
}
// *****************************************************************
void RmTmpFile(void)
{
 if(szTFile[0])
   {int i=remove(szTFile);
    if(i)perror("ERROR Fit: Removing tmp-file");
   }
 } 		
// ******************************************************************
int ReadPars(LineString &FP) {
 char pFitFn[MAXPATH+1];
 int pnCurves, pmaxIt, pnPars, pncPars=0, fn; // fn: funtion index in FIT_FUNCS
 double *plfP, *plfL, *plfStp, *plfC;

 pFitFn[0]=0;
 pnCurves=0;
 pmaxIt=0;
 int i,j,icnt=0;
 if(FP.GetPar("StdPars","pFitFn","%S",pFitFn)) { icnt++;}
 if(FP.GetPar("StdPars","pnCurves","%d",&pnCurves)) { icnt++;}
 if(FP.GetPar("StdPars","pmaxIt","%d",&pmaxIt)) { icnt++;}
 if (icnt<3) {
  fprintf(stderr,"ERROR Fit : Parameters missing \n");
  exit(EXIT_FAILURE);
 }
 for (i=0;i<MAX_FITFN;i++)  {
  if (strcmp(FF[i]->GetFnName(),pFitFn)==0) {fn=i; break; }
 }
 if (i==MAX_FITFN) {
  fprintf(stderr,"ERROR Fit : %s : Unknown Fit Function\n", pFitFn);
  exit(EXIT_FAILURE);
 }
 const int iSP=FF[fn]->GetNSumP();
 const int iFP=FF[fn]->GetNFixP();
 const int iCP=FF[fn]->GetNConstP();
 //*********************
 pnPars=pnCurves*iSP+iFP;
 //*********************

 plfP = new double [pnPars];
 plfL = new double [pnPars+1];
 plfStp = new double [pnPars];
 if(iCP>0)plfC = new double [iCP];

 // FUNCTION PARAMETERS
 int k=0;
 icnt=0;
 for (j=0;j<pnCurves;j++) {
  for (i=0;i<iSP;i++) {
      String B;
      B.Set(FF[fn]->GetParName(i));
      B.Addf( "[%d]",j+1);    
      if(FP.GetPar("SumPars",(const char *)B,"%lf",&(plfP[k++])))
        { icnt++; }
      else { printf("%s missing\n",(const char *)B); }
  }
 }
  for (j=i;j<iSP+iFP;j++) {
  if(FP.GetPar("FixPars",FF[fn]->GetParName(j),"%lf",&(plfP[k++])))  
              { icnt++; }
       else { printf("%s missing\n",(const char *)FF[fn]->GetParName(j)); }
  }
 if (icnt<pnPars) {
  fprintf(stderr,"ERROR Fit : Function Parameters missing (%d/%d)\n",icnt,pnPars);
  exit(EXIT_FAILURE);
 }
  // STEPS
 for (j=0,k=0;j<pnCurves;j++) {
  for (i=0;i<iSP;i++,k++) {
      String B;
      B.Set(FF[fn]->GetParName(i));
      B.Addf( "Stp[%d]",j+1);    
      if(FP.GetPar("SumParsStp",(const char *)B,"%lf",&plfStp[k])) { icnt++; }
      else { plfStp[k]=plfP[k]/10; } // default = 10%
  }
 }
  for (j=i;j<iFP+iSP;j++,k++) {
      String B;
      B.Set(FF[fn]->GetParName(j));
      B.Add("Stp");
      if(FP.GetPar("FixParsStp",(const char *)B,"%lf",&plfStp[k])) { icnt++;}
      else { plfStp[j]=plfP[k]/10; }
  }
  // LIMITS
 for (j=0,k=0;j<pnCurves;j++) {
  for (i=0;i<iSP;i++,k++) {
      String B;
      B.Set(FF[fn]->GetParName(i));
      B.Addf("L[%d]",j+1);    
   if(FP.GetPar("SumParsStp",(const char *)B,"%lf",&plfL[k])) { icnt++;}
   else { plfL[k]=1e-06; } // default = 1e-06
  }
 }
  for (j=i;j<iFP+iSP;j++,k++) {
      String B;
      B.Set(FF[fn]->GetParName(j));
      B.Add("L");    
   if(FP.GetPar("FixParsL",(const char *)B,"%lf",&plfL[k])) { icnt++;}
   else { plfL[j]=1e-06; }
  }
 if (iCP!=0) {
  for (i=0;i<iCP;i++) {
   if(FP.GetPar("ConstPars",FF[fn]->GetParName(j+i),"%lf",&plfC[i])) { pncPars++; }
  }
  if(strcmp(pFitFn,"Lorentz Alfa 12")==0 && pncPars==0 )
   {pncPars=1;plfC[0]=1.79278/1.78892;}
  else {
   if (pncPars<iCP) {
    fprintf(stderr,"ERROR Fit : %s : Constant Parameters missing\n", pFitFn);
    exit(EXIT_FAILURE);
   }
  }	
 }

 PIterate(pmaxIt, plfStp, plfL, pnPars, plfP, pnCurves, fn, pncPars, plfC);
 return 1;
}
// ******************************************************************
void PIterate(int nMaxI, double *plfStp, double *plfL, int pnPars, double *plfP, int pnCurves, int fn, int ncPars, double *lfC) {
int i,it=0;
double lfH[pnPars];
double lfV=0, lfM=0, lfE;

Simplex *SI;
SI = new Simplex(DF->GetCRData(), nMaxI, plfStp, plfL, pnPars, plfP, pnCurves, FF[fn], ncPars, lfC);
CHECK_POINTER_RET(SI)
if(SI->GetStatus()!=STATUS_OK) {PRINT_DEBUG("Error creating Simplex\n") return; }

//********************************
while (it==0) {it=SI->Iterate(); }
//********************************
FILE *ef;
ef=fopen(szEFile,"w");
fprintf(ef, "Fit function %s: End at iteration %d of %d\n",SI->GetFnName(),SI->GetNIter(), nMaxI);
for(i=0; i<pnPars; i++)
{lfH[i]=SI->GetFitPar(i);
 fprintf(ef,"%15.15s: %-14.7g Error: %-12.5g Limit: %-12.5g\n",
 SI->GetPName(i),lfH[i],
 SI->GetErr(i),SI->GetLimit(i));
}

double lfSoS=SI->GetFitPar(i);
fprintf(ef,"\nSquaresum of deviations: %-15.5g   ",lfSoS);

if(DF->SaveData(DF->GetInfoText(),szAFile)==0)
  {PRINT_DEBUG("Can not create ausw.fit !\n")
   return;
  }
AsciiFile *OF;
OF = new AsciiFile(szAFile);
 CHECK_POINTER_RET(OF)
 if(OF->ReadData()>0)
   {PRINT_DEBUG("Can not read ausw.fit !\n")
    return;
   }
AsciiFile *EF;
EF = new AsciiFile(szAFile);
 CHECK_POINTER_RET(EF)
 if(EF->ReadData()>0)
   {PRINT_DEBUG("Can not read ausw.fit !\n")
    return;
   }
CRData *ED, *TD;
ED=EF->GetCRData();
TD=OF->GetCRData();
CHECK_POINTER_RET(TD)
for(i=0;i<TD->GetSteps(); i++)
{TD->SetPointY(FF[fn]->GetFFStruct().FitFunc(TD->GetPointX(i),pnCurves,lfH, ncPars, lfC),i);
 lfM+=TD->GetPointY(i);
}
lfM/=TD->GetSteps();

for(i=0;i<ED->GetSteps(); i++)
{lfE=(FF[fn]->GetFFStruct().FitFunc(TD->GetPointX(i),pnCurves,lfH, ncPars, lfC) -
 DF->GetCRData()->GetPointY(i))/TD->GetPointY(i);
 ED->SetPointY(lfE,i);
 lfV+=SQUARE(TD->GetPointY(i)-lfM);
}

fprintf(ef, "r-value: %-15.5g\n", sqrt(fabs(1-lfSoS/lfV)) );
fprintf(ef, "Mean deviation: %-12.5g Error: %-12.5g Limit: %-12.5g\n",
       sqrt(SI->GetFitPar(pnPars))/DF->GetCRData()->GetSteps(),
       SI->GetErr(pnPars),SI->GetLimit(pnPars));
fclose(ef);
EF->SaveData(szAFile);
OF->SaveData(szLFile);
}
