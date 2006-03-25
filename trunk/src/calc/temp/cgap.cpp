// File: cgap.cpp
// $Log:$

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "dfile.h"
#include "formula.h"
#include "stdfunc.h"
#include "thefunc.h"

#define PARAMETER_FILE "ausw.conf"
#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id:$"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet=1,iPrint=0;
 char szPFile[MAXPATH+1]="";
 const char *szManPath=getenv(MAN_PATH);
 const char *szConfPath=getenv("AUSW_CONF_PATH");

 if(iArgC<1)exit(EXIT_FAILURE);

 //if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "F:t:vVh");
   if(iO=='F'){strncpy(szPFile,optarg,MAXPATH); szPFile[MAXPATH]=0;} 
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
    
   if(iO=='v')iPrint=1;
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

  //if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);
 
  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "CGap", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }
    	    
  DataFile *DF;     //=new DataFile(szArgV[optind]);

  LineString *FP=0,*TH=0;

  if(!szPFile[0])AppendPath(szConfPath,PARAMETER_FILE,szPFile,MAXPATH);
  szPFile[MAXPATH]=0;
  
  int iType,nSets;
 enum EndLine  LT;
  iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"CGap: ERROR-exit\n");
     exit(EXIT_FAILURE);
    }

       switch ((FileGroup)(int)(iType/10)*10)
        {case TheGROUP:
              DF=new TheFile(szFile);
              nSets=DF->GetNSets();
              FP=new LineString();
              if(FP==0 || !FP->Read(szPFile) || FP->GetNLines()==0)
                {FP=0;
                 PRINT_DEBUG("Error GAP: Cannot open parameter file %s\a\n",szPFile)
                 exit(EXIT_FAILURE);
                }
              TH=new LineString();
              if(TH==0 || !TH->Read(szFile," ",LineString::S_NUMBER) ||
              TH->GetNLines()==0)
	       {TH=0;
                PRINT_DEBUG("Error GAP: Cannot read file header %s\a\n",szFile)
                exit(EXIT_FAILURE);
               }
              break;
       default:fprintf(stderr,"ERROR CGap: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR CGap: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR CGap: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR CGap: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
    
  CAP_CELL CC;
  if(!GetCellPars(&CC,*FP,*TH))
    {fprintf(stderr,"ERROR CGap: Invalid parameter (header) in %s\n",szFile);
     exit(EXIT_FAILURE);
    }
  if(! (   (CC.uFlag & CC_CELL_DIAM) && (CC.uFlag & CC_CAL1_FILE) &&
           (CC.uFlag & CC_CELL_GAP)  && (CC.uFlag & CC_CELL_PIVO)) )
   {fprintf(stderr,"ERROR CGap: Invalid parameter in %s\n",szFile);
     exit(EXIT_FAILURE);
   }

   
   if(iPrint)
   fprintf(stderr,"Ag file:%s\n",CC.szCal);

   TheFile *AgF=new TheFile(CC.szCal);

  if( AgF->ReadData(1)>0 || AgF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Cgap: Error reading data of %s\n",CC.szCal);
     exit(EXIT_FAILURE);
    }

   CRData *Ag=AgF->GetCRData();
   if(AgF->GetError() || !Ag || Ag->GetSteps()==0)
     {fprintf(stderr,"ERROR Cgap: %s not found\n",CC.szCal);
      exit(EXIT_FAILURE);
     }

  Ag->SetColX(0);
  Ag->SetColY(1);
   
   double ri=0.5*CC.fCGap;
   double ro=0.5*CC.fCDiam;

   TiltPlate *TP=new TiltPlate(ro,ri,CC.fCPivot,0.8, CC.K0,Ag,0);
   if(!TP->IsValid())
     {fprintf(stderr,"ERROR Cgap: Invalid parameter\n");
      exit(EXIT_FAILURE);
     }

   if(iPrint)
   {String S;
    fprintf(stderr,"Calculating gap of File:%s Parameter file: %s ",
                    szFile,szPFile);
    DF->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

    TP->CalcGap(DF->GetCRData(),FUNC_COMPENS);

    delete TP;
    delete Ag;
    delete FP;
 
 DF->SetColID("Gap [m]",1);
   
 if(LineT!=NOEnd)DF->SetLineType(LineT);
 if(!DF->SaveData(stdout))
   {fprintf(stderr,"ERROR CGap: Can not write to stdout\n");
    exit(EXIT_FAILURE);
   } 

 if(iPrint)
   {String S;
    DF->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

 delete DF;
 
 exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR CGap: Removing tmp-file");
   }
 
} 		
// *******************************************************************
