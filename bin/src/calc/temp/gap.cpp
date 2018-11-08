// File: gap.cpp
// $Log: gap.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// -V added
//

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "dfile.h"
#include "formula.h"
#include "stdfunc.h"
#include "thefunc.h"
#include "strings.h"

#define PARAMETER_FILE "ausw.conf"
#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: gap.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet=1,iPrint=0;
 int iC=2;
 double lfG=-10,lfD=0;
 int iDHere=0;
 const char *szManPath=getenv(MAN_PATH);
 const char *szConfPath=getenv("AUSW_CONF_PATH");

 if(iArgC<1)exit(EXIT_FAILURE);

 //if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "D:G:s:t:c:vVh");
   if(iO=='s'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               else iSet=strtol(optarg,(char **)NULL, 10);
              }
   if(iO=='D'){if(!IsNumString(optarg,FLOAT_NUM))Usage(szManPath,szArgV[0],stderr);
 	       lfD=strtod(optarg,(char **)NULL);
               iDHere=1;
              } 
   if(iO=='G'){if(!IsNumString(optarg,FLOAT_NUM))Usage(szManPath,szArgV[0],stderr);
 	       lfG=strtod(optarg,(char **)NULL);
              } 
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
   if(iO=='c'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               iC=strtol(optarg,(char **)NULL, 10);
              } 
              
   if(iO=='v')iPrint=1;
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

  //if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);
 
  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Gap", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }
    	    
  DataFile *DF;     //=new DataFile(szArgV[optind]);

  LineString *FP=0,*TH=0;


  int iType,nSets;
 enum EndLine  LT;
  iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Gap: ERROR-exit\n");
     exit(EXIT_FAILURE);
    }

       switch ((FileGroup)(int)(iType/10)*10)
        {case TheGROUP:
              DF=new TheFile(szFile);
              nSets=DF->GetNSets();
              if(!iDHere)
                {char szB[MAXPATH+1];
                 AppendPath(szConfPath,PARAMETER_FILE,szB,MAXPATH);
                 FP=new LineString();
                 if(FP==0 || !FP->Read(szB) || FP->GetNLines()==0)
	           {FP=0;
                    PRINT_DEBUG("Error GAP: Cannot open parameter file %s\a\n",PARAMETER_FILE)
                    exit(EXIT_FAILURE);
                   }
                 TH=new LineString();
                 if(TH==0 || !TH->Read(szFile," ",LineString::S_NUMBER) ||
                    TH->GetNLines()==0)
	           {TH=0;
                    PRINT_DEBUG("Error GAP: Cannot read file header %s\a\n",szFile)
                    exit(EXIT_FAILURE);
                   }
		   
		}   
              break;
         case SxSGROUP:
              DF=new SxSFile(szFile);
	      nSets=DF->GetNSets();
              break;
        case XDifGROUP:
             DF=new XDifFile(szFile);
	     nSets=DF->GetNSets();
             break;
        case SplineGROUP:
             DF=new SplineFile(szFile);
	     nSets=DF->GetNSets();
             break;
        case NoGROUP:
             DF=new AsciiFile(szFile);
	     nSets=DF->GetNSets();
             break;
       default:fprintf(stderr,"ERROR Gap: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Gap: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Gap: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Gap: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
    
  if(FP)
    { CAP_CELL CC;
     if(!GetCellPars(&CC,*FP,*TH))
       {fprintf(stderr,"ERROR Gap: Invalid parameter (header) in %s\n",szFile);
        exit(EXIT_FAILURE);
       }
     if(! (CC.uFlag & CC_CELL_DIAM) )
       {fprintf(stderr,"ERROR Gap: Invalid parameter (Diam) in %s\n",szFile);
        exit(EXIT_FAILURE);
       }
     if(! (CC.uFlag & CC_CELL_GAP) )CC.fCGap=-10.;
     lfG=CC.fCGap;
     lfD=CC.fCDiam;
    }
   if(iPrint)
   {String S;
    fprintf(stderr,"Calculating gap of File:%s\nColumn:%d D=%f [mm] G=%f [mm]\n",
                    szFile,iC,lfD,lfG);
    DF->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

 if(!(  (*DF->GetCRData())[iC-1].CalculateGap(lfD*0.1,lfG*0.1) )  ) 
   {fprintf(stderr,"ERROR Gap: Calculation failed\n");
    exit(EXIT_FAILURE);
   }
 
  DF->SetColID("Gap [mm]",1);
    
 if(LineT!=NOEnd)DF->SetLineType(LineT);
 if(!DF->SaveData(stdout))
   {fprintf(stderr,"ERROR Gap: Can not write to stdout\n");
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
{fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Gap: Removing tmp-file");
   }
 
} 		
// *******************************************************************
