// File: zero.cpp
// $Log: zero.cpp,v $
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

#define MAN_PATH "AUSW_MANPATH"

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

const char* szRCSID={"$Id: zero.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet=1,iPrint=0;
 int ix=1,  iy=2;
 double x0=0;
 const char *szManPath=getenv(MAN_PATH);

 if(iArgC<1)exit(EXIT_FAILURE);

 //if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "f:s:t:x:y:vVh");
   if(iO=='s'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               else iSet=strtol(optarg,(char **)NULL, 10);
              }
   if(iO=='f'){if(!IsNumString(optarg,FLOAT_NUM))Usage(szManPath,szArgV[0],stderr);
 	       x0=strtod(optarg,(char **)NULL);
              } 
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
   if(iO=='x'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               ix=strtol(optarg,(char **)NULL, 10);
              } 
              
   if(iO=='y'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               iy=strtol(optarg,(char **)NULL, 10);
              }

   if(iO=='v')iPrint=1;
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

  //if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);
 
  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Zero", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }
    	    
  DataFile *DF;     //=new DataFile(szArgV[optind]);

  int iType,nSets;
 enum EndLine  LT;
  iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Zero: ERROR-exit\n");
     exit(EXIT_FAILURE);
    }

       switch ((FileGroup)(int)(iType/10)*10)
        {case TheGROUP:
              DF=new TheFile(szFile);
              nSets=DF->GetNSets();
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
       default:fprintf(stderr,"ERROR Zero: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Zero: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Zero: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Zero: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
    

  if(iPrint)
   {String S;
    fprintf(stderr,"Shifting zero of File:%s x-Column:%d x-Column:%d\n",szFile,ix,iy);
    DF->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

 
 if(!DF->GetCRData()->ZeroShift(x0,ix-1,iy-1))
   {fprintf(stderr,"ERROR Zero: Shift failed\n");
    exit(EXIT_FAILURE);
   }
    
 if(LineT!=NOEnd)DF->SetLineType(LineT);
 if(!DF->SaveData(stdout))
   {fprintf(stderr,"ERROR Zero: Can not write to stdout\n");
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
    if(i)perror("ERROR Zero: Removing tmp-file");
   }
 
} 		
// *******************************************************************
