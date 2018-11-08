// File: remove.cpp
// $Log: remove.cpp,v $
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

const char *szRCSID={"$Id: remove.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet=1,iPrint=0;
 
 int iCs=0,iCe=0;
 int iRs=0,iRe=0;
 double lfS=0,lfE=0;
 int ifHere=0,irHere=0;
 const char *szManPath=getenv(MAN_PATH);

 if(iArgC<1)exit(EXIT_FAILURE);

 if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "s:t:c:r:f:vVh");
   if(iO=='s'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               else iSet=strtol(optarg,(char **)NULL, 10);
              }
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
   if(iO=='c'){if(strchr(optarg,':'))
                 {if(GetIntRange(optarg,iCs,iCe))Usage(szManPath,szArgV[0],stderr);
                 }
                else
		{if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
                 iCs=strtol(optarg,(char **)NULL,10);
		 iCe=iCs;
	        }
	       }	              
   
   if(iO=='r'){if(strchr(optarg,':'))
                 {if(GetIntRange(optarg,iRs,iRe))Usage(szManPath,szArgV[0],stderr);
                 }
                else
		{if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
                 iRs=strtol(optarg,(char **)NULL,10);
		 iRe=iRs;
	        }
               irHere=1;
	       }	              
   if(iO=='f'){if(strchr(optarg,':'))
                 {if(GetFloatRange(optarg,lfS,lfE))Usage(szManPath,szArgV[0],stderr);
                 }
                else
		{if(!IsNumString(optarg,FLOAT_NUM))Usage(szManPath,szArgV[0],stderr);
                 lfS=strtod(optarg,(char **)NULL);
		 lfE=lfS;
	        }
                ifHere=1;
	       }	              
       
   if(iO=='v')iPrint=1;
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

  //if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);
  if( ifHere)// float range specified
    { if(irHere)
        {fprintf(stderr,"ERROR Remove: Cannot use -f and -r together\n");
         exit(EXIT_FAILURE);
        }
      if(iCs==0 && iCe==0)
        {fprintf(stderr,"ERROR Remove: No column (-c) specified\n");
         exit(EXIT_FAILURE);
        }
      if(iCs!=iCe) 
       {fprintf(stderr,"ERROR Remove: -c Illegal column range\n");
         exit(EXIT_FAILURE);
        }

      if(lfS>lfE)
        {fprintf(stderr,"ERROR Remove: -f Illegal range\n");
         exit(EXIT_FAILURE);
        }
    }  

  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Remove", TMPFILE_EXIT);
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
    {fprintf(stderr,"Remove: ERROR-exit\n");
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
       default:fprintf(stderr,"ERROR Remove: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Remove: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Remove: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Remove: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
    
  fprintf(stderr,"Remove Col: %d:%d Row: %d:%d\n",iCs,iCe,iRs,iRe);

  int ir;
   if(ifHere)// float range specified
    {ir=DF->GetCRData()->DelVal(iCs-1,lfS,lfE);}	    

  else 
    {ir=DF->GetCRData()->DelRange(iCs-1,iCe-1,iRs-1,iRe-1);}

  if(ir<=0)
    {fprintf(stderr,"ERROR Remove: %s\n",szFile);
     exit(EXIT_FAILURE);
    }
		  
     	  	  
  if(iPrint)
   {String S;
    fprintf(stderr,"Removing file:%s Column: %d:%d Row: %d:%d\n",szFile,
		                                            iCs,iCe,iRs,iRe);
    DF->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

//  DF->GetCRData()->RemoveData(iC-1);
 
   
 if(LineT!=NOEnd)DF->SetLineType(LineT);
 if(!DF->SaveData(stdout))
   {fprintf(stderr,"ERROR Remove: Can not write to stdout\n");
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
    if(i)perror("ERROR Remove: Removing tmp-file");
   }
 
} 		
// *******************************************************************
