// File: tfile1.cpp
// $Log: tfile.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#ifndef DFILE_H
#include "dfile.h"
#endif

#ifndef FORMULA_H
#include "formula.h"
#endif

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: tfile.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet=1,iPrint=0;
 
 int iCs=-1,iCe=-1;
 int iRs=-1,iRe=-1;
 double lfS=0,lfE=0;
 int ifHere=0, irHere=0;
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
   if(iO=='r'){if(GetIntRange(optarg,iRs,iRe))Usage(szManPath,szArgV[0],stderr);
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
        {fprintf(stderr,"ERROR Select: Cannot use -f and -r together\n");
         exit(EXIT_FAILURE);
        }
      if(iCs==0 && iCe==0)
        {fprintf(stderr,"ERROR Select: No column (-c) specified\n");
         exit(EXIT_FAILURE);
        }
      if(iCs!=iCe)
        {fprintf(stderr,"ERROR Select: -c Illegal column range\n");
         exit(EXIT_FAILURE);
        }

      if(lfS>lfE)
        {fprintf(stderr,"ERROR Select: -f Illegal range\n");
         exit(EXIT_FAILURE);
        }
    }  
 
  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Select", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }
    	    
  DataFile *DF;     //=new DataFile(szArgV[optind]);

  int iType,nSets;
 enum EndLine LT;
   iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Select: ERROR-exit\n");
     exit(EXIT_FAILURE);
    }

       switch ((FileGroup) (iType/10)*10)
        {case TheGROUP:
              DF=new TheFile(szFile);
              nSets=DF->GetNSets();
              break;
         case SxSGROUP:
              DF=new SxSFile(szFile);
	      nSets=DF->GetNSets();
              break;
//         case SplineGROUP:
//              DF=new SplineFile(szFile);
// 	     nSets=DF->GetNSets();
//              break;
        case NoGROUP:
             DF=new AsciiFile(szFile);
	     nSets=DF->GetNSets();
             break;
       default:fprintf(stderr,"ERROR Select: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Select: Error reading file %s\n",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Select: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Select: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
    

   if(ifHere)// float range specified
    {int ics=(iCs==-1 ? 1 : iCs);
     if(DF->GetCRData()->SelVal(ics-1,lfS,lfE)<=0)
              fprintf(stderr,"ERROR Select: %s\n",szFile);
    }	    
   else    
    {int ics=(iCs==-1 ? 0 : iCs);
     int ice=(iCe==-1 ? DF->GetCRData()->GetCols() : iCe);
     fprintf(stderr,"Col: %d:%d  %d:%d\n",ics,ice,iCs,iCe);
     if(!CHECK_INDEX(ics, DF->GetCRData()->GetCols()) ||
        !CHECK_INDEX(ics, DF->GetCRData()->GetCols()))exit(EXIT_FAILURE);
		  
     int irs=(iRs==-1 ? 0 : iRs);
     int ire=(iRe==-1 ? DF->GetCRData()->GetSteps(): iRe);

     fprintf(stderr,"Row: %d:%d  %d:%d\n",irs,ire,iRs,iRe);
     if(!CHECK_INDEX(ics, DF->GetCRData()->GetSteps()) ||
        !CHECK_INDEX(ics, DF->GetCRData()->GetSteps()) )exit(EXIT_FAILURE);

  
     if(iCs!=-1 && iCe!=-1)
       {if(!DF->GetCRData()->NewColRange(ics-1,ice-1) || !DF->SetColRange(ics-1,ice-1))
          {fprintf(stderr,"ERROR Select: %s Illegal column range (%d:%d)\n",szFile,ics,ice);
           exit(EXIT_FAILURE);
          }
        }			  

     if(iRs!=-1 && iRe!=-1)
       {if(!DF->GetCRData()->NewRowRange(irs-1,ire-1))
          {fprintf(stderr,"ERROR Select: %s Illegal row range (%d:%d)\n",szFile,irs,ire);
           exit(EXIT_FAILURE);
          }
        }			  
     } 	  	  

 if(iPrint)
   {String B;
    fprintf(stderr,"Selecting file:%s Column: %d:%d Row: %d:%d\n",szFile,
                                            iCs,iCe,iRs,iRe);
    DF->SPrintInfo(B,PR_ALL);
    fprintf(stderr,"%s",(const char *)B);    
   }
   
 if(LineT!=NOEnd)DF->SetLineType(LineT);
 if(!DF->SaveData(stdout))
   {fprintf(stderr,"ERROR Select: Can not write to stdout\n");
    exit(EXIT_FAILURE);
   } 

 if(iPrint)
   {String B;
    DF->SPrintInfo(B,PR_ALL);
    fprintf(stderr,"%s",(const char *)B);    
   }

 delete DF;
 
 exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Select: Removing tmp-file");
   }
 
} 		
// *******************************************************************
