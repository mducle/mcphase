// File: tfile.cpp
// $Log: tfile1.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#include <stdio.h>
#include <unistd.h>

#ifndef NDFILE_H
#include "ndfile.h"
#endif

#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: tfile1.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet=1,iF=1, iL=0;
 int iLHere=0, iFHere=0;
 const char *szManPath=getenv(MAN_PATH);
  
 if(iArgC<1)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "s:f:n:hVt");
   if(iO=='s')iSet=strtol(optarg,(char **)NULL, 10);
   if(iO=='f'){iF=strtol(optarg,(char **)NULL, 10); iFHere=1;}
   if(iO=='n'){iL=strtol(optarg,(char **)NULL, 10); iLHere=1;}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
  }	   
 while (iO!=EOF);

  //if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);

  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Head", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);}

  DataFile *DF; 
  int iType,nSets;
 enum EndLine  LT;
  iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"ERROR tfile: file %s has unsupported type %d\n",szFile,iType);
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
    {fprintf(stderr,"ERROR tfile: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR tfile: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR tfile:Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
     
 String B;
 if(iFHere || (!iFHere && ! iLHere) )
   {DF->SPrintInfo(B,iF);
    fprintf(stdout,"%s",(const char *)B);    
   }

 if(iLHere || (!iFHere && !iLHere))
   {int i;
    fprintf(stdout,"File has %d text line(s)\n",DF->GetNoTxtL());
    if(iL==-1)
      {if(DF->GetNoTxtL()>1)
         {for(i=0; i<DF->GetNoTxtL(); i++)
             fprintf(stdout," Text %d: %s\n",i+1,DF->GetInfoText(i));
         }
       
       else fprintf(stdout," Text: %s\n",DF->GetInfoText());
      }
    else
      {if(DF->GetNoTxtL()>1 && iL>0)
             fprintf(stdout," Text %d: %s\n",iL,DF->GetInfoText(iL-1));

       else fprintf(stdout," Text: %s\n",DF->GetInfoText());
      }
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
    if(i)perror("ERROR Calc: Removing tmp-file");
   }
} 		
// *******************************************************************
