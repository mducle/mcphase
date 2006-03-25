// File: head.cpp
// $Log: head.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.4  1999/07/12 09:35:07  herbie
// New dfile routines
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// -V added
//

#include <stdio.h>
#include <unistd.h>

#include "dfile.h"

#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: head.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

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
    {fprintf(stderr,"Head: ERROR-exit\n");
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
     default:fprintf(stderr,"ERROR Head: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
     }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Head: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Head: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Head:Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
     
 String S;
 if(iFHere || (!iFHere && ! iLHere) )
   {DF->SPrintInfo(S,iF);
    fprintf(stdout,"%s",(const char *)S);    
   }

 if(iLHere || (!iFHere && ! iLHere))
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
      {if(DF->GetNoTxtL()>1)
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
