//File: addspec.cpp
// $Log: addspec.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.6  1999/07/09 11:27:45  herbie
// new dfile routines
//
// Revision 1.5  1999/03/15 09:08:37  herbie
// -V added
//
// Revision 1.3  1999/02/15 16:13:44  herbie
// switch -V corrected
//
// * Revision 1.2  1999/02/15 16:08:01  herbie
// * switch -V added
// *

#include <stdio.h>
#include <unistd.h>

#include "dfile.h"

#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: addspec.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

int main(int iArgC, char ** szArgV)
{int iO,iPrint=0;
 enum EndLine LineT=NOEnd;
 
 const char *szManPath=getenv(MAN_PATH);
		 
 if(iArgC<1)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "t:vhV");
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
      }
   if(iO=='v')iPrint=1;
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_SUCCESS);}
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);


  //fprintf(stderr,"Opt: %c: %x: OptInd: %c %x OptErr: %c %x OptOpt: %c,%x\n",
  //             iO,iO, optind, optind, opterr, opterr, optopt, optopt);

  // fprintf(stderr,"OptArg:%s\n", optarg);
  //fprintf(stderr,"LineEnd: %d InFile:%s  OutFile:%s\n",
  //          (int)LineT,szArgV[optind],szArgV[optind+1]);

  //fprintf(stderr,"ArgC: %d    optind: %d\n", iArgC,optind);

  //if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);

  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Addspec", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);}
 
  DataFile *DF;
  DF=new SxSFile(szFile);
  if(DF->GetFType()!=FtSXS_NEW)
    {fprintf(stderr,"File %s is not a SXS-data-file (NEW)\n",szFile);
     exit(EXIT_FAILURE); 
    }	    

  int nSets=DF->GetNSets();
  if(nSets==1)
    {fprintf(stderr,"Nothing to do. File %s has only one data set\n",szFile);
     exit(EXIT_FAILURE); 
    }	    

 if(DF==0 || DF->ReadData(1)>0 || DF->GetCRData()->GetSteps()==0 )
   {fprintf(stderr,"Error reading InputFile %s\n",szFile);
    exit(EXIT_FAILURE);
   } 

 int ir=DF->AddSets();
 if(ir < 2)
   {fprintf(stderr,"Error adding sets from %s\n",szFile);
    exit(EXIT_FAILURE);
   } 

 if(LineT!=NOEnd)DF->SetLineType(LineT);
 if(!DF->SaveData(stdout))
   {fprintf(stderr,"Error writing output to stdout\n");
    exit(EXIT_FAILURE);
   } 

 if(iPrint)
   {String S;
    DF->SPrintInfo(S, PR_TYPE | PR_COLS | PR_TEXT);
    fprintf(stderr,"%s",(const char *)S);    
   }
   
 fprintf(stderr,"\n%d of %d spectra in %s added\n",
                 ir,nSets,szFile);
 delete DF;
 exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Addspec: Removing tmp-file");
   }
} 		
// *******************************************************************
