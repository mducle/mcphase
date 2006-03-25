// File: merge.cpp
// $Log: merge.cpp,v $
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

const char *szRCSID={"$Id: merge.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet[2]={1,1},iPrint=0;
 int iC=0, iR=0;
 const char *szManPath=getenv(MAN_PATH);

 if(iArgC<1)exit(EXIT_FAILURE);

 if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "s:t:crvVh");
   if(iO=='s')
     {if(*optarg==',')iSet[1]=strtol(optarg+1,(char **)NULL, 10);
      else 
       {if(strchr(optarg,','))
	   {if(GetIntRange(optarg,iSet[0],iSet[1],','))Usage(szManPath,szArgV[0],stderr);
	   }
	  
        else iSet[0]=strtol(optarg,(char **)NULL, 10);
       }
      }
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
   if(iO=='c')iC=1;
   if(iO=='r')iR=1;
   if(iO=='v')iPrint=1;
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

  if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);

  if(iC && iR)
    {fprintf(stderr,"Can not merge columns and rows. Use either -c or -r\n");
     exit(EXIT_FAILURE);
    } 

  if( !iC && !iR)
    {fprintf(stderr,"Specify either -c or -r\n");
     exit(EXIT_FAILURE);
    } 
   
  char szFile[2][MAXPATH+1];

  if(!szArgV[optind+1]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Merge", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     
     strncpy(szFile[0],szTFile,MAXPATH);
     strncpy(szFile[1],szArgV[optind],MAXPATH);
    }
  else
    {strncpy(szFile[0],szArgV[optind],MAXPATH);
     strncpy(szFile[1],szArgV[optind+1],MAXPATH);
     }
    	    
  DataFile *DF[2];     //=new DataFile(szArgV[optind]);

  int iType,nSets[2],i;
 enum EndLine  LT;
  for(i=0;i<2;i++)
     {iType=CheckFileType(szFile[i],LT);

      if( (FileType) iType == FtNONE)
        {fprintf(stderr,"Merge: ERROR-exit\n");
         exit(EXIT_FAILURE);
        }

       switch ((FileGroup)(int)(iType/10)*10)
        {case TheGROUP:
              DF[i]=new TheFile(szFile[i]);
              nSets[i]=DF[i]->GetNSets();
              break;
         case SxSGROUP:
              DF[i]=new SxSFile(szFile[i]);
	      nSets[i]=DF[i]->GetNSets();
              break;
        case XDifGROUP:
             DF[i]=new XDifFile(szFile[i]);
	     nSets[i]=DF[i]->GetNSets();
             break;
        case SplineGROUP:
             DF[i]=new SplineFile(szFile[i]);
	     nSets[i]=DF[i]->GetNSets();
             break;
        case NoGROUP:
             DF[i]=new AsciiFile(szFile[i]);
	     nSets[i]=DF[i]->GetNSets();
             break;
       default:fprintf(stderr,"ERROR Merge: file %s has unsupported type %d",
                           szFile[i],iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF[i]->GetFType()==FtERROR || DF[i]->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Merge: Error reading file %s",szFile[i]);
     exit(EXIT_FAILURE);
    }	    

  if(iSet[i]>nSets[i])   
    {fprintf(stderr,"ERROR Merge: Illegal # of set %d (1...%d)\n",iSet[i],nSets[i]);
     exit(EXIT_FAILURE);
    }		
  if( DF[i]->ReadData(iSet[i])>0 || DF[i]->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Merge: Error reading data of %s\n",szFile[i]);
     exit(EXIT_FAILURE);
    }
  }//for     

  if(iPrint)
   {String S;
    fprintf(stderr,"Merging: %s of %s and %s \n",(iC ? "Columns" : "Rows")
		                                 ,szFile[0],szFile[1]);
    DF[0]->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }
 if(iR)
   {if(! DF[0]->GetCRData()->MergeRows(*(DF[1]->GetCRData())))
      {fprintf(stderr,"ERROR Merge: Merge failed\n");
       exit(EXIT_FAILURE);
      }
   }  
 else
  {fprintf(stderr," Merging columns: not implemented now!!\n");
   exit(EXIT_FAILURE);
  } 	 
 if(LineT!=NOEnd)DF[0]->SetLineType(LineT);
 if(!DF[0]->SaveData(stdout))
   {fprintf(stderr,"ERROR Merge: Can not write to stdout\n");
    exit(EXIT_FAILURE);
   } 

 if(iPrint)
   {String S;
    DF[0]->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

 delete DF[0];
 delete DF[1];
 
 exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Merge: Removing tmp-file");
   }
 
} 		
// *******************************************************************
