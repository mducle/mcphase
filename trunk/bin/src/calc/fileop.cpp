// File: fileop.cpp
// $Log: fileop.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.4  1999/04/29 09:41:45  herbie
// skip-bug removed
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

const char *szRCSID={"$Id: fileop.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet[2]={1,1},iPrint=0;
 int ix1=1, ix2=1, iy1=2, iy2=2;
 char cOp=0;
 const char *szManPath="Usage: Fileop -f op [-x #[,#]] [-y #[,#]] [-s #,#] [-t] [-v] [-h]  (File1) File2\n"
"        -s #,#: (int) number of data set to use for the calculation\n"
"	        Two numbers separated by commas can be specified, refering \n"
"		to File1 or File2, respectively   \n"
"                (Only valid if file type supports multiple data sets)  \n"
"         -f op: A single character defining operation\n"
"        -x #,#:             \n"
"        -y #,#: (int) number of x-, y-columns used for the calculation.\n"
"	        Two numbers separated by commas can be specified, refering \n"
"		to File1 or File2, respectively   \n"
"	        If ommited the default values x:1, y:2 are assumed.\n"
"       -t dos : OutFile in DOS <lf><cr> format\n"
"       -t unix: OutFile in UNIX <lf> format\n"
"            -t: ommited: OutFile same as IputFile format\n"
"            -v: verify -> print header before and after operation (stderr)          \n"
"            -h: Print this help message \n"
"     InputFile: Input data file\n"
"RESULT:\n"
"The operation defined by the character is performed:\n"
"y1(x1) = y1(x1) <op> y2(x1). \n"
"x(File1) and x(File2) need not be same, linear interpolation is used to get\n"
"the value y2(x1) from y2(x2)\n"
"Columns not involved in the operation are unchanged.\n"
"Values MUST be sorted by x-column\n"
"File1 can be piped (>).\n"
"Output is written to stdout.\n"
"$Id: Fileop.man,v 1.4 1999/04/29 09:44:09 herbie Exp herbie $\n";

 if(iArgC<1)exit(EXIT_FAILURE);

 if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "f:s:t:x:y:vVh");
   if(iO=='s')
     {if(*optarg==',')iSet[1]=strtol(optarg+1,(char **)NULL, 10);
      else 
       {if(strchr(optarg,','))
	   {if(GetIntRange(optarg,iSet[0],iSet[1],','))Usage(szManPath,szArgV[0],stderr);
	   }
        else iSet[0]=strtol(optarg,(char **)NULL, 10);
       }  
     }
   if(iO=='f')cOp=optarg[0];
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
   if(iO=='x')
     {if(*optarg==',')ix2=strtol(optarg+1,(char **)NULL, 10);
      else 
       {if(strchr(optarg,','))
	   {if(GetIntRange(optarg,ix1,ix2,','))Usage(szManPath,szArgV[0],stderr);
	   }
        else ix1=strtol(optarg,(char **)NULL, 10);
	}  
     }
   if(iO=='y')
     {if(*optarg==',')iy2=strtol(optarg+1,(char **)NULL, 10);
      else 
       {if(strchr(optarg,','))
	   {if(GetIntRange(optarg,iy1,iy2,','))Usage(szManPath,szArgV[0],stderr);
	   }
        else iy1=strtol(optarg,(char **)NULL, 10);
	}  
     }
   if(iO=='v')iPrint=1;
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

  if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);
 
  char szFile[2][MAXPATH+1];

  if(!szArgV[optind+1]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Fileop", TMPFILE_EXIT);
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
        {fprintf(stderr,"Fileop: ERROR-exit\n");
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
       default:fprintf(stderr,"ERROR Fileop: file %s has unsupported type %d",
                           szFile[i],iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF[i]->GetFType()==FtERROR || DF[i]->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Fileop: Error reading file %s",szFile[i]);
     exit(EXIT_FAILURE);
    }	    

  if(iSet[i]>nSets[i])   
    {fprintf(stderr,"ERROR Fileop: Illegal # of set %d (1...%d)\n",iSet[i],nSets[i]);
     exit(EXIT_FAILURE);
    }		
  if( DF[i]->ReadData(iSet[i])>0 || DF[i]->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Fileop: Error reading data of %s\n",szFile[i]);
     exit(EXIT_FAILURE);
    }
  }//for     

  if(iPrint)
   {String S;
    fprintf(stderr,"Calcultaing: File:%s %c File: %s\n",szFile[0],cOp,szFile[1]);
    DF[0]->SPrintInfo(S,PR_TEXT | PR_COLS | PR_TYPE);
    fprintf(stderr,"%s",(const char *)S);    
   }

   DF[0]->GetCRData()->SetColX(ix1-1);
   DF[0]->GetCRData()->SetColY(iy1-1);

   int ir=DF[0]->GetCRData()->MathOper(DF[1]->GetCRData(),cOp,ix2-1,iy2-1);
 
 if(ir<0)
   {fprintf(stderr,"ERROR Fileop: Illegal range or operation\n");
    exit(EXIT_FAILURE);
   }

 if(ir>0)
   fprintf(stderr,"WARNING Fileop: %d data points skipped\n",ir);
    
 if(LineT!=NOEnd)DF[0]->SetLineType(LineT);
 if(!DF[0]->SaveData(stdout))
   {fprintf(stderr,"ERROR Fileop: Can not write to stdout\n");
    exit(EXIT_FAILURE);
   } 

 if(iPrint)
   {String S;
    DF[0]->SPrintInfo(S,PR_TEXT | PR_COLS | PR_TYPE);
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
    if(i)perror("ERROR Fileop: Removing tmp-file");
   }
 
} 		
// *******************************************************************
