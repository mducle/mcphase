// File: calc.cpp
// $Log: calc.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.4  1999/07/12 10:50:46  herbie
// new dfile outines
// -
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// -V added
//

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>

#include "dfile.h"
#include "formula.h"
#include "stdfunc.h"
#include "strings.h"

#define MAN_PATH "$MCPHASE_DIR/src/calc"

const char *szRCSID={"$Id: calc.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iSet=1,iPrint=0;
 
 char szF[MAX_LINELENGTH+1];
 const char *szManPath="Usage: Calc [-s #] [-f string] [-t type] [-v] [-h] InputFile\n"
                       "        -s #: (int) number of data set to perform calculation\n"
                       "                (Only valid if file type supports multiple data sets)\n"
                       "   -f string: Formula defining calculation\n"
                       "               Allowed formulas  (#: Column number : 1 - n)\n"
                       "                c# = 1.23 <+-*/^> c#\n"
                       "                c# = c# <+-*/^> 1.234 \n"
                       "                c# = c# <+-*/^> c# \n"
                       "                c# = FUN(c#) \n"
                       "                c# = FUN(1.234) \n" 
                       "                FUN: ABS, COS, EXP, LOG, SIN,  TAN, SEQ,\n"
                       " 	             ACOS, ASIN, ATAN, SQRT,\n"
                       "                Not case sensitive; blancs are ignored\n"
                       "    -t dos : OutFile in DOS <lf><cr> format\n"
                       "    -t unix: OutFile in UNIX <lf> format\n"
                       "         -t: ommited: OutFile same as InputFile format\n"
                       "         -v: verify: print header before and after operation (stderr)\n"
                       "         -h: Print this help message \n"
                       "  InputFile: Input data file\n"
                       "RESULT:\n"
                       "The operation defined by the formula is performed.\n"
                       "Columns not involved in the formula are unchanged.\n"
                       "InputFile can be piped (|). \n"
                       "Output is written to stdout.\n"
                       "$Id: Calc.man,v 1.4 1999/07/12 10:51:10 herbie Exp herbie $\n";
 
 if(iArgC<1)exit(EXIT_FAILURE);

 if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "s:f:t:vVh");
   if(iO=='s')iSet=strtol(optarg,(char **)NULL, 10);
   if(iO=='f')strncpy(szF,optarg,MAX_LINELENGTH);
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

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Calc", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);}
    	    
  DataFile *DF;     //=new DataFile(szArgV[optind]);
  Formula F(szF);

  int iType,nSets;
 enum EndLine  LT;
  iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Calc: ERROR-exit\n");
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
     default:fprintf(stderr,"ERROR Calc: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
     }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Calc: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Calc: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Calc: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
     
  if(iPrint)
   {String S;
    fprintf(stderr,"Calculating: %s of %s\n",F.GetFormula(),szFile);
    DF->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

 int ir=DF->GetCRData()->Calculate(F);
 
 if(ir<=0)
   {fprintf(stderr,"ERROR %d Calc: %s of %s\n",ir,F.GetFormula(),szFile);
    exit(EXIT_FAILURE);
   }
    
 if(LineT!=NOEnd)DF->SetLineType(LineT);
 if(!DF->SaveData(stdout))
   {fprintf(stderr,"ERROR Calc: Can not write to stdout\n");
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
    if(i)perror("ERROR Calc: Removing tmp-file");
   }
} 		
// *******************************************************************
