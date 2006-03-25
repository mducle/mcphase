// File: xport.cpp
// $Log: xport.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// -V added
//

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <qlist.h>
#include <qapplication.h>
#include <qpoint.h>

#include "dfile.h"
#include "stdfunc.h"
#include "grafic.h"
#include "WidGen.h"

#define MAN_PATH "AUSW_MANPATH"

const char* szRCSID={"$Id: xport.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

 QPoint SCREEN=QPoint(0,0);

int main(int iArgC, char ** szArgV)
{         int  iO;
 enum EndLine  LineT=NOEnd;
          int  iPrint=0;
          int  iSet =  1;
   const char *szManPath=getenv(MAN_PATH);
 
  
 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "s:t:vVh");
   if(iO=='s'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               else iSet=strtol(optarg,(char **)NULL, 10);
              }
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

 char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Xport", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }


  DataFile *DF;     //=new DataFile(szArgV[optind]);

  int iType;
 enum EndLine  LT;
  iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Xport: ERROR-exit\n");
     exit(EXIT_FAILURE);
    }

       switch ((FileGroup)(int)(iType/10)*10)
        {case TheGROUP:
              DF=new TheFile(szFile);
              break;
         case SxSGROUP:
              DF=new SxSFile(szFile);
              break;
       case XDifGROUP:
             DF=new XDifFile(szFile);
             break;
        case SplineGROUP:
             DF=new SplineFile(szFile);
             break;
        case NoGROUP:
             DF=new AsciiFile(szFile);
             break;
       default:fprintf(stderr,"ERROR Xport: file %s has unsupported type %d\n",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  int nSets=DF->GetNSets();

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Xport: Error reading file %s\n",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Xport: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Xport: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }

  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( iArgC, szArgV );
  QApplication::setFont( QFont("Helvetica") );
  
 QWidget *d=QApplication::desktop();
 SCREEN.setX(d->width());
 SCREEN.setY(d->height());

 SaveFile *sf=new SaveFile(szFile,DF);
 
 a.setMainWidget( sf );
 if(sf->exec())
   {if(LineT!=NOEnd)DF->SetLineType(LineT);
    if(!DF->SaveData(stdout))
      {fprintf(stderr,"ERROR Fileop: Can not write to stdout\n");
       exit(EXIT_FAILURE);
      } 
    if(iPrint)fprintf(stderr,"Data exported\n");
   }
 else
   { if(iPrint)fprintf(stderr,"Xport cancelled\n");

   }
 if(iPrint)
  {String S;
   DF->SPrintInfo(S,PR_TEXT | PR_TYPE);
   fprintf(stderr,"%s",(const char *)S);    
  }
 delete sf; 
 a.quit();

  //exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Xport: Removing tmp-file");
   }
 
} 		
// *******************************************************************
