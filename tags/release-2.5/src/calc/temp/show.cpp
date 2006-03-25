// File: show.cpp
// $Log: show.cpp,v $
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

#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: show.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

 char szAnswer [2048]; 

 QPoint SCREEN=QPoint(0,0);

int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iPrint=0;
 int ix=1, iy=2;
 int iSet=1; 
 const char *szManPath=getenv(MAN_PATH);
  
 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "s:x:y:vVh");
   if(iO=='s'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               else iSet=strtol(optarg,(char **)NULL, 10);
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

 char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Show", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }


  DataFile *DF;     //=new DataFile(szArgV[optind]);

  int iType,nSets;
      iType=CheckFileType(szFile,LineT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Show: ERROR-exit\n");
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
       default:fprintf(stderr,"ERROR Show: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Show: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Deriv: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		

  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Show: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
  if(!DF->GetCRData()->SetColX(ix-1))
    {fprintf(stderr,"ERROR Show: Illegal x column %d\n",ix);
     exit(EXIT_FAILURE);
    }
   if(!DF->GetCRData()->SetColY(iy-1))
    {fprintf(stderr,"ERROR Show : Illegal y column %d\n",iy);
     exit(EXIT_FAILURE);
    }

  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( iArgC, szArgV );
  QApplication::setFont( QFont("Helvetica") );
  
 QWidget *d=QApplication::desktop();
 SCREEN.setX(d->width());
 SCREEN.setY(d->height());

  ShowWid *sw=new ShowWid(DF);
  sw->SetEndExit();
  sw->show();
  a.setMainWidget( sw );
  return a.exec();

//  exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Show: Removing tmp-file");
   }
 
} 		
// *******************************************************************
