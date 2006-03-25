// File: theview.cpp
// $Log: theview.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.1  1999/05/06 10:22:55  herbie
// Initial revision
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

const char *szRCSID={"$Id: theview.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

 char szAnswer [2048]; 

 QPoint SCREEN=QPoint(0,0);

int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 const char *szManPath=getenv(MAN_PATH);
  
 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "Vh");
  if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

 char szDFile[MAXPATH+1];
 char szAFile[MAXPATH+1];

 strncpy(szDFile,szArgV[optind],MAXPATH);
 szDFile[MAXPATH]=0;
 char *p=strrchr(szDFile,'.');
 if(!p)
   {fprintf(stderr,"Theview: Bad filename %s\n",szDFile);
    exit(EXIT_FAILURE);
   }

 strncpy(szAFile,szDFile, p-szDFile+1);
 szAFile[p-szDFile+1]=0;
 strcat( szAFile,"aux");
 
  TheFile *DF,*AF;     //=new DataFile(szArgV[optind]);

  int iType=CheckFileType(szDFile,LineT);

  if( (FileType) iType != FtTHE_CAPV1)
    {fprintf(stderr,"Theview: Unsupported file type: %d\n",iType);
     exit(EXIT_FAILURE);
    }

  iType=CheckFileType(szAFile,LineT);
  if( (FileType) iType != FtTHE_CAPV1)
    {fprintf(stderr,"Theview: Unsupported file type: %d\n",iType);
     exit(EXIT_FAILURE);
    }

   DF=new TheFile(szDFile);

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Theview: Error reading file %s",szDFile);
     exit(EXIT_FAILURE);
    }	    

  if( DF->ReadData(1)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Theview: Error reading data of %s\n",szDFile);
     exit(EXIT_FAILURE);
    }

   AF=new TheFile(szAFile);

  CHECK_POINTER_EXIT(AF,EXIT_FAILURE)
 		    
  if(AF->GetFType()==FtERROR || AF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Theview: Error reading file %s",szAFile);
     exit(EXIT_FAILURE);
    }	    

  if( AF->ReadData(1)>0 || AF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Theview: Error reading data of %s\n",szAFile);
     exit(EXIT_FAILURE);
    }
    
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( iArgC, szArgV );
  QApplication::setFont( QFont("Helvetica") );
  
 QWidget *d=QApplication::desktop();
 SCREEN.setX(d->width());
 SCREEN.setY(d->height());

  TheView *sw=new TheView(DF,AF);
  a.setMainWidget( sw );
  sw->SetEndExit();
  sw->show();
  return a.exec();

 exit(EXIT_SUCCESS);
}
// ******************************************************************
