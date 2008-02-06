// File: xhelp.cpp
// $Log: xhelp.cpp,v $
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

const char* szRCSID={"$Id: xhelp.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

 char szAnswer [2048]; 

 QPoint SCREEN=QPoint(0,0);

 extern char *optarg;
 extern int optind, opterr, optopt;

int main(int iArgC, char ** szArgV)
{int iO;
 const char *szManPath=getenv(MAN_PATH);
  
 if(iArgC<1)exit(EXIT_FAILURE);

 char szFile[MAXPATH+1];
 do
  {iO=getopt(iArgC, szArgV, "vVh");
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Xhelp", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }

  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( iArgC, szArgV );
  QApplication::setFont( QFont("Helvetica") );
  
 QWidget *d=QApplication::desktop();
 SCREEN.setX(d->width());
 SCREEN.setY(d->height());

 ShowTxt *st=new ShowTxt(szFile);
 st->SetEndExit();
  if(!st->FileValid())
   {fprintf(stderr,"ERROR Xhelp: Can not open %s\n",szFile);
    a.quit();
   }
 
 a.setMainWidget( st );
 return a.exec();

  //exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Xhelp: Removing tmp-file");
   }
 
} 		
// *******************************************************************
