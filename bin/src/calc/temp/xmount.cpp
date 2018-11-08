// File: xmount.cpp
// $Log: xmount.cpp,v $
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
#include "browser.h"

#define MAN_PATH "AUSW_MANPATH"

const char* szRCSID={"$Id: xmount.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

 char szAnswer [2048]; 

 QPoint SCREEN=QPoint(0,0);

int main(int iArgC, char ** szArgV)
{         int  iO;
          int  iPrint=0;
   const char *szManPath=getenv(MAN_PATH);
      
 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "hvV");
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='v')iPrint=1;
  }	   
 while (iO!=EOF);

  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( iArgC, szArgV );
  QApplication::setFont( QFont("Helvetica") );
  
 QWidget *d=QApplication::desktop();
 SCREEN.setX(d->width());
 SCREEN.setY(d->height());

 
 MountShare *ms=new MountShare();
 a.setMainWidget( ms );
 return a.exec();

  //exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Xmount: Removing tmp-file");
   }
 
} 		
// *******************************************************************
