// File: xfun.cpp
// $Log: xfun.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// *** empty log message ***
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
#include "sfit.h"

#define MAN_PATH "AUSW_MANPATH"

const char* szRCSID={"$Id: xfun.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szAnswer [2048]; 

 QPoint SCREEN=QPoint(0,0);

int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine LineT=NOEnd;
 int iPrint=0;
 const char *szManPath=getenv(MAN_PATH);
 
 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "t:vVh");
  if(iO=='v')iPrint=1;
  if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
     }	   
  }	   
 while (iO!=EOF);

  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( iArgC, szArgV );
  QApplication::setFont( QFont("Helvetica") );
  
 QWidget *d=QApplication::desktop();
 SCREEN.setX(d->width());
 SCREEN.setY(d->height());

  FunWid *fw=new FunWid;
  
  a.setMainWidget( fw );
  return a.exec();

//  exit(EXIT_SUCCESS);
}
// ******************************************************************
