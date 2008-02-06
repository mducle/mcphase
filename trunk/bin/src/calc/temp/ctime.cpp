// File: cdate.cpp
// $Log: ctime.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.1  1999/03/17 14:17:25  herbie
// Initial revision
//

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "stdfunc.h"

#define MAN_PATH "AUSW_MANPATH"

const char* szRCSID={"$Id: ctime.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

int main(int iArgC, char ** szArgV)
{int iO;
 const char *szManPath=getenv(MAN_PATH);
  
 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "vVh");
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);

 long int cT=atol(szArgV[optind]);
 fprintf(stdout,"%ld > %s",cT,ctime(&cT));
 exit(EXIT_SUCCESS);
}
// ******************************************************************
