// File: convtrd.cpp
// $Log: convtrd.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.2  1999/03/15 09:08:37  herbie
// -V added
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

#define BUFLEN 1024

#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: convtrd.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

int main(int iAcnt, char *szArgV[])
{

 const char *szManPath=getenv(MAN_PATH);
 if(iAcnt!=3){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}

 fprintf(stderr,"Converting *.TRD to *.UTH \n");
 fprintf(stderr,"Sensor File:%s\n",szArgV[1]);


 if(iAcnt!=3){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}

 char szFile[MAXPATH+1];

 if(!strcmp(szArgV[1],"-")) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Deriv", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[1],MAXPATH);
    }

 FILE *f=fopen(szFile,"r");
 CHECK_FILE_POINTER_EXIT(f,szFile,EXIT_FAILURE);

 char szB[BUFLEN+2];
 while(1)
      {if(fgets(szB,BUFLEN,f)==0)
         {fprintf(stderr,"Unexpected End of File in %s\n",szArgV[1]);
          fclose(f);
          exit(EXIT_FAILURE);}

       if(strchr(szB,'}')!=0)break;
      }

 FILE *fo=stdout;
 char *pc;
 strcpy(szB,szArgV[1]);
 pc=strchr(szB,'.');
 pc[1]='u';pc[2]='t';pc[3]='h';pc[4]=0;

 fprintf(fo,"FIELDTABLE\n");
 int i=0,iCnt=0;
 double t,r[10],dr[10],fi=atof(szArgV[2]);
 char szN[80];
 while(fgets(szB,BUFLEN,f)!=NULL)
	  {iCnt++;
       pc=strtok(szB," ,");
       t=atof(pc);
	   for(i=0;i<10;i++)
		 {if((pc=strtok(NULL," ,"))==NULL)break;
		  r[i]=atof(pc);

          if((pc=strtok(NULL," ,"))==NULL)break;
		  dr[i]=atof(pc);
	     }//for i
       memset(szB,0,BUFLEN);
       sprintf(szN,"%15.8g ",t);
       strcat(szB,szN);
       for(i=0;i<10;i++)
          {sprintf(szN,"%15.8g  ",r[i]*fi);
           strcat(szB,szN);
          }//for i
           fprintf(fo,"%s\n",szB);
       }//while
 fprintf(stderr,"File %s converted. %d values I=%12.5g\n",szArgV[1],iCnt,fi);
 fclose (fo);
 fclose(f);
 exit (EXIT_SUCCESS);

}
// *****************************************
void RmTmpFile(void)
{fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Deriv: Removing tmp-file");
   }
 
} 		
// *****************************************
