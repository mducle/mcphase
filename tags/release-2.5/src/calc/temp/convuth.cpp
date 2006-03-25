// File: convuth.cpp
// $Log: convuth.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.2  1999/03/15 09:08:37  herbie
// -V added
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spline.h"
#include "stdinc.h"

#define BUFLEN 1024

const char *szRCSID={"$Id: convuth.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

int main(int iAcnt, char *szArgV[])
{

 printf("Converting *.UTH to *.TTH \n");
 printf("Sensor File:%s\n",szArgV[1]);


 if(iAcnt!=2){printf("Wrong arguments!\nUsage CONVUTH filename current\n");
              exit(0);}

 FILE *f;
 if( (f=fopen(szArgV[1],"r"))==NULL)
    {printf("Bad filename or file not found (%s)\n",szArgV[1]);
      exit(1);}

 char szB[BUFLEN+2];

 fgets(szB,BUFLEN,f);
 if(strstr(szB,"FIELDTABLE")==0)
   {printf("%s not a valid data file\n",szArgV[1]);exit(1);}
 int iRow,iCol;

 fscanf(f,"%d,%d\n",&iRow,&iCol);
 printf("Converting %s %d rows, %d columns (%u)\n",szB,iRow,iCol,(unsigned)iRow*iCol*sizeof(double));
 char *pc;
 double *lfField=new double[iCol];
 int i;
 if(lfField==NULL){printf("Out of memory\n");exit(1);}

 fgets(szB,BUFLEN,f);
 pc=strtok(szB," ,");
 lfField[0]=atof(pc);
 for(i=1;i<iCol;i++)
	 {if((pc=strtok(NULL," ,"))==NULL)break;
	 lfField[i]=atof(pc);
     }


 double **lfU,*lfT;

 lfT=new double [iRow];

 lfU=new double * [iCol];

 if(lfU==NULL || lfT==NULL){printf("Out of memory\n");exit(1);}

 for(i=0;i<iCol;i++)
    { if( !(lfU[i]=new double [iRow]) ){printf("Out of memory\n");exit(1);} }

 int iCnt;
 for(iCnt=0;iCnt<iRow;iCnt++)
    {if(fgets(szB,BUFLEN,f)==NULL){iCnt=0;break;}
	 pc=strtok(szB," ,");
     lfT[iCnt]=atof(pc);
	 for(i=0;i<iCol;i++)
	    {if((pc=strtok(NULL," ,"))==NULL)break;
	     lfU[i][iCnt]=atof(pc);
	     }//for i
     }//for iCnt
 if(iCnt==0)printf("Unexpected end of file\n");
   else
   {SplineTab *ST=new SplineTab(lfU[0],lfT,iRow);
    if(ST==NULL)printf("Out of memory\n");
      else
       {FILE *fo;
        double t;
         strcpy(szB,szArgV[1]);
         pc=strchr(szB,'.');
         pc[1]='t';pc[2]='t';pc[3]='h';pc[4]=0;

        if( (fo=fopen(szB,"w"))==NULL)printf("Error opening %s\n",szB);
        else
         { fprintf(fo,"FIELDTABLE\n");
           fprintf(fo,"%d,%d\n",iRow,iCol-1);
           szB[0]=0;

           for(i=1;i<iCol;i++)fprintf(fo,"%12.5g,",lfField[i]);
           fprintf(fo,"\n");

           for(iCnt=0;iCnt<iRow;iCnt++)
              {fprintf(fo,"%15.7g %15.7g ",lfT[iCnt],lfU[0][iCnt]);
                for(i=1;i<iCol;i++)
                   {ST->Intpol(lfU[i][iCnt],t);
                    fprintf(fo,"%15.7g",t);
                   }
               fprintf(fo,"\n");
              }//for iCnt
           fclose(fo);
          }//fo!=NULL
        }//ST!==0
      if(ST!=NULL)delete ST;
   }//iCnt!=0

   delete lfField;
   delete lfT;
   for(i=0;i<iCol;i++)delete lfU[i];
   delete lfU;

 printf("File %s converted.\n",szArgV[1]);
 fclose(f);
 exit (0);

}
