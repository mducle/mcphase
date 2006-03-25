// File: polyfit.cpp
// $Log: polyfit.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// -V added
//

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "dfile.h"
#include "formula.h"
#include "stdfunc.h"

#define MAN_PATH "AUSW_MANPATH"
const char *szRCSID={"$Id: polyfit.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);
 
int main(int iArgC, char ** szArgV)
{int iO;
 enum EndLine  LineT = NOEnd;
          int   iSet =  1;
          int iPrint =  0;
          int     ix =  1;
	  int     iy =  2;
	  int     nP =  0;
          int   nOut = -1;
   const char *szManPath=getenv(MAN_PATH);

 if(iArgC<1)exit(EXIT_FAILURE);

 //if(iArgC<2)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "n:s:t:x:y:N:vVh");
   if(iO=='s'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
               else iSet=strtol(optarg,(char **)NULL, 10);
              }
   if(iO=='n'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
 	       nP=strtol(optarg,(char **)NULL, 10);
              } 
   if(iO=='N'){if(!IsNumString(optarg,UINT_NUM))Usage(szManPath,szArgV[0],stderr);
 	       nOut=strtol(optarg,(char **)NULL, 10);
              } 
   if(iO=='t')
     {if(strcmp(optarg,"dos")==0)LineT=DOSEnd;
      else {if(strcmp(optarg,"unix")==0)LineT=UNIXEnd;
            else Usage(szManPath,szArgV[0],stderr);
           } 
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

  //if(iArgC<optind+1)Usage(szManPath,szArgV[0],stderr);
 
  char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file1 name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Fitpoly", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);
    }
    	    
  DataFile *DF;     //=new DataFile(szArgV[optind]);

  int iType,nSets;
 enum EndLine  LT;
  iType=CheckFileType(szFile,LT);

  if( (FileType) iType == FtNONE)
    {fprintf(stderr,"Polyfit: ERROR-exit\n");
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
       default:fprintf(stderr,"ERROR Fitpoly: file %s has unsupported type %d",
                           szFile,iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
 if(nP<0 || nP>10)
    {fprintf(stderr,"ERROR Fitpoly: Illegal order %d (0 ...10)\n",nP);
     exit(EXIT_FAILURE);
    }	    

 if(DF->GetFType()==FtERROR || DF->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Fitpoly: Error reading file %s",szFile);
     exit(EXIT_FAILURE);
    }	    

  if(iSet>nSets)   
    {fprintf(stderr,"ERROR Fitpoly: Illegal # of set %d (1...%d)\n",iSet,nSets);
     exit(EXIT_FAILURE);
    }		
  if( DF->ReadData(iSet)>0 || DF->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Fitpoly: Error reading data of %s\n",szFile);
     exit(EXIT_FAILURE);
    }
    

  if(iPrint)
   {String S;
    fprintf(stderr,"Fit polynom of File:%s\n",szFile);
    fprintf(stderr,"x-Column:%d x-Column:%d Order: %d\n",ix,iy,nP);
    DF->SPrintInfo(S,PR_ALL);
    fprintf(stderr,"%s",(const char *)S);    
   }

   CRData *TD=DF->GetCRData();
   TD->SetColX(ix-1);
   TD->SetColY(iy-1);

 
   Polynom P(nP);   

   if(TD->PolyFit(P,nP,ix-1,iy-1)<=0)
     {fprintf(stderr,"ERROR Fitpoly: Fit failed\n");
      exit(EXIT_FAILURE);
     }

   fprintf(stderr,"Result of polynomal fit:\n");

   double lfM;
      int i;
   for(i=0, lfM=0 ;i<TD->GetSteps(); i++)lfM+=TD->GetPointY(i);
   lfM/=TD->GetSteps();
     
   double lfSoS,lfV;

   for(i=0, lfV=0, lfSoS=0; i<TD->GetSteps(); i++)
      {lfSoS+=SQUARE( P(TD->GetPointX(i)) - TD->GetPointY(i) );
       lfV+=SQUARE(TD->GetPointY(i)-lfM);
      }

   fprintf(stderr,"Squaresum of deviations: %-15.5g\n",lfSoS);
   fprintf(stderr,"r-value: %-15.5g  ", sqrt(fabs(1-lfSoS/lfV)) );
   fprintf(stderr,"Mean deviation: %-12.5g\n\n",lfSoS/TD->GetSteps());

   for(i=0; i<=P.Degree(); i++)
      {fprintf(stderr,"%-14.7g * x^%d\n",P[i],i);}

  
   if(nOut==-1)nOut=TD->GetSteps();
   if(nOut<=0)
     {fprintf(stderr,"ERROR Fitpoly: Illegal steps (-N) %d; cannot make output file\n",nOut);
      exit(EXIT_FAILURE);
     }
     
   TD->NewMinMax();
   
   double lfS = TD->GetMinX();
   double lfE = TD->GetMaxX();
   String S;
   S.Setf("Data file: %s",szFile);
   if(!P.Save(stdout,lfS,lfE,nOut,(const char *)S,LineT))
     fprintf(stderr,"WARNING Fitpoly: Cannot save polynom\n");
	    
 delete DF;
 
 exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i = remove(szTFile);
    fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Fitpoly: Removing tmp-file");
   }
 
} 		
// *******************************************************************
