// File: clpar.cpp
// $Log: clpar.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "stdinc.h"
#include "stdfunc.h"
#include "strings.h"

#define MAN_PATH "AUSW_MANPATH"

const char *szRCSID={"$Id: clpar.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

 extern char *optarg;
 extern int optind, opterr, optopt;

 char szTFile[MAXPATH+1]="";

 void RmTmpFile(void);

// ***************************************************
// Data file to calculate lattice pars
// [GenPars]
// System=ORTHO
// # Other  values are HEXAG, TETRA, CUBIC
// WaveL=1.792
// # Wavelength in Angstroem
// nRef=4
// # Number of reflexes which are following
// [ReflPars]
// R1=1,1,1,35.798,1
// R2=2,2,0,77.988,1
// ...
// R4=4,4,2,120.345,0
// # h,k,l,2theta,validity
// # h,k,l: miller indizes of the reflex
// # 2theta: two theta value of the reflex
// # validity: 0 or 1 1: reflex is used for calculation 0: not used
// *****************************************************
#define MAX_CRY 4
#define MAX_SYS_LEN 6
char szCryStr[MAX_CRY][MAX_SYS_LEN]={"ORTHO","HEXAG","TETRA","CUBIC"};

struct Refl {int h,k,l;
             double th2;
             int v;
            };
// *****************************************************
int main(int iArgC, char ** szArgV)
{int iO;
 const char *szManPath=getenv(MAN_PATH);
  
 if(iArgC<1)Usage(szManPath,szArgV[0],stderr);

 do
  {iO=getopt(iArgC, szArgV, "hV");
   if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
   if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);


 char szFile[MAXPATH+1];

  if(!szArgV[optind]) // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Clpar", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     strncpy(szFile,szTFile,MAXPATH);
    }
  else
    {strncpy(szFile,szArgV[optind],MAXPATH);}

 LineString LS;  
 if(LS.Read(szFile)<=0)
   {fprintf(stderr,"Can not read input file %s\n",szFile);
    exit(EXIT_FAILURE);
   }

 int iErr=0;
 char szCSyst[MAX_SYS_LEN]
 if(!GetPar("[GenPars]","System", "%5s",szCSyst))
   {iErr=1;goto err_ret;}

 double lfWl;
 if(!GetPar("[GenPars]","WaveL", "%lf",&lfWl) || lfWl<=1)
   {iErr=2;goto err_ret;}

 int nRef;
 if(!GetPar("[GenPars]","nRef", "%d",&nRef) || nRef<=1)
   {iErr=3;goto err_ret;}

 int i,iC;

 int h,k,l,v;
 double th2;

 String S;
 struct Refl * R;
 R=new struct Refl [nRef];
 CHECK_POINTER_EXIT(R,EXIT_FAILURE)

 for(i=0,iC=0;i<nRef;i++)
    {S.Setf("R%d",i+1);
     if(!GetPar("[ReflPars]",(const char *)S, "%d,%d,%d,%lf,%d",
                            &(R->h),&(R->k),&(R->l),&(R->th2),&(R->v)))
        {iErr=4+i;goto err_ret;}
     if(R->v)iC++;
    }

 
 LS.Write(stdout);

 exit(EXIT_SUCCESS);

err_ret:
 fprintf(stderr,"Error reading parameter %d from file\n",iErr);
 exit(EXIT_FAILURE);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Clpar: Removing tmp-file");
   }
} 		
// *******************************************************************
