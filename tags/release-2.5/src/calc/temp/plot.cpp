// File: plot.cpp
// $Log: plot.cpp,v $
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

const char *szRCSID={"$Id: plot.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

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

 char szB[MAX_LINELENGTH+1];
 char szX[MAX_LINELENGTH+1]=""; 
 char szY[MAX_LINELENGTH+1]=""; 
 const char *szManPath=getenv(MAN_PATH);

//   QList<DataFile> DataList;


 if(iArgC<1)exit(EXIT_FAILURE);

 do
  {iO=getopt(iArgC, szArgV, "x:y:vVh");
  if(iO=='x')strncpy(szX,optarg,MAX_LINELENGTH); 
  if(iO=='y')strncpy(szY,optarg,MAX_LINELENGTH); 
  if(iO=='v')iPrint=1;
  if(iO=='V'){PrintRCS(szRCSID,szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  if(iO=='h'){Usage(szManPath,szArgV[0],stderr); exit(EXIT_FAILURE);}
  }	   
 while (iO!=EOF);


  int iFCount;
  if(!szArgV[optind] || strcmp(szArgV[optind],"-")==0)
             // no file name; is it a pipe ?  
    {CreateTmpFile(szTFile, "Plot", TMPFILE_EXIT);
     if(atexit(RmTmpFile))
	fprintf(stderr,"Internal error: cannot register atexit()\n");
     if(!szArgV[optind])iFCount=1;     
     else {iFCount=iArgC-optind;}
    }
  else
    {iFCount=iArgC-optind;}

  int *ix=0, *iy=0;
  int i;

  ix=new int[iFCount];
  CHECK_POINTER_EXIT(ix,EXIT_FAILURE)

  iy=new int[iFCount];
  CHECK_POINTER_EXIT(iy,EXIT_FAILURE)

  for(i=0;i<iFCount;i++){ix[i]=1;iy[i]=2;}

  for(i=0;i<iFCount && StrTok(szX,',',i,szB);i++)
     {if(*szB)ix[i]=strtol(szB,(char **)NULL,10);}

  for(i=0;i<iFCount && StrTok(szY,',',i,szB);i++)
     {if(*szB)iy[i]=strtol(szB,(char **)NULL,10);}

  char **szFile = new char * [iFCount];
  CHECK_POINTER_EXIT(szFile,EXIT_FAILURE)
 
  for(i=0;i<iFCount;i++)
      {szFile[i]=new char [MAX_LINELENGTH+1];
       CHECK_POINTER_EXIT(szFile[i],EXIT_FAILURE)
      }
  if(*szTFile)strncpy(szFile[0],szTFile,MAX_LINELENGTH);
  //else
    {for(i=(*szTFile ? 1:0);i<iFCount;i++)
            strncpy(szFile[i],szArgV[optind+i],MAX_LINELENGTH); 
    }

  if(iPrint)
    {fprintf(stderr,"Files: %d\n",iFCount);
     for(i=0;i<iFCount;i++)
         fprintf(stderr,"File %d: %s x:%d y:%d\n",i,szFile[i],ix[i],iy[i]);
    }   

  DataFile **DF = new DataFile * [iFCount];     //=new DataFile(szArgV[optind]);
  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)


  int iType;
  for(i=0;i<iFCount;i++)
     {iType=CheckFileType(szFile[i],LineT);

      if( (FileType) iType == FtNONE)
        {fprintf(stderr,"Plot: ERROR-exit\n");
         exit(EXIT_FAILURE);
        }

       switch ((FileGroup)(int)(iType/10)*10)
        {case TheGROUP:
              DF[i]=new TheFile(szFile[i]);
              break;
         case SxSGROUP:
              DF[i]=new SxSFile(szFile[i]);
              break;
        case XDifGROUP:
             DF[i]=new XDifFile(szFile[i]);
             break;
        case SplineGROUP:
             DF[i]=new SplineFile(szFile[i]);
             break;
        case NoGROUP:
             DF[i]=new AsciiFile(szFile[i]);
             break;
       default:fprintf(stderr,"ERROR Plot: file %s has unsupported type %d",
                           szFile[i],iType);
          exit(EXIT_FAILURE);
      }//switch

  CHECK_POINTER_EXIT(DF,EXIT_FAILURE)
 		    
  if(DF[i]->GetFType()==FtERROR || DF[i]->GetFType()==FtILLG)
    {fprintf(stderr,"ERROR Plot: Error reading file %s\n",szFile[i]);
     exit(EXIT_FAILURE);
    }	    

  if( DF[i]->ReadData()>0 || DF[i]->GetCRData()->GetSteps()==0 )
    {fprintf(stderr,"ERROR Plot: Error reading data of %s\n",szFile[i]);
     exit(EXIT_FAILURE);
    }
  if(!DF[i]->GetCRData()->SetColX(ix[i]-1))
    {fprintf(stderr,"ERROR Plot %d: x column %d\n",i,ix[i]);
     exit(EXIT_FAILURE);
    }
   if(!DF[i]->GetCRData()->SetColY(iy[i]-1))
    {fprintf(stderr,"ERROR Plot %d: y column %d\n",i,iy[i]);
     exit(EXIT_FAILURE);
    }

  }//for     


  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( iArgC, szArgV );
  QApplication::setFont( QFont("Helvetica") );
  
 QWidget *d=QApplication::desktop();
 SCREEN.setX(d->width());
 SCREEN.setY(d->height());

  DisplayData *dd=new DisplayData();
  dd->SetEndExit();
  for(i=0;i<iFCount; i++) dd->AddData(DF[i]);	
  a.setMainWidget( dd );
  dd->InitPlot();
  //dd->show();
  return a.exec();

//  delete ix;
//  delete iy;
//  if(szFile){for(i=0;i<iFCount;i++)delete szFile[i];}
//  delete szFile;
//  exit(EXIT_SUCCESS);
}
// ******************************************************************
void RmTmpFile(void)
{//fprintf(stderr,"temp file:%s\n",szTFile);
 if(szTFile[0])
   {int i=remove(szTFile);
    //fprintf(stderr,"Removing tmp-file: %s\n",szTFile);
    if(i)perror("ERROR Fileop: Removing tmp-file");
   }
 
} 		
// *******************************************************************
