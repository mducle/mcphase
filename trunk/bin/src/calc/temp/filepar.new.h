#ifndef FILEPAR_H
#define FILEPAR_H

// $Id: filepar.new.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: filepar.new.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.2  1999/06/28 15:37:30  herbie
// Drastic change from char to String
//

#include <stdio.h>

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

#define PFILE_NFOUND 1
#define TOPIC_NFOUND 2
#define  NAME_NFOUND 3
#define       NO_MEM 4
#define    BAD_INPUT 5
#define INVALID_FILE 6
#define        FOUND 0

#define DECIMAL 1
#define HEXADEC 2

#define  F_CLOSE 1
#define F_NCLOSE 2

 // *******************************************************
class FilePar
{

protected:
   String Pars;
   String DeLim;
      int nLines;
      int iStatus;
  static char *szMsg[];

  void ReadFile(FILE *f, const int iClose);

public:
   FilePar(const char *szFileName, const int iClose=F_CLOSE);
   FilePar(FILE *f, const int iClose=F_CLOSE);
   ~FilePar(void);
   int GetPar(const char *szTopic, const char *szName, char *szFormat, ...);
   const char * FindPar(const char *szTopic, const char *szName);
   void Print(FILE *f) const;
   void Update(FILE *f) {ReadFile(f,F_NCLOSE);}

   int GetStatus(void) const {return iStatus;}
   const char *GetMsg(void){return szMsg[iStatus];}

};
//*************************************************************

#endif //SFILEPAR_H



