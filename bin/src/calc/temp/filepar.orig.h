#if !defined(FILEPAR_H)
#define FILEPAR_H

// $Id: filepar.orig.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: filepar.orig.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#define MAX_PAR_LEN 82

#define PFILE_NFOUND -1
#define TOPIC_NFOUND -2
#define  NAME_NFOUND -3
#define       NO_MEM  0
#define    BAD_INPUT -4
#define        FOUND  1

#define DECIMAL 1
#define HEXADEC 2

#define  F_CLOSE 1
#define F_NCLOSE 2

//#define FPAR_DEBUG

#if defined (FPAR_DEBUG)
  void InitDebug(void);
  void CloseDebug(void);
#endif


typedef struct fpar_elem
  {char *szLine; struct fpar_elem *Next;} FLINE;

#include <stdio.h>
 // *******************************************************
class FilePar
  {public:

   FLINE *Lines;
     int iLen;
     int iStatus;

   FilePar(const char *szFileName, const int iClose=F_CLOSE);
   FilePar(FILE *f, const int iClose=F_CLOSE);
   ~FilePar(void);
   inline  int GetStatus(void){return iStatus;}
           int GetPar(const char *szTopic, const char *szName, char *szFormat, ...);
   FLINE * FindPar(const char *szTopic, const char *szName);
      void Print(void);
      void PrintTree(void);
      void Update(FILE *f) {DeleteLines();
                            ReadFile(f,F_NCLOSE);}
   private:
   void ReadFile(FILE *f, const int iClose);
    int Insert(const char *NewE);
   void DeleteLines(void);

  };
//*************************************************************

#endif //FILEPAR_H



