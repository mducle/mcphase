// File: strings.h
// $Id: strings.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: strings.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.1  1999/06/29 16:30:56  herbie
// Initial revision
//
#ifndef STRINGS_H
#define STRINGS_H 1

#include "stdinc.h"
// ******************************************************************
class String
{
protected:
 char *b;
   int iSize;

  static char Null;

 int Copy(const String &);
 int Copy(const char *);

 public:
 enum Dir {FROM_LEFT,FROM_RIGHT};

              String(const char * c=NULL){Copy(c);}
              String(const String &s){Copy(s);}

          int Add(const char *c, const int iLen=0 );
          int Addf(const char *szFormat, ...);
          int Set(const char *c, const int iLen=0 );
          int Setf(const char *szFormat, ...);
          int Cut(const unsigned iBehind);
          int Select(const char *szD, const unsigned iStart, int iEnd=-1);
          int Count(const char *szD);
          int GetSize(void) const {return iSize;}
          int IsEmpty(void) const {return (iSize==0 || b[0]==0);}
          int IsZero(void)  const {return (iSize==0 || b==0);}

          int DelSpace(void);
          int Strupr(void);
         void RemoveCR_LF(void);
          int Read(const char *szFile,const int iAppend=0);
          int Read(FILE *f, const int iAppend=0);

          int Write(const char *szFile) const ;
          int Write(FILE *f) const ;

          int Contains(const char *s) const;
  const char * GetBuf(void) const {return b;}
        char  operator[](const int i);

        String & operator=(const String &s)
	                   {if(&s!=this)Copy(s);
			    return *this;
			   }
        String & operator=(const char *c)
	                   {Copy(c);
			    return *this;
			   }
        friend int operator==(const String &s1, const String &s2);
        friend int operator>=(const String &s1, const String &s2);
        friend int operator<=(const String &s1, const String &s2);
	      
          operator const char *() const {return b;}
              ~String();
    const char * operator()(const int i, const enum String::Dir d=FROM_LEFT) const;
    const char * operator()(const char *s) const;
    const char * Find(const char *s, const int iNext=0) const;

};

// ****************************************************
// ****************************************************
int SetDelim(const char *s, String &Delim);

inline int operator==(const String &s1, const String &s2)
                   {return strcmp(s1,s2)==0;}
inline int operator>=(const String &s1, const String &s2)
                   {return strcmp(s1,s2)>0;}
inline int operator<=(const String &s1, const String &s2)
                   {return strcmp(s1,s2)<0;}

// *****************************************************
// *****************************************************
class LineString
{
protected:
 String ** Lines;
    String DeLim;
       int nLines;
       int iPlceFor;
 static String Result;

 void Init(const char *s);
  int Alloc(const int iL);
 void Dealloc(void); 
  int Realloc(const int iNewS);
  int CountLines(const char *s);

public:
  enum Search {S_EQUAL, S_CONTAINS, S_NUMBER, S_ISCOL, S_ISCOL_BREAK};

  LineString(void){Lines=0; nLines=iPlceFor=0; DeLim.Set("");}
  LineString(const char *s){Init(s);}
  LineString(const String &Ls){Init((const char *)Ls);}

           int GetNLines(void) const {return nLines;}
  const char * Line(const int i, const int iWithDel=0);
  const char * GetDelim(void) const {return DeLim.GetBuf();}

           int InsertLine(const char * szL, int iBehind=-1);
           int ReplaceLine(const char * szL, const int iL);
           int RemoveLine(const int iL);

           void Set(const String &S){Dealloc(); Init((const char *)S);}

           int Read(const char *szFile, const char *szEnd=0,
                    const enum Search es=S_CONTAINS);
           int Read(FILE *f, const char *szEnd=0, 
                    const enum Search es=S_CONTAINS);

          int Write(const char *szFile) const;
          int Write(FILE *f) const;
          int Write(String &S, const int iAdd) const;

   const char * Find(int &iL, const char *s,
                     const unsigned iStart=0, const unsigned iNext=1);

  int GetPar(const char *szTopic, const char *szName, char *szFormat, ...);
   const char *FindPar(const char *szTopic, const char *szName, const int iWithDel=0);

};
// *********************************************************
#endif //STRINGS_H
