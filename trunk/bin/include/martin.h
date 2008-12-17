// martins functions
// to get on better in c


#ifndef MARTIN_FOPEN_ERRCHK
#define MARTIN_FOPEN_ERRCHK

#if __GNUC__ > 2
#include <unistd.h>
#endif

#include <stdio.h>
// extract parameter 'parameter'  from string instr (z.B. "blabla dmin=0.2 blabla") -
// output: var ... value of parameter
// returns 1 on error and 0 if successful
extern  int extract(char * instr,const char * parameter,double & var);

// extract a variable named [parmeter] into var out of a string [instr]
extern   int extract(char * instr,const char * parameter,int & var);
extern   int extract(char * instr,const char * parameter,float & var);
extern   int extract(char * instr,const char * parameter,char * var, size_t n);


// open file: like fopen but with error check 
extern  FILE * fopen_errchk (const char * filename, const char * mode);
// get string instr: like fgets but with error check
extern  char * fgets_errchk (char * instr,int size, FILE * file);


// input line and split into numbers nn[1 ... ]
// on input: nn[0]... size of array nn[] (not changed in routine)
//           fin_coq..filepointer
// returns nof numbers read and
// 0 if end of file or no numbers have been read 
extern  int inputline (FILE * fin_coq, float *nn);

// function to input a line of numbers separated by delimiters
// example:
// parname=3 23 542 23
// returns:0 .... it is a comment line (starting with #) or parameter parname not found
//         n .... number of numbers read
extern int inputparline (const char * parname, FILE * fin_coq, float *nn)
;

// return random number between 0 and z
extern  float rnd(float z);

// return integer of floating number (like basic integer function)
extern  float integer (float s);

//return rounded integer of floating number
extern int cint (float s);

// return threej symbol 
extern float threej (float AJ1,float AJ2,float AJ3,float AM1,float AM2,float AM3);


#ifndef __linux__
// return rounded integer (not needed any more in MINGW 3.1.3)
//extern double rint(double value);

//extern double copysign(double a,double b);
extern int sleep(int a);
#endif


#endif

