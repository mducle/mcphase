// martins functions
// to get on better in c

#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<martin.h>
#include<vector.h>

#ifndef __linux__
#include<cfloat>
#endif

// extract parameter 'parameter'  from string instr (z.B. "blabla dmin=0.2 blabla") -
// output: var ... value of parameter
// returns 1 on error and 0 if successful
  int extract(char * instr,const char * parameter,double & var)
{ //const char delimiters[] = " =:\n";
  char *token;
  
// check if line is comment line -> if yes return 0
//  if (instr[0]=='#') return 1; //removed 26.5.02 in order to be able to place parameters in comment lines
 
// strip /r (dos line feed) from line if necessary
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

  if ((token = strstr (instr, parameter))==NULL) return 1; // parameter string not found
  token+=strlen(parameter);
  if (strstr (token, "=")==NULL) return 1;  // no '=' found after parameter string
  while(strstr(token," ")==token||strstr(token,"\t")==token)++token;
  if (strstr(token,"=")!=token) return 1; // there are other characters than tab or spaces between parameter and =
  ++token;
  var = strtod (token, NULL);
  return 0;
}
// same for in and double 
// extract a variable named [parmeter] into var out of a string [instr]
  int extract(char * instr,const char * parameter,int & var)
      {double dd;if(0==extract(instr,parameter,dd)){var=(int)dd;return 0;}else{return 1;}}
  int extract(char * instr,const char * parameter,float & var)
      {double dd;if(0==extract(instr,parameter,dd)){var=(float)dd;return 0;}else{return 1;}}
//the same for a string ... maximal n characters are copied
  int extract(char * instr,const char * parameter,char * var,size_t n)
{ const char delimiters[] = " \n";
  char *token;
  
// check if line is comment line -> if yes return 1
//  if (instr[0]=='#') return 1; //removed 26.5.02 in order to be able to place parameters in comment lines
 
// strip /r (dos line feed) from line if necessary
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

  if ((token = strstr (instr, parameter))==NULL) return 1;
  token+=strlen(parameter);
  token = strstr (token, "=")+1;
  while (*token==' '){++token;} // remove starting spaces
  strncpy (var,token, n);
  //remove from string var all characters after delimiters
  strtok(var,delimiters);
  return 0;
}


//open a file: similar fopen but with error check 
FILE * fopen_errchk (const char * filename,const char * mode)
{ FILE *file;
 
  errno = 0;
  file = fopen (filename,mode);
  if (file == NULL)
    { 
      fprintf (stderr, "Couldn't open file %s: %s\n",filename, strerror (errno));
      exit (EXIT_FAILURE);
     }
  return file;
} 

//input file from string with error check
char * fgets_errchk (char * instr,int size, FILE * file)
{char * s;
 errno=0;
 s=fgets(instr,size,file);
 if (s==NULL)
 {fprintf (stderr,"ERROR fgets_errchk: end of file\n"); 
  exit(EXIT_FAILURE);}
 return s;
}

// **************************************************************************
// input line and split into numbers
// on input: nn[0]... size of array nn[] (not changed in routine)
//           fin_coq..filepointer
// returns nof numbers read and
// 0 if no numbers have been read (e.g. empty line or line is comment line starting with # or EOF)

#define maxnofcharinline  1024


char * mystrtok (char * s, char * delimiters)
{char * pointer;
 pointer=s+strcspn(s,delimiters);
 pointer+=strspn(pointer,delimiters);
 if ((*pointer)=='\0')pointer=NULL;
return pointer;
}

// function to input a line of numbers separated by delimiters
// example:
// 3 23 542 23
// returns:0 .... it is a comment line (starting with #) or empty line
//         n .... number of numbers read
int inputline (FILE * fin_coq, float *nn)
{
  char instr[maxnofcharinline];
  char delimiters[] = " \n\t";
  char *token;
  int i;
  errno=0;

  if (fgets (instr, sizeof (instr), fin_coq) == NULL)
    { return 0;}

 
  if (feof(fin_coq)==0&&strchr(instr,'\n')==NULL)
    { fprintf (stderr, "Error in function inputline: input string too long");
      exit (EXIT_FAILURE);
     }

// strip /r (dos line feed) from line if necessary
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

// check if line is comment line -> if yes return 0
  if (instr[strspn(instr," \t")]=='#') return 0;
 

  //initialize token to first nonspace character
  token = instr+strspn(instr," \t");
  i=0;
  if (*token!='\n'){
  for (i = 1; token != NULL ; ++i)
    {
if(i>=(int)nn[0])
        { fprintf (stderr, "Error in function inputline: maximum value of numbers in line exceeded,more numbers in line (>%i) to be read.\n",i);
           exit (EXIT_FAILURE);
        }
    
      nn[i] = strtod (token, NULL);
      token = mystrtok (token, delimiters);
 // printf("i=%i token=%g\n",i,nn[i]);
    }
  }        
//  printf("i=%i token=%g\n",i,nn[i]);

  if (i<1) {return 0;}
  return i-1;
}

// function to input a line of numbers separated by delimiters
// example:
// parname=3 23 542 23
// returns:0 .... it is a comment line (starting with #) or parameter parname not found
//         n .... number of numbers read
int inputparline (const char * parname, FILE * fin_coq, float *nn)
{
  char instr[maxnofcharinline];
  char delimiters[] = " \n";
  char *token;
  int i;
  errno=0;

  if (fgets (instr, sizeof (instr), fin_coq) == NULL)
    { return 0;}
  
  if (feof(fin_coq)==0&&strchr(instr,'\n')==NULL)
    { fprintf (stderr, "Error in function inputparline: input string too long");
      exit (EXIT_FAILURE);
     }
// check if line is comment line -> if yes return 0
  if (instr[strspn(instr," \t")]=='#') return 0;

// strip /r (dos line feed) from line if necessary
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
 
// check line contains parname -> if no return 0
  if ((token = strstr (instr, parname))==NULL) return 0;

  token+=strlen(parname);
  token = strstr (token, "=")+1;
  while (*token==' '){++token;} // remove starting spaces
  //+ initialize token to first nonspace character
  i=0;
  if (*token!='\n'){
  for (i = 1; token != NULL ; ++i)
    {
if(i>=(int)nn[0])
        { fprintf (stderr, "Error in function inputparline: maximum value of numbers in line exceeded,more numbers in line (>%i) to be read.\n",i);
           exit (EXIT_FAILURE);
        }
    
      nn[i] = strtod (token, NULL);
      token = mystrtok (token, delimiters);
 // printf("i=%i token=%g\n",i,nn[i]);
    }
  }        
//  printf("i=%i token=%g\n",i,nn[i]);

  if (i<1) {return 0;}
  return i-1;
}

 // *************************************************************************

// return random number between 0 and z
float rnd(float z){return  z*rand()/RAND_MAX;};

// return integer of floating number (like basic integer function)
float integer (float s){  double result;modf(s,&result);  return result;};

//return rounded integer of floating number
int cint (float s){  double result;if(modf(s,&result)<0.5){return (int)result;}{return (int)result+1;}};


// factorial of an integer number
int factorial(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}

// return threej symbol 
float threej (float AJ1,float  AJ2,float  AJ3,float AM1,float AM2,float AM3)
{//  THIS IS A PROGRAM TO CALCULATE THE 3_J SYMBOL
//    (AJ1,AJ2,AJ3;AM1,AM2,AM3)
//    F(I)=I! IS THE FACTORIAL FUNCTION
//    I = 0 TO 30 IS USUSALLY ENOUGH
//    PLEASE NOTE THAT I STARTS FROM ZERO
//
//    "CRYSTAL FIELD HANDBOOK", EDITED BY D.J. NEWMAN AND BETTY NG
//    CAMBRIDGE UNIVERSITY PRESS 
//    ... and transfered to c and improved by M. Rotter 2008

        int I,J,M,N,J7;
	double R0,R4,R5,R6,R7,R8,R9;
	double F[31];      

//  DEFINE THE FACTORIAL FUNCTION
        F[0] = 1;for(I=1;I<=30;++I){F[I]=F[I-1]*I;}

//   SELECTION RULES
     if (AM1 + AM2 + AM3 != 0) { return 0;}
     if (AJ1 + AJ2 - AJ3 < 0) { return 0;}
     if (AJ3 + AJ1 - AJ2 < 0) { return 0;}
     if (AJ3 + AJ2 - AJ1 < 0) { return 0;}
     if (AJ1 + AJ2 + AJ3 + 1 < 0) { return 0;}
     if (AJ1 - AM1 < 0) { return 0;}
     if (AJ2 - AM2 < 0) { return 0;}
     if (AJ3 - AM3 < 0) { return 0;}
     if (AJ1 + AM1 < 0) { return 0;}
     if (AJ2 + AM2 < 0) { return 0;}
     if (AJ3 + AM3 < 0) { return 0;}

   R4 = F[int(AJ1 + AJ2 - AJ3)] * F[int(AJ1 - AM1)] * F[int(AJ2 - AM2)] * F[int(AJ3 - AM3)] * F[int(AJ3 + AM3)];
   R5 = F[int(AJ1 + AJ2 + AJ3 + 1)] * F[int(AJ3 + AJ1 - AJ2)] * F[int(AJ3 + AJ2 - AJ1)] * F[int(AJ1 + AM1)] * F[int(AJ2 + AM2)];
   R6 = 0;
   for (J7 = 0;J7<=25;++J7)
     { if (AJ1 + AM1 + J7 >= 0 && AJ2 + AJ3 - AM1 - J7 >= 0 && AJ3 + AM3 - J7 >= 0 && AJ1 - AM1 - J7 >= 0
           && AJ2 - AJ3 + AM1 + J7 >= 0) 
       {R8 = F[int(AJ1 + AM1 + J7)] * F[int(AJ2 + AJ3 - AM1 - J7)];
        R8 *= odd(int(AJ1 - AM1 - J7)) ? -1 : 1; // .... equiv to  * (-1) ^ (int(AJ1 - AM1 - J7));
        R9 = F[J7] * F[int(AJ3 + AM3 - J7)] * F[int(AJ1 - AM1 - J7)] * F[int(AJ2 - AJ3 + AM1 + J7)];
        R6 +=  R8 / R9;
       }    
     }
     R0 = sqrt(R4 / R5) * R6 ;
     R0 *= odd(int(AJ1 - AJ2 - AM3)) ? -1 : 1; //.... equiv to * (-1) ^ int((AJ1 - AJ2 - AM3))
  return R0;
};



#ifndef __linux__
// /* not needed any more in mingw3.1
// return rounded value
//double rint(double value)
//{if (fabs(floor(value)-value)<fabs(ceil(value)-value))
//  {return floor(value);}
// else
//  {return ceil(value);}
//}


//copysign function
double copysign(double a,double b)
{//if fabs(a)
 //return a*a/fabs(a)*b/fabs(b);
#if defined CYGWIN || defined __CYGWIN__ 
 return copysign(a,b);
#else
 return _copysign(a,b);
#endif
}
// */
//sleep function
int sleep(int a)
{
#if defined CYGWIN || defined __CYGWIN__ 
#include<unistd.h>
  return usleep(a*1000);
#else
  _sleep(a*1000);
  return true;
#endif
}

#endif

