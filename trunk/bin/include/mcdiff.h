//libraries, classes and  functions in mcdiff




#define MAXNOFREFLECTIONS 6000
#define SMALL 1e-8
#include "../../version"
#include <mpspecfunp.h>
#include <martin.h>
#include <myev.h>
#include <jjjpar.hpp>
#include <complex>
#include<cstddef>
#include<spincf.hpp>

// get intensity of one reflection
int getint(jjjpar ** jjjpars,int hi,int ki,int li,float thetamax,Vector rez1,Vector rez2, Vector rez3,float scale,double T,float lambda,float ovalltemp,int lorenz,int & n,int * J,float & d,float & Theta,float & Imag,float & Imagdip,float & inuc,float & SF,float & lorentzf,complex <double> & mqx,complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2);

// calculate pattern
void neutint(jjjpar ** jjjpars,int code,double T,float lambda, float thetamax, float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int & n,int * J,int & m,Vector *  hkl,float * D,float * theta,float * intmag,float * intmagdip,float * ikern,float * sf,float * lpg,complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2);

// mcdiff - output of results
void printeln(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T, Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int * J,int m,Vector * hkl,float * ikern,float * intmag,float * intmagdip,float * D,float * theta,float * sf,float * lpg,complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,float a,float b,float c);

// output to mcdiff.sps
void print_sps(const char * filename,int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H);
