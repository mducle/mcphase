//  class inipar ... initial parameters for program mcphas
//
#ifndef INIMCDIS
#define INIMCDIS


#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<martin.h>
#include<vector.h>
#include<mfcf.hpp>

class inimcdis
{ private:
  char * savfilename;
    
  public:

  char * info;
  double ** hkls;
  int hkllist; 
  int nofatoms; //nofatoms in primitive cryst unit cell
  int nofcomponents; //number of components of mean field (including magnetic, quadrupolar fields ...
  double T;
  double Ha;
  double Hb;
  double Hc;
  double emax;
  double emin; // energy boundary for dispersion (used for calc. of sta - see manual)
  double ki;
  double kf; // constant ki/kf
  // qvectors to be considered
  Vector qmin,qmax,deltaq;
  mfcf mf;
  
   void inimcdis::errexit();
  inimcdis (char * file,char * spinfile); //constructor
  inimcdis (const inimcdis & p);//kopier-konstruktor
 ~inimcdis ();//destruktor
};

#endif

