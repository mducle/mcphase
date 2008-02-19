// class of cf and exchange parameters for a given 
// crystal: corresponds to paramters given in mcphas.j
#ifndef PAR
#define PAR

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<ctime>
#include<martin.h>
#include<vector.h>
#include "jjjpar.hpp"

#define MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL 4000

class par
{ 
  private:
   
   
  public:
  
  char *rems[MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL];
  jjjpar **jjj; // pointer of field of exchange parameter sets
  Vector gJ;
     
  //lattice
  float a,b,c,alpha,beta,gamma;
  Matrix r,rez;
  int nofatoms;
  int nofcomponents;

  //jjjpar 
   
par (const char *filejjj);	//konstruktor
   par (const par & pars);	// kopier-konstruktor
   
~par ();		//destruktor

int newatom(jjjpar * p); //creates new atom from an existing and returns its index
void add(par & b); // add exchange parameters
void save(FILE * fout); // save lattice, atoms and exchange parameters to file
void savelattice(FILE *fout);// save lattice to file
void saveatoms(FILE *fout);// save atom positions and properties  to file
};

#endif