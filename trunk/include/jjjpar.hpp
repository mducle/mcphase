// class of cf and exchange parameters for one specific atom
// in crystallographic unit cell
#ifndef JJJPAR
#define JJJPAR

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<ctime>
#include<martin.h>
#include<vector.h>
#include<ionpars.hpp>

#ifdef __linux__
#include<dlfcn.h>
#endif


class jjjpar
{
  private:

  void (*m)(Vector*,double*,Vector*,double*,Vector*,double*,double*);  
  int  (*dm)(int*,double*,Vector*,double*,Vector*,ComplexMatrix*,float*);
  void *handle;

  Vector & kramer (double & T,Vector & H, double & Z,double & U);
  int  kramerdm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);

  ionpars * iops;
  int intern_mcalc;
  char * cffilename;
  Vector ABC;         
  Vector magFF; // magnetic formfactor numbers
  double DWF; // DebeyWallerFactor 


  public:

  Vector xyz,J;
  int paranz;
  int nofcomponents; // number of moments (components of moment vector)
  double gJ;
  Matrix *jij;
  Vector *dn; // exchange, coordinates 
  int *sublattice; // sublattice of neighbour
  int diagonalexchange;  // switches 1=exchange is diagonal

// subroutine to calculate momentum <J> from effective field gjmbH [meV]
   Vector &  mcalc (double & T,Vector &  gjmbH, double & Z,double & U);
   Vector &  tetan(); //returns stevens parameters if possible
   int  dmcalc (double & T,Vector &  gjmbheff, ComplexMatrix & mat,float & delta);
   int transitionnumber; // the transition associated with the ion (important if there are more in the single ion spectrum)

   jjjpar (FILE * fin); //konstruktor with filehandle
   jjjpar (int n=1,int diag=0,int nofmom=3); // konstructor without file
   jjjpar (const jjjpar & jjjpars);	// kopier-konstruktor
   
~jjjpar ();		//destruktor
   void add(jjjpar & b, Vector & abc); // add parameters b to this
   void addpars (int number, jjjpar & addjjj); // enlarge the set of parameters by
                                                        // inserting a new exchange parameters addjjj
							// into field at position number
   void save (FILE *file); // to save the parameters to a filehandle
   void saveatom (FILE *file); // to save the atom coordinates and properties to a filehandle

//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
// 'magnetic formfactors
//  j0 = ff(1) * EXP(-ff(2) * s * s) + ff(3) * EXP(-ff(4) * s * s)
//  j0 = j0 + ff(5) * EXP(-ff(6) * s * s) + ff(7)
//  j2 = ff(8) * s * s * EXP(-ff(9) * s * s) + ff(10) * s * s * EXP(-ff(11) * s * s)
//  j2 = j2 + ff(12) * s * s * EXP(-ff(13) * s * s) + s * s * ff(14)
//  F = (j0 + j2 * (2 / gJ - 1))  formfactor F(Q)
   double F(double & Q);

//   debeywallerfactor = EXP(-2 * DWF *s*s)
   double debeywallerfactor(double & Q);

//    jjjpar & operator= (const jjjpar & op2); // zuweisung

};

#endif
