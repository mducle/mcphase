// jjjpar is a class to store all parameters associated
// with a specific ion, for example CEF parameters
// and exchange parameters 
// it is used by many programs in the package
// moreover, it loads also the user defined single ion module functions (linux only)
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
#else
#include <windows.h>
#endif


class jjjpar
{

  public:
   
   // subroutine to calculate momentum <J> from effective field gjmbH [meV]
   Vector &  mcalc (double & T,Vector &  gjmbH, double & Z,double & U,ComplexMatrix & ests);

   // returns transition element matrix M  and transition energy delta (to calculate chi0 in mcdisp,see manual)
   int  dmcalc (double & T,Vector &  gjmbheff, ComplexMatrix & mat,float & delta,ComplexMatrix & ests);
   int transitionnumber; // the transition associated with the ion (important if there are more in the single ion spectrum)

   // returns transition element matrix N(Q) in order to be able to go beyond 
   // dipolar approximation in mcdisp - it requires a call to eigenstates first
   int dncalc(Vector & Qvec, double & T, ComplexMatrix & nat, ComplexMatrix & ests);

   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est, requires a call to eigenstates first
   ComplexVector & MQ(Vector & Qvec);
   ComplexVector Mq;

   ComplexMatrix est; // eigenstates 
   // returns eigenvalues and eigenstates matrix parameters of ion (if possible)
   ComplexMatrix & eigenstates (Vector & gjmbheff, double & T);

  char * cffilename; // single ion parameter filename
  double SLR,SLI; // scattering length
  double DWF; // DebeyWallerFactor [A^2] 
  Vector magFFj0; // magnetic formfactor numbers
  Vector magFFj2; // magnetic formfactor numbers
  Vector magFFj4; // magnetic formfactor numbers
  Vector magFFj6; // magnetic formfactor numbers
  Vector Zc;      // Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)

//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
// 'magnetic formfactors
//  j0 = ff(1) * EXP(-ff(2) * s * s) + ff(3) * EXP(-ff(4) * s * s)
//  j0 = j0 + ff(5) * EXP(-ff(6) * s * s) + ff(7)
//  j2 = ff(8) * s * s * EXP(-ff(9) * s * s) + ff(10) * s * s * EXP(-ff(11) * s * s)
//  j2 = j2 + ff(12) * s * s * EXP(-ff(13) * s * s) + s * s * ff(14)
//  F = (j0 + j2 * (2 / gJ - 1))  formfactor F(Q)
//  RETURN TOTAL FORMFACTOR, 
//    however if gJ=0 and Q>0 return spin form factor FS(Q)=<j0(Q)>
//            if gJ=0 and Q<0 return angular  form factor FL(Q)=<j0(Q)>+<j2(Q)>
   double F(double Q);

//   debyewallerfactor = EXP(-2 * DWF *s*s)
   double debyewallerfactor(double & Q);

   double J(); // returns total angular momentum if possible
   Vector &  tetan(); //returns stevens parameters if possible

  Vector xyz,mom; // atom position, moment
  int paranz;   // number of exchange parameters
  int nofcomponents; // number of moments (components of moment vector)
  double gJ;
  Matrix *jij; // exchange constants 
  Vector *dn; // exchange - coordinates of neighbors 
  int *sublattice; // sublattice of neighbours
  int diagonalexchange;  // switch 1=exchange is diagonal, 0=exchange is not diagonal
   void increase_nofcomponents(int n); // increase nofcomponents by n
   void add(jjjpar & b, Vector & abc); // add parameters b to this
   void addpars (int number, jjjpar & addjjj); // enlarge the set of parameters by
                                                        // inserting a new exchange parameters addjjj
							// into field at position number

   void save (FILE *file); // to save the parameters to a filehandle
   void saveatom (FILE *file); // to save the atom coordinates and properties to a filehandle


   jjjpar (FILE * fin); //konstruktor with filehandle of mcphas.j file
   jjjpar (double x, double y, double z,char * sipffile); // constructor with filename of single ion parameter file
               // constructor with positions scattering length dwf
   jjjpar(double x,double y,double z, double slr,double sli, double dwf);
   jjjpar (int n=1,int diag=0,int nofmom=3); // konstructor without file
   jjjpar (const jjjpar & jjjpars);	// kopier-konstruktor
   
  ~jjjpar ();		//destruktor
  
  private:

  // integer to tell which module is loaded 
  int intern_mcalc;

  // external module functions, intern_mcalc=0
  void (*m)(Vector*,double*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*);  
  int  (*dm)(int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*,float*,ComplexMatrix*);
  int  (*ddnn)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexMatrix*);

  void (*mq)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*);
  
  void (*estates)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**);
  
#ifdef __linux__
  void *handle;
#else
//  HANDLE handle;
  HINSTANCE__* handle;
#endif

  // kramers internal module functions, intern_mcalc=1
  Vector & kramer (double & T,Vector & H, double & Z,double & U);
  int  kramerdm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);

  // realisation of class iops - cfield internal module functions, intern_mcalc=2
  // the class iops calls for some functionality the program cfield (e.g. for
  // getting stevens factors and other parameters, for the matrices Olm etc.)
  ionpars * iops;

  // brillouin internal module functions, intern_mcalc=3
  Vector & brillouin (double & T,Vector & H, double & Z,double & U);
  int  brillouindm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);

    
  Vector ABC;   // storage for single ion module paramters
  void getpolar(double x,double y, double z, double & r, double & th, double & ph);// calculates polar coordinates from Vector X(1..3)

  void get_parameters_from_sipfile(char * cffilename); // function to read single ion parameter files



};

#endif
