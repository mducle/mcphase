#ifndef IONPARS
#define IONPARS

#include <vector.h>
#include <cstdio>
#include <mpspecfunp.h>

// ionpars: class to 
//          - read single ion parameterfile for cfield module
//          - load and store matrices for internal module cfield
//          - diagonalize cfproblem and calculate moment and transition matrix elements

class ionpars  
{private: 
   double rk_from_radial_wavefunction(int k); // needed for public radial wave function <r^n> calculation
   // calculates scattering operator 
   void MQM(ComplexMatrix & MQXM,ComplexMatrix & MQYM,ComplexMatrix & MQZM, double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc);

 public:
   char * iontype; // description string
   double J;// momentum quantum number
   double gJ; // Lande factor
   double alpha;
   double beta;  // stevens factors
   double gamma;
   double r2;
   double r4;  // radial wave function exp values
   double r6;
  
   Matrix Ja; Matrix Jb; Matrix Jc; Matrix Hcf;
   ComplexMatrix Jaa;
   ComplexMatrix Jbb;
   ComplexMatrix Jcc;
 

   Matrix **Olm; ComplexMatrix **OOlm; // array of matrices
 
   Vector Blm; // Cf parameters  
   Vector Llm; // Cf parameters  

   Vector Np,Xip,Cp; // radial wave function parameters

   // evaluate radial wave function // r given in Angstroems, returns R(r) in units of 1/A^1.5
   double radial_wavefunction(double r);
   void save_radial_wavefunction(char * filename);

   //functions to calculate radial matrix elements <r^n> from radial wave function
   int r2_from_radial_wavefunction();
   int r4_from_radial_wavefunction();
   int r6_from_radial_wavefunction();
 
   // functions needed to calculate thermal expectation value of moment  
   Vector & cfield (double & T,Vector & H, double & Z,double & U);
   ComplexMatrix & cfeigenstates (Vector & H, double & T);
   // and transition matrix elements
   int  cfielddm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);
   int cfielddn(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexMatrix & nat);
   // calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
   // according to stored eigenstate matrix est
   // calculates the scattering operator given the polar angles th, ph (with respect to the CEF coordinate 
   // system xyz and the <jl(qr)> and the eigenstate matrix with eigenstates and thermal population numbers
   ComplexVector & MQ(double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc,ComplexMatrix & est);


   void savBlm(FILE * file); // saving Blm to file 
   void savLlm(FILE * file); // saving Blm to file 
 
   void save_radial_wavefunction(const char * filename);

   ionpars(int dimj);
   ionpars(FILE * cf_file);
   ionpars (char * iontype); // constructor from iontype (mind:no matrices filled with values !)
   ~ionpars();
   ionpars(const ionpars & p);
};
#endif

