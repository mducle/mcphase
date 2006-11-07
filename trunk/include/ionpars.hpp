#ifndef IONPARS
#define IONPARS

#include <vector.h>
#include <cstdio>

// ionpars: class to 
//          - read single ion parameterfile for cfield module
//          - load and store matrices for internal module cfield
//          - diagonalize cfproblem and calculate moment and transition matrix elements

class ionpars  
{private: 
 public:
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

   // functions needed to calculate thermal expectation value of moment  
   Vector & cfield (double & T,Vector & H, double & Z,double & U);
   // and transition matrix elements
   int  cfielddm (int & tn,double & T,Vector &  heff, ComplexMatrix & mat,float & delta);
   void savBlm(FILE * file); // saving Blm to file 

   ionpars(int dimj);
   ionpars(FILE * cf_file);
   ionpars (char * iontype); // constructor from iontype (mind:no matrices filled with values !)
   ~ionpars();
   ionpars(const ionpars & p);
};
#endif

