#ifndef IONPARS
#define IONPARS

#include <vector.h>
#include <cstdio>

// ionpars: class to load and store matrices for internal module cfield

class ionpars  
{private: 
 public:
   double J;// momentum quantum number
   double alpha;
   double beta;  // stevens factors
   double gamma;
   
   Matrix Ja; Matrix Jb; Matrix Jc; Matrix Hcf;
   ComplexMatrix Jaa;
   ComplexMatrix Jbb;
   ComplexMatrix Jcc;
 

   Matrix **Olm; ComplexMatrix **OOlm; // array of matrices

   ionpars(int dimj);
   ionpars(FILE * cf_file);
   ~ionpars();
   ionpars(const ionpars & p);
};
#endif

