// klasse testmdcf zum speichern der single ion - crystal field uebergangsmatrixelemente

#include<cstdlib>
#include<cerrno>
#include<vector.h>
#include<complex>

#ifndef MDCF_H
#define MDCF_H

class mdcf 
{
  private:
 // number of spins  
   int nofa,nofb,nofc;
   int mxa,mxb,mxc;
   ComplexMatrix ** m; //matrix to store M
   ComplexMatrix ** s; //matrix to store S
   ComplexMatrix ** l; //matrix to store eigenvalues sqrt_gamma
   Vector ** d;
   IntVector * nt; // vector to store number of transitions for each atom
   int iv[4];
       
  public:
 // array of spins 
    int nofatoms,nofcomponents;
    int in(int i,int j, int k); 
   
    ComplexMatrix & M(int i,int j,int k); // returns pointer to  matrix M(ijk) 
    ComplexMatrix & Mi(int in); // returns pointer to matrix M(i)
    ComplexMatrix & U(int i,int j,int k); // returns pointer to eigenvector matrix (ijk) 
    ComplexMatrix & Ui(int in); // returns pointer to eigenvector matrix i
    ComplexMatrix & sqrt_gamma(int i,int j,int k); // returns pointer to eigenvaluematrix (ijk) 
    ComplexMatrix & sqrt_gammai(int in); // returns pointer to eigenvalue matrix i

    Vector & delta(int i,int j,int k); // returns pointer to matrix (ijk) 
    Vector & deltai(int in); // returns pointer to mean field i
    int * ijk(int in);  // returns mf indizes (ijk)(in): in=0,...,n(=na*nb*nc)
    
    int n(); // returns total number of spins
    int na(); // returns number of spins
    int nb(); // returns number of spins
    int nc(); // returns number of spins
    int baseindex (int i, int j, int k, int l, int tn); // returns base index for atom l and transition tn
                                                       // which is needed for  setting up matrix U,M, sqrt_gamma and Vector D
    int baseindex_max(int i, int j, int k);
    int noft(int i, int j, int k, int l); // returns number of transitions of ion l in cryst unit ijk
//    mdcf & operator= (const mdcf & op2); // zuweisung

   mdcf (int n1,int n2,int n3,int n,int nc);	//konstruktor

// initialisierung 
    void set_noftransitions (int i, int j, int k, IntVector & notr);

//   mdcf (const mdcf & spins);	// kopier-konstruktor
   ~mdcf ();		//destruktor

};

#endif
