 // *************************************************************************
 // ************************ jq *************************************
 // *************************************************************************
// methods for class jq
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <martin.h>
#include <vector.h>
#include <complex>
#include "jq.hpp"

#define MAXNOFSPINS  200

// returns md of spin [i=na,j=nb,k=nc] 
ComplexMatrix & jq::mat(int na, int nb, int nc,int ma,int mb,int mc)
{ return jj[iin(in(na,nb,nc),in(ma,mb,mc))];
}
// the same but for spin number "i,j"
ComplexMatrix & jq::mati(int i,int j)
{ return jj[iin(i,j)];
}
// get index ijk=iv(1-3)  of spinconfiguration number in 
int * jq::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of spin from indizes i,j,k
int jq::in(int i, int j, int k)
{return ((i*mxb+j)*mxc+k);}

int jq::iin(int i, int j)
{return (i*mx+j);}


// return number of spins
int jq::n()
{return (nofa*nofb*nofc);
}
int jq::na()
{return nofa;
}
int jq::nb()
{return nofb;
}
int jq::nc()
{return nofc;
}


/**************************************************************************/

//constructors
jq::jq (int n1,int n2,int n3,int n, int nc)
{  nofa=n1;nofb=n2;nofc=n3;
   nofatoms=n;nofcomponents=nc;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
   mx=mxa*mxb*mxc+1;
   
//dimension arrays
  jj = new ComplexMatrix[mx*mx+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (jj == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 

}

//kopier-konstruktor
jq::jq (const jq & p)
{ int i,j;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  mx=p.mx;
  
//dimension arrays
  jj = new ComplexMatrix[mx*mx+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (jj == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 

 for (i=1;i<=in(nofa,nofb,nofc);++i)
  {for (j=1;j<=in(nofa,nofb,nofc);++j)
     {jj[iin(i,j)]=p.jj[iin(i,j)];
    }
  }           

}


//destruktor
jq::~jq ()
{
  delete []jj;
}


