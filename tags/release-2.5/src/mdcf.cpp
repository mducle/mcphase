 // *************************************************************************
 // ************************ spincf *************************************
 // *************************************************************************
// methods for class spincf 
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <martin.h>
#include <vector.h>
#include <complex>
#include "mdcf.hpp"

#define MAXNOFSPINS  200

// returns md of spin [i=na,j=nb,k=nc] 
ComplexMatrix & mdcf::U(int na, int nb, int nc)
{ return s[in(na,nb,nc)];
}
ComplexMatrix & mdcf::M(int na, int nb, int nc)
{ return m[in(na,nb,nc)];
}
ComplexMatrix & mdcf::lambda(int na, int nb, int nc)
{ return l[in(na,nb,nc)];
}
Vector & mdcf::delta(int na, int nb, int nc)
{ return d[in(na,nb,nc)];
}
// the same but for spin number "i"
ComplexMatrix & mdcf::Ui(int i)
{ return s[i];
}
ComplexMatrix & mdcf::Mi(int i)
{ return m[i];
}
ComplexMatrix & mdcf::lambdai(int i)
{ return l[i];
}
Vector  & mdcf::deltai(int i)
{ return d[i];
}
// get index ijk=iv(1-3)  of spinconfiguration number in 
int * mdcf::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of spin from indizes i,j,k
int mdcf::in(int i, int j, int k)
{return ((i*mxb+j)*mxc+k);}


// return number of spins
int mdcf::n()
{return (nofa*nofb*nofc);
}
int mdcf::na()
{return nofa;
}
int mdcf::nb()
{return nofb;
}
int mdcf::nc()
{return nofc;
}

/**************************************************************************/

//zuweisung
mdcf & mdcf::operator= (const mdcf & op2)
{int i,j,k;
 nofa=op2.nofa; nofb=op2.nofb; nofc=op2.nofc;
 mxa=op2.mxa; mxb=op2.mxb; mxc=op2.mxc;
 nofatoms=op2.nofatoms;
 nofcomponents=op2.nofcomponents;
  
  delete[]s;delete []m;delete []d;delete []l;
//dimension arrays
  s = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  m = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (m == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  l = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (l == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  d = new Vector[mxa*mxb*mxc+1](1,nofatoms);
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {s[in(i,j,k)]=op2.s[in(i,j,k)];
      m[in(i,j,k)]=op2.m[in(i,j,k)];
      l[in(i,j,k)]=op2.l[in(i,j,k)];
      d[in(i,j,k)]=op2.d[in(i,j,k)];
     } 
    }
  }           
  return *this;
}



//constructors
mdcf::mdcf (int n1,int n2,int n3,int n,int nc)
{  nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
   nofatoms=n;nofcomponents=nc;
//dimension arrays
  s = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  m = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (m == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  l = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (l == NULL){ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  d = new Vector[mxa*mxb*mxc+1](1,nofatoms);
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

}

//kopier-konstruktor
mdcf::mdcf (const mdcf & p)
{ int i,j,k;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  nofatoms=p.nofatoms;nofcomponents=p.nofcomponents;
//dimension arrays
  s = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (s == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  m = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (m == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  l = new ComplexMatrix[mxa*mxb*mxc+1](1,nofcomponents*nofatoms,1,nofcomponents*nofatoms);
  if (l == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 
  d = new Vector[mxa*mxb*mxc+1](1,nofatoms);
  if (d == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

 for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {s[in(i,j,k)]=p.s[in(i,j,k)];
      m[in(i,j,k)]=p.m[in(i,j,k)];
      l[in(i,j,k)]=p.l[in(i,j,k)];
      d[in(i,j,k)]=p.d[in(i,j,k)];
     } 
    }
  }           

}


//destruktor
mdcf::~mdcf ()
{
 delete []s;delete []m;delete []d;delete []l;
}


