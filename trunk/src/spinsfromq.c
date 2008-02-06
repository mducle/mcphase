/***********************************************************************
 *
 * spins.c - program to create spinsconfiguration from q vector
 *
 ***********************************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100
#include "spincf.hpp"
#include "martin.h"
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>
#include<par.hpp>

/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ 

// check command line
  if (argc < 6)
    { printf (" program spinsfromq - create spinconfiguration from q vector\n \
                use as: spinsfromq  n1 n2 n3 h k l\n \
		n1 n2 n3 .... periodicity of supercell\n \
		h k l ....... components of qvector\n");
      exit (1);
    }

  par inputpars("./mcphas.j");
  
 
  int i,n1,n2,n3;
  double lnz,u;
  float d;
  double T;
  Vector h(1,inputpars.nofcomponents);
  Vector nettom(1,inputpars.nofcomponents);
  Vector qvector (1,3);
  Vector phi(1,inputpars.nofcomponents);
  phi=0;
  Vector momentq0(1,inputpars.nofcomponents);
  momentq0=0;
     
   T=1;
   h=0;for(i=1;i<=inputpars.nofcomponents;++i)h(i)=0.1;
  for(i=1;i<=inputpars.nofatoms;++i)
  {nettom=(*inputpars.jjj[i]).mcalc(T,h,lnz,u); }

  n1=strtol(argv[1],NULL,10);  
  n2=strtol(argv[2],NULL,10);  
  n3=strtol(argv[3],NULL,10);
  qvector(1)=strtod(argv[4],NULL);
  qvector(2)=strtod(argv[5],NULL);
  qvector(3)=strtod(argv[6],NULL);
   
  spincf savspin (n1,n2,n3,inputpars.nofcomponents,inputpars.nofatoms);
  spincf savspin1 (n1,n2,n3,inputpars.nofcomponents,inputpars.nofatoms);
  spincf savspin2 (n1,n2,n3,inputpars.nofcomponents,inputpars.nofatoms);
  // n1 n2 n3 .. periodicity of supercell
// nettom .... saturation moment (positive)
// qvector ... wave vector in units of reciprocal lattice
// momentq0 .. ferromagnetic component (between 0 and 1)
// phi ....... phase (for each component)
  savspin.spinfromq (n1,n2,n3,qvector,nettom, momentq0, phi);
 if (argc>7)
  {double s,t;
   printf("double q option\n q=%g %g %g:\n",qvector(1),qvector(2),qvector(3));
   s=nettom(1);t=nettom(2); nettom(2)=0;nettom(1)=s;
   savspin1.spinfromq (n1,n2,n3,qvector,nettom, momentq0, phi);
   savspin1.print(stdout);

   nettom=0;nettom(2)=t;
   qvector(1)=strtod(argv[7],NULL);
   qvector(2)=strtod(argv[8],NULL);
   qvector(3)=strtod(argv[9],NULL);
   printf(" q=%g %g %g:\n",qvector(1),qvector(2),qvector(3));
   savspin2.spinfromq (n1,n2,n3,qvector,nettom, momentq0, phi);
   savspin2.print(stdout);
   savspin=savspin1+savspin2;
   printf("double q structure:\n");   
  }
 
  savspin.print(stdout);
  return 0;
}


