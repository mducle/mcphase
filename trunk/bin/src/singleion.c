/***********************************************************************
 *
 * singleion.c - program to display single ion momentum  at given htpoint
 *
 ***********************************************************************/

#define MAXNOFCHARINLINE 1000
#include<cstdio>
#include "martin.h"
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>
#include<par.hpp>

/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ int i,j,nt;
  double lnz,u;
  float d;
  double T;


  par inputpars("./mcphas.j");
// check command line
  if (argc < 3||argc-2>inputpars.nofcomponents)
    { printf (" program single ion  - display single ion properties at HT point\n \
                use as: single ion  T[K] gjmbHa gmbHb gjmbHc[meV] ...\n");
      exit (1);
    }

  Vector h(1,inputpars.nofcomponents);
  Vector m(1,inputpars.nofcomponents);
  // transition matrix Mij
  ComplexMatrix Mijkl(1,inputpars.nofcomponents,1,inputpars.nofcomponents);

      T=strtod(argv[1],NULL);
   h=0;for(i=1;i<=argc-2;++i)h(i)=strtod(argv[i+1],NULL);

printf("#\n#atom-number T[K] ");
for(j=1;j<=inputpars.nofcomponents;++j)printf("gjmbH%c(meV) ",'a'-1+j);
for(j=1;j<=inputpars.nofcomponents;++j)printf(" <J%c> ",'a'-1+j);
  printf("transition-energies(meV)...\n");

  for(i=1;i<=inputpars.nofatoms;++i)
{


            m=(*inputpars.jjj[i]).mcalc(T,h,lnz,u,(*inputpars.jjj[i]).eigenstates(h,T));
	    (*inputpars.jjj[i]).transitionnumber=-1;
            nt=(*inputpars.jjj[i]).dmcalc(T,h,Mijkl,d,(*inputpars.jjj[i]).est); 

  printf("%3i %13g ",i,T);
  for(j=1;j<=inputpars.nofcomponents;++j)printf(" %10g ",h(j));
//  printf("        ");
  for(j=1;j<=inputpars.nofcomponents;++j)printf(" %4g ",m(j));
  for(j=1;j<=nt;++j)
  {(*inputpars.jjj[i]).transitionnumber=-j;
  (*inputpars.jjj[i]).dmcalc(T,h,Mijkl,d,(*inputpars.jjj[i]).est);
  printf(" %4g ",d);
  }
  printf("\n");
}




}
