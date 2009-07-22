/***********************************************************
 *
 * chrgplt - program to display charge density of 
 *           an ion given its CF pars, T and effective H
 * Author: M. Rotter
 **********************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100
#define MU_B  5.788378E-02 // Bohrmagneton in meV/tesla

#include "../../version"
#include "chargedensity.hpp"
#include "spincf.hpp"
#include "jjjpar.hpp"
#include "martin.h"
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>


/**********************************************************************/
// main program
int main (int argc, char **argv)
{
printf("***********************************************************\n");
printf("*\n");
printf("* chrgplt - program to display charge density of \n");
printf("*           an ion given its CF pars, T and effective H\n");
printf("* Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("* %s\n",MCPHASVERSION);
printf("***********************************************************\n");
// check command line
  if (argc < 6)
    { printf ("\nprogram chrgplt - calculate chargdensity of a single ion\n\n\
                use as: chrgplt T Ha Hb Hc mcphas.cf \n\n\
                given is temperature T[K] and magnetic effective field H[T]\n \
		and crystal field  parameters Blm should be read from a \n\
	        standard mcphas single ion property file mcphas.cf \n\
                options: if T<0 then no thermal boltzmann distribution is taken\n\
		the statistical probability of each CF state has to be entered by hand.\n\
		\n");
      printf ("\n ... to view chargeplots type: javaview chrgplt.jvx\n\n");
      printf ("mind: 1)the 4f wavefunction R4f(r) is taken the same for all RE.\n");
      printf ("      2) the coordinate system is xyz||cab\n\n\n");
      exit (1);
    }
    else
    { printf ("#mind: 1) the 4f wavefunction R4f(r) is taken the same for all RE.\n");
      printf ("#      2) the coordinate system is xyz||cab, i.e. Ha=Hy=%sT  Hb=Hz=%sT  Hc=Hx=%sT \n",argv[2],argv[3],argv[4]);
    }

  double T,ha,hb,hc;
  T=strtod(argv[1],NULL);
  ha=strtod(argv[2],NULL);
  hb=strtod(argv[3],NULL);
  hc=strtod(argv[4],NULL);

  double dtheta=0.1; //stepwidth to step surface
  double dfi=0.1;
  double ccc=0.05; //surface value of chargedensity

 FILE * fout, * cf_file;

 // read cf-parameters into class object jjjpar
 jjjpar jjjps(0.0,0.0,0.0,argv[5]);
  int dim;
  double lnz,u;
  if(jjjps.module_type==0) 
  {dim=51;} // external module
  else  
  {if(jjjps.module_type==2) {dim=48;} // cfield
   else {fprintf(stderr,"ERROR chrgplt: calculation not possible for this single ion module\n");exit(EXIT_FAILURE);}
  }
  Vector h(1,dim);
  Vector moments(1,dim);
  h=0;

  if(jjjps.module_type==0){ h(1)=2.0*MU_B*ha;h(3)=2.0*MU_B*hb;h(5)=2.0*MU_B*hc;h(2)=MU_B*ha;h(4)=MU_B*hb;h(6)=MU_B*hc;} 
  if(jjjps.module_type==2){ h(1)=jjjps.gJ*MU_B*ha;h(2)=jjjps.gJ*MU_B*hb;h(3)=jjjps.gJ*MU_B*hc;}
  //int dj=(int)(2.0*(*iops).J+1);
 // ComplexMatrix ests(0,dj,1,dj);
  jjjps.eigenstates(h,T);
  moments=jjjps.mcalc(T,h,lnz,u,jjjps.est);
//  cfield  has to be used to calculate all the <Olm>.
//  printf("Stevens factors: alpha beta gamma = %4g %4g %4g \n",(*iops).alpha,(*iops).beta,(*iops).gamma);
//  printf("Lande Factor: gJ = %4g\n",(*iops).gJ);

 char text[1000];
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT in javaview xyz=cab</title>\n", T, ha, hb, hc);


 char * cffilenames[MAXNOFATOMS];
 cffilenames[1]=new char[MAXNOFCHARINLINE];
 float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS],gJ[MAXNOFATOMS];

 strcpy(cffilenames[1],argv[5]);
 x[1]=0;y[1]=0;z[1]=0;

// read pointcharge-parameters 
 cf_file = fopen_errchk (argv[5], "rb");
 float par[100];par[0]=99;
int pchere;int i,nofpc=0;
while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){;}
while(pchere>0)
{printf("pointcharge %g |e| at xyz=%g %g %g mind: xyz=cab\n",par[1],par[2],par[3],par[4]);
 ++nofpc;if(nofpc>MAXNOFATOMS){fprintf(stderr,"Error chrgplt - too many pointcharges");exit(1);}
  cffilenames[1+nofpc]=new char[MAXNOFCHARINLINE];
  sprintf(cffilenames[1+nofpc],"pointcharge radius=%g",par[1]);
  x[nofpc+1]=par[2]/5.0;
  y[nofpc+1]=par[3]/5.0;
  z[nofpc+1]=par[4]/5.0;
while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){}
}
fclose(cf_file);

  spincf s(1,1,1,nofpc+1,dim);
  
  Vector abc(1,3);abc(1)=5.0;abc(2)=5.0;abc(3)=5.0;
  
  Matrix r(1,3,1,3);r=0;r(1,1)=1.0;r(2,2)=1.0;r(3,3)=1.0;
  Vector gJJ(1,nofpc+1);gJJ=0;gJJ(1)=jjjps.gJ;
  Vector hkl(1,3);hkl=0;s=s*0;
  spincf ev_real(s),ev_imag(s);

  for(i=1;i<=48;++i)s.m(1,1,1)(i)=moments(i);
  double show_atoms=1;
  double show_spindensity=1;
  fout = fopen_errchk ("results/chrgplt.jvx", "w");
   s.jvx_cd(fout,text,abc,r,x,y,z,gJJ,0,0,0,show_atoms,1.0,1.0,1.0,
            1,0.0,ev_real,ev_imag,0.0,hkl,0.0,0.0,
            cffilenames,1.0,show_spindensity);
  fclose (fout);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program chrgplt\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/chrgplt.jvx\n");
fprintf(stderr,"# ************************************************************************\n");

  int dj;
  for(dj=1;dj<=nofpc+1;++dj){delete cffilenames[dj];}
  return 0;
}


