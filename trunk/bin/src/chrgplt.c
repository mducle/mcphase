/***********************************************************************
 *
 * chrgplt.c - new program substituting basic program to display 
 *             charge density of an ion given its CF pars, T and effective H
 *
 ***********************************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100
#define MU_B  5.788378E-02 // Bohrmagneton in meV/tesla

#include "chargedensity.hpp"
#include "spincf.hpp"
#include "ionpars.hpp"
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
      printf ("\n ... to view chargeplots type: java javaview chrgplt.jvx\n\n");
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

 FILE * cf_file, * fout;

// read cf-parameters into class object ionpars
 ionpars * iops;
 cf_file = fopen_errchk (argv[5], "rb");
 char instr[MAXNOFCHARINLINE];
 fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
  // strip /r (dos line feed) from line if necessary
  char *token;  
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}  
  // determine the file type
  if(strncmp(instr,"#!",2)!=0)
    {fprintf(stderr,"Error: single ion property file %s does not start with '#!'\n",argv[5]);
     exit(EXIT_FAILURE);}

  if(strncmp(instr,"#!cfield ",9)==0||strncmp(instr,"#!cfield\n",9)==0)
    {iops=new ionpars(cf_file);  
     // get 1ion parameters - operator matrices
    }else {fprintf(stderr,"Error program chrgplt: single ion property file %s does not start with '#!cfield'\nother modules not supported yet\n",argv[5]);
     exit(EXIT_FAILURE);}
     fclose(cf_file);


  double lnz,u;
  Vector h(1,48);
  Vector moments(1,48);
  h=0; h(1)=(*iops).gJ*MU_B*ha;h(2)=(*iops).gJ*MU_B*hb;h(3)=(*iops).gJ*MU_B*hc;
  int dj=(int)(2.0*(*iops).J+1);
  ComplexMatrix ests(0,dj,1,dj);
  ests=(*iops).cfeigenstates(h,T);
//  cfield  has to be used to calculate all the <Olm>.
  moments=(*iops).cfield(T,h,lnz,u,ests);
  printf("Stevens factors: alpha beta gamma = %4g %4g %4g \n",(*iops).alpha,(*iops).beta,(*iops).gamma);
  printf("Lande Factor: gJ = %4g\n",(*iops).gJ);

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

  spincf s(1,1,1,nofpc+1,48);
  
  Vector abc(1,3);abc(1)=5.0;abc(2)=5.0;abc(3)=5.0;
  
  Matrix r(1,3,1,3);r(1,1)=1.0;r(2,2)=1.0;r(3,3)=1.0;
  Vector gJJ(1,nofpc+1);gJJ=0;gJJ(1)=(*iops).gJ;
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
 
  for(dj=1;dj<=nofpc+1;++dj){delete cffilenames[dj];}
  return 0;
}


