/***********************************************************
 *
 * chrgplt - program to display charge density of 
 *           an ion given its CF pars, T and effective H
 * Author: M. Rotter
 **********************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100
#define MU_B  5.788378E-02 // Bohrmagneton in meV/tesla
#define PI 3.1415926535

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
  if (argc < 7)
    { printf ("\nprogram chrgplt - calculate chargdensity of a single ion\n\n\
                use as: chrgplt threshhold T Ha Hb Hc  mcphas.sipf \n\n\
                given is temperature T[K] and magnetic effective field H[T]\n \
		and crystal field  parameters Blm should be read from a \n\
	        standard mcphas single ion property file mcphas.cf \n\
                options: if T<0 then no thermal boltzmann distribution is taken\n\
		the statistical probability of each CF state has to be entered by hand.\n\
		\n");
      printf ("\n ... to view chargeplots type: javaview chrgplt.jvx\n\n");
      exit (1);
    }

  double T,ha,hb,hc, threshhold;
  threshhold=strtod(argv[1],NULL);
  T=strtod(argv[2],NULL);
  ha=strtod(argv[3],NULL);
  hb=strtod(argv[4],NULL);
  hc=strtod(argv[5],NULL);

  double dtheta=0.1; //stepwidth to step surface
  double dfi=0.1;
  double ccc=0.05; //surface value of chargedensity

 FILE * fout, * cf_file;

 // read cf-parameters into class object jjjpar
 jjjpar jjjps(0.0,0.0,0.0,argv[6]);
  int dim;
  double lnz,u;
  if(jjjps.module_type==0) 
  {dim=51;} // external module
  else  
  {if(jjjps.module_type==2||jjjps.module_type==4) {dim=48;} // cfield
   else {fprintf(stderr,"ERROR chrgplt: calculation not possible for this single ion module\n");exit(EXIT_FAILURE);}
  }
  Vector h(1,dim);
  Vector moments(1,dim);
  h=0;

  if(jjjps.module_type==0){ h(1)=2.0*MU_B*ha;h(3)=2.0*MU_B*hb;h(5)=2.0*MU_B*hc;h(2)=MU_B*ha;h(4)=MU_B*hb;h(6)=MU_B*hc;} 
  if(jjjps.module_type==2||jjjps.module_type==4){ h(1)=jjjps.gJ*MU_B*ha;h(2)=jjjps.gJ*MU_B*hb;h(3)=jjjps.gJ*MU_B*hc;}
  //int dj=(int)(2.0*(*iops).J+1);

  jjjps.mcalc_parameter_storage_init(h,T);
  printf("calculating expectation values ....\n");
  jjjps.mcalc(moments,T,h,lnz,u,jjjps.mcalc_parstorage);
//  cfield  has to be used to calculate all the <Olm>.
//  printf("Stevens factors: alpha beta gamma = %4g %4g %4g \n",(*iops).alpha,(*iops).beta,(*iops).gamma);
//  printf("Lande Factor: gJ = %4g\n",(*iops).gJ);

int pchere;int i,nofpc=0;
const char lm[]=             "Jy  Jz  Jx  O22SO21SO20 O21 O22 O33SO32SO31SO30 O31 O32 O33 O44SO43SO42SO41SO40 O41 O42 O43 O44 O55SO54SO53SO52SO51SO50 O51 O52 O53 O54 O55 O66SO65SO64SO63SO62SO61SO60 O61 O62 O63 O64 O65 O66 ";
const char lme[]="Sx  Lx  Sy  Ly  Sz  Lz  O22SO21SO20 O21 O22 O33SO32SO31SO30 O31 O32 O33 O44SO43SO42SO41SO40 O41 O42 O43 O44 O55SO54SO53SO52SO51SO50 O51 O52 O53 O54 O55 O66SO65SO64SO63SO62SO61SO60 O61 O62 O63 O64 O65 O66 ";
char lm4[5];lm4[4]='\0';
printf("#chargedensity is expanded in tesseral harmonics as\n#   ro(r)= -|e||R4f(r)|^2 sum_lm clm Zlm(Omega)\n#\n");
   for(i=1;i<=dim;++i){int l,m;double factor;
                       if (dim==48) strncpy(lm4,lm+(i-1)*4,4);
                       if (dim==51) strncpy(lm4,lme+(i-1)*4,4);
                       factor=0;
                       if(jjjps.module_type==2||jjjps.module_type==4){if(i>3){l=lm4[1]-48;m=lm4[2]-48;
                                                        factor=jjjps.tetan()(l)*jjjps.cnst(l,m);
                                                        }
                                               } // these are prefactors in case of module cfield (stevens parameters tetan and zlm prefactors)
                                      else     {if(i>6){l=lm4[1]-48;m=lm4[2]-48;
                                                       if(m!=0){factor=sqrt((2.0*l+1)/8/PI);}
                                                           else{factor=sqrt((2.0*l+1)/4/PI);}
                                                       }
                                               } 
                       printf(" <J%c> = <%s> =%12.6f   clm=%12.6f\n",'a'-1+i,lm4,moments(i),moments(i)*factor);}
printf("\n");

 char text[1000];
 if(jjjps.module_type==0||jjjps.module_type==4){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=abc</title>\n", T, ha, hb, hc);}
 if(jjjps.module_type==2){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=cab</title>\n", T, ha, hb, hc);}

 char * cffilenames[MAXNOFATOMS];
 cffilenames[1]=new char[MAXNOFCHARINLINE];
 float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS],gJ[MAXNOFATOMS];

Vector abc(1,3);abc(1)=6.0;abc(2)=6.0;abc(3)=6.0;
Matrix r(1,3,1,3);r=0;r(1,1)=1.0;r(2,2)=1.0;r(3,3)=1.0;

 strcpy(cffilenames[1],argv[6]);
 x[1]=0.5;y[1]=0.5;z[1]=0.5; // put atom in middle of cell
//printf("hello\n");

if(jjjps.calcmagdensity==0)
{
// read pointcharge-parameters 
 cf_file = fopen_errchk (argv[6], "rb");
 float par[100];par[0]=99;
while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){;}
while(pchere>3)
{if(jjjps.module_type==0||jjjps.module_type==4){printf("pointcharge %g |e| at xyz=%g %g %g mind: xyz=abc\n",par[1],par[2],par[3],par[4]);}
 if(jjjps.module_type==2){printf("pointcharge %g |e| at xyz=%g %g %g mind: xyz=cab\n",par[1],par[2],par[3],par[4]);}
 ++nofpc;if(nofpc>MAXNOFATOMS){fprintf(stderr,"Error chrgplt - too many pointcharges");exit(1);}
  cffilenames[1+nofpc]=new char[MAXNOFCHARINLINE];
  sprintf(cffilenames[1+nofpc],"pointcharge radius=%g",0.529177*copysign(1.0,par[1])*pow(fabs(par[1]),0.3333));
  if(jjjps.module_type==0||jjjps.module_type==4){
  x[nofpc+1]=par[2]/abc(1)+x[1];// these are the positions in Angstroem
  y[nofpc+1]=par[3]/abc(2)+y[1];// however in order to be in line with the cfield xyz=cab
  z[nofpc+1]=par[4]/abc(3)+z[1];// and ic1ion xyz=abc we have to set these parameters
                          }
  if(jjjps.module_type==2){
  x[nofpc+1]=par[3]/abc(1)+x[1];// these are the positions in Angstroem (we set a=b=c=1A below)
  y[nofpc+1]=par[4]/abc(2)+y[1];// however in order to be in line with the cfield xyz=cab
  z[nofpc+1]=par[2]/abc(3)+z[1];// and ic1ion xyz=abc we have to set these parameters
                          }
while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){}
}
fclose(cf_file);
}
  spincf s(1,1,1,nofpc+1,dim);
  
  
  Vector gJJ(1,nofpc+1);gJJ=0;gJJ(1)=jjjps.gJ;
  Vector hkl(1,3);hkl=0;s=s*0;
  spincf ev_real(s),ev_imag(s);

  for(i=1;i<=dim;++i)s.m(1,1,1)(i)=moments(i);
  double show_atoms=1;
  double spins_scale_static_moment=1;
  fout = fopen_errchk ("results/chrgplt.jvx", "w");
   s.jvx_cd(fout,text,abc,r,x,y,z,gJJ,0,0,0,show_atoms,1.0,1.0,1.0,
            1,0.0,ev_real,ev_imag,0.0,hkl,0.0,spins_scale_static_moment,
            cffilenames,1.0,0.0,threshhold);
  fclose (fout);
  fout = fopen_errchk ("results/chrgplti.grid", "w");
  s.cd(fout,abc,r,x,y,z,cffilenames,1,1,200,200,1.0,1.0,1.0,ev_real,ev_imag,0.0,0.0,hkl);
  fclose (fout);
  fout = fopen_errchk ("results/chrgpltj.grid", "w");
  s.cd(fout,abc,r,x,y,z,cffilenames,1,200,1,200,1.0,1.0,1.0,ev_real,ev_imag,0.0,0.0,hkl);
  fclose (fout);
  fout = fopen_errchk ("results/chrgpltk.grid", "w");
  s.cd(fout,abc,r,x,y,z,cffilenames,1,200,200,1,1.0,1.0,1.0,ev_real,ev_imag,0.0,0.0,hkl);
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


