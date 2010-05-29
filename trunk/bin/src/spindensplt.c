/***********************************************************
 *
 * spindensplt - program to display spin density of
 *           an ion given its CF pars, T and effective H
 * Author: M. Rotter
 **********************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100


#include "../../version"
#include "density.hpp"
#include "spincf.hpp"
#include "jjjpar.hpp"
#include "martin.h"
#include "plt_func.c"



/**********************************************************************/
// main program
int main (int argc, char **argv)
{
printf("***********************************************************\n");
printf("*\n");
printf("* spindensplt - program to display spin density of \n");
printf("*               an ion given its CF pars, T and effective H\n");
printf("* Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("* %s\n",MCPHASVERSION);
printf("***********************************************************\n");
// check command line
  if (argc < 10)
    { printf ("\nprogram spindensplt - calculate spindensity of a single ion\n\n\
                use as: spindensplt threshhold T Ha Hb Hc i j k mcphas.sipf \n\n\
                - given is temperature T[K] and magnetic effective field H[T]\n\
		- the spindensity vector component along direction (i,j,k) is\n\
                  calculated\n\
                - crystal field  parameters Blm should be read from a \n\
	          standard mcphas single ion property file mcphas.sipf \n\
                options: if T<0 then no thermal boltzmann distribution is taken\n\
		the statistical probability of each CF state has to be entered \n\
                by hand.\n\
		\n");
      printf ("\n\
...to view spindensityplots type: javaview results/spindensplt.jvx\n\
                                  displaycontour 2 3 4 results/spindenplti.grid\n\
                                  displaycontour 1 3 4 results/spindenpltj.grid\n\
                                  displaycontour 1 2 4 results/spindenpltk.grid\n\
                \n");
      exit (1);
    }

 graphic_parameters gp;
gp.show_abc_unitcell=0;
gp.show_primitive_crystal_unitcell=0;
gp.show_magnetic_unitcell=0;
gp.threshhold=strtod(argv[1],NULL);
gp.scale_pointcharges=1;
gp.show_pointcharges=1;

  double T,ha,hb,hc,xx,yy,zz;
  T=strtod(argv[2],NULL);
  ha=strtod(argv[3],NULL);
  hb=strtod(argv[4],NULL);
  hc=strtod(argv[5],NULL);
  xx=strtod(argv[6],NULL);
  yy=strtod(argv[7],NULL);
  zz=strtod(argv[8],NULL);
  double rr;
  // normalize direction vector
  rr=sqrt(xx*xx+yy*yy+zz*zz);
  xx/=rr;yy/=rr;zz/=rr;

 FILE * fout;

 // read cf-parameters into class object jjjpar
 jjjpar jjjps(0.0,0.0,0.0,argv[9]);
 if(jjjps.module_type!=0){fprintf(stderr,"ERROR chrgplt: calculation not possible for this single ion module\n");exit(EXIT_FAILURE);}
  
  int dim=49;
// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

  Vector h(1,dim);
  Vector moments(1,dim);
  Vector momentsx(1,dim);
  Vector momentsy(1,dim);
  Vector momentsz(1,dim);
  h=0;

  if(jjjps.module_type==0){ h(1)=2.0*MU_B*ha;h(3)=2.0*MU_B*hb;h(5)=2.0*MU_B*hc;h(2)=MU_B*ha;h(4)=MU_B*hb;h(6)=MU_B*hc;} 
  if(jjjps.module_type==2||jjjps.module_type==4){ h(1)=jjjps.gJ*MU_B*ha;h(2)=jjjps.gJ*MU_B*hb;h(3)=jjjps.gJ*MU_B*hc;}
  //int dj=(int)(2.0*(*iops).J+1);

  jjjps.mcalc_parameter_storage_init(h,T);
  printf("calculating expectation values ....\n");
double lnz,u;
  jjjps.mcalc(moments,T,h,lnz,u,jjjps.mcalc_parstorage);
if(xx!=0)jjjps.spindensity_mcalc (momentsx,1, T, h, jjjps.mcalc_parstorage);
if(yy!=0)jjjps.spindensity_mcalc (momentsy,2, T, h, jjjps.mcalc_parstorage);
if(zz!=0)jjjps.spindensity_mcalc (momentsz,3, T, h, jjjps.mcalc_parstorage);
 moments=xx*momentsx+yy*momentsy+zz*momentsz;

int i,nofpc=0;
printf("#spindensity is expanded in tesseral harmonics Zlm\n\
#   Ms(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
   for(i=1;i<=dim;++i){int l,m;double factor;
                       printf(" aS(%i,%i) =%12.6f\n",k[i],q[i],moments(i));}
printf("\n");

 char text[1000];
 if(jjjps.module_type==0||jjjps.module_type==4){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=abc, spindensity Ms(r).(%g,%g,%g)</title>\n", T, ha, hb, hc,xx,yy,zz);}
 if(jjjps.module_type==2){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=cab, spindensity Ms(r).(%g,%g,%g)</title>\n", T, ha, hb, hc,xx,yy,zz);}

 char * cffilenames[MAXNOFATOMS];
 cffilenames[1]=new char[MAXNOFCHARINLINE];
 float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS],gJ[MAXNOFATOMS];

Vector abc(1,3);abc(1)=6.0;abc(2)=6.0;abc(3)=6.0;
Matrix r(1,3,1,3);r=0;r(1,1)=1.0;r(2,2)=1.0;r(3,3)=1.0;

 strcpy(cffilenames[1],argv[9]);
 x[1]=0.5;y[1]=0.5;z[1]=0.5; // put atom in middle of cell


// read pointcharge-parameters
if(gp.show_pointcharges>0) nofpc=read_pointcharge_parameters(gp.scale_pointcharges,cffilenames,argv[9],x,y,z,jjjps,abc);

  spincf s(1,1,1,nofpc+1,dim);
  Vector gJJ(1,nofpc+1);gJJ=0;gJJ(1)=jjjps.gJ;
  Vector hkl(1,3);hkl=0;s=s*0;
  spincf ev_real(s),ev_imag(s);

  for(i=1;i<=dim;++i)s.m(1,1,1)(i)=moments(i);
  gp.show_atoms=1;gp.showprim=1;
  gp.spins_scale_static_moment=0;
  gp.spins_show_static_moment_direction=0;
  sprintf(gp.title,"spindensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  fout = fopen_errchk ("results/spindensplt.jvx", "w");
   s.jvx_cd(fout,text,abc,r,x,y,z,gJJ,gp,0.0,ev_real,ev_imag,hkl,cffilenames);
  fclose (fout);
  fout = fopen_errchk ("results/spindensplti.grid", "w");
  gp.showprim=1;
  s.cd(fout,abc,r,x,y,z,cffilenames,gp,1,200,200,ev_real,ev_imag,0.0,hkl);
  fclose (fout);
  fout = fopen_errchk ("results/spindenspltj.grid", "w");
  s.cd(fout,abc,r,x,y,z,cffilenames,gp,200,1,200,ev_real,ev_imag,0.0,hkl);
  fclose (fout);
  fout = fopen_errchk ("results/spindenspltk.grid", "w");
  s.cd(fout,abc,r,x,y,z,cffilenames,gp,200,200,1,ev_real,ev_imag,0.0,hkl);
  fclose (fout);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program spindensplt\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/spindensplt.jvx\n");
fprintf(stderr,"# * displaycontour 2 3 4 results/spindensplti.grid\n");
fprintf(stderr,"# * displaycontour 1 3 4 results/spindenspltj.grid\n");
fprintf(stderr,"# * displaycontour 1 2 4 results/spindenspltk.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  int dj;
  for(dj=1;dj<=nofpc+1;++dj){delete cffilenames[dj];}
  return 0;
}


