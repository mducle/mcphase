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
#include "ionpars.hpp"
#include "martin.h"
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>


// sub for calculation of charge density given a radiu R and polar angles teta, fi and expansion coeff. alm

double rocalc (double & teta,double & fi,double & R, Matrix & a)
{double ro,ct,ct2,st,st2,sfi,cfi,rs,rr;
if (R>3.0||R<0){ro = 1e+10;}else{
ct = cos(teta);                      //z
ct2 = ct * ct;
st = sin(teta);
st2 = st * st;
sfi = sin(fi);
cfi = cos(fi);

//  we have to find the 4f wavefunction R4f(r) for each single ion and the Zlm, cfield has nothing: so we have
//     to take this from chrgplt.bas - a little problem: how do we get the correct R4f(r) ? for a first attempt
//     we could just take the same for all RE.

//k = 11 / 10 * 11 / 9 * 11 / 8 * 11 / 7 * 11 / 6 * 11 / 5 * 11 / 4 * 11 / 3 * 11 / 2 * 11 / 1 * 11
rs = R * exp(-R);
rr = 78624.0 * rs * rs * rs * rs;
rr = rr * rs * rs * rs * rs * exp(-3.0 * R);

ro = a(0, 0) / sqrt(4.0 * 3.1415);
ro = ro + a(2, -2)  * 2 * st2 * sfi * cfi;
ro = ro + a(2, -1)  * st * sfi * ct;
ro = ro + a(2, 0)  * (3 * ct2 - 1);
ro = ro + a(2, 1)  * st * cfi * ct;
ro = ro + a(2, 2)  * st2 * (cfi * cfi - sfi * sfi);

ro = ro + a(4, -4) * st2 * st2 * 4 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
ro = ro + a(4, -3) * ct * st * st2 * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
ro = ro + a(4, -2) * (7 * ct2 - 1) * 2 * st2 * cfi * sfi;
ro = ro + a(4, -1) * st * sfi * ct * (7 * ct2 - 3);
ro = ro + a(4, 0) * (35 * ct2 * ct2 - 30 * ct2 + 3);
ro = ro + a(4, 1)  * st * cfi * ct * (7 * ct2 - 3);
ro = ro + a(4, 2)  * (7 * ct2 - 1) * st2 * (cfi * cfi - sfi * sfi);
ro = ro + a(4, 3)  * ct * st * st2 * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
ro = ro + a(4, 4)  * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);

ro = ro + a(6, -6) * st2 * st2 * st2 * (6 * cfi * cfi * cfi * cfi * cfi * sfi - 20 * cfi * cfi * cfi * sfi * sfi * sfi + 6 * cfi * sfi * sfi * sfi * sfi * sfi);
ro = ro + a(6, -5) * ct * st * st2 * st2 * (5 * cfi * cfi * cfi * cfi * sfi - 10 * cfi * cfi * sfi * sfi * sfi + sfi * sfi * sfi * sfi * sfi);
ro = ro + a(6, -4) * (11 * ct2 - 1) * 4 * st2 * st2 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
ro = ro + a(6, -3) * (11 * ct * ct2 - 3 * ct) * st2 * st * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
ro = ro + a(6, -2) * 2 * st2 * sfi * cfi * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2);
ro = ro + a(6, -1) * ct * st * sfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
ro = ro + a(6, 0)  * (231 * ct2 * ct2 * ct2 - 315 * ct2 * ct2 + 105 * ct2 - 5);
ro = ro + a(6, 1)  * ct * st * cfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
ro = ro + a(6, 2)  * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2) * st2 * (cfi * cfi - sfi * sfi);
ro = ro + a(6, 3)  * (11 * ct * ct2 - 3 * ct) * st2 * st * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
ro = ro + a(6, 4)  * (11 * ct2 - 1) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);
ro = ro + a(6, 5)  * ct * st * st2 * st2 * (cfi * cfi * cfi * cfi * cfi - 10 * cfi * cfi * cfi * sfi * sfi + 5 * cfi * sfi * sfi * sfi * sfi);
ro = ro + a(6, 6) * st2 * st2 * st2 * (cfi * cfi * cfi * cfi * cfi * cfi - 15 * cfi * cfi * cfi * cfi * sfi * sfi + 15 * cfi * cfi * sfi * sfi * sfi * sfi - sfi * sfi * sfi * sfi * sfi * sfi);
ro = ro * rr;
}
return ro;
}

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
 fout = fopen_errchk ("./chrgplt.jvx", "w");

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

  Vector tetan(1,6);
  double gJ;
// returns stevens parameters and landefactor of ion
  tetan(2)=(*iops).alpha;tetan(4)=(*iops).beta;tetan(6)=(*iops).gamma;
  gJ=(*iops).gJ;
  
  chargedensity cd(dtheta,dfi);


  double lnz,u;
  Vector h(1,48);
  Vector moments(1,48);
  h=0; h(1)=gJ*MU_B*ha;h(2)=gJ*MU_B*hb;h(3)=gJ*MU_B*hc;
  int dj=(int)(2.0*(*iops).J+1);
  ComplexMatrix ests(0,dj,1,dj);
  ests=(*iops).cfeigenstates(h,T);
//  cfield  has to be used to calculate all the <Olm>.
  moments=(*iops).cfield(T,h,lnz,u,ests);
  printf("Stevens factors: alpha beta gamma = %4g %4g %4g \n",tetan(2),tetan(4),tetan(6));
  printf("Lande Factor: gJ = %4g\n",gJ);


 int l,m;   
// cnst is the Zlm constants - put them into the matrix
Matrix cnst(0,6,-6,6);

cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;
cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;
cnst(6,0)=  0.06357014;
cnst(6,1)=  1.032669;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}


Matrix a(0,6,-6,6);

a(0, 0) = 1 / sqrt(4 * 3.1415);

a(2,-2)=moments(4);
a(2,-1)=moments(5);
a(2,0)=moments(6);
a(2,1)=moments(7);
a(2,2)=moments(8);

a(4,-4)=moments(16);
a(4,-3)=moments(17);
a(4,-2)=moments(18);
a(4,-1)=moments(19);
a(4, 0)=moments(20);
a(4, 1)=moments(21);
a(4, 2)=moments(22);
a(4, 3)=moments(23);
a(4, 4)=moments(24);

a(6,-6)=moments(36);
a(6,-5)=moments(37);
a(6,-4)=moments(38);
a(6,-3)=moments(39);
a(6,-2)=moments(40);
a(6,-1)=moments(41);
a(6,-0)=moments(42);
a(6, 1)=moments(43);
a(6, 2)=moments(44);
a(6, 3)=moments(45);
a(6, 4)=moments(46);
a(6, 5)=moments(47);
a(6, 6)=moments(48);

for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m)a(l,m)*=tetan(l)*cnst(l,m)*cnst(l,m);}


/*FOR l = 2 TO 6 STEP 2: FOR m = -l TO l
PRINT USING "coefficient## ## ##.###^^^^ "; l, m, a(l, m)/cnst(l,m);
PRINT #9, USING "coefficient## ## ##.###^^^^ "; l, m, a(l, m)/cnst(l,m);
a(l, m) = a(l, m) * 1
NEXT m: PRINT : PRINT #9, : NEXT l*/

// here the set of points for this ion should be created corresponding to its
// charge density: 


//REM dieses programm berechnet aus alm die ladungsdichteverteilung
//REM und zeichnet diese auf

Vector rrttff(1,3,1,3);
int iii,iv,nt;
int anzahl = 0;
int imax = 3;
Vector rp(1,imax);
double rmax = 0;
double rstp = .1;
double max = .01 * ccc;  //end of intervalschachtelung to find r(ro=ccc)

double theta,fi,R,rin,ro,deltaa,rstpp,delta1,dd;
int tt,ff;
for(tt=0;tt<=3.1415/dtheta;++tt){for(ff=0;ff<=2*3.1415/dfi;++ff){
   rp=0;theta=(double)tt*dtheta;fi=(double)ff*dfi;
   nt=0;
   for(rin=0.1;rin<=2.2;rin+=0.2){
   R=rin;ro=rocalc(theta,fi,R,a);
   deltaa=fabs(ro-ccc);rstpp=rstp;delta1=1e4;
   for(iii=1;iii<=100&delta1>=max;++iii)
     {R+=rstpp;ro=rocalc(theta,fi,R,a);
      delta1=fabs(ro-ccc);
      if(delta1>=max){
                      if(delta1<deltaa){deltaa=delta1;}
		      else{R-=rstpp;rstpp*=-0.5;}
                     }
     }
     if(delta1<max)
     {for(iv=1;iv<=nt&fabs(rp(iv)-R)>=0.1;++iv);
      if(fabs(rp(iv)-R)>=0.1)
       {++nt;rp(nt)=R;
        if (R>rmax)rmax=R;
       }
     }
   }//next rin
   if(nt==0) {rp(1)=0.05;rp(2)=0.06;}
   if(nt==1) {rp(2)=rp(1)+0.0001;}
   // select most outsinde value of rp
   R=0;
   for(iv=1;iv<=imax;++iv){if(rp(iv)>R)R=rp(iv);}
	++anzahl;
	rrttff(1)=R;rrttff(2)=theta;rrttff(3)=fi;
   cd.rtf(anzahl)=rrttff;
   if(tt==0){ff=(int)(2*3.1415/dfi+1);}
}}


  
// 3.: copy from chrgplot.bas the calculation module of the charge density given the <Olm> - creating both a
//     postscript output and a jvx file. Call this for every spin subsequently. Maybe for the postscript
//     files look a bit on the .eps3d function in how to create such a plot


  

 fprintf(fout,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n");
 fprintf(fout,"<!DOCTYPE jvx-model SYSTEM \"http://www.javaview.de/rsrc/jvx.dtd\">\n");
 fprintf(fout,"<jvx-model>	<meta generator=\"JavaView v.2.00.008\"/>\n");
 fprintf(fout,"<meta date=\"Wed Mar 07 22:30:58 GMT+01:00 2001\"/>\n");
 fprintf(fout,"<version type=\"dump\">0.02</version>\n");

fprintf(fout,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT in javaview xyz=cab</title>\n", T, ha, hb, hc);
fprintf(fout,"<geometries>\n");


fprintf(fout,"<geometry name=\"");
fprintf(fout,"T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT in javaview xyz=cab\">\n",T,ha,hb,hc);
fprintf(fout,"<pointSet color=\"hide\" point=\"show\" dim=\"1\">\n");
fprintf(fout,"<points num=\"%i\">\n",cd.nofpoints());
int ii;
 for(ii=1;ii<=cd.nofpoints();++ii)
    {double dx,dy,dz,R,fi,theta;
     R=cd.rtf(ii)(1);
     theta=cd.rtf(ii)(2);
     fi=cd.rtf(ii)(3);
     // mind abc||yzx
     dx=R*sin(theta)*cos(fi);
     dy=R*sin(theta)*sin(fi);
     dz=R*cos(theta);
     fprintf(fout,"<p>%4g %4g %4g</p>\n",dx,dy,dz);
     }

    fprintf(fout,"<thickness>0.0</thickness><color type=\"rgb\">255 0 0</color><colorTag type=\"rgb\">255 0 255</colorTag>\n");
    fprintf(fout,"</points>			</pointSet>\n");
    fprintf(fout,"<faceSet face=\"show\" edge=\"show\">\n");
    fprintf(fout,"<faces num=\"%i\">\n",cd.nofpoints());
    int i,offset=0;

int ntt,nff,pointnr,ffnr,p1,p2,p3,p4;
ntt=(int)(3.1415/dtheta);
nff=(int)(2*3.1415/dfi);
pointnr=ntt*(nff+1);
ffnr=nff+1;
for(tt=1;tt<=ntt;++tt){for(ff=0;ff<=nff;++ff){
   p1 = ff + 1 + (tt - 2) * ffnr+offset;
   p2 = ff + 2 + (tt - 2) * ffnr+offset;
   p3 = ff + 2 + (tt - 1) * ffnr+offset;
   p4 = ff + 1 + (tt - 1) * ffnr+offset;
   if (ff==nff){p3 = p3 - ffnr; p2 = p2 - ffnr;}
   if (tt==1) {p1 = offset; p2 = offset;}
   fprintf(fout,"<f> %i %i %i %i </f>\n",p1,p2,p3,p4);
  }}
offset+=pointnr+1;
     
fprintf(fout,"<color type=\"rgb\">100 230 255</color>\n");
fprintf(fout,"<colorTag type=\"rgb\">255 0 255</colorTag>\n");
fprintf(fout,"</faces></faceSet></geometry>\n");

dtheta*=2;dfi*=2;

  chargedensity cdp(dtheta,dfi);
anzahl=0;

for(tt=0;tt<=3.1415/dtheta;++tt){for(ff=0;ff<=2*3.1415/dfi;++ff){
   theta=(double)tt*dtheta;fi=(double)ff*dfi;
	++anzahl;
	rrttff(1)=0;rrttff(2)=theta;rrttff(3)=fi;
   cdp.rtf(anzahl)=rrttff;
   if(tt==0){ff=(int)(2*3.1415/dfi+1);}
}}

// read pointcharge-parameters 
 cf_file = fopen_errchk (argv[5], "rb");
 float par[100];par[0]=99;
int pchere;
while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){;}
while(pchere>0)
{
printf("pointcharge %g |e| at xyz=%g %g %g mind: xyz=cab\">\n",par[1],par[2],par[3],par[4]);

fprintf(fout,"<geometry name=\"");
fprintf(fout,"pointcharge %g at xyz=%g %g %g in javaview xyz=cab\">\n",par[1],par[2],par[3],par[4]);
fprintf(fout,"<pointSet color=\"hide\" point=\"show\" dim=\"1\">\n");
fprintf(fout,"<points num=\"%i\">\n",cdp.nofpoints());
int ii;
 for(ii=1;ii<=cdp.nofpoints();++ii)
    {double dx,dy,dz,R,fi,theta;
//     R=0.5*pow(abs(par[1]),0.33333333);
     R=0.5*cbrt(abs(par[1]));
     theta=cdp.rtf(ii)(2);
     fi=cdp.rtf(ii)(3);
     // mind abc||yzx
     dx=R*sin(theta)*cos(fi)+par[2];
     dy=R*sin(theta)*sin(fi)+par[3];
     dz=R*cos(theta)+par[4];
     fprintf(fout,"<p>%4g %4g %4g</p>\n",dx,dy,dz);
     }

    fprintf(fout,"<thickness>0.0</thickness><color type=\"rgb\">255 0 0</color><colorTag type=\"rgb\">255 0 255</colorTag>\n");
    fprintf(fout,"</points>			</pointSet>\n");
    fprintf(fout,"<faceSet face=\"show\" edge=\"show\">\n");
    fprintf(fout,"<faces num=\"%i\">\n",cd.nofpoints());
    int i,offset=0;

int ntt,nff,pointnr,ffnr,p1,p2,p3,p4;
ntt=(int)(3.1415/dtheta);
nff=(int)(2*3.1415/dfi);
pointnr=ntt*(nff+1);
ffnr=nff+1;
for(tt=1;tt<=ntt;++tt){for(ff=0;ff<=nff;++ff){
   p1 = ff + 1 + (tt - 2) * ffnr+offset;
   p2 = ff + 2 + (tt - 2) * ffnr+offset;
   p3 = ff + 2 + (tt - 1) * ffnr+offset;
   p4 = ff + 1 + (tt - 1) * ffnr+offset;
   if (ff==nff){p3 = p3 - ffnr; p2 = p2 - ffnr;}
   if (tt==1) {p1 = offset; p2 = offset;}
   fprintf(fout,"<f> %i %i %i %i </f>\n",p1,p2,p3,p4);
  }}
offset+=pointnr+1;
if (par[1]>0){fprintf(fout,"<color type=\"rgb\"> 255 0 0</color>\n");}
             else
             {fprintf(fout,"<color type=\"rgb\">0  0 255</color>\n");}
fprintf(fout,"<colorTag type=\"rgb\">255 0 255</colorTag>\n");
fprintf(fout,"</faces></faceSet></geometry>\n");
while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){}
}
fclose(cf_file);


fprintf(fout,"</geometries></jvx-model>\n");
  fclose (fout);

  return 0;
}


