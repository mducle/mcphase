/***********************************************************************
 *
 * charges.c - program to display charges at given htpoint
 *
 ***********************************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100
#define MUB  5.788378E-02 // Bohrmagneton in meV/tesla

#include "chargedensity.hpp"
#include "spincf.hpp"
#include "martin.h"
#include "myev.h"
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>
#include<par.hpp>


// sub for calculation of charge density given a radiu R and polar angles teta, 
// fi and expansion coeff. alm

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
// hauptprogramm
int main (int argc, char **argv)
{ spincf savmf;
 FILE * fin_coq, * fout;
 float delta,dd,ddT,ddHa,ddHb,ddHc,alpha,beta,gamma;
 int n=0,nofatoms=0,nofcomponents=3;
 long int pos=0,j;
 float numbers[11];numbers[9]=1;numbers[10]=3;
 numbers[0]=11;
 char instr[MAXNOFCHARINLINE];
 char outstr[MAXNOFCHARINLINE];
 float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS];
 char * cffilenames[MAXNOFATOMS];
// ComplexMatrix * eigenstates[MAXNOFATOMS];
  Matrix r(1,3,1,3);
  Vector abc(1,3);

// check command line
  if (argc < 5)
    { printf (" program charges - display charges at HT point\n\
                use as: charges T Ha Hb Hc [file.mf]\n");
      printf ("\n ... to view chargeplots type: java javaview charges.jvx\n\n");
      printf ("mind: 1)the 4f wavefunction R4f(r) is taken the same for all RE.\n");
      printf ("      2) the coordinate system is xyz||cab\n");
      exit (1);
    }
    else
    { printf ("mind: the 4f wavefunction R4f(r) is taken the same for all RE.\n");}

 if (argc <6) 
 { fin_coq = fopen_errchk ("./mcphas.mf", "rb");}
 else
 { fin_coq = fopen_errchk (argv[5], "rb");}
    
 fout = fopen_errchk ("./charges.out", "w");



abc=0;char *token;
 // input file header ------------------------------------------------------------------
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
   if (pos==-1) 
       {fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq);
   // strip /r (dos line feed) from line if necessary
    while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

   if (instr[strspn(instr," \t")]=='#'){fprintf(fout,"%s",instr);}
   if(abc[1]==0){extract(instr,"a",abc[1]);extract(instr,"b",abc[2]); extract(instr,"c",abc[3]); 
                 extract(instr,"alpha",alpha);  extract(instr,"beta",beta);extract(instr,"gamma",gamma); 
   }
   extract(instr,"r1x",r[1][1]);extract(instr,"r2x",r[1][2]); extract(instr,"r3x",r[1][3]); 
   extract(instr,"r1y",r[2][1]); extract(instr,"r2y",r[2][2]); extract(instr,"r3y",r[2][3]);
   extract(instr,"r1z",r[3][1]); extract(instr,"r2z",r[3][2]); extract(instr,"r3z",r[3][3]);
   extract(instr,"r1a",r[1][1]);extract(instr,"r2a",r[1][2]); extract(instr,"r3a",r[1][3]); 
   extract(instr,"r1b",r[2][1]); extract(instr,"r2b",r[2][2]); extract(instr,"r3b",r[2][3]);
   extract(instr,"r1c",r[3][1]); extract(instr,"r2c",r[3][2]); extract(instr,"r3c",r[3][3]);
   extract(instr,"nofatoms",nofatoms);    extract(instr,"nofcomponents",nofcomponents); 
   if (nofatoms>0&&(extract(instr,"x",x[n+1])+
                   extract(instr,"y",y[n+1])+
  		       extract(instr,"z",z[n+1])==0)||
		       (extract(instr,"da",x[n+1])+
                   extract(instr,"db",y[n+1])+
		       extract(instr,"dc",z[n+1])==0))
		  {++n;if(n>nofatoms||nofatoms>MAXNOFATOMS) 
                    {fprintf(stderr,"ERROR charges.c reading file:maximum number of atoms in unit cell exceeded\n");exit(EXIT_FAILURE);}
                   cffilenames[n]=new char[MAXNOFCHARINLINE];
                   extract(instr,"cffilename",cffilenames[n],(size_t)MAXNOFCHARINLINE);
//		   printf("%s\n",cffilenames[n]);
                  }
  }
  if (alpha!=90||beta!=90||gamma!=90)
  {fprintf(stderr,"ERROR: non orthogonal lattice not supported yet\n");exit(EXIT_FAILURE);}
  
  
// load mfconfigurations and check which one is nearest -------------------------------   
  double T,ha,hb,hc;

   j=fseek(fin_coq,pos,SEEK_SET); 
    if (j!=0){fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
   
 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      spins.load(fin_coq);
      ddT=strtod(argv[1],NULL)-numbers[3];ddT*=ddT;
      ddHa=strtod(argv[2],NULL)-numbers[5];ddHa*=ddHa;
      ddHb=strtod(argv[3],NULL)-numbers[6];ddHb*=ddHb;
      ddHc=strtod(argv[4],NULL)-numbers[7];ddHc*=ddHc;
      dd=sqrt(ddT+ddHa+ddHb+ddHc+0.000001);
      if (dd<delta)
       {delta=dd;
        sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g n=%g spins nofatoms=%i in primitive basis nofcomponents=%i",numbers[3],numbers[5],numbers[6],numbers[7],numbers[8],(int)numbers[9],(int)numbers[10]);
        savmf=spins;
	T=numbers[3];ha=numbers[5];hb=numbers[6];hc=numbers[7];
       }
    }
  fclose (fin_coq);
  
// create plot of spin+chargeconfiguration -----------------------------------------------------------
  int i,ii,nt,k,l,m;
  double lnz,u;
  float d;


Matrix a(0,6,-6,6);
  a(0, 0) = 1 / sqrt(4 * 3.1415);
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

  Vector tetan(1,6);
  par inputpars("./mcphas.j");
  double dtheta=0.2; //stepwidth to step surface
  double dfi=0.2;
  double ccc=0.05; //surface value of chargedensity
  Vector rrttff(1,3,1,3);
int tt,ff,iii,iv;
  
  Vector h(1,48);
  Vector moments(1,48);
  Vector hh(1,savmf.nofcomponents*savmf.nofatoms);
  spincf extendedspincf(savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,48);
  chargedensity * cd[savmf.na()*savmf.nb()*savmf.nc()*savmf.nofatoms+1];
  for(i=1;i<=savmf.na()*savmf.nb()*savmf.nc()*savmf.nofatoms;++i)cd[i]=new chargedensity(dtheta,dfi);

          // the following is for the printout of charges.out ...........................
           fprintf(fout,"#T=%g K Ha=%g T Hb= %g T Hc= %g T: nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",T,ha,hb,hc,savmf.na(),savmf.nb(),savmf.nc(),inputpars.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
            fprintf(fout,"#J=value {atom-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3]  <Ja> <Jb> <Jc> ...\n");
            fprintf(fout,"# Eigenvalues [meV] and eigenvectors [as columns]\n");
	  // determine primitive magnetic unit cell
           Vector nofabc(1,3),dd3(1,3),pa(1,3),pb(1,3),pc(1,3);
           Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3);
           nofabc(1)=savmf.na();nofabc(2)=savmf.nb();nofabc(3)=savmf.nc();
             for (i=1;i<=3;++i){for(j=1;j<=3;++j) {dd3(j)=nofabc(j)*r(i,j)*abc(i);p(i,j)=dd3(j);}}
             pa=p.Column(1);  //primitive magnetic unit cell
             pb=p.Column(2);
             pc=p.Column(3);
        // .............................................................................                                
	       
//  1. from the meanfieldconfiguration (savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: mcalc has to be used to calculate all the <Olm>.

 for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
 {
    hh=savmf.m(i,j,k);
  for(ii=1;ii<=inputpars.nofatoms;++ii)
 {
    h=0;
   for(nt=1;nt<=savmf.nofcomponents;++nt){h(nt)=hh(nt+savmf.nofcomponents*(ii-1));}

            moments=(*inputpars.jjj[ii]).mcalc(T,h,lnz,u); // here we trigger single ion 
                                                           // module to calculate all 48
                                                           // higher order moments 
 

          // output atoms and moments in primitive unit cell to stdout
         dd3(1)=x[ii]*abc(1);
         dd3(2)=y[ii]*abc(2);
         dd3(3)=z[ii]*abc(3);
         dd3+=pa*(double)(i-1)/nofabc(1)+pb*(double)(j-1)/nofabc(2)+pc*(double)(k-1)/nofabc(3);
         dd0=p.Inverse()*dd3;dd0(1)*=savmf.na();dd0(2)*=savmf.nb();dd0(3)*=savmf.nc();
              fprintf(fout,"J=%4.1f {%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
	              (*inputpars.jjj[ii]).J(),cffilenames[ii],dd3(1)/abc(1),dd3(2)/abc(2),dd3(3)/abc(3),dd0(1),dd0(2),dd0(3));
                     for(nt=1;nt<=48;++nt)
		        {extendedspincf.m(i,j,k)(nt+48*(ii-1))=moments(nt);
                         fprintf(fout," %4.4f",extendedspincf.m(i,j,k)(nt+48*(ii-1)));}
                         fprintf(fout,"\n");
                             
	                 myPrintComplexMatrix(fout,(*inputpars.jjj[ii]).eigenstates(h));      
							   // ... and the eigenvalues + eigenvectors !


// mind  if there is kramers or another si module, not all the Olm vectors are available - then we plot a sphere

// how do we get the stevens factors - look into cfield if they are there !!!
 //  tetan(2) = alpha: tetan(4) = beta: tetan(6) = gamma
  tetan=(*inputpars.jjj[ii]).tetan();
  printf("alpha beta gamma  %4g %4g %4g ",tetan(2),tetan(4),tetan(6));
  printf("\n");

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
int atomnr=ii+(i-1)*savmf.nofatoms+(j-1)*savmf.nofatoms*savmf.na()+(k-1)*savmf.nofatoms*savmf.na()*savmf.nb();
int anzahl = 0;
int imax = 3;
Vector rp(1,imax);
double rmax = 0;
double rstp = .1;
double max = .01 * ccc;  //end of intervalschachtelung to find r(ro=ccc)
double theta,fi,R,rin,ro,deltaa,rstpp,delta1,dd;

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
   (*cd[atomnr]).rtf(anzahl)=rrttff;
   if(tt==0){ff=(int)(2*3.1415/dfi+1);}
}}


  
// 3.: copy from chrgplot.bas the calculation module of the charge density given the <Olm> - creating both a
//     postscript output and a jvx file. Call this for every spin subsequently. Maybe for the postscript
//     files look a bit on the .eps3d function in how to create such a plot


  }}}}
//  extendedspincf.printall(fout,abc,r,x,y,z,cffilenames);
  fclose(fout);





//print out the long vector of moments 1-48
  printf("%s - spin configuration <Olm>(i)\n",outstr);
  extendedspincf.print(stdout);

  fout = fopen_errchk ("./charges.jvx", "w");

 fprintf(fout,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n");
 fprintf(fout,"<!DOCTYPE jvx-model SYSTEM \"http://www.javaview.de/rsrc/jvx.dtd\">\n");
 fprintf(fout,"<jvx-model>	<meta generator=\"JavaView v.2.00.008\"/>\n");
 fprintf(fout,"<meta date=\"Wed Mar 07 22:30:58 GMT+01:00 2001\"/>\n");
 fprintf(fout,"<version type=\"dump\">0.02</version>\n");

fprintf(fout,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT in javaview xyz=cab</title>\n", T, ha, hb, hc);
fprintf(fout,"<geometries> 		<geometry name=\"");
fprintf(fout,"T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT in javaview xyz=cab\">\n",T,ha,hb,hc);
fprintf(fout,"<pointSet color=\"hide\" point=\"show\" dim=\"1\">\n");
fprintf(fout,"<points num=\"%i\">\n",savmf.na()*savmf.nb()*savmf.nc()*savmf.nofatoms*(*cd[1]).nofpoints());
 for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k){for(l=1;l<=savmf.nofatoms;++l)
{Vector dadbdc(1,3);
 dadbdc=savmf.pos(i,j,k,l, abc, r,x,y,z);
 int ic=l+(i-1)*savmf.nofatoms+(j-1)*savmf.nofatoms*savmf.na()+(k-1)*savmf.nofatoms*savmf.na()*savmf.nb();
 for(ii=1;ii<=(*cd[i]).nofpoints();++ii)
    {double dx,dy,dz,R,fi,theta;
     R=(*cd[ic]).rtf(ii)(1);
     theta=(*cd[ic]).rtf(ii)(2);
     fi=(*cd[ic]).rtf(ii)(3);
     // mind abc||yzx
     dx=R*sin(theta)*cos(fi)+dadbdc(3);
     dy=R*sin(theta)*sin(fi)+dadbdc(1);
     dz=R*cos(theta)+dadbdc(2);
     fprintf(fout,"<p>%4g %4g %4g</p>\n",dx,dy,dz);
     }
}}}}

    fprintf(fout,"<thickness>0.0</thickness><color type=\"rgb\">255 0 0</color><colorTag type=\"rgb\">255 0 255</colorTag>\n");
    fprintf(fout,"</points>			</pointSet>\n");
    fprintf(fout,"<faceSet face=\"show\" edge=\"show\">\n");
    fprintf(fout,"<faces num=\"%i\">\n",savmf.na()*savmf.nb()*savmf.nc()*savmf.nofatoms*(*cd[1]).nofpoints());
    int offset=0;
for(i=1;i<=savmf.na()*savmf.nb()*savmf.nc()*savmf.nofatoms;++i)
{int ntt,nff,pointnr,ffnr,p1,p2,p3,p4;
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
}     
fprintf(fout,"<color type=\"rgb\">100 230 255</color>\n");
fprintf(fout,"<colorTag type=\"rgb\">255 0 255</colorTag>\n");
fprintf(fout,"</faces></faceSet></geometry></geometries></jvx-model>\n");
  fclose (fout);



  for(i=1;i<=savmf.na()*savmf.nb()*savmf.nc()*savmf.nofatoms;++i)delete cd[i];
  for(i=1;i<=nofatoms;++i){  delete cffilenames[i];}

  return 0;
}


