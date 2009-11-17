// *************************************************************************
// ************************ class jjjpar     *******************************
// *************************************************************************
// jjjpar is a class to store all parameters associated
// with a specific ion, for example CEF parameters
// and exchange parameters 
// it is used by many programs in the package
// moreover, it loads also the user defined single ion module functions (linux only)
#include "jjjpar.hpp"
#include "../../version"

#define MAXNOFNUMBERSINLINE 200
#define MAXNOFCHARINLINE 1024

#define MU_B 0.05788
#define K_B  0.0862
#define SMALL 1e-6   //!!! must match SMALL in mcdisp.c and ionpars.cpp !!!
                     // because it is used to decide wether for small transition
		     // energy the matrix Mijkl contains wn-wn' or wn/kT
#define PI 3.1415926535

#include "jjjpar_intmod_kramer.cpp"   // some functions for module_type=1
#include "jjjpar_intmod_brillouin.cpp"// some functions for module_type=3



/****************************************************************************/
// function to calculate magnetisation M from effective field H
// this is the heart of the meanfield algorithm an it is necessary to
// keep this routine as efficient as possible
// at the moment we do only groundstate doublet
/****************************************************************************/
void jjjpar::mcalc (Vector &mom, double & T, Vector &  gjmbH, double & lnZ,double & U,ComplexMatrix & ests)
{switch (module_type)
  {case 1: kramer(mom,T,gjmbH,lnZ,U);break;
   case 2: (*iops).cfieldJJ(mom,T,gjmbH,lnZ,U,ests);break;
   case 3: brillouin(mom,T,gjmbH,lnZ,U);break;
   default: (*m)(&mom,&T,&gjmbH,&gJ,&ABC,&cffilename,&lnZ,&U,&ests);
  }
}

/****************************************************************************/
// this function returns n (the number of transitions in the single ion susceptibility)
// the transition matrix mat corresponding to jjjpar.transitionnumber and delta
// for effective field heff and temperature given on input
/****************************************************************************/
int jjjpar::dmcalc(double & T,Vector & gjmbheff,ComplexMatrix & mat,float & delta,ComplexMatrix & ests)
{ switch (module_type)
  {case 0: if (dm!=NULL){return (*dm)(&transitionnumber,&T,&gjmbheff,&gJ,&ABC,&cffilename,&mat,&delta,&ests);}
           else return 0;
           break;
   case 1: return kramerdm(transitionnumber,T,gjmbheff,mat,delta);break;
   case 2: return (*iops).cfielddm(transitionnumber,T,gjmbheff,mat,delta,ests);break;
   case 3: return brillouindm(transitionnumber,T,gjmbheff,mat,delta);break;
   default: return 0;
  }
}


/****************************************************************************/
// returns eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::eigenstates (Vector & gjmbheff,double & T)
{switch (module_type)
  {case 0:  if(estates!=NULL){(*estates)(&est,&gjmbheff,&gJ,&T,&ABC,&cffilename);}
            return est;break;
   case 2:  est=(*iops).cfeigenstates(gjmbheff,T);return est;break;
   default: est=0;return est;
  }
}


/****************************************************************************/
// returns transition element matrix N(Q) in order to be able to go beyond 
//
// dipolar approximation in mcdisp - it requires a call to eigenstates first
//
//on input
//    transitionnumber has to be set correctly to that one which is to be computed 
//    sign(transitionnumber)... 1... without printout, -1 with extensive printout
//    est		matrix with eigenstates, eigenvalues [meV], population numbers
//    T                 temperature
//     Q                 components of Qvector in euclidian coordinates 123=abc
//  on output    
//    int   	total number of transitions
//    N(i,j)	<-|Q|+><+|Q|-> (n+-n-),  n+,n- population numbers 
//    with Q the scattering operator according to Lovesey 11.4, p 222, eq 6.87b
//     (note that  <M(Q)>=-2x<Q>_TH in units of mb)
//    .... occupation number of states (- to + transition chosen according to transitionnumber)
//   
/****************************************************************************/
int jjjpar::dncalc(Vector & Qvec,double & T, ComplexMatrix & nat,ComplexMatrix & ests)

{double J0,J2,J4,J6;
 double Q,d,s,th,ph;
            Q = Norm(Qvec); //dspacing
            d = 2.0 * PI / Q; s=0.5 / d; 
      J0=magFFj0(1)*exp(-magFFj0(2)*s*s)+magFFj0(3)*exp(-magFFj0(4)*s*s)+magFFj0(5)*exp(-magFFj0(6)*s*s)+magFFj0(7);
      J2=magFFj2(1)*exp(-magFFj2(2)*s*s)+magFFj2(3)*exp(-magFFj2(4)*s*s)+magFFj2(5)*exp(-magFFj2(6)*s*s)+magFFj2(7);
      J2*=s*s;
      J4=magFFj4(1)*exp(-magFFj4(2)*s*s)+magFFj4(3)*exp(-magFFj4(4)*s*s)+magFFj4(5)*exp(-magFFj4(6)*s*s)+magFFj4(7);
      J4*=s*s;
      J6=magFFj6(1)*exp(-magFFj6(2)*s*s)+magFFj6(3)*exp(-magFFj6(4)*s*s)+magFFj6(5)*exp(-magFFj6(6)*s*s)+magFFj6(7);
      J6*=s*s;
	 // calculate th and ph (polar angles of Q with respect to xyz of CEF)
 switch (module_type)
  {static int washere=0;
   
   case 0:if (ddnn!=NULL){getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph);
                          return (*ddnn)(&transitionnumber,&th,&ph,&J0,&J2,&J4,&J6,&ests,&T,&nat);break;}
          else {return 0;}
   case 2: getpolar(Qvec(3),Qvec(1),Qvec(2),Q,th,ph); // for internal module cfield xyz||cba and we have to give cfielddn polar angles with respect to xyz
           return (*iops).cfielddn(transitionnumber,th,ph,J0,J2,J4,J6,Zc,ests,T,nat);break;
   default: if(washere==0){fprintf(stderr,"Warning in scattering operator function dncalc - for ion %s \ngoing beyond dipolar approximation is not implemented\n",cffilename);
                           washere=1;}
            return 0;
  }

}




/****************************************************************************/
// calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
// according to stored eigenstate matrix est
// input: Qvec ..... Q Vector components 123=xyz=cab
/****************************************************************************/
ComplexVector & jjjpar::MQ(Vector & Qvec)
{double J0,J2,J4,J6;
 double Q,d,s,th,ph;
            Q = Norm(Qvec); //dspacing
            d = 2.0 * PI / Q; s=0.5 / d; 
      J0=magFFj0(1)*exp(-magFFj0(2)*s*s)+magFFj0(3)*exp(-magFFj0(4)*s*s)+magFFj0(5)*exp(-magFFj0(6)*s*s)+magFFj0(7);
      J2=magFFj2(1)*exp(-magFFj2(2)*s*s)+magFFj2(3)*exp(-magFFj2(4)*s*s)+magFFj2(5)*exp(-magFFj2(6)*s*s)+magFFj2(7);
      J2*=s*s;
      J4=magFFj4(1)*exp(-magFFj4(2)*s*s)+magFFj4(3)*exp(-magFFj4(4)*s*s)+magFFj4(5)*exp(-magFFj4(6)*s*s)+magFFj4(7);
      J4*=s*s;
      J6=magFFj6(1)*exp(-magFFj6(2)*s*s)+magFFj6(3)*exp(-magFFj6(4)*s*s)+magFFj6(5)*exp(-magFFj6(6)*s*s)+magFFj6(7);
      J6*=s*s;
            complex<double>dummy;
switch (module_type)
  {case 0:  getpolar(Qvec(2),Qvec(3),Qvec(1),Q,th,ph); // for external module we must provide th and ph with respect 
                                                       // to abc coordinate system
            (*mq)(&Mq,&th,&ph,&J0,&J2,&J4,&J6,&est);
             // external module provide Mq(123)=Mq(abc)
             // we must transform this to mcdiff internal xyz||cab coordinate system
            dummy=Mq(3);Mq(3)=Mq(2);Mq(2)=Mq(1);Mq(1)=dummy;
            return Mq;break;
   case 2:  getpolar(Qvec(1),Qvec(2),Qvec(3),Q,th,ph); // internal module cfield does not need transformation
            return (*iops).MQ(th,ph,J0,J2,J4,J6,Zc,est);break;
   default: fprintf(stderr,"ERROR in scattering operator function M(Q) for ion %s \nM(Q) is currently only implemented for internal module cfield:\n",cffilename);exit(EXIT_FAILURE);
  }
}



/************************************************************************************/
//  RETURN TOTAL FORMFACTOR, 
//    however if gJ=0 and Q>0 return spin form factor FS(Q)=<j0(Q)>
//            if gJ=0 and Q<0 return angular  form factor FL(Q)=<j0(Q)>+<j2(Q)>
//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
/************************************************************************************/
   double jjjpar::F(double Q)
   {double s,j0,j2;    
    s=Q/4/PI;
    j0=magFFj0(1)*exp(-magFFj0(2)*s*s)+magFFj0(3)*exp(-magFFj0(4)*s*s)+magFFj0(5)*exp(-magFFj0(6)*s*s)+magFFj0(7);
    if(gJ==0&&Q>0){return j0;} // in case of intermediate coupling return spin form factor 
    j2=magFFj2(1)*exp(-magFFj2(2)*s*s)+magFFj2(3)*exp(-magFFj2(4)*s*s)+magFFj2(5)*exp(-magFFj2(6)*s*s)+magFFj2(7);
    j2*=s*s;
    if(gJ==0&&Q<0){return j0+j2;} // in case of intermediate coupling return angular form factor 
   return (j0 + j2 * (2 / gJ - 1)); // formfactor F(Q) for rare earth 

   }

/************************************************************************************/
//   debyewallerfactor = exp(-2 * DWF *s*s)      (sf ~ exp(-2 DWF sin^2(theta) / lambda^2)=EXP (-W),  (2*DWF=B=8 pi^2 <u^2>)
/************************************************************************************/
   double jjjpar::debyewallerfactor(double & Q)
   {double s;
    s=Q/4/PI;
    return exp(-2*DWF*s*s);
   }

/************************************************************************************/
// calculates polar coordinates from Vector X(1..3)
/************************************************************************************/
void jjjpar::getpolar(double x,double y, double z, double & r, double & th, double & ph)
{	 r=sqrt(x*x+y*y+z*z);
         th=acos(z/r);
	 if(sin(th)>=SMALL){
	                    if(x>0) {ph=acos(x/(r*sin(th))-SMALL);}
			    else    {ph=acos(x/(r*sin(th))+SMALL);}
			   }
			 else{ph=0;}
	 if (y<0){ph=2*PI-ph;} 
}


/************************************************************************************/
// returns total angular momentum quantum number J
/************************************************************************************/
double jjjpar::J()
{
 switch (module_type)
  {
   case 2: return (*iops).J;break;
   case 3:  return ABC[1];break;
   default: fprintf (stderr, "error class jjjpar: single ion module does not allow to calculate stevens parameters alpha beta gamma \n"); 
            exit (EXIT_FAILURE);
  }
}

/************************************************************************************/
// returns stevens parameters of ion
/************************************************************************************/
Vector & jjjpar::tetan ()
{static Vector tt(1,6);
 tt=0;
 switch (module_type)
  {
   case 2: tt(2)=(*iops).alpha;tt(4)=(*iops).beta;tt(6)=(*iops).gamma;
            return tt;break;
   default: fprintf (stderr, "error class jjjpar: single ion module does not allow to calculate stevens parameters alpha beta gamma \n"); 
            exit (EXIT_FAILURE);
  }
}


/************************************************************************************/
// evaluate radial wave function
/************************************************************************************/
double jjjpar::radial_wavefunction(double rr) // rr given in Angstroems, returns R(r) in units of 1/A^1.5
   {//printf("%g ",rr);

    double R=0;int p;double a0=0.5292;
    int ok=0;
    double r=rr/a0;// r is the distance in units of a0
    for(p=1;p<=9;++p){if(Np(p)!=0){ok=1;
                                   R+=exp(-Xip(p)*r)*pow(r,Np(p)-1)*Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of radial wave function R(r=%g) failed due to Xi%i<=0 - continuing with R(r=%g)=0\n",r,p,r);return 0;}
                     }            }    
    // now we have R in units of 1/a0^1.5
    R/=sqrt(a0*a0*a0);
    // now we have R in units of 1/A^1.5
//printf("%g ",R);
    if (ok==1) return R;


//  we have to find the 4f wavefunction R4f(r) for each single ion and the Zlm, cfield has nothing: so we have
//     to take this from chrgplt.bas - a little problem: how do we get the correct R4f(r) ? for a first attempt
//     we could just take the same for all RE.

    static int washere=0;    
    if(washere==0){washere=1;fprintf (stderr,"Warning: radial wave function parameters not found, will use 4f hydrogen radial wave function\n");}
double rs;
//k^2 = 11 / 10 * 11 / 9 * 11 / 8 * 11 / 7 * 11 / 6 * 11 / 5 * 11 / 4 * 11 / 3 * 11 / 2 * 11 / 1 * 11
rs = rr * exp(-rr);
R = 280.4 * rs * rs * rs * rs  * exp(-1.5 * rr);
//printf("R4f(%g)=%g\n ",rr,R);
return R;
   }

   //functions to calculate radial matrix elements <r^n> from radial wave function in units of a0=0.5292 A
   double jjjpar::rk_from_radial_wavefunction(int k)
   {int p,q, pmax=0;
    Vector coeff(1,9);
    
    for(p=1;p<=9;++p){if(Np(p)!=0){pmax=p;
                                   coeff(p)=Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of <r^%i> failed due to Xi%i<=0 - continuing with <r^%i>=0\n",k,p,k);return 0;}
                     }            }
    if(pmax==0){fprintf (stderr,"Warning: calculation of <r^%i> failed - continuing with <r^%i>=0\n",k,k);return 0;}
    double rk=0;
    for(p=1;p<=pmax;++p){
    for(q=1;q<=pmax;++q){
                         rk+=coeff(p)*coeff(q)*factorial((int)Np(p)+(int)Np(q)+k)/pow(Xip(p)+Xip(q),Np(p)+Np(q)+k+1);
    }}
   return rk;
   }

   int jjjpar::r2_from_radial_wavefunction() {r2=rk_from_radial_wavefunction(2);if(module_type==2){(*iops).r2=r2;}}
   int jjjpar::r4_from_radial_wavefunction() {r4=rk_from_radial_wavefunction(4);if(module_type==2){(*iops).r4=r4;}}
   int jjjpar::r6_from_radial_wavefunction() {r6=rk_from_radial_wavefunction(6);if(module_type==2){(*iops).r6=r6;}}

void jjjpar::save_radial_wavefunction(const char * filename)
   {double r=0.1;
    FILE * fout;
    if (radial_wavefunction(r)==0){fprintf(stderr,"Warning: save_radial_wavefunction not possible\n");return;}
    fout=fopen_errchk(filename,"w");
    fprintf(fout,"# radial wave function for %s\n",cffilename);
    fprintf(fout,"# the radial wave function is expanded as \n");
    fprintf(fout,"# R(r)=sum_p C_p R_Np,XIp(r)\n");
    fprintf(fout,"# R_Np,XIp(r)=r^(Np-1).exp(-XIp * r).(2 * XIp)^(Np+0.5)/sqrt(2Np!)\n");
    fprintf(fout,"# radial wave function parameters Np XIp Cp values are\n");
    fprintf(fout,"# tabulated in clementi & roetti Atomic data and \n");
    fprintf(fout,"# nuclear data tables 14 (1974) 177-478 for the transition metals\n");
    fprintf(fout,"# for rare earth parameters can be found in Freeman and Watson PR 127 (1962) 2058\n");
    fprintf(fout,"# and O. Sovers, J. Phys. Chem. Solids Vol 28 (1966) 1073\n");
    fprintf(fout,"# the parameters used are: \n");
    int p;    
    for(p=1;p<=9;++p){if(Np(p)!=0){fprintf(fout,"#! N%i=%g XI%i=%g C%i=%g\n",p,Np(p),p,Xip(p),p,Cp(p));}}
    fprintf(fout,"# r[A]  vs R(r)[1/A^1.5]\n");
    for(r=0.01;r<=10;r*=1.05){fprintf(fout,"%8.8g  %8.8g\n",r,radial_wavefunction(r));}
    fclose(fout);
   }


// sub for calculation of charge density given a radiu R and polar angles teta, 
// fi and expansion coeff. alm
double jjjpar::rocalc (double & teta,double & fi,double & R, Vector & moments)
{double ro,ct,ct2,st,st2,sfi,cfi,rs,rr;

if (R>4.0||R<0){ro=0;}else{
ct = cos(teta);                      //z
ct2 = ct * ct;
st = sin(teta);
st2 = st * st;
sfi = sin(fi);
cfi = cos(fi);
 int l,m;
  Vector tetan(1,6);
 int offset=0; // for cfield module
 if(module_type==0){offset=3;} // for external modules this is to treat ic1ion correctly
 if(module_type==2){
                    tetan(2)=(*iops).alpha;// Stevens factors
                    tetan(4)=(*iops).beta;
                    tetan(6)=(*iops).gamma;
                   }

Matrix a(0,6,-6,6);
a(0, 0) = 1 / sqrt(4.0 * 3.1415);
if(calcmagdensity>0)a(0, 0) = moments(calcmagdensity) / sqrt(4.0 * 3.1415);

a(2,-2)=moments(offset+4);
a(2,-1)=moments(offset+5);
a(2,0)=moments(offset+6);
a(2,1)=moments(offset+7);
a(2,2)=moments(offset+8);

a(4,-4)=moments(offset+16);
a(4,-3)=moments(offset+17);
a(4,-2)=moments(offset+18);
a(4,-1)=moments(offset+19);
a(4, 0)=moments(offset+20);
a(4, 1)=moments(offset+21);
a(4, 2)=moments(offset+22);
a(4, 3)=moments(offset+23);
a(4, 4)=moments(offset+24);

a(6,-6)=moments(offset+36);
a(6,-5)=moments(offset+37);
a(6,-4)=moments(offset+38);
a(6,-3)=moments(offset+39);
a(6,-2)=moments(offset+40);
a(6,-1)=moments(offset+41);
a(6,-0)=moments(offset+42);
a(6, 1)=moments(offset+43);
a(6, 2)=moments(offset+44);
a(6, 3)=moments(offset+45);
a(6, 4)=moments(offset+46);
a(6, 5)=moments(offset+47);
a(6, 6)=moments(offset+48);


// r given in Angstroems, returns R(r) in units of 1/A^1.5
rr=radial_wavefunction(R);
rr=rr*rr;// then the chargedensity will be in units of 1/A^3

if(module_type==2){for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m){a(l,m)*=tetan(l)*cnst(l,m)*cnst(l,m);}}
         } // these are prefactors in case of module cfield (stevens parameters tetan and zlm prefactors)
else     {for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m){if(m!=0){a(l,m)*=cnst(l,m)*sqrt((2.0*l+1)/8/PI);}else{a(l,m)*=cnst(l,m)*sqrt((2.0*l+1)/4/PI);}}}
         } // in case
           // of module ic1ion we just take the prefactors of the Zlm ... ??? what should we take here ???

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



void jjjpar::set_zlm_constants()
{// cnst is the Zlm constants - put them into the matrix ... (same code is reused in ionpars.cpp)
 cnst= Matrix(0,6,-6,6);
 
cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;
cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;
cnst(6,0)=  0.06357014;
cnst(6,1)=  0.582621;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
int l,m;
for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);} 
}

// additional
// function to look if string s2 lies in string s1, checking the first n characters of s2
int strncomp(const char * s1,const char * s2, size_t n)
{size_t i;
 if (strlen(s1)>=n)
 {for (i=0;i<=strlen(s1)-n;++i)
  {if (strncmp(&s1[i],s2,n)==0){return 0;}
  }
 }
 return strncmp(s1,s2,n);
}












void jjjpar::increase_nofcomponents(int n) // increase nofcomponents by n
{int i,j,k,nold;
  nold=nofcomponents;
  nofcomponents+=n;
  mom.Resize(1,nofcomponents); 

  Matrix * jijstore;
  jijstore = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jijstore[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jijstore == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  for (i=1;i<=paranz;++i)
   {jijstore[i]=0;
    for (j=1;j<=nold;++j)
    {for (k=1;k<=nold;++k)
     {jijstore[i](j,k)=jij[i](j,k);
   }}}
 

 delete []jij;
  jij = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jij[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  for (i=1;i<=paranz;++i)
  {jij[i]=jijstore[i];}

  delete[] jijstore;
  fprintf(stderr,"Warning: increasing nofcomponents not tested  yet ... addition of parameter sets may be erroneous\n");
//  exit(0);

}


void jjjpar::add(jjjpar & b,Vector & abc) // add set b to this (abc: lattice constants)
{int i,j; 
 if(diagonalexchange==1&&b.diagonalexchange==0)diagonalexchange=0;
  if (nofcomponents!=b.nofcomponents)
  { fprintf (stderr, "class jjjpar: function add - nofcomponents does not match (check number of columns)\n"); 
    exit (EXIT_FAILURE);}

 for(i=1;i<=b.paranz;++i)
 {int found=0;
  for(j=1;j<=paranz&&found==0;++j)
  {if(Norm(dn[j]-b.dn[i])<SMALL)
    {//parameter found in list
     jij[j]+=b.jij[i];found=1;
    }
  }
  if (found==0){ // parameter not found in list !!!
                 jjjpar c(1,diagonalexchange,nofcomponents);
                 for(j=1;j<=paranz&&abc(1)*abc(1)*dn[j](1)*dn[j](1)+
		                    abc(2)*abc(2)*dn[j](2)*dn[j](2)+
				    abc(3)*abc(3)*dn[j](3)*dn[j](3)
				    <
				    abc(1)*abc(1)*b.dn[i](1)*b.dn[i](1)+
		                    abc(2)*abc(2)*b.dn[i](2)*b.dn[i](2)+
				    abc(3)*abc(3)*b.dn[i](3)*b.dn[i](3)
				    ;++j); //look for matching distance
		c.dn[1]=b.dn[i];
		c.jij[1]=b.jij[i];
		addpars(j,c);
               }
 }
}

// enlarge the set of parameters 
// inserting a set of exchange parameters
// into field at position number
void jjjpar::addpars (int number, jjjpar & addjjj)
{ Matrix * jijn;
  Vector * dnn;
  int i;
  jijn = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jijn[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  dnn = new Vector[paranz+1];for(i=0;i<=paranz;++i){dnn[i]=Vector(1,3);}
  
  for (i=1;i<=paranz;++i)
  {jijn[i]=jij[i];
   dnn[i]=dn[i];
  }
  
  if (diagonalexchange!=addjjj.diagonalexchange)
  { fprintf (stderr, "class jjjpar: function addpar - diagonalexchange does not match\n"); 
    exit (EXIT_FAILURE);}
  if (nofcomponents!=addjjj.nofcomponents)
  { fprintf (stderr, "class jjjpar: function addpar - nofcomponents does not match (check number of columns)\n"); 
    exit (EXIT_FAILURE);}
      
    paranz+=addjjj.paranz;  // increase parameters   
  
  delete []jij;
  delete []dn;
  delete []sublattice;
  dn = new Vector[paranz+1];for(i=0;i<=paranz;++i){dn[i]=Vector(1,3);}
  if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  jij = new Matrix[paranz+1];for(i=0;i<=paranz;++i){jij[i]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

// setup new field jij, dn
  for (i=1;i<number;++i)
  {jij[i]=jijn[i];dn[i]=dnn[i];}
  
  for (i=number;i<number+addjjj.paranz;++i)
  {jij[i]=addjjj.jij[i-number+1];dn[i]=addjjj.dn[i-number+1];}
  
  for (i=number+addjjj.paranz;i<=paranz;++i)
  {jij[i]=jijn[i-addjjj.paranz];dn[i]=dnn[i-addjjj.paranz];}
  delete []jijn;
  delete []dnn;
}

/************************************************************************************/
// save/get parameters 
/************************************************************************************/

//saving parameters to file
void jjjpar::save(FILE * file) 
{ int i,i1,j1;
  saveatom(file);
  fprintf(file,"#da[a]   db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
// save the exchange parameters to file (exactly paranz parameters!)
  for  (i=1;i<=paranz;++i)
  {fprintf(file,"%-+8.6g %-+8.6g %-+8.6g  ",dn[i](1),dn[i](2),dn[i](3));
    // format of matrix 
  // 11 22 33 12 21 13 31 23 32 (3x3 matrix)
  // 11 22 33 44 12 21 13 31 14 41 23 32 24 42 34 43 (4x4 matrix)
  // 11 22 33 44 55 12 21 13 31 14 41 15 51 23 32 24 42 25 52 34 43 35 53 45 54 (5x5 matrix)
  // etc ...
  //save diagonal components of exchange matrix
  for(i1=1;i1<=nofcomponents;++i1){fprintf(file,"%-+8.6e ",jij[i](i1,i1));}
  //save off-diagonal components of exchange matrix (if required)
  if (diagonalexchange==0){for(i1=1;i1<=nofcomponents-1;++i1)
                              {for(j1=i1+1;j1<=nofcomponents;++j1)
                               {fprintf(file,"%-+8.6e %-+8.6e ",jij[i](i1,j1),jij[i](j1,i1));
			       }
			      }
                          }
   fprintf(file,"\n");  
  }
}

void jjjpar::saveatom(FILE * file) 
{   fprintf(file,"#! da=%4.6g [a] db=%4.6g [b] dc=%4.6g [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g cffilename=%s\n",xyz(1),xyz(2),xyz(3),paranz,diagonalexchange,gJ,cffilename);
}

//save single ion parameter file filename to path*
void jjjpar::save_sipf(const char * path)
{char  instr[MAXNOFCHARINLINE];
 char * savfilename;
 int i;
 savfilename= new char[strlen(cffilename)+strlen(path)+2];
 strcpy(savfilename,path);
 strcpy(savfilename+strlen(path),cffilename);
 FILE * fout; FILE * cfin;
 fout = fopen_errchk (savfilename, "w");

 switch (module_type)
  {case 1: fprintf(fout,"#!MODULE=kramer\n#<!--mcphase.sipf-->\n");
           fprintf(fout,"#***************************************************************\n");
           fprintf(fout,"# Single Ion Parameter File for Module Kramer for  \n");
           fprintf(fout,"# %s\n",MCPHASVERSION);
           fprintf(fout,"# - program to calculate static magnetic properties\n");
           fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           fprintf(fout,"# %s\n",MCDISPVERSION);
           fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           fprintf(fout,"# %s\n",MCDIFFVERSION);
           fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           fprintf(fout,"#***************************************************************\n#\n");
           fprintf(fout,"#\n# this is a crystal field ground state doublet\n");
           fprintf(fout,"# module, parameters are the following 3 matrix\n# elements\n#\n");
           fprintf(fout,"# A=|<+-|Ja|-+>| B=|<+-|Jb|-+>| C=|<+-|Jc|+->|\n");
           fprintf(fout,"A=%10f \n B=%10f \n C=%10f\n\n",ABC(1),ABC(2),ABC(3));
            
          break;
   case 2: fprintf(fout,"#!MODULE=cfield\n#<!--mcphase.sipf-->\n");
           fprintf(fout,"#***************************************************************\n");
           fprintf(fout,"# Single Ion Parameter File for Module Cfield for  \n");
           fprintf(fout,"# %s\n",MCPHASVERSION);
           fprintf(fout,"# - program to calculate static magnetic properties\n");
           fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           fprintf(fout,"# %s\n",MCDISPVERSION);
           fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           fprintf(fout,"# %s\n",MCDIFFVERSION);
           fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           fprintf(fout,"#***************************************************************\n#\n");
           fprintf(fout,"#\n# crystal field paramerized in Stevens formalism\n#\n");
           (*iops).save(fout);
          break;
   case 3: fprintf(fout,"#!MODULE=brillouin\n#<!--mcphase.sipf-->\n");
           fprintf(fout,"#***************************************************************\n");
           fprintf(fout,"# Single Ion Parameter File for Module Brillouin for\n");
           fprintf(fout,"# %s\n",MCPHASVERSION);
           fprintf(fout,"# - program to calculate static magnetic properties\n");
           fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
           fprintf(fout,"# %s\n",MCDISPVERSION);
           fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
           fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
           fprintf(fout,"# %s\n",MCDIFFVERSION);
           fprintf(fout,"# - program to calculate neutron and magnetic xray diffraction\n");
           fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
           fprintf(fout,"#****************************************************************\n#\n");
           fprintf(fout,"#\n# single ion parameterized by Brillouin function\n");
           fprintf(fout,"# BJ(x) with angular momentum number J=S,\n# no crystal field\n#\n");
           fprintf(fout,"J = %g\n\n",ABC(1));

          break;
   default: // in case of external single ion module just save a copy of the input file 
           char *token;
           cfin=fopen_errchk(cffilename,"rb");
           while(feof(cfin)==false){fgets(instr, MAXNOFCHARINLINE, cfin);
                      // strip /r (dos line feed) from line if necessary
                      while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
                      fprintf(fout,"%s",instr);
                                    }
           fclose(cfin);
   }

  if(module_type>0) // in case of internal modules save common information
   {fprintf(fout,"#----------------\n# Lande factor gJ\n#----------------\nGJ=%g\n\n",gJ);
    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"# Neutron Scattering Length (10^-12 cm) (can be complex)\n");
    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"SCATTERINGLENGTHREAL=%g\nSCATTERINGLENGTHIMAG=%g\n",SLR,SLI);
    fprintf(fout,"#  ... note: - if an occupancy other than 1.0 is needed, just reduce \n");
    fprintf(fout,"#              the scattering length linear accordingly\n\n");

    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"# Debye-Waller Factor: sqr(Intensity)~|sf|~EXP(-2 * DWF *s*s)=EXP (-W)\n");
    fprintf(fout,"#                      with s=sin(theta)/lambda=Q/4pi\n");
    fprintf(fout,"# relation to other notations: 2*DWF=Biso=8 pi^2 <u^2>\n");
    fprintf(fout,"# unit of DWF is [A^2]\n");
    fprintf(fout,"#-------------------------------------------------------\n");
    fprintf(fout,"DWF=%g\n",DWF);

    fprintf(fout,"#--------------------------------------------------------------------------------------\n");
    fprintf(fout,"# Neutron Magnetic Form Factor coefficients - thanks to J Brown\n");
    fprintf(fout,"#   d = 2*pi/Q      \n");
    fprintf(fout,"#   s = 1/2/d = Q/4/pi   \n");   
    fprintf(fout,"#   sin(theta) = lambda * s\n");
    fprintf(fout,"#    s2= s*s = Q*Q/16/pi/pi\n");
    fprintf(fout,"#\n");
    fprintf(fout,"#   <j0(Q)>=   FFj0A*EXP(-FFj0a*s2) + FFj0B*EXP(-FFj0b*s2) + FFj0C*EXP(-FFj0c*s2) + FFj0D\n");
    fprintf(fout,"#   <j2(Q)>=s2*(FFj2A*EXP(-FFj2a*s2) + FFj2B*EXP(-FFj2b*s2) + FFj2C*EXP(-FFj2c*s2) + FFj2D\n");
    fprintf(fout,"#   <j4(Q)>=s2*(FFj4A*EXP(-FFj4a*s2) + FFj4B*EXP(-FFj4b*s2) + FFj4C*EXP(-FFj4c*s2) + FFj4D\n");
    fprintf(fout,"#   <j6(Q)>=s2*(FFj6A*EXP(-FFj6a*s2) + FFj6B*EXP(-FFj6b*s2) + FFj6C*EXP(-FFj6c*s2) + FFj6D\n");
    fprintf(fout,"#\n");
    fprintf(fout,"#   Dipole Approximation for Neutron Magnetic Formfactor:\n");
    fprintf(fout,"#        -Spin Form Factor       FS(Q)=<j0(Q)>\n");
    fprintf(fout,"#        -Angular Form Factor    FL(Q)=<j0(Q)>+<j2(Q)>\n");
    fprintf(fout,"#        -Rare Earth Form Factor F(Q) =<j0(Q)>+<j2(Q)>*(2/gJ-1)\n\n");
    fprintf(fout,"#--------------------------------------------------------------------------------------\n");
    fprintf(fout,"FFj0A=%+7.4f FFj0a=%+7.4f FFj0B=%+7.4f FFj0b=%+7.4f FFj0C=%+7.4f FFj0c=%+7.4f FFj0D=%+7.4f\n",magFFj0[1],magFFj0[2],magFFj0[3],magFFj0[4],magFFj0[5],magFFj0[6],magFFj0[7]);
    fprintf(fout,"FFj2A=%+7.4f FFj2a=%+7.4f FFj2B=%+7.4f FFj2b=%+7.4f FFj2C=%+7.4f FFj2c=%+7.4f FFj2D=%+7.4f\n",magFFj2[1],magFFj2[2],magFFj2[3],magFFj2[4],magFFj2[5],magFFj2[6],magFFj2[7]);
    fprintf(fout,"FFj4A=%+7.4f FFj4a=%+7.4f FFj4B=%+7.4f FFj4b=%+7.4f FFj4C=%+7.4f FFj4c=%+7.4f FFj4D=%+7.4f\n",magFFj4[1],magFFj4[2],magFFj4[3],magFFj4[4],magFFj4[5],magFFj4[6],magFFj4[7]);
    fprintf(fout,"FFj6A=%+7.4f FFj6a=%+7.4f FFj6B=%+7.4f FFj6b=%+7.4f FFj6C=%+7.4f FFj6c=%+7.4f FFj6D=%+7.4f\n",magFFj6[1],magFFj6[2],magFFj6[3],magFFj6[4],magFFj6[5],magFFj6[6],magFFj6[7]);
    fprintf(fout,"\n\n");

  if(abs(Zc)>1e-10){
    fprintf(fout,"#----------------------------------------------------------------------\n");
    fprintf(fout,"# coefficients of Z(K') according to Lovesey (Neutron Scattering) vol.2\n");
    fprintf(fout,"# chapter 11.6.1 page 233: Z(K)= ZKcK-1 * <jK-1(Q)> + ZKcK+1 * <jK+1(Q)>\n");
    fprintf(fout,"#  ... these coefficients are needed to go beyond dipolar approx.\n");
    fprintf(fout,"#      for the neutron magnetic formfactor in rare earth ions\n");
    fprintf(fout,"#----------------------------------------------------------------------\n");
    fprintf(fout,"Z1c0=%+10.8f  Z1c2=%+10.8f\n",Zc(1),Zc(2));
    fprintf(fout,"		  Z3c2=%+10.8f  Z3c4=%+10.8f\n",Zc(3),Zc(4)); 
    fprintf(fout,"				    Z5c4=%+10.8f  Z5c6=%+10.8f\n",Zc(5),Zc(6)); 
    fprintf(fout,"						      Z7c6=%+10.8f\n\n",Zc(7));
                    }

  if(abs(Np)>1e-10){fprintf(fout,"#---------------------------------------------------------------------------------------------------\n");
                    fprintf(fout,"# radial wave function parameters, for transition metal ions the the values are tabulated in\n");
                    fprintf(fout,"# Clementi & Roetti Atomic data and nuclear data tables 14 (1974) 177-478, the radial wave\n");
                    fprintf(fout,"# function is expanded as R(r)=sum_p Cp r^(Np-1) . exp(-XIp r) . (2 XIp)^(Np+0.5) / sqrt(2Np!)\n");
                    fprintf(fout,"# for rare earth ions see Freeman & Watson PR 127(1962)2058, Sovers J. Phys. Chem. Sol. 28(1966)1073\n");
                    fprintf(fout,"#---------------------------------------------------------------------------------------------------\n");
                    for(i=Np.Lo();i<=Np.Hi();++i){if(Np(i)!=0){fprintf(fout,"N%i=%i XI%i=%g C%i=%g\n",i,(int)Np(i),i,Xip(i),i,Cp(i));}
                                                 }
                   fprintf(fout,"\n");
                   }

   }

 fclose(fout);
 delete []savfilename;
}

void jjjpar::get_parameters_from_sipfile(char * sipffilename)
{FILE * cf_file;
 int i,j;
 float nn[MAXNOFNUMBERSINLINE];
 nn[0]=MAXNOFNUMBERSINLINE;
  char modulefilename[MAXNOFCHARINLINE];

 char instr[MAXNOFCHARINLINE];
  cf_file = fopen_errchk (sipffilename, "rb");
  fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
  if(extract(instr,"MODULE=",modulefilename,(size_t)MAXNOFCHARINLINE))
   {if(extract(instr,"#!",modulefilename,(size_t)MAXNOFCHARINLINE))
    {fprintf(stderr,"Error: single ion property file %s does not start with '#!' or 'MODULE='\n",sipffilename);
     exit(EXIT_FAILURE);}
   }

  fprintf (stderr,"#parsing single ion property file: %s - loading module %s\n",sipffilename,modulefilename);

  if(strcmp(modulefilename,"kramer")==0)
    {module_type=1;fprintf (stderr,"[internal]\n");
      ABC=Vector(1,3);i=3;
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
                                           i+=extract(instr,"A",ABC(1))-1; 
                                           i+=extract(instr,"B",ABC(2))-1; 
                                           i+=extract(instr,"C",ABC(3))-1;   
                                          }
      }
      // input all  lines starting with comments
      if(i!=0){fprintf(stderr,"Error reading |<+-|Ja|-+>|,|<+-|Jb|-+>|,|<+-|Jc|+->| from file %s\ncorrect file format is:\n",sipffilename);
              fprintf(stderr,"\nMODULE=kramer\n#comment lines ..\n#matrix elements A=|<+-|Ja|-+>| B=|<+-|Jb|-+>| C=|<+-|Jc|+->|\nA=2 \nB=3 \nC=1\n\n");exit(EXIT_FAILURE);}
      // now we have the numbers corresponding to the vector ABC() in nn[]
      fprintf(stderr," ... kramers doublet with A=<+|Ja|->=%g B=<+-|Jb|+->=+-%g C=<+|Jc|->/i=%g\n",ABC(1),ABC(2),ABC(3));
      est=ComplexMatrix(0,2,1,2); // not used, just initialize to prevent errors
    }
  else
    {if(strcmp(modulefilename,"brillouin")==0)
     {module_type=3;fprintf (stderr,"[internal]\n");
      ABC=Vector(1,1);i=1;
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
                                           i+=extract(instr,"J",ABC(1))-1; 
                                          }
      }// input all  lines starting with comments
      if(i!=0){fprintf(stderr,"Error reading spin quantum number J=S from file %s\ncorrect file format is:\n",sipffilename);
              fprintf(stderr,"\n#!brillouin\n#comment lines ..\n# Quantum number  J\nJ=3.5\n\n");exit(EXIT_FAILURE);}
      // now we have the numbers corresponding to the vector ABC() in nn[]
      fprintf(stderr," ... Brillouin function with J=S=%g\n",ABC(1));
      est=ComplexMatrix(0,2,1,2);// not used, just initialize to prevent errors
      est=0;
     }
     else 
     {if(strcmp(modulefilename,"cfield")==0)
     {module_type=2;fprintf (stderr,"#[internal]\n");
      fclose(cf_file);cf_file = fopen_errchk (sipffilename, "rb"); // reopen file
      iops=new ionpars(cf_file);  
      int dj;dj=(int)(2*J()+1);
      est=ComplexMatrix(0,dj,1,dj);
     // get 1ion parameters - operator matrices
     
     }
     else
     {fprintf (stderr,"#[external]\n");
      i=0;
      while(feof(cf_file)==false)
      {fgets(instr, MAXNOFCHARINLINE, cf_file);
       if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
                                           i-=extract(instr,"MODPAR1",nn[1])-1; 
                                           i-=extract(instr,"MODPAR2",nn[2])-1; 
                                           i-=extract(instr,"MODPAR3",nn[3])-1; 
                                           i-=extract(instr,"MODPAR4",nn[4])-1; 
                                           i-=extract(instr,"MODPAR5",nn[5])-1; 
                                           i-=extract(instr,"MODPAR6",nn[6])-1; 
                                           i-=extract(instr,"MODPAR7",nn[7])-1; 
                                           i-=extract(instr,"MODPAR8",nn[8])-1; 
                                           i-=extract(instr,"MODPAR9",nn[9])-1; 
                                          }
      }
       // input all  lines starting with comments
    //while((i=inputparline ("params",cf_file, nn))==0&&feof(cf_file)==false);
    // now we have the numbers corresponding to vector ABC() in nn[] - these are the module parameters !
    fprintf(stderr,"#parameters: ");
    if(i>0){
             ABC=Vector(1,i);for(j=1;j<=i;++j){ABC(j)=nn[j];fprintf(stderr,"%g ",nn[j]);} 
            }else{
             ABC=Vector(1,1);
	    } 
    fprintf(stderr,"\n");

  char * error;module_type=0;
#ifdef __linux__
  handle=dlopen (modulefilename,RTLD_NOW | RTLD_GLOBAL);
  if (!handle){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n");
               if ((error=dlerror())!=NULL) 
	         {fprintf (stderr,"%s\n",error);}
	       exit (EXIT_FAILURE);
	      }
  m=(void(*)(Vector*,double*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*))dlsym(handle,"mcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s\n",error);exit (EXIT_FAILURE);}
  dm=(int(*)(int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*,float*,ComplexMatrix*))dlsym(handle,"dmcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dm=NULL;}
  mq=(void(*)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*))dlsym(handle,"mq");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);mq=NULL;}
  estates=(void(*)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**))dlsym(handle,"estates");
  
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);estates=NULL;
                                est=ComplexMatrix(0,2,1,2);// not used, just initialize to prevent errors
                                est=0;
                               }

  ddnn=(int(*)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexMatrix*))dlsym(handle,"dncalc");
  
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);ddnn=NULL;}


#else
  handle=LoadLibrary(modulefilename);
  if ((int)handle<= HINSTANCE_ERROR){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n");
	       exit (EXIT_FAILURE);
	      }
  
    m=(void(*)(Vector*,double*,Vector*,double*,Vector*,char**,double*,double*,ComplexMatrix*))GetProcAddress(handle,"mcalc");
     if (m==NULL) {fprintf (stderr,"jjjpar::jjjpar error %d  module %s loading function mcalc not possible\n",GetLastError(),modulefilename);exit (EXIT_FAILURE);}
    dm=(int(*)(int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*,float*,ComplexMatrix*))GetProcAddress(handle,"dmcalc");
     if (dm==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d module %s loading function dmcalc not possible - continuing\n",GetLastError(),modulefilename);}
    mq=(void(*)(ComplexVector*,double*,double*,double*,double*,double*,double*,ComplexMatrix*))GetProcAddress(handle,"mq");
     if (mq==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function mq not possible - continuing\n",GetLastError(),modulefilename);}
    estates=(void(*)(ComplexMatrix*,Vector*,double*,double*,Vector*,char**))GetProcAddress(handle,"estates");
     if (estates==NULL) {fprintf (stderr,"jjjpar::jjjpar warning %d  module %s loading function estates not possible - continuing\n",GetLastError(),modulefilename);
                                est=ComplexMatrix(0,2,1,2);// not used, just initialize to prevent errors
                                est=0;
                               }

  ddnn=(int(*)(int*,double*,double*,double*,double*,double*,double*,ComplexMatrix*,double*,ComplexMatrix*))GetProcAddress(handle,"dncalc");
     if (ddnn==NULL) {fprintf (stderr,"jjjpar::jjjpar warning  %d  module %s loading function dncalc not possible - continuing\n",GetLastError(),modulefilename);}
  
#endif

    }
   }
  }
  fclose(cf_file);

  magFFj0=Vector(1,7);magFFj0=0;  magFFj0[1]=1;
  magFFj2=Vector(1,7);magFFj2=0;
  magFFj4=Vector(1,7);magFFj4=0;
  magFFj6=Vector(1,7);magFFj6=0;
  Zc=Vector(1,7);Zc=0;

   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;

  DWF=0;  gJ=0;
   calcmagdensity=0;

  cf_file = fopen_errchk (sipffilename, "rb");

  while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
    extract(instr,"SCATTERINGLENGTHREAL",SLR);  
    extract(instr,"SCATTERINGLENGTHIMAG",SLI);  
    extract(instr,"GJ",gJ);  
    extract(instr,"gJ",gJ);  
    extract(instr,"calcmagdensity",calcmagdensity);
    // read formfactor if given
    extract(instr,"FFj0A",magFFj0[1]);
    extract(instr,"FFj0a",magFFj0[2]);
    extract(instr,"FFj0B",magFFj0[3]);
    extract(instr,"FFj0b",magFFj0[4]);
    extract(instr,"FFj0C",magFFj0[5]);
    extract(instr,"FFj0c",magFFj0[6]);
    extract(instr,"FFj0D",magFFj0[7]);
    extract(instr,"FFj2A",magFFj2[1]);
    extract(instr,"FFj2a",magFFj2[2]);
    extract(instr,"FFj2B",magFFj2[3]);
    extract(instr,"FFj2b",magFFj2[4]);
    extract(instr,"FFj2C",magFFj2[5]);
    extract(instr,"FFj2c",magFFj2[6]);
    extract(instr,"FFj2D",magFFj2[7]);
    extract(instr,"FFj4A",magFFj4[1]);
    extract(instr,"FFj4a",magFFj4[2]);
    extract(instr,"FFj4B",magFFj4[3]);
    extract(instr,"FFj4b",magFFj4[4]);
    extract(instr,"FFj4C",magFFj4[5]);
    extract(instr,"FFj4c",magFFj4[6]);
    extract(instr,"FFj4D",magFFj4[7]);
    extract(instr,"FFj6A",magFFj6[1]);
    extract(instr,"FFj6a",magFFj6[2]);
    extract(instr,"FFj6B",magFFj6[3]);
    extract(instr,"FFj6b",magFFj6[4]);
    extract(instr,"FFj6C",magFFj6[5]);
    extract(instr,"FFj6c",magFFj6[6]);
    extract(instr,"FFj6D",magFFj6[7]);
   // coefficients of Z(K') according to Lovesey chapter 11.6.1 page 233
    extract(instr,"Z1c0",Zc(1));  
    extract(instr,"Z1c2",Zc(2));  
    extract(instr,"Z3c2",Zc(3));  
    extract(instr,"Z3c4",Zc(4));  
    extract(instr,"Z5c4",Zc(5));  
    extract(instr,"Z5c6",Zc(6));  
    extract(instr,"Z7c6",Zc(7));  
   // read debeywallerfactor if given    
    extract(instr,"DWF",DWF);
   // read radial wavefunction parameters
        extract(instr,"N1",Np(1));extract(instr,"XI1",Xip(1));extract(instr,"C1",Cp(1));
        extract(instr,"N2",Np(2));extract(instr,"XI2",Xip(2));extract(instr,"C2",Cp(2));
        extract(instr,"N3",Np(3));extract(instr,"XI3",Xip(3));extract(instr,"C3",Cp(3));
        extract(instr,"N4",Np(4));extract(instr,"XI4",Xip(4));extract(instr,"C4",Cp(4));
        extract(instr,"N5",Np(5));extract(instr,"XI5",Xip(5));extract(instr,"C5",Cp(5));
        extract(instr,"N6",Np(6));extract(instr,"XI6",Xip(6));extract(instr,"C6",Cp(6));
        extract(instr,"N7",Np(7));extract(instr,"XI7",Xip(7));extract(instr,"C7",Cp(7));
        extract(instr,"N8",Np(8));extract(instr,"XI8",Xip(8));extract(instr,"C8",Cp(8));
        extract(instr,"N9",Np(9));extract(instr,"XI9",Xip(9));extract(instr,"C9",Cp(9));
  }    
 }
 
 fclose (cf_file); 
// check gJ
if(module_type==2&&fabs(gJ-(*iops).gJ)>0.00001)
{fprintf(stderr,"Error internal module cfield : Lande Factor read from %s (gJ=%g) does not conform to internal module value gJ=%g\n",sipffilename,gJ,(*iops).gJ);exit(EXIT_FAILURE);}
if (gJ==0){printf("# reading gJ=0 in single ion property file %s -> entering intermediate coupling mode by assigning Ja=Sa Jb=La Jc=Sb Jd=Lb Je=Sc Jf=Lc (S... Spin, L... angular momentum)\n",sipffilename);
           if (module_type==1){fprintf(stderr,"Error internal module kramers: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
           if (module_type==2){fprintf(stderr,"Error internal module cfield : intermediate coupling not supported\n");exit(EXIT_FAILURE);}
           if (module_type==3){fprintf(stderr,"Error internal module brillouin: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
          }

}

/*****************************************************************************************/
//constructor with file handle of mcphas.j
jjjpar::jjjpar(FILE * file,int nofcomps) 
{ FILE * cf_file;   
  char instr[MAXNOFCHARINLINE];
  cffilename= new char [MAXNOFCHARINLINE];
  int i,j,i1,j1,k1,l;
  double gjcheck;
  float nn[MAXNOFNUMBERSINLINE];
  nn[0]=MAXNOFNUMBERSINLINE;
  xyz=Vector(1,3);
  set_zlm_constants();
  i=7;
  while(i>0){fgets_errchk (instr, MAXNOFCHARINLINE, file);
             if(instr[strspn(instr," \t")]!='#'){fprintf (stderr, "Error reading mcphas.j - exchangeparameters start before all variables (da,db,dc,gJ,nofneighbors,diagonalexchange and cffilename) have been given\n");
                                                 exit (EXIT_FAILURE);}
             i+=extract(instr,"x",xyz[1])-1;
             i+=extract(instr,"y",xyz[2])-1;
             i+=extract(instr,"z",xyz[3])-1;
             i+=extract(instr,"da",xyz[1])-1;
             i+=extract(instr,"db",xyz[2])-1;
             i+=extract(instr,"dc",xyz[3])-1;
             i+=extract(instr,"nofneighbours",paranz)-1;
             i+=extract(instr,"diagonalexchange",diagonalexchange)-1;
             i+=extract(instr,"gJ",gjcheck)-1;
             i+=extract(instr,"cffilename",cffilename,(size_t)MAXNOFCHARINLINE)-1;
            }

  fgets_errchk (instr, MAXNOFCHARINLINE, file);

// read single ion parameter file and see which type it is (internal module or loadable)
  transitionnumber=1;
  
  //start reading again at the beginning of the file to get formfactors, debye waller factor
  get_parameters_from_sipfile(cffilename);
  if (gJ!=gjcheck){fprintf (stderr, "Error: Lande factor gJ in file mcphas.j and %s are not the same\n",cffilename);
                   exit (EXIT_FAILURE);}
  Mq=ComplexVector(1,3);

  nofcomponents=nofcomps; // default value for nofcomponents - (important in case nofparameters=0)
// read the exchange parameters from file (exactly paranz parameters!)
  for  (i=1;i<=paranz;++i)
  {while((j=inputline(file, nn))==0&&feof(file)==0){}; // returns 0 if comment line or eof, exits with error, if input string too long
   if(feof(file)!=0){ fprintf (stderr, "Error in jjjpar.cpp: input jjj parameters - \n");
                      fprintf(stderr," end of file reached while reading exchange parameter %i(%i)",i,paranz);
                      exit (EXIT_FAILURE);
                    }
    if(i==1){// determine nofcomponents from number of parameters read in first line of mcphas.j
             if(diagonalexchange==1){nofcomponents=j-3;}else{nofcomponents=(int)sqrt((double)(j-3));}
             if(module_type==1)
	     {// check dimensions of vector if internal kramers is used
              if(nofcomponents!=3)
              {fprintf(stderr,"Error reading mcphas.j: number of dimensions (not equal 3) not compatible with internal single ion module kramer - check number of columns in file mcphas.j\n");
               exit(EXIT_FAILURE);}
             }
             // dimension arrays
             dn = new Vector[paranz+1];for(i1=0;i1<=paranz;++i1){dn[i1]=Vector(1,3);}
             if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
             sublattice = new int[paranz+1];if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
             jij = new Matrix[paranz+1];for(i1=0;i1<=paranz;++i1){jij[i1]=Matrix(1,nofcomponents,1,nofcomponents);}
             if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
            }
   //check if correct number of columns has been read	        
    if((diagonalexchange==1&&nofcomponents!=j-3)||(diagonalexchange==0&&nofcomponents!=(int)sqrt((double)(j-3))))
              {fprintf(stderr,"Error reading mcphas.j line %i: check number of columns\n",i);
               exit(EXIT_FAILURE);}

  mom=Vector(1,nofcomponents); 

   //(1-3) give the absolute coordinates of the neighbour and are transformed here to relative
   // coordinates !!
   dn[i](1) = nn[1];dn[i](2) = nn[2];dn[i](3) = nn[3];
   jij[i]=0;

  // format of matrix 
  // 11 22 33 12 21 13 31 23 32 (3x3 matrix)
  // 11 22 33 44 12 21 13 31 14 41 23 32 24 42 34 43 (4x4 matrix)
  // 11 22 33 44 55 12 21 13 31 14 41 15 51 23 32 24 42 25 52 34 43 35 53 45 54 (5x5 matrix)
  // etc ...
  //read diagonal components of exchange matrix
  for(i1=1;i1<=nofcomponents;++i1){jij[i](i1,i1)= nn[i1+3];}
  //read off-diagonal components of exchange matrix (if required)
  if (diagonalexchange==0){k1=3+nofcomponents;
                           for(i1=1;i1<=nofcomponents-1;++i1)
                              {for(j1=i1+1;j1<=nofcomponents;++j1)
                               {++k1;jij[i](i1,j1)= nn[k1];
			        ++k1;jij[i](j1,i1)= nn[k1];
			       }
			      }
                          }
  }
}

// constructor with filename of singleion parameter  used by mcdiff and charges-chargeplot
jjjpar::jjjpar(double x,double y,double z, char * sipffile)
{xyz=Vector(1,3);xyz(1)=x;xyz(2)=y;xyz(3)=z;
  mom=Vector(1,9); mom=0; 
  Mq=ComplexVector(1,3);
  cffilename= new char [MAXNOFCHARINLINE];
  strcpy(cffilename,sipffile);
  get_parameters_from_sipfile(cffilename);
  set_zlm_constants();

}

// constructor with positions scattering length dwf
jjjpar::jjjpar(double x,double y,double z, double slr,double sli, double dwf)
{xyz=Vector(1,3);xyz(1)=x;xyz(2)=y;xyz(3)=z;
 mom=Vector(1,9); mom=0; 
 DWF=dwf;SLR=slr;SLI=sli;
  magFFj0=Vector(1,7);magFFj0=0;  magFFj0[1]=1;
  magFFj2=Vector(1,7);magFFj2=0;
  magFFj4=Vector(1,7);magFFj4=0;
  magFFj6=Vector(1,7);magFFj6=0;
  Zc=Vector(1,7);Zc=0;
  set_zlm_constants();
   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;
   r2=0;r4=0;r6=0;
  calcmagdensity=0;

}
//constructor without file
jjjpar::jjjpar(int n,int diag,int nofmom) 
{ cffilename= new char [MAXNOFCHARINLINE];
  diagonalexchange=diag;
  paranz=n;xyz=Vector(1,3);xyz=0;
  set_zlm_constants();
  int i1;r2=0;r4=0;r6=0;
  module_type=1;ABC=Vector(1,3);ABC=0;
  transitionnumber=1;
  nofcomponents=nofmom;
  mom=Vector(1,nofcomponents);
  mom=0;
  calcmagdensity=0;

  dn = new Vector[n+1];for(i1=0;i1<=n;++i1){dn[i1]=Vector(1,3);}
  if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  jij = new Matrix[n+1];for(i1=0;i1<=n;++i1){jij[i1]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  magFFj0=Vector(1,7);magFFj0=0;  magFFj0[1]=1;
  magFFj2=Vector(1,7);magFFj2=0;
  magFFj4=Vector(1,7);magFFj4=0;
  magFFj6=Vector(1,7);magFFj6=0;
  Zc=Vector(1,7);Zc=0;
   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;
  DWF=0;gJ=0;

}

//copy constructor
jjjpar::jjjpar (const jjjpar & p)
{ int i;
  xyz=Vector(1,3);
  nofcomponents=p.nofcomponents;
  mom=Vector(1,nofcomponents); 
  xyz=p.xyz;paranz=p.paranz;
  set_zlm_constants();
  SLR=p.SLR;SLI=p.SLI;

  diagonalexchange=p.diagonalexchange;
  gJ=p.gJ;module_type=p.module_type;
  Mq=ComplexVector(1,3);
  Mq=p.Mq;
   calcmagdensity=p.calcmagdensity;

  Np=p.Np; Xip=p.Xip;Cp=p.Cp;
  r2=p.r2;r4=p.r4;r6=p.r6;

  transitionnumber=p.transitionnumber;
  cffilename= new char [strlen(p.cffilename)+1];
  strcpy(cffilename,p.cffilename);
  if (p.module_type==1||p.module_type==0)  ABC=p.ABC;
  if (p.module_type==2)  {iops=new ionpars((int)(2*(*p.iops).J+1));iops=p.iops;
                           int dj;dj=(int)(2*J()+1);
                           est=ComplexMatrix(0,dj,1,dj);est=p.est;
                           }
//  if (module_type==2)  iops=new ionpars(4);iops=p.iops;
//  if (module_type==2)  iops=p.iops;
  
#ifdef __linux__
/*  if (module_type==0)
  {char * error;
   handle=dlopen (cffilename,RTLD_NOW | RTLD_GLOBAL);
   if (!handle){fprintf (stderr, "Could not load dynamic library\n");
               if ((error=dlerror())!=NULL) 
	         {fprintf (stderr,"%s\n",error);}
	       exit (EXIT_FAILURE);
	      }
*/
   m=p.m;
   dm=p.dm;
   mq=p.mq;
   ddnn=p.ddnn;
   estates=p.estates;
/*  }*/
#endif
  magFFj0=Vector(1,7);magFFj0=p.magFFj0;
  magFFj2=Vector(1,7);magFFj2=p.magFFj2;
  magFFj4=Vector(1,7);magFFj4=p.magFFj4;
  magFFj6=Vector(1,7);magFFj6=p.magFFj6;
  Zc=Vector(1,7);Zc=p.Zc;
  DWF=p.DWF;  
int i1;
//dimension arrays
  jij = new Matrix[paranz+1];for(i1=0;i1<=paranz;++i1){jij[i1]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  dn = new Vector[paranz+1];for(i1=0;i1<=paranz;++i1){dn[i1]=Vector(1,3);}
  if (dn == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  for (i=1;i<=paranz;++i)
  {jij[i]=p.jij[i];dn[i]=p.dn[i];sublattice[i]=p.sublattice[i];}
}


//destruktor
jjjpar::~jjjpar ()
{ //delete []jij; //will not work in linux 
  //delete []dn;  // will not work in linux
  //delete []sublattice;
  //delete []cffilename;// will not work in linux
  if (module_type==2) delete iops;
#ifdef __linux__
   if (module_type==0)dlclose(handle);
#endif
}


