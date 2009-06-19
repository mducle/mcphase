/***********************************************************************
 *
 * mcdiff.c - program to calculate neutron diffraction pattern
 *
 ***********************************************************************/

#define PI 3.141592654
#define KB 0.0862     // Boltzmanns constant in mev/K
#define MAXNOFCHARINLINE 1000
#define MAXNOFREFLECTIONS 10000
#define SMALL 1e-8
#include <mpspecfunp.h>
#include <martin.h>
#include <myev.h>
#include <jjjpar.hpp>
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<ctime>
#include<cerrno>
#include<cstring>
#include<cmath>
#include<vector.h>
#include <complex>
#include<cstddef>


/**********************************************************************/
void xproduct(Vector & result,Vector a, Vector b)
{
 result(1)=a(2)*b(3)-a(3)*b(2);
 result(2)=-a(1)*b(3)+a(3)*b(1);
 result(3)=a(1)*b(2)-a(2)*b(1);
 
 return ;
}

void rezcalc(Vector r1,Vector  r2,Vector  r3,Vector  rez1,Vector  rez2,Vector  rez3)
{// calculate reciprocal lattice rezi from real lattice ri
 float vol;
 xproduct(rez1,r2,r3); vol=rez1*r1; rez1*=2.0*PI/vol;
 xproduct(rez2,r1,r3); vol=rez2*r2; rez2*=2.0*PI/vol;
 xproduct(rez3,r1,r2); vol=rez3*r3; rez3*=2.0*PI/vol;
return;}

//double ZZ(int K, float J0, float J2, float J4, float J6, Vector Zc)
//{// calculate Z(K)
// if (K==1) return Zc(1)*J0+Zc(2)*J2;
// if (K==3) return Zc(3)*J2+Zc(4)*J4;
// if (K==5) return Zc(5)*J4+Zc(6)*J6;
// if (K==7) return Zc(7)*J6;
// 
// return 0;
//}

int getint(jjjpar ** jjjpars,int hi,int ki,int li,float thetamax,Vector rez1,Vector rez2, Vector rez3,float scale,double T,float lambda,float ovalltemp,int lorenz,int & n,int * J,float & d,float & Theta,float & Imag,float & Imagdip,float & inuc,float & SF,float & lorentzf,complex <double> & mqx,complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2)
{
// this routine calculates the intensity of elastic neutrons for a reflection (hi ki li)
//
// input:
// (*jjjpar[1...n]).xyz(1..3)         atomic positional parameters dr1 dr2 dr3
//                                        '(with respect to primitive lattice)
// (*jjjpars[1...n]).DWF              debye waller factors [A^2]
// (*jjjpars[1...n]).SLR,SLI          nuclear scattering length[10^-12cm]
// (*jjjpars[1...n]).mom(1..3)(45)(67)(89)        atomic magnetic moment Ma Mb Mc [mb] and (if input) Sa La Sb Lb Sc Lc
//                                       ' (with respect to coordinates 1,2,3=yzx)
// (*jjjpars[1...n]).gj		      Lande factor
// J[1..n] // code for indicating if ion is nonmagnetic (J=1), 
            //rare earth with dipole approx (J=-1), 
            //rare earth beyond dipole approx, but with given nonzero gJ (stevens-balcar formalism) (J=0),
            //gJ=0,general L and S moments given, use dipole approximation and separate formfactor for spin and orbital moment (J=-2)
            //intermediate coupling (gJ=0), go beyond dipole approximation (J=-3)
// (*jjjpars[1...n]).magFFj0(1..7)         formfactor j0 for atom 1...n <j0(kr)>-terms A,a,B,b,C,c,D
// (*jjjpars[1...n]).magFFj2(1..7)         formfactor j2 for atom 1...n <j2(kr)>-terms A,a,B,b,C,c,D
//     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
//     where U(r) is the Radial wave function for the unpaired electrons in the atom
// (*jjjpars[1...n]).magFFj4(1..7)         formfactor j4 for atom 1...n  (needed to go beyond dipole approx)
// (*jjjpars[1...n]).magFFj6(1..7)         formfactor j6 for atom 1...n  (needed to go beyond dipole approx)
// (*jjjpars[1...n]).Zc		         Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
// (*jjjpars[1...n]).eigenstates(1..2J+1,1..2J+1)   CF+MF eigenstates (needed to go beyond dipole approx)
// thetamax  			      maximum theta value, if theta larger, routine returns false
// rez1,rez2,rez3                     vectors of reciprocal lattice 
// scale                              scaling factor for intensity
// T				         temperature [K] (needed to go beyond dipole approx)
// lambda                             wavelength[A]
// ovalltemp                          overall temperature factor [A^2]
// lorenz                             code for lorentzfactor to be used
// n                                  number of atoms per unit cell
//
// output:
// d                                  d spacing in A
// theta                              scattering angle
// Imag, Imagdip, inuc                scattering intensities
// SF                                 structure factor
// lorentzf                           Lorentz Factor
// mx,my,mz,mxmy,mxmz,mymz,mx2my2mz2[].. fouriertransform of momentunitvectors (for mag xray scattering)


            double s,Q,FQ,FQL,sintheta,qr,sin2theta,ovallt,mux,muy,muz;
            int i,j;
            Vector Qvec(1,3);
            //calculate d spacing and intensity for h,k,l triple (d,intmag,ikern)************
            Qvec=rez1*(double)(hi) + rez2*(double)(ki)  + rez3*(double)(li) ;
            Q = Norm(Qvec); //dspacing
            d = 2.0 * PI / Q;
            s=0.5 / d; 
	    sintheta = lambda * s;
            if (sintheta >= sin(thetamax / 180 * PI)) return false;
               Theta = 180 / PI * atan(sintheta / sqrt(1 - sintheta * sintheta));
               //nuclear(|nsfr+i nsfc|^2) and magnetic structure factor(msf) calculation
               complex <double> nsf=0;
               complex <double> msfx=0,msfdipx=0;
               complex <double> msfy=0,msfdipy=0;
               complex <double> msfz=0,msfdipz=0;
               complex <double> im(0,1);

               for(i=1;i<=n;++i){
                                 complex <double> scl((*jjjpars[i]).SLR,(*jjjpars[i]).SLI); 
                                 qr=hi*(*jjjpars[i]).xyz(1)+ki*(*jjjpars[i]).xyz(2)+li*(*jjjpars[i]).xyz(3);

                                 //nuclear structure factor nsfr,nsfc
                                 nsf+=scl*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);

                                //magnetic structure factors
                                if(J[i]<=0){   // i.e. atom is magnetic
                                             // formfactor F(Q)
                                             FQ = (*jjjpars[i]).F(Q); //rare earth

                                             if(J[i]==0){ // go beyond dipole approximation for rare earth
                                                         ComplexVector MQ(1,3);MQ=(*jjjpars[i]).MQ(Qvec);
					               msfx+=0.5*MQ(1)*exp(-2*PI*qr*im);//MQ(123)=MQ(xyz)
					               msfy+=0.5*MQ(2)*exp(-2*PI*qr*im);
					               msfz+=0.5*MQ(3)*exp(-2*PI*qr*im);
                                                       msfdipx+=(*jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im);// mom(123)=mom(abc)=mom(yzx)
					               msfdipy+=(*jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im);
					               msfdipz+=(*jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im);
					                  }
 					      if(J[i]==-1){// dipole approximation - use magnetic moments and rare earth formfactor
                                                           //                        for transition metals always set gJ=2 (spin only moment)
                                                        msfx+=(*jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im);
					                msfy+=(*jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im);
					                msfz+=(*jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im);
                                                        msfdipx+=(*jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im);
					                msfdipy+=(*jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im);
					                msfdipz+=(*jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im);
					               }
					      if(J[i]==-2){// dipole approximation - use S and L moments (only if gJ=0)
                                                        FQL = (*jjjpars[i]).F(-Q); // orbital formfactor
                                                        msfx+=(*jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im); // spin FF
					                msfy+=(*jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im);
					                msfz+=(*jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im);
					                msfx+=(*jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im); // orbital FF
					                msfy+=(*jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im);
					                msfz+=(*jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im);
                                                        msfdipx+=(*jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im); // spin FF
					                msfdipy+=(*jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im);
					                msfdipz+=(*jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im);
					                msfdipx+=(*jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im); // orbital FF
					                msfdipy+=(*jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im);
					                msfdipz+=(*jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im);
					               }
                                     if(J[i]==-3){ // go beyond dipole approximation for gJ=0 (intermediate coupling)
                                                       ComplexVector MQ(1,3);MQ=(*jjjpars[i]).MQ(Qvec);
//                                             printf("MQxyz=(%g %+g i, %g %+g i,%g %+g i)",real(MQ(1)),imag(MQ(1)),real(MQ(2)),imag(MQ(2)),real(MQ(3)),imag(MQ(3)));
					               msfx+=0.5*MQ(1)*exp(-2*PI*qr*im);//MQ(123)=MQ(xyz)
					               msfy+=0.5*MQ(2)*exp(-2*PI*qr*im);
					               msfz+=0.5*MQ(3)*exp(-2*PI*qr*im);
                                             FQL = (*jjjpars[i]).F(-Q); // orbital formfactor
                                              msfdipx+=(*jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im); // spin FF
					                msfdipy+=(*jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im);
					                msfdipz+=(*jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im);
					                msfdipx+=(*jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im); // orbital FF
					                msfdipy+=(*jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im);
					                msfdipz+=(*jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im);
					               }
// myPrintVector(stdout,(*jjjpars[i]).mom);//equivalent to moment ...
                                          
                                                         //mux=(*jjjpars[i]).mom(3); // this should be done in future to implement
                                                         //muy=(*jjjpars[i]).mom(1); // also nonortholattices in corr functions and res mag scatt
                                                         //muz=(*jjjpars[i]).mom(2);
                                                         mux=(*jjjpars[i]).mom(1); // this is still here because correlation functions are calculated
                                                         muy=(*jjjpars[i]).mom(2); // only for orhtogonal lattices (see printeln sub) and so we take 
                                                         muz=(*jjjpars[i]).mom(3); // the old convention of the mcdiff program (a||x,b||y,c||z)
                                                         mqx+=mux*exp(-2*PI*qr*im);
                                                         mqy+=muy*exp(-2*PI*qr*im);
                                                         mqz+=muz*exp(-2*PI*qr*im);
                                                         mqx2+=mux*mux*exp(-2*PI*qr*im);
                                                         mqy2+=muy*muy*exp(-2*PI*qr*im);
                                                         mqz2+=muz*muz*exp(-2*PI*qr*im);
                                                         mqxy+=mux*muy*exp(-2*PI*qr*im);
                                                         mqxz+=mux*muz*exp(-2*PI*qr*im);
                                                         mqyz+=muy*muz*exp(-2*PI*qr*im);
                             }
                                }

             //magnetic structure factors + polarisation factor===>msf
             float msf,msfdip;
             msf = norm(msfx)+norm(msfy)+norm(msfz);
             msfdip = norm(msfdipx)+norm(msfdipy)+norm(msfdipz);

            msf -=  2 * Qvec(1) * Qvec(2) / Q / Q * (real(msfx) * real(msfy) + imag(msfx) * imag(msfy));
            msf -=  2 * Qvec(1) * Qvec(3) / Q / Q * (real(msfx) * real(msfz) + imag(msfx) * imag(msfz));
            msf -=  2 * Qvec(2) * Qvec(3) / Q / Q * (real(msfy) * real(msfz) + imag(msfy) * imag(msfz));

            msf -=  Qvec(1) * Qvec(1) / Q / Q * norm(msfx);
            msf -=  Qvec(2) * Qvec(2) / Q / Q * norm(msfy);
            msf -=  Qvec(3) * Qvec(3) / Q / Q * norm(msfz);

            msfdip -=  2 * Qvec(1) * Qvec(2) / Q / Q * (real(msfdipx) * real(msfdipy) + imag(msfdipx) * imag(msfdipy));
            msfdip -=  2 * Qvec(1) * Qvec(3) / Q / Q * (real(msfdipx) * real(msfdipz) + imag(msfdipx) * imag(msfdipz));
            msfdip -=  2 * Qvec(2) * Qvec(3) / Q / Q * (real(msfdipy) * real(msfdipz) + imag(msfdipy) * imag(msfdipz));

            msfdip -=  Qvec(1) * Qvec(1) / Q / Q * norm(msfdipx);
            msfdip -=  Qvec(2) * Qvec(2) / Q / Q * norm(msfdipy);
            msfdip -=  Qvec(3) * Qvec(3) / Q / Q * norm(msfdipz);


            //lorentzfactor*************************************************************
            sin2theta = 2.0 * sintheta * sqrt(1.0 - sintheta * sintheta);
            if(lorenz == 0){lorentzf = 100;} // no lorentzfactor
            if(lorenz == 1){lorentzf = 1.0 / sin2theta / sin2theta;} // powder flat sample
            if(lorenz == 2){lorentzf = 1.0 / sin2theta / sintheta;}  // powder cyl. sample
            if(lorenz == 3){lorentzf = 1.0 / sin2theta;}             //single crystal
            if(lorenz == 4){lorentzf = d * d * d;}      //TOF powder cyl sample... log scaled d-pattern
            if(lorenz == 5){lorentzf = d * d * d * d;}  //TOF powder cyl sample... d-pattern

             //overall temperature factor*************************************************
             ovallt = exp(-2 * ovalltemp * (sintheta * sintheta / lambda / lambda));
             //***************************************************************************

             //A)nuclear intenisty
            SF = abs(nsf);
            inuc = SF * SF * lorentzf * scale * ovallt;

             //B)magnetic intensity
            Imag = msf * 3.65 / 4 / PI * lorentzf * scale * ovallt;
            Imagdip = msfdip * 3.65 / 4 / PI * lorentzf * scale * ovallt;
return true;
}

void neutint(jjjpar ** jjjpars,int code,double T,float lambda, float thetamax, float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int & n,int * J,int & m,Vector *  hkl,float * D,float * theta,float * intmag,float * intmagdip,float * ikern,float * sf,float * lpg,complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2)
{//****************************************************************************
// this routine calculates the intensity of elastic neutrons
// for a given magnetic unit cell (crystal axis orthogonal)
// the magnetic scattering is treated in the dipole approximation
// input:
// code                               governs if a list of hkl given in hkl[] or all hkls should be generated
// lambda                             wavelength[A]
// ovalltemp                          overall temperature factor [A^2]
// r1(1..3),r2(),r3()                 vectors of primitive unit cell[A]
// n                                  number of atoms per unit cell
// (*jjjpar[1...n]).xyz(1..3)         atomic positional parameters dr1 dr2 dr3
//                                        '(with respect to primitive lattice)
// (*jjjpars[1...n]).DWF              debye waller factors [A^2]
// (*jjjpars[1...n]).SLR,SLI          nuclear scattering length[10^-12cm]
// (*jjjpars[1...n]).mom(1..3)(45)(67)(89)        atomic magnetic moment Ma Mb Mc [mb] and (if input) Sa La Sb Lb Sc Lc
//                                       ' (with respect to coordinates 1,2,3=yzx)
// (*jjjpars[1...n]).gj		      Lande factor
// J[1..n] // code for indicating if ion is nonmagnetic (J=1), 
            //rare earth with dipole approx (J=-1), 
            //rare earth beyond dipole approx, but with given nonzero gJ (stevens-balcar formalism) (J=0),
            //gJ=0,general L and S moments given, use dipole approximation and separate formfactor for spin and orbital moment (J=-2)
            //intermediate coupling (gJ=0), go beyond dipole approximation (J=-3)
// (*jjjpars[1...n]).magFFj0(1..7)         formfactor j0 for atom 1...n <j0(kr)>-terms A,a,B,b,C,c,D
// (*jjjpars[1...n]).magFFj2(1..7)         formfactor j2 for atom 1...n <j2(kr)>-terms A,a,B,b,C,c,D
//     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
//     where U(r) is the Radial wave function for the unpaired electrons in the atom
// (*jjjpars[1...n]).magFFj4(1..7)         formfactor j4 for atom 1...n  (needed to go beyond dipole approx)
// (*jjjpars[1...n]).magFFj6(1..7)         formfactor j6 for atom 1...n  (needed to go beyond dipole approx)
// T				         temperature [K] (needed to go beyond dipole approx)
// (*jjjpars[1...n]).Zc		         Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
// (*jjjpars[1...n]).eigenstates(1..2J+1,1..2J+1)   CF+MF eigenstates (needed to go beyond dipole approx)

// output
// m                                  number of calculated reflections
// hkl[1...m](1..3)                   hkl values
// D[1...m]                           d spacing
// theta[1...m]                       scattering angle theta
// intmag[1...m]                        magnetic intensity
// intmagdip[1...m]                     magnetic intensity in dipole approx
// ikern[1...m]                       nuclear intensity
// sf[1...m]                          nuclear structurfactor |sf|
// lpg[1...m]                         lorentzfactor
// mx,my,mz,mxmy,mxmz,mymz,mx2my2mz2[].. fouriertransform of momentunitvectors (for mag xray scattering)

//****experimental parameters*************************************************
float scale,inuc;
float SF,Imag,Imagdip,Theta,lorentzf,d;

scale = 1 /(double)(n) /(double)(n); // scalingfactor of intensities
//***************************************************************************


D[0] = 10000;
int i;
//calculate reciprocal lattice vectors from r1,r2,r3
  Vector rez1(1,3),rez2(1,3),rez3(1,3);
  rezcalc(r1, r2, r3, rez1, rez2, rez3);

if(code==0){ m = 0;// reset m
 double qmax,rr;
 int hmax,kmax,lmax,ahi,aki,ali,sh,sk,sl,htrue,ktrue,ltrue,msort,hi,ki,li;
 qmax = 4.0 * PI * sin(thetamax / 180 * PI) / lambda;
 rr=r1*r1; hmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);
 rr=r2*r2; kmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);
 rr=r3*r3; lmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);
 for(ahi=0;ahi<=hmax;++ahi){
  for(aki=0;aki<=kmax;++aki){
   for(ali=0;ali<=lmax;++ali){
    for(sh=-1;sh<=1;sh+=2){
     for(sk=-1;sk<=1;sk+=2){
      for(sl=-1;sl<=1;sl+=2){
        if(ahi==0){sh=1;hi=0;}else{hi=sh*ahi;}
        if(aki==0){sk=1;ki=0;}else{ki=sk*aki;}
        if(ali==0){sl=1;li=0;}else{li=sl*ali;}

        if(hi==0&&li==0&&ki==0){htrue=1;ktrue=1;ltrue=1;} //goto 30
        else {  complex <double> mqx=0,mqx2=0,mqxy=0;
                complex <double> mqy=0,mqy2=0,mqxz=0;
                complex <double> mqz=0,mqz2=0,mqyz=0;  

          if(getint(jjjpars,hi,ki,li,thetamax,rez1,rez2,rez3,scale,T,lambda,ovalltemp,lorenz,n,J,d,Theta,Imag,Imagdip,inuc,SF,lorentzf,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2))
          {// reflection was found below thetamax....

            //sort according to descending d spacing
             if((Imag + inuc) > .0001||Imagdip > .0001||abs(mqx)*sqrt(scale)>0.0001
              ||abs(mqy)*sqrt(scale)>0.0001||abs(mqz)*sqrt(scale)>0.0001
              ||abs(mqx2)*sqrt(scale)>0.0001||abs(mqy2)*sqrt(scale)>0.0001||abs(mqz2)*sqrt(scale)>0.0001
              ||abs(mqxy)*sqrt(scale)>0.0001||abs(mqxz)*sqrt(scale)>0.0001||abs(mqyz)*sqrt(scale)>0.0001){
 
               ++m; if(m > MAXNOFREFLECTIONS){fprintf(stderr,"ERROR mcdiff: out of memory - too many reflections - chose smaller thetamax or recompile program with larger MAXNOFREFLECTIONS\n");exit(EXIT_FAILURE);}
               msort = m;
               while(D[msort-1]<=d){
                D[msort] = D[msort - 1];
                theta[msort]= theta[msort - 1];
                hkl[msort] = hkl[msort - 1];
                intmag[msort] = intmag[msort - 1];
                intmagdip[msort] = intmagdip[msort - 1];
                ikern[msort] = ikern[msort - 1];
                sf[msort] = sf[msort - 1];
                lpg[msort] = lpg[msort - 1];
 
                mx[msort] = mx[msort - 1];
                my[msort] = my[msort - 1];
                mz[msort] = mz[msort - 1];
                mxmy[msort] = mxmy[msort - 1];
                mxmz[msort] = mxmz[msort - 1];
                mymz[msort] = mymz[msort - 1];
                mx2[msort] = mx2[msort - 1];
                my2[msort] = my2[msort - 1];
                mz2[msort] = mz2[msort - 1];
                  --msort; 
               }
               hkl[msort](1) = hi;hkl[msort](2) = ki; hkl[msort](3) = li; 
               D[msort] = d; theta[msort] = Theta;
               intmag[msort] = Imag;intmagdip[msort] = Imagdip;  ikern[msort] = inuc;
               sf[msort] = SF; lpg[msort] = lorentzf;
               mx[msort]=(double)sqrt(scale)*mqx;
               my[msort]=(double)sqrt(scale)*mqy;
               mz[msort]=(double)sqrt(scale)*mqz;
               mxmy[msort]=(double)sqrt(scale)*mqxy;
               mxmz[msort]=(double)sqrt(scale)*mqxz;
               mymz[msort]=(double)sqrt(scale)*mqyz;
               mx2[msort]=(double)sqrt(scale)*mqx2;
               my2[msort]=(double)sqrt(scale)*mqy2;
               mz2[msort]=(double)sqrt(scale)*mqz2;
              }
            }
          }
     }}} //30 NEXT sl: NEXT sk: NEXT sh
   }}// NEXT ali NEXT aki
   printf("%i ", ahi);
  }
 printf("\n");
 }
else
 {for(i=1;i<=m;++i){int ii;
                complex <double> mqx=0,mqx2=0,mqxy=0;
                complex <double> mqy=0,mqy2=0,mqxz=0;
                complex <double> mqz=0,mqz2=0,mqyz=0;  
//               printf("%g %g %g\n",hkl[i](1),hkl[i](2),hkl[i](3));
               hkl[i](1)=rint(hkl[i](1));
               hkl[i](2)=rint(hkl[i](2));
               hkl[i](3)=rint(hkl[i](3));
          if(!getint(jjjpars,(int)hkl[i](1),(int)hkl[i](2),(int)hkl[i](3),thetamax,rez1,rez2,rez3,scale,T,lambda,ovalltemp,lorenz,n,J,d,Theta,Imag,Imagdip,inuc,SF,lorentzf,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2))
                {fprintf(stderr,"ERROR mcdiff: theta for reflection number %i above thetamax=%g\n",i,thetamax);exit(1);}
               D[i] = d; theta[i] = Theta;
               intmag[i] = Imag;intmagdip[i] = Imagdip;  ikern[i] = inuc;
               sf[i] = SF; lpg[i] = lorentzf;
               if(code==1){
               mx[i]=(double)sqrt(scale)*mqx;
               my[i]=(double)sqrt(scale)*mqy;
               mz[i]=(double)sqrt(scale)*mqz;
               mxmy[i]=(double)sqrt(scale)*mqxy;
               mxmz[i]=(double)sqrt(scale)*mqxz;
               mymz[i]=(double)sqrt(scale)*mqyz;
               mx2[i]=(double)sqrt(scale)*mqx2;
               my2[i]=(double)sqrt(scale)*mqy2;
               mz2[i]=(double)sqrt(scale)*mqz2;
                          }
                  }
 }

return;}


void printeln(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T, float Ha, float Hb, float Hc,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int * J,int m,Vector * hkl,float * ikern,float * intmag,float * intmagdip,float * D,float * theta,float * sf,float * lpg,complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,float a,float b,float c)
{// ausgabe auf file filename
 FILE * fout;char l[MAXNOFCHARINLINE];
 int i,j,chinr=0,ortho=1;
 double isave=0,isavedip=0,alpha,beta,gamma;
   extract(unitcell, " alpha", alpha); extract(unitcell, " beta", beta); extract(unitcell, " gamma", gamma); 
   if(alpha!=90||beta!=90||gamma!=90){ortho=0;}
 time_t curtime;
 struct tm * loctime;
  fout = fopen_errchk (filename, "w");
 fprintf(fout, "#{%s input file: %s ",filename,infile);
 curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
 fprintf(fout,"# unit cell:%s",unitcell);
 fprintf(fout,"#                   / %6.3f A \\     / %6.3f A \\     / %6.3f A \\ \n", r1(1), r2(1), r3(1));
 fprintf(fout,"#                r1=| %6.3f A |  r2=| %6.3f A |  r3=| %6.3f A |\n", r1(2), r2(2), r3(2));
 fprintf(fout,"#                   \\ %6.3f A /     \\ %6.3f A /     \\ %6.3f A /\n", r1(3), r2(3), r3(3));
 fprintf(fout, "# Wavelength=%g A   number of atoms: %i\n",lambda, n);
 fprintf(fout, "# T= %g K Ha= %g T Hb= %g T Hc= %g T\n",T,Ha,Hb,Hc);
 fprintf(fout, "# Overall temperature factor B=%g A^2: Intensity is proportional to exp(-2*B*(sin(theta)/lambda)^2)\n",ovalltemp);

 if(lorenz == 0){sprintf(l,"100 no lorentz factor calculated");}
 if(lorenz == 1){sprintf(l,"1 / sin^2(2theta)   neutron powder flat sample");}
 if(lorenz == 2){sprintf(l,"1 / sin(2theta) / sin(theta)    neutron powder cyl. sample");}
 if(lorenz == 3){sprintf(l,"1 / sin(2theta)     neutron single crystal");}
 if(lorenz == 4){sprintf(l,"d^3  neutron TOF powder cyl sample... log scaled d-pattern");}
 if(lorenz == 5){sprintf(l,"d^4  neutron TOF powder cyl sample... d-pattern");}
 fprintf(fout, "# Lorentz Factor: %s\n#\n",l);
 if(code<2)
 {fprintf(fout, "# Lorentz Factor not considered for resonant magnetic xray scattering - F1 and F2 transition intensities calculated\n");
  fprintf(fout, "# according to fRMXS as given in equation (2) of Longfield et al. PRB 66 054417 (2002) and maximized with respect to azimuth.\n#\n");
 }
 fprintf(fout, "# List of atomic positions dr1 dr2 dr3, moments m scattering lengths sl,\n# Debye Waller factor (sf ~ exp(-2 DWF sin^2(theta) / lambda^2)=EXP (-W),  (2*DWF=B=8 pi^2 <u^2>)\n");
 fprintf(fout, "#  and  Lande factors total angular momentum J (=0 if dipole approximation is used) <j0> and <j2> formfactor\n# coefficients\n");
 if (ortho==1)
 {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
 } else
 {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]my[MuB]mz[MuB]mx[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
  fprintf(fout, "#                         ...with x||(a x b), z||b and y normal to x and z\n");
 }
 for (i = 1;i<=n;++i)
 {if((double)(i)/50==(double)(i/50))
  {
   if (ortho==1)
   {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
   } else
   {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]my[MuB]mz[MuB]mx[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
    fprintf(fout, "#                         ...with x||a x b, z||b and y normal to x and z\n");
   }
  }
  fprintf(fout, "# %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f%+6.3fi %6.3f %6.3f ",(*jjjpars[i]).xyz(1),(*jjjpars[i]).xyz(2),(*jjjpars[i]).xyz(3),(*jjjpars[i]).mom(1),(*jjjpars[i]).mom(2),(*jjjpars[i]).mom(3),(*jjjpars[i]).SLR,(*jjjpars[i]).SLI,(*jjjpars[i]).DWF,(*jjjpars[i]).gJ);
  if(J[i]==0||J[i]==-3){fprintf(fout,"F(Q) beyond dip.approx.");}
  if(J[i]==-1){fprintf(fout,"F(Q)=j0-(1-2/gJ)j2 formfactor for rare earth/transition metals with gJ=2");}
  if(J[i]==-2){fprintf(fout,"FL(Q)=(j0+j2)/2 and FS(Q)=j0 formfactors separate for spin and orb. moments");}

  for (j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",(*jjjpars[i]).magFFj0(j));}
  for (j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",(*jjjpars[i]).magFFj2(j));}
  fprintf(fout, "\n");
 }
 fprintf(fout, "#}\n");
 hkl[0] = 0; D[0] = 100; theta[0] = 0; ikern[0] = 0; intmag[0] = 0;intmagdip[0] = 0; sf[0] = 0; lpg[0] = 0;mx[0]=0;my[0]=0;mz[0]=0;mx2[0]=0;my2[0]=0;mz2[0]=0;mxmy[0]=0;mxmz[0]=0;mymz[0]=0;
 double rpvalue=0,chisquared=0,total=0,rpvaluedip=0,chisquareddip=0;
 int imin=1;
 if(code==0)imin=0;
 for(i = imin;i<=m;++i)
 {if(code<2)
   {
    if((double)(i-imin)/50==(double)((i-imin)/50))
    {if (ortho==1)
     {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t) Itot(2t) |sf|     LF   Imag_dip(2t) F1:max-Isigpi azim Ipisig azim Ipipig azim F2:max-Isigpi azim Ipisig azim Ipipig azim  |^ma_q| |^mb_q| |^mc_q| |^ma^2_q||^mb^2_q||^mc^2_q||(^ma*^mb)_q||(^ma*^mc)_q||(^mb*^mc)_q|}\n");}
     else
     {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t) Itot(2t) |sf|     LF   Imag_dip(2t) \n");}
    }
    // calculate alpha_i delta_i for reflection hkl[i](1..3)  [currently ok only for ortholattices !!!]
    double alpha1,alpha2,alpha3,delta1,delta2,delta3,Q,sqr1,sqr2;
    alpha1=acos(-0.999999*hkl[i](1)*D[i]/a);   // the following lines should be extended to non ortho lattices !!!
    alpha2=acos(-0.999999*hkl[i](2)*D[i]/b);   // mind: in this section still the old convention is used: a||x,b||y,c||z ... this should be changed for nonortholattices
    alpha3=acos(-0.999999*hkl[i](3)*D[i]/c);
    delta1=acos(-1.0);
    sqr1=sqrt((1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)*(1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)+(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b)*(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b)+(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c)*(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c));
    sqr2=sqrt((1.0/b-hkl[i](2)*hkl[i](2)*D[i]*D[i]/b/b/b)*(1.0/b-hkl[i](2)*hkl[i](2)*D[i]*D[i]/b/b/b)+(hkl[i](2)*hkl[i](1)*D[i]*D[i]/b/b/a)*(hkl[i](2)*hkl[i](1)*D[i]*D[i]/b/b/a)+(hkl[i](2)*hkl[i](3)*D[i]*D[i]/b/b/c)*(hkl[i](2)*hkl[i](3)*D[i]*D[i]/b/b/c));
    delta2=acos(0.99999*hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b/b/sqr1/sqr2);
    // mind that delta2 is larger than pi if l is positive
    if(hkl[i](3)>0){delta2*=-1.0;}

    sqr1=sqrt((1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)*(1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)+(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c)*(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c)+(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b)*(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b));
    sqr2=sqrt((1.0/c-hkl[i](3)*hkl[i](3)*D[i]*D[i]/c/c/c)*(1.0/c-hkl[i](3)*hkl[i](3)*D[i]*D[i]/c/c/c)+(hkl[i](3)*hkl[i](1)*D[i]*D[i]/c/c/a)*(hkl[i](3)*hkl[i](1)*D[i]*D[i]/c/c/a)+(hkl[i](3)*hkl[i](2)*D[i]*D[i]/c/c/b)*(hkl[i](3)*hkl[i](2)*D[i]*D[i]/c/c/b));
    delta3=acos(0.99999*hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c/c/sqr1/sqr2);
    // mind that delta3 is larger than pi if k is negative
    if(hkl[i](2)<0){delta3*=-1.0;}
    //printf("%g %g %g\n",delta1,delta2,delta3);
    // maximize IspF1 IppF1 IpsF1  IspF2 IppF2 IpsF2  and remember corresponding azimuth 
    double IspF1=0,IppF1=0,IpsF1=0, IspF2=0, IppF2=0, IpsF2=0 ;
    double IspF1a=0,IppF1a=0,IpsF1a=0, IspF2a=0, IppF2a=0, IpsF2a=0;
    double azimuth;
    //printf("%g %g %g %g %g\n",hkl[i](1),hkl[i](2),hkl[i](3),delta3,hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c/c/sqr1/sqr2);
    
    for(azimuth=0.0;azimuth<=2*PI;azimuth+=PI/90)
     {complex <double> z1,z2,z3,z1z2,z2z3,z12,z32;
      double f,f1ps,f1sp,f1pp,f2ps,f2sp,f2pp;
      double st,ct,s2t;

      Matrix ang(1,3,1,3);
      ang(1,1)=sin(alpha1)*cos(azimuth+delta1);
      ang(1,2)=sin(alpha2)*cos(azimuth+delta2);
      ang(1,3)=sin(alpha3)*cos(azimuth+delta3);
      ang(2,1)=sin(alpha1)*sin(azimuth+delta1);
      ang(2,2)=sin(alpha2)*sin(azimuth+delta2);
      ang(2,3)=sin(alpha3)*sin(azimuth+delta3);
      ang(3,1)=cos(alpha1);
      ang(3,2)=cos(alpha2);
      ang(3,3)=cos(alpha3);
   
      z1=mx[i]*ang(1,1)+my[i]*ang(1,2)+mz[i]*ang(1,3);
      z2=mx[i]*ang(2,1)+my[i]*ang(2,2)+mz[i]*ang(2,3);
      z3=mx[i]*ang(3,1)+my[i]*ang(3,2)+mz[i]*ang(3,3);
   
      z1z2=mx2[i]*ang(1,1)*ang(2,1);
      z1z2+=my2[i]*ang(1,2)*ang(2,2);
      z1z2+=mz2[i]*ang(1,3)*ang(2,3);
      z1z2+=mxmy[i]*(ang(1,1)*ang(2,2)+ang(1,2)*ang(2,1));
      z1z2+=mxmz[i]*(ang(1,1)*ang(2,3)+ang(1,3)*ang(2,1));
      z1z2+=mymz[i]*(ang(1,2)*ang(2,3)+ang(1,3)*ang(2,2));
   
      z2z3=mx2[i]*ang(2,1)*ang(3,1);
      z2z3+=my2[i]*ang(2,2)*ang(3,2);
      z2z3+=mz2[i]*ang(2,3)*ang(3,3);
      z2z3+=mxmy[i]*(ang(2,1)*ang(3,2)+ang(2,2)*ang(3,1));
      z2z3+=mxmz[i]*(ang(2,1)*ang(3,3)+ang(2,3)*ang(3,1));
      z2z3+=mymz[i]*(ang(2,2)*ang(3,3)+ang(2,3)*ang(3,2));
   
      z12=mx2[i]*ang(1,1)*ang(1,1);
      z12+=my2[i]*ang(1,2)*ang(1,2);
      z12+=mz2[i]*ang(1,3)*ang(1,3);
      z12+=mxmy[i]*(ang(1,1)*ang(1,2)+ang(1,2)*ang(1,1));
      z12+=mxmz[i]*(ang(1,1)*ang(1,3)+ang(1,3)*ang(1,1));
      z12+=mymz[i]*(ang(1,2)*ang(1,3)+ang(1,3)*ang(1,2));
   
      z32=mx2[i]*ang(3,1)*ang(3,1);
      z32+=my2[i]*ang(3,2)*ang(3,2);
      z32+=mz2[i]*ang(3,3)*ang(3,3);
      z32+=mxmy[i]*(ang(3,1)*ang(3,2)+ang(3,2)*ang(3,1));
      z32+=mxmz[i]*(ang(3,1)*ang(3,3)+ang(3,3)*ang(3,1));
      z32+=mymz[i]*(ang(3,2)*ang(3,3)+ang(3,3)*ang(3,2));
      st=sin(theta[i]*PI/180);
      ct=cos(theta[i]*PI/180);
      s2t=sin(2.0*theta[i]*PI/180);
   
      f1ps=abs(z1*ct+z3*st); if(f1ps*f1ps>IpsF1){IpsF1=f1ps*f1ps;IpsF1a=azimuth*180/PI;}
      f1sp=abs(z3*st-z1*ct); if(f1sp*f1sp>IspF1){IspF1=f1sp*f1sp;IspF1a=azimuth*180/PI;}
      f1pp=-abs(z2*s2t); if(f1pp*f1pp>IppF1){IppF1=f1pp*f1pp;IppF1a=azimuth*180/PI;}
   
      f2ps=abs(-z1z2*st+z2z3*ct); if(f2ps*f2ps>IpsF2){IpsF2=f2ps*f2ps;IpsF2a=azimuth*180/PI;}
      f2sp=abs(z1z2*st+z2z3*ct); if(f2sp*f2sp>IspF2){IspF2=f2sp*f2sp;IspF2a=azimuth*180/PI;}
      f2pp=abs(-ct*ct*(z12*st*st/ct/ct+z32)); if(f2pp*f2pp>IppF2){IppF2=f2pp*f2pp;IppF2a=azimuth*180/PI;}
      if(code==1&&ortho==1)
       {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f        %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
       hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],
       f1ps*f1ps,azimuth*180/PI,
       f1sp*f1sp,azimuth*180/PI,
       f1pp*f1pp,azimuth*180/PI,
       f2ps*f2ps,azimuth*180/PI,
       f2sp*f2sp,azimuth*180/PI,
       f2pp*f2pp,azimuth*180/PI,
       abs(mx[i]),abs(my[i]),abs(mz[i]),abs(mx2[i]),abs(my2[i]),abs(mz2[i]),abs(mxmy[i]),abs(mxmz[i]),abs(mymz[i]));}
     }

    if(IspF1+IpsF1+IppF1+IspF2+IpsF2+IppF2+ikern[i]+intmag[i]>0.0001)
      {if(ortho==1)
       {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f        %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
        hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],
        IspF1,IspF1a,IpsF1,IpsF1a,IppF1,IppF1a,IspF2,IspF2a,IpsF2,IpsF2a,IppF2,IppF2a,
        abs(mx[i]),abs(my[i]),abs(mz[i]),abs(mx2[i]),abs(my2[i]),abs(mz2[i]),abs(mxmy[i]),abs(mxmz[i]),abs(mymz[i]));}
       else
       {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
        hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i]);}
      }
    if(code==1&&ortho==1){fprintf(fout,"#\n");}
   }
   if(code==2)//calculate rpvalue and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t) Itot(2t) |sf|     LF   Imag_dip(2t) Iobs\n");}
      fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
      hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],real(mx[i]));
     if(real(mx[i])>=0){
                      rpvalue+=abs(isave+ikern[i]+intmag[i]-abs(mx[i])); total+=abs(mx[i]);
                      rpvaluedip+=abs(isavedip+ikern[i]+intmagdip[i]-abs(mx[i]));
                      isave=0;isavedip=0;
                      }
     else {isave+=ikern[i]+intmag[i];isavedip+=ikern[i]+intmagdip[i];
          }
   }
   if(code==3)//calculate also rpvalue and chisquared and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t) Itot(2t) |sf|     LF   Imag_dip(2t) Iobs error\n");}
      fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
      hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],real(mx[i]),abs(my[i]));
     if(real(mx[i])>=0){
      rpvalue+=abs(isave+ikern[i]+intmag[i]-abs(mx[i])); total+=abs(mx[i]);
      rpvaluedip+=abs(isavedip+ikern[i]+intmagdip[i]-abs(mx[i]));
      chisquared+=(isave+ikern[i]+intmag[i]-abs(mx[i]))*(isave+ikern[i]+intmag[i]-abs(mx[i]))/abs(my[i])/abs(my[i]);
      chisquareddip+=(isavedip+ikern[i]+intmagdip[i]-abs(mx[i]))*(isavedip+ikern[i]+intmagdip[i]-abs(mx[i]))/abs(my[i])/abs(my[i]);
      isave=0;isavedip=0;++chinr;
                      }
     else {isave+=ikern[i]+intmag[i];isavedip+=ikern[i]+intmagdip[i];
          }
   } 

 }
if (code>=2){rpvalue*=100.0/total;fprintf(fout,"#rpvalue=%6.2f\n",rpvalue);
             rpvaluedip*=100.0/total;fprintf(fout,"#rpvaluedip=%6.2f\n",rpvaluedip);}
if (code==3){chisquared*=1.0/(double)chinr;fprintf(fout,"#chisquared=%6.4f\n",chisquared);
             chisquareddip*=1.0/(double)chinr;fprintf(fout,"#chisquareddip=%6.4f\n",chisquareddip);}

fclose(fout);
return;}


// hauptprogramm
int main (int argc, char **argv)
{ FILE * fin, * fin_coq, * fout;
  float ovalltemp,thetamax,lambda,a=0,b=0,c=0,Ha=0,Hb=0,Hc=0,alpha=0,beta=0,gamma=0;
  double T=0;
  int i,j,k,n,lorenz,nat, nofatoms,nr1=0,nr2=0,nr3=0,natmagnetic;
  long int pos=0;
  char instr[MAXNOFCHARINLINE+1];
  char cffilename[MAXNOFCHARINLINE+1];
  char unitcellstr[MAXNOFCHARINLINE+1];
  float numbers[60];numbers[0]=60;
  Vector r1(1,3),r2(1,3),r3(1,3);
  Vector rez1(1,3),rez2(1,3),rez3(1,3);

  // check command line
   if (argc > 1)
    { if (strcmp(argv[1],"-h")==0)
      {printf (" program mcdiff - calculate neutron diffraction pattern\n \
                use as: mcdiff [hkllistfilename]\n \
		- for format of input file mcdiff.in see mcphase manual\n \
                - optional an hkl list can be given in file hkllistfilename-in\n \
                  this case the program computes reflections in this list.\n \
                  if it is a 3column list,the azimuth dependence of\n \
                  magnetic scattering is calculated.If neutron intensities\n \
                  are listed in column 4 the program computes\n \
                  rpvalue=100 * sum_i |Icalc_i-Iexp_i|/sum_i |Iobs_i|\n \
                  and does not output the azimuth dependence \n \
                  if in addition experimental errors are given in column 5,\n \
                  chisquared=1/N *sum_i (Icalc_i - Iexp_i)^2/err_i^2 is\n \
                  calculated also.\n \
                - results are saved in results/mcdiff.out\n");
        exit (1);}
    }

// check if directory results exists and can be written to ...
   fout = fopen_errchk ("./results/mcdiff.out", "w");fclose(fout);

// test to test threej function
/* if (argc>6) {printf ("cint(%g)=%i\n", strtod(argv[6],NULL),cint(strtod(argv[6],NULL)));
		// test matpack routine  
                    int n,ndim=20; double thrcof[20];int errflag;
		    double min,max;
		     ThreeJSymbolM	(strtod(argv[1],NULL),
	             strtod(argv[2],NULL),
	             strtod(argv[3],NULL),
	             strtod(argv[4],NULL),
	             min,max, thrcof, ndim, 
			 errflag);
		     n=(int)(strtod(argv[5],NULL)-min);	 

              printf ("threej symbol=%g=%g\n",
              threej(strtod(argv[1],NULL),
	             strtod(argv[2],NULL),
	             strtod(argv[3],NULL),
	             strtod(argv[4],NULL),
	             strtod(argv[5],NULL),
	             strtod(argv[6],NULL)
	             ),thrcof[n]);
		     return 0;}
*/
// test spherical harmonics
/*if (argc>3) {
  int l,m; double theta,phi;
  l=(int)strtod(argv[1],NULL);
  m=(int)strtod(argv[2],NULL);
  theta=strtod(argv[3],NULL);
  phi=strtod(argv[4],NULL);
  
  printf("Y%i%i(%g,%g)=%g %+g i\n",l,m,theta,phi,real(SphericalHarmonicY (l,m,theta,phi)),imag(SphericalHarmonicY (l,m,theta,phi)));
 return 0;}
*/


fin_coq = fopen_errchk ("./mcdiff.in", "rb");
 fprintf(stdout,"\n reading file mcdiff.in\n\n");
 pos=ftell(fin_coq);
        char *token;
        fout = fopen_errchk ("./results/mcdiff.in", "w"); //copy input file to results
        while(feof(fin_coq)==false){if(fgets(instr, MAXNOFCHARINLINE, fin_coq))
                                     {// strip /r (dos line feed) from line if necessary
                                      while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
                                      fprintf(fout,"%s",instr);}
                                    }
        fclose(fout);
 fseek(fin_coq,pos,SEEK_SET); 
// input section 1 *******************************************************

  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&strstr (instr, "%SECTION 2%")==NULL) // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
    if (pos==-1) 
       {fprintf(stderr,"Error mcdiff: wrong sps file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   extract(instr,"nofatoms",nofatoms);  
   extract(instr,"lambda", lambda);
   extract(instr, "thetamax", thetamax);
   extract(instr, "nat", nat);
   extract(instr, "ovalltemp", ovalltemp);
   extract(instr, "lorentz", lorenz);
  }
  fseek(fin_coq,pos,SEEK_SET); 

if (lorenz == 0){fprintf(stderr,"Warning mcdiff: read lorentz=0, will calculate no Lorentzfactor.\n");}
if (lambda == 0){fprintf(stderr,"ERROR mcdiff: no wavelength lambda given or line does not start with # in section 1\n");exit(EXIT_FAILURE);}
if (thetamax == 0){fprintf(stderr,"ERROR mcdiff: no thetamax given or line does not start with # in section 1\n");exit(EXIT_FAILURE);}
printf("     section 1 - lambda=%g A thetamax= %g deg\n",lambda, thetamax);
printf("                 ovalltemp=%g A^2 lorentz-type=%i\n",ovalltemp,lorenz);


// input section 2 *********************************************************

  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&a==0&&b==0&&c==0) // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
    if (pos==-1) 
       {fprintf(stderr,"Error mcdiff: wrong sps file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   extract(instr, "nat", nat);
   extract(instr, " a", a);
   extract(instr, " b", b);
   extract(instr, " c", c);
   extract(instr, " alpha", alpha);
   extract(instr, " beta", beta);
   extract(instr, " gamma", gamma);
  }
  fseek(fin_coq,pos,SEEK_SET); 
  printf ("     section 2 - nat=%i\n",nat);
  float x1[nat+1],y1[nat+1],z1[nat+1];
  float sl1r[nat+1],sl1i[nat+1],dwf1[nat+1];

  if (nat!=0){ for(i=1;i<=nat;++i) { pos=ftell(fin_coq); 
                                     n=inputline(fin_coq,numbers);
                                     if (n==0) {if(feof(fin_coq)==0){fprintf(stderr,"Error mcdiff: end of input file in section 2\n");exit (EXIT_FAILURE);}
                                                fseek(fin_coq,pos,SEEK_SET); 
                                                fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                                                if(strstr (instr, "%SECTION 3%")!=NULL){fprintf (stderr,"ERROR mcdiff: Section 3 started before all nat=%i atoms of crystallographic unit cell were listed !\n",nat);exit (EXIT_FAILURE);}
                                               }
                                     else      {if (n<9) {fprintf (stderr,"ERROR mcdiff: Section 2 - Nonmagnetic Atoms: too few positional parameters for atom %i!\n",i);exit (EXIT_FAILURE);}
                                                sl1r[i]=numbers[1];sl1i[i]=numbers[2]; x1[i] = numbers[6]; y1[i] = numbers[7]; z1[i] = numbers[8];dwf1[i]=numbers[9];
                                                printf("                 sl=%g%+gi 10^-12cm at %g*r1%+g*r2%+g*r3 DWF=%g A^2\n",sl1r[i],sl1i[i],x1[i],y1[i],z1[i],dwf1[i]);
                                               }
                                    }
              }

// input section 3 *********************************************************
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&nr1*nr2*nr3==0) 
  { pos=ftell(fin_coq); 
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   if(a==0)extract(instr, " a", a);
   if(b==0)extract(instr, " b", b);
   if(c==0)extract(instr, " c", c);
   if(alpha==0)extract(instr, " alpha", alpha);
   if(beta==0)extract(instr, " beta", beta);
   if(gamma==0)extract(instr, " gamma", gamma);
    extract(instr, "r1x", r1(1));
    extract(instr, "r1y", r1(2));
    extract(instr, "r1z", r1(3));
    extract(instr, "r2x", r2(1));
    extract(instr, "r2y", r2(2));
    extract(instr, "r2z", r2(3));
    extract(instr, "r3x", r3(1));
    extract(instr, "r3y", r3(2));
    extract(instr, "r3z", r3(3));
    extract(instr, "r1a", r1(1));
    extract(instr, "r1b", r1(2));
    extract(instr, "r1c", r1(3));
    extract(instr, "r2a", r2(1));
    extract(instr, "r2b", r2(2));
    extract(instr, "r2c", r2(3));
    extract(instr, "r3a", r3(1));
    extract(instr, "r3b", r3(2));
    extract(instr, "r3c", r3(3));
    extract(instr, "nr1", nr1);
    extract(instr, "nr2", nr2);
    extract(instr, "nr3", nr3);
    extract(instr, "nat", natmagnetic);
    extract(instr, "T", T);
    extract(instr, "Ha", Ha);
    extract(instr, "Hb", Hb);
    extract(instr, "Hc", Hc);
  }
  fseek(fin_coq,pos,SEEK_SET); 
if (a == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant a given in section 3 or line does not start with # or nat too small: \n%s\n",instr);exit (EXIT_FAILURE);}
if (b == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant b given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
if (c == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant c given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
printf("     section 3 - a=%g A  b=%g A c=%g A alpha=%g  beta=%g gamma=%g\n",a,b,c,alpha,beta,gamma);
sprintf(unitcellstr," a= %g A  b= %g A c= %g A  alpha=%g  beta=%g gamma=%g\n",a,b,c,alpha,beta,gamma);
printf("                 r1= %5.3ga + %5.3gb + %5.3gc\n", r1(1), r1(2), r1(3));
printf("                 r2= %5.3ga + %5.3gb + %5.3gc\n", r2(1), r2(2), r2(3));
printf("                 r3= %5.3ga + %5.3gb + %5.3gc\n", r3(1), r3(2), r3(3));

//printf("                    / %5.3ga \\     / %5.3ga \\     / %5.3ga \\    x||c \n", r1(1), r2(1), r3(1));
//printf("                 r1=| %5.3gb |  r2=| %5.3gb |  r3=| %5.3gb |    y||a\n", r1(2), r2(2), r3(2));
//printf("                    \\ %5.3gc /     \\ %5.3gc /     \\ %5.3gc /    z||b\n", r1(3), r2(3), r3(3));
Matrix rtoxyz(1,3,1,3); // define transformation matrix to calculate components of
                           // r1 r2 and r3 with respect to the xyz coordinate system
rtoxyz(1,1)=0;
rtoxyz(2,1)=a*sin(gamma*PI/180);
rtoxyz(3,1)=a*cos(gamma*PI/180);

rtoxyz(1,2)=0;
rtoxyz(2,2)=0;
rtoxyz(3,2)=b;

rtoxyz(3,3)=c*cos(alpha*PI/180);
rtoxyz(2,3)=(a*c*cos(beta*PI/180)-rtoxyz(3,3)*rtoxyz(3,1))/rtoxyz(2,1);
if (fabs(rtoxyz(2,3))>c){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoxyz(1,3)=c*c-rtoxyz(2,3)*rtoxyz(2,3)-rtoxyz(3,3)*rtoxyz(3,3);
if (rtoxyz(1,3)<=0){fprintf(stderr,"ERROR mcdiff: alpha beta and gamma geometrically inconsistent\n");exit(EXIT_FAILURE);}
rtoxyz(1,3)=sqrt(rtoxyz(1,3));

r1=(rtoxyz*r1)*(double)nr1;
r2=(rtoxyz*r2)*(double)nr2;
r3=(rtoxyz*r3)*(double)nr3;

//r1(1) = a * r1(1) * nr1;
//r2(1) = a * r2(1) * nr2;
//r3(1) = a * r3(1) * nr3; 
//r1(2) = b * r1(2) * nr1;
//r2(2) = b * r2(2) * nr2;
//r3(2) = b * r3(2) * nr3;
//r1(3) = c * r1(3) * nr1;
//r2(3) = c * r2(3) * nr2;
//r3(3) = c * r3(3) * nr3;

// input section 4 *********************************************************

if (nr1 == 0){fprintf(stderr,"ERROR mcdiff: nr1 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
if (nr2 == 0){fprintf(stderr,"ERROR mcdiff: nr2 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
if (nr3 == 0){fprintf(stderr,"ERROR mcdiff: nr3 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
printf ("     section 4 - nr1=%i nr2=%i nr3=%i\n",nr1,nr2,nr3);
printf ("                 nat=%i magnetic atoms\n",natmagnetic);

n = nr1 * nr2 * nr3 * nat + natmagnetic; //atoms in der magnetic unit cell


int J[n+1]; // code for indicating if ion is nonmagnetic (J=1), 
            // go beyond dipole approx for rare earth (J=0)
            // use magnetic moment with dipole approx (J=-1)
            // use L and S values with dipole approx (J=-2) 
            // go beyond dipole approx for gJ=0 (L and S separately)
jjjpar ** jjjpars = new jjjpar * [n+1];
double lnZ,U;
printf("                 reading magnetic atoms and moments ...\n");

for(i=1;i<=natmagnetic;++i){ 
                            instr[0]='#';J[i]=-1;
                            while(instr[strspn(instr," \t")]=='#'){pos=ftell(fin_coq);
                                                                   if(feof(fin_coq)==1){fprintf(stderr,"mcdiff Error: end of file before all magnetic atoms could be read\n");exit(EXIT_FAILURE);}
                                                                  fgets(instr,MAXNOFCHARINLINE,fin_coq);
                                                                  }
			     // get cffilename out of "{filename}   ..."

                            if(instr[strspn(instr," \t")]!='{'){fprintf(stderr,"ERROR mcdiff: magnetic atom line has to start with '{'\n");exit (EXIT_FAILURE);}
			    if (strchr(instr,'}')==NULL){fprintf(stderr,"ERROR mcdiff: no '}' found after filename for magnetic atom %s\n",instr);exit (EXIT_FAILURE);}

                            instr[strspn(instr," \t")]='=';
			    extract(instr,"",cffilename,(size_t)MAXNOFCHARINLINE);
                            if(strchr(cffilename,'}')!=NULL){*strchr(cffilename,'}')='\0';}
                            if(strchr(cffilename,' ')!=NULL){*strchr(cffilename,' ')='\0';}
                            if(strchr(cffilename,'\t')!=NULL){*strchr(cffilename,'\t')='\0';}
                            //printf("%s\n",cffilename);

                             // read the rest of the line and split into numbers
                            fseek(fin_coq,pos+strchr(instr,'}')-instr+1,SEEK_SET); 
                            j=inputline(fin_coq,numbers);
                            if (j<9) {fprintf(stderr,"ERROR mcdiff: too few parameters for magnetic atom %i: %s\n",i,instr);exit(EXIT_FAILURE);}
                             jjjpars[i]=new jjjpar(numbers[4] / nr1,numbers[5] / nr2,numbers[6] / nr3, cffilename);
                              // store moment and components of S and L (if given)
                              for(k=7;k<=j&&k<=15;++k){(*jjjpars[i]).mom(k-6) = numbers[k];}
                              if((*jjjpars[i]).gJ==0){if(j>=15){J[i]=-2; // do not use input moment but spin and angular momentum for calculation
                                                                // do some consistency checks
                                                                if (fabs((*jjjpars[i]).mom(1)-2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5))>0.001){fprintf(stderr,"Warning mcdiff: a-component magnetic moment and <La>+2<Sa> not consistent for atom %i - setting moment=<L>+2<S> \n",i);}
                                                                if (fabs((*jjjpars[i]).mom(2)-2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7))>0.001){fprintf(stderr,"Warning mcdiff: b-component magnetic moment and <Lb>+2<Sb> not consistent for atom %i - setting moment=<L>+2<S>\n",i);}
                                                                if (fabs((*jjjpars[i]).mom(3)-2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9))>0.001){fprintf(stderr,"Warning mcdiff: c-component magnetic moment and <Lc>+2<Sc> not consistent for atom %i - setting moment=<L>+2<S>\n",i);}
                                                                (*jjjpars[i]).mom(1)=2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5);
                                                                (*jjjpars[i]).mom(2)=2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7);
                                                                (*jjjpars[i]).mom(3)=2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9);
                                                                }
                                                           else {J[i]=-1;(*jjjpars[i]).gJ=2;} // just use spin formfactor
                                                      }
                            instr[0]='#';
                            while(instr[strspn(instr," \t")]=='#'&&feof(fin_coq)==0){pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}

                            if (strchr(instr,'>')==NULL)
                             {fseek(fin_coq,pos,SEEK_SET);} // no ">" found --> do dipole approx
                             else          
                             {J[i]=0; // J=0 tells that full calculation should be done for this ion
                              fseek(fin_coq,pos+strchr(instr,'>')-instr+1,SEEK_SET); 
                              j=inputline(fin_coq,numbers);printf("dimension of mf = %i\n",j);
                              Vector heff(1,j);for(k=1;k<=j;++k){heff(k)=numbers[k];}
                              if ((*jjjpars[i]).gJ==0)
 			      {J[i]=-3;fprintf(stderr,"mcdiff: gJ=0 - going beyond dipolar approximation for intermediate coupling");
   			             (*jjjpars[i]).eigenstates(heff,T); // calculate eigenstates
                               // do some consistency checks
                               ComplexMatrix est((*jjjpars[i]).est.Rlo(),(*jjjpars[i]).est.Rhi(),(*jjjpars[i]).est.Clo(),(*jjjpars[i]).est.Chi());
                                             est=(*jjjpars[i]).est;
                               Vector moment(1,j);moment=(*jjjpars[i]).mcalc(T,heff,lnZ,U,est);
                               for(k=1;k<=j;++k){if (fabs((*jjjpars[i]).mom(k+3)-moment(k))>0.001){fprintf(stderr,"Warning mcdiff: meanfields and <J> read from input file not consistent for atom %i - using values calculated from meanfield\n",i);}
                                                  (*jjjpars[i]).mom(3+k)=moment(k);
                                                 }
                               (*jjjpars[i]).mom(1)=2*(*jjjpars[i]).mom(4)+(*jjjpars[i]).mom(5);
                               (*jjjpars[i]).mom(2)=2*(*jjjpars[i]).mom(6)+(*jjjpars[i]).mom(7);
                               (*jjjpars[i]).mom(3)=2*(*jjjpars[i]).mom(8)+(*jjjpars[i]).mom(9);
			      }
			      else
			      {// beyond formalism for rare earth		    
                               (*jjjpars[i]).eigenstates(heff,T); //calculate some eigenstates
                               // do some consistency checks
                               ComplexMatrix est((*jjjpars[i]).est.Rlo(),(*jjjpars[i]).est.Rhi(),(*jjjpars[i]).est.Clo(),(*jjjpars[i]).est.Chi());
                                             est=(*jjjpars[i]).est;
                               Vector moment(1,j);moment=(*jjjpars[i]).mcalc(T,heff,lnZ,U,est);
                               if (fabs((*jjjpars[i]).mom(1)-(*jjjpars[i]).gJ*moment(1))>0.001){fprintf(stderr,"Warning mcdiff: a-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                               if (fabs((*jjjpars[i]).mom(2)-(*jjjpars[i]).gJ*moment(2))>0.001){fprintf(stderr,"Warning mcdiff: b-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                               if (fabs((*jjjpars[i]).mom(3)-(*jjjpars[i]).gJ*moment(3))>0.001){fprintf(stderr,"Warning mcdiff: c-component meanfields and moments not consistent for atom %i - using values calculated from meanfield\n",i);}
                               (*jjjpars[i]).mom(1)=(*jjjpars[i]).gJ*moment(1);
                               (*jjjpars[i]).mom(2)=(*jjjpars[i]).gJ*moment(2);
                               (*jjjpars[i]).mom(3)=(*jjjpars[i]).gJ*moment(3);
                            
                     
                               if(Norm((*jjjpars[i]).Zc)==0){fprintf(stderr,"WARNING mcdiff: Z(K) coefficients not found or zero in file %s\n",cffilename);}
                               }
                               if(Norm((*jjjpars[i]).magFFj4)==0){fprintf(stderr,"WARNING mcdiff: <j4(Q)> coefficients not found or zero in file %s\n",cffilename);}
                               if(Norm((*jjjpars[i]).magFFj6)==0){fprintf(stderr,"WARNING mcdiff: <j6(Q)> coefficients not found or zero in file %s\n",cffilename);}
 			      }

                             if((*jjjpars[i]).SLR==0){fprintf(stderr,"WARNING mcdiff: SCATTERINGLENGTHREAL not found or zero in file %s\n",cffilename);}
//                             if((*jjjpars[i]).gJ==0){fprintf(stderr,"WARNING mcdiff: GJ not found or zero in file %s - gJ=0 means Ja=Sa Jb=La Jc=Sb Jd=Lb Je=Sc Jf=Lc !\n",cffilename);}
                             if(Norm((*jjjpars[i]).magFFj0)==0){fprintf(stderr,"WARNING mcdiff: <j0(Q)> coefficients not found or zero in file %s\n",cffilename);}
                             if(Norm((*jjjpars[i]).magFFj2)==0){fprintf(stderr,"WARNING mcdiff: <j2(Q)> coefficients not found or zero in file %s\n",cffilename);}
                           }
  fclose(fin_coq);
printf ("calculating ...\n");  

//now insert also nonmagnetic elements into the unit cell
int ncryst,na,nb,nc;
ncryst = natmagnetic;
for(na = 1;na<=nr1;++na){
 for(nb = 1;nb<=nr2;++nb){
  for(nc = 1;nc<=nr3;++nc){
   if(nat!=0){
    for(i=1;i<=nat;++i){
      ++ncryst;
      J[ncryst]=1;
      jjjpars[ncryst]=new jjjpar((na + x1[i] - 1) / nr1,(nb + y1[i] - 1) / nr2,(nc + z1[i] - 1) / nr3,sl1r[i],sl1i[i],dwf1[i]);
      (*jjjpars[ncryst]).mom=0;
      (*jjjpars[ncryst]).gJ=0;
      }
    }
}}}

int m=0;
Vector * hkl = new Vector[MAXNOFREFLECTIONS+1];for(i=0;i<=MAXNOFREFLECTIONS;++i){hkl[i]=Vector(1,3);}
float D[MAXNOFREFLECTIONS+1];
float theta[MAXNOFREFLECTIONS+1];
float intmag[MAXNOFREFLECTIONS+1];
float intmagdip[MAXNOFREFLECTIONS+1];
float ikern[MAXNOFREFLECTIONS+1];
float sf[MAXNOFREFLECTIONS+1];
float lpg[MAXNOFREFLECTIONS+1];
complex <double> mx[MAXNOFREFLECTIONS+1];
complex <double> my[MAXNOFREFLECTIONS+1];
complex <double> mz[MAXNOFREFLECTIONS+1];
complex <double> mxmy[MAXNOFREFLECTIONS+1];
complex <double> mxmz[MAXNOFREFLECTIONS+1];
complex <double> mymz[MAXNOFREFLECTIONS+1];
complex <double> mx2[MAXNOFREFLECTIONS+1];
complex <double> my2[MAXNOFREFLECTIONS+1];
complex <double> mz2[MAXNOFREFLECTIONS+1];

rezcalc (r1, r2, r3, rez1, rez2, rez3);
Vector hhkkll(1,3);
int code=0;
// if hkllist is given, read the file and put hkls to hkl[i], m is number of reflections to be considered
if (argc>1){int nr;
      float nn[20];nn[0]=19;
     // open hkllist file
     fprintf(stdout,"reading hkl list from file %s\n",argv[1]);
       fin = fopen_errchk (argv[1], "rb");
       while(feof(fin)==false){nr=inputline(fin,nn);
                               if(nr>2)
                               {hhkkll(1)=nn[1];hhkkll(2)=nn[2];hhkkll(3)=nn[3];++m;
//printf("%g %g %g\n", hhkkll(1),hhkkll(2),hhkkll(3));
                                code=1;                                
                               // transformieren der millerindizes auf magnetische einheitszelle
                                  hkl[m](1)=hhkkll*(rtoxyz.Inverse()*r1);
                                  hkl[m](2)=hhkkll*(rtoxyz.Inverse()*r2);
                                  hkl[m](3)=hhkkll*(rtoxyz.Inverse()*r3);
//printf("%g %g %g\n\n", hkl[m](1), hkl[m](2), hkl[m](3));                                
                                if(nr>3){mx[m]=complex <double> (nn[4],0);code=2;}// intensities given                                
                                if(nr>4){my[m]=complex <double> (nn[5],0);code=3;}// errors given                                
                              }}
       fclose(fin);      
           }

// transformieren der millerindizes auf kristallographische einheitszelle


neutint(jjjpars,code,T,lambda, thetamax, ovalltemp, lorenz, r1, r2, r3, n,  J, m, hkl, D, theta, intmag,intmagdip, ikern, sf, lpg,mx,my,mz,mxmy,mxmz,mymz,mx2,my2,mz2);



// transformieren der millerindizes auf kristallographische einheitszelle

for(i=1;i<=m;++i){hhkkll=hkl[i];
                  hkl[i]=hhkkll(1)*rez1+hhkkll(2)*rez2+hhkkll(3)*rez3;
                  hkl[i]/=2.0*PI;
                  hhkkll=hkl[i];
                  hkl[i](1)=hhkkll*rtoxyz.Column(1);
                  hkl[i](2)=hhkkll*rtoxyz.Column(2);
                  hkl[i](3)=hhkkll*rtoxyz.Column(3);
                 }


printeln(jjjpars,code,"./results/mcdiff.out","mcdiff.in", unitcellstr,T,Ha,Hb,Hc, lambda, ovalltemp, lorenz, r1, r2, r3, n,  J, m, hkl, ikern, intmag,intmagdip, D, theta, sf, lpg,mx,my,mz,mxmy,mxmz,mymz,mx2,my2,mz2,a,b,c);

fprintf (stderr,"...results written to ./results/mcdiff.out\n");

//  for (i=1;i<=n;++i){delete jjjpars[i];}
//  delete []jjjpars;
  delete []hkl;
 return 0;
}


