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

double Z(int K, float J0, float J2, float J4, float J6, Vector Zc)
{// calculate Z(K)
 if (K==1) return Zc(1)*J0+Zc(2)*J2;
 if (K==3) return Zc(3)*J2+Zc(4)*J4;
 if (K==5) return Zc(5)*J4+Zc(6)*J6;
 if (K==7) return Zc(7)*J6;
 
 return 0;
}

int getint(int hi,int ki,int li,float thetamax,Vector rez1,Vector rez2, Vector rez3,float scale,float T,float lambda,float ovalltemp,int lorenz,int & n,Vector * xyz,float * slr,float * sli,float * dwf,Vector * mom,float * gj,float * J,float & d,Vector * j0,Vector * j2,Vector * j4,Vector * j6,Vector * Zc,ComplexMatrix ** eigenstates,float & Theta,float & Imag,float & Imagdip,float & inuc,float & SF,float & lorentzf,complex <double> & mqx,complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2)
{
            double s,Q,sintheta,qr,FQ,sin2theta,ovallt,mux,muy,muz;
            int i;
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
                                 complex <double> scl(slr[i],sli[i]); 
                                 qr=hi*xyz[i](1)+ki*xyz[i](2)+li*xyz[i](3);

                                 //nuclear structure factor nsfr,nsfc
                                 nsf+=scl*exp(-2*PI*qr*im)*exp(-2 * dwf[i] *s*s);;

                                //magnetic structure factors
                                if(gj[i]!=0){
                                             float J0,J2,J4,J6;
                                             J0=j0[i](1)*exp(-j0[i](2)*s*s)+j0[i](3)*exp(-j0[i](4)*s*s)+j0[i](5)*exp(-j0[i](6)*s*s)+j0[i](7);
                                             J2=j2[i](1)*exp(-j2[i](2)*s*s)+j2[i](3)*exp(-j2[i](4)*s*s)+j2[i](5)*exp(-j2[i](6)*s*s)+j2[i](7);
					     J2*=s*s;
                                             FQ = (J0 + J2 * (2 / gj[i] - 1)); // formfactor F(Q)
                                             if(J[i]>=0){ // go beyond dipole approximation
                                                         J4=j4[i](1)*exp(-j4[i](2)*s*s)+j4[i](3)*exp(-j4[i](4)*s*s)+j4[i](5)*exp(-j4[i](6)*s*s)+j4[i](7);
	                				     J4*=s*s;
                                                         J6=j6[i](1)*exp(-j6[i](2)*s*s)+j6[i](3)*exp(-j6[i](4)*s*s)+j6[i](5)*exp(-j6[i](6)*s*s)+j6[i](7);
				                	     J6*=s*s;
//printf("sintheta/lambda=%g J0=%g J2=%g J4=%g J6=%g\n",s,J0,J2,J4,J6);
							     int dj=(int)(2*J[i]+1);
							 // calculate th and ph (polar angles of Q with respect to xyz of CEF)
                                                         double th,ph,Qx,Qy,Qz;							 
							 Qx=Qvec(3);Qy=Qvec(1);Qz=Qvec(2);
							 th=acos(Qz/Q);
							 if(sin(th)>=SMALL){
							                    if(Qx>0){ph=acos(Qx/(Q*sin(th))-SMALL);}
									    else    {ph=acos(Qx/(Q*sin(th))+SMALL);}
									   }
							 else{ph=0;}
							 if (Qy<0){ph=2*PI-ph;} 

							 ComplexMatrix MQXM(1,dj,1,dj),MQYM(1,dj,1,dj),MQZM(1,dj,1,dj);
                                                             // .... calculate these matrices ...(formula 11.141-143 in lovesey)
							    int K,Qd,M,Md;double MJ,MJd,PKdQd,thj,factor;
							    complex <double>bracketx,brackety,bracketz;
							    MQXM=0;MQYM=0;MQZM=0;
							    for(K=1;K<=7;K+=2){
							     factor=sqrt(4.0*PI)*Z(K,J0,J2,J4,J6,Zc[i])/K;
                               if (factor!=0){
							     thj=threej((float)K,J[i],J[i],0,J[i],-J[i]);
							     for(Qd=-K;Qd<=K;Qd+=1){
                                                 bracketx=0;brackety=0;
							       if(K-1>=Qd+1&&K-1>=-Qd-1)
							       {bracketx+=SphericalHarmonicY (K-1,Qd+1,th,ph)*sqrt((double)(K-Qd)*(K-Qd-1));
							        brackety+=SphericalHarmonicY (K-1,Qd+1,th,ph)*sqrt((double)(K-Qd)*(K-Qd-1));
							       }
							       if(K-1>=Qd-1&&K-1>=-Qd+1)
                                                 {bracketx-=SphericalHarmonicY (K-1,Qd-1,th,ph)*sqrt((double)(K+Qd)*(K+Qd-1));
							        brackety+=SphericalHarmonicY (K-1,Qd-1,th,ph)*sqrt((double)(K+Qd)*(K+Qd-1));
							       }
							       if(K-1>=Qd&&K-1>=-Qd)
                                                 {bracketz=SphericalHarmonicY (K-1,Qd,th,ph)*sqrt((double)(K-Qd)*(K+Qd));
                                                 }else {bracketz=0;}

//.(1)..USE !		     ThreeJSymbolM	(J1,J2,J3,M1,&M2min,&M2max,*thrcof,ndim,errflag);
                                                 double thrj[30];int ndim=30; double MJdmin,MJdmax; int errflag;
                                                                      
                                                 ThreeJSymbolM ((float)K,J[i],J[i],-(float)Qd,MJdmin,MJdmax,thrj,ndim,errflag);
                                                 if (errflag!=0){fprintf(stderr,"ERROR mcdiff: threejsymbol error %i\n",errflag);exit(EXIT_FAILURE);}           
                                                 for (Md=int(MJdmin+1+J[i]);Md<=int(MJdmax+1+J[i]);++Md){
						                 MJd=(float)Md-1-J[i];
								 MJ=-Qd+MJd;M=int(MJ+1+J[i]);
							         PKdQd=thrj[Md-int(MJdmin+1+J[i])]/thj; 
							         PKdQd*=odd(int(J[i]-MJd)) ? -1 : 1;
							         MQXM(M,Md)+=0.5*factor*PKdQd*bracketx;
							         MQYM(M,Md)+=-im*brackety*0.5*factor*PKdQd;
							         MQZM(M,Md)+=factor*PKdQd*bracketz;
                                                                 }

/*
//.(2)..                     3jsymb=threej(J1,J2,J3,M1,M2,M3) 
//							       for(M=1;M<=dj;++M){MJ=(float)M-1-J[i];  
//							        for(Md=1;Md<=dj;++Md){MJd=(float)Md-1-J[i]; 
//							         // according to 11.140 lovesey book        
//							         PKdQd=threej((float)K,J[i],J[i],-(float)Qd,MJd,-MJ)/thj; 
//							         PKdQd*=odd(int(J[i]-MJd)) ? -1 : 1;
//							         MQXM(M,Md)+=im*0.5*factor*PKdQd*bracketx;
//							         MQYM(M,Md)+=-0.5*factor*PKdQd*brackety;
//							         MQZM(M,Md)+=factor*PKdQd*bracketz;
//							        }
//							       }
*/
                                             }
							      }
                                                             }
							     // ... calculate thermal expectation values
							     // using the eigenstates and T 
							     // mom[i](1) = ....
          						     //
							     ComplexVector mm(1,3); 
                                                             mm=0;
	                                                     for(K=1;K<=dj;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){
                                                                mm(3)+=(*eigenstates[i])(0,K)*conj((*eigenstates[i])(M,K))*MQXM(M,Md)*(*eigenstates[i])(Md,K); 
                                                                mm(1)+=(*eigenstates[i])(0,K)*conj((*eigenstates[i])(M,K))*MQYM(M,Md)*(*eigenstates[i])(Md,K); 
                                                                mm(2)+=(*eigenstates[i])(0,K)*conj((*eigenstates[i])(M,K))*MQZM(M,Md)*(*eigenstates[i])(Md,K); 
                                                                }}}
// myPrintComplexMatrix(stdout,MQXM);
// myPrintComplexMatrix(stdout,MQYM);
// myPrintComplexMatrix(stdout,MQZM);
// 					       myPrintComplexMatrix(stdout,(*eigenstates[i]));}
					       
                                                         mm*=2;
// myPrintComplexVector(stdout,mm);//equivalent to moment ...
							 msfx+=0.5*mm(1)*exp(-2*PI*qr*im);
							 msfy+=0.5*mm(2)*exp(-2*PI*qr*im);
							 msfz+=0.5*mm(3)*exp(-2*PI*qr*im);
					                }else
							{msfx+=mom[i](1)*FQ/2*exp(-2*PI*qr*im);
							 msfy+=mom[i](2)*FQ/2*exp(-2*PI*qr*im);
							 msfz+=mom[i](3)*FQ/2*exp(-2*PI*qr*im);
							}
// myPrintVector(stdout,mom[i]);//equivalent to moment ...
                                                         msfdipx+=mom[i](1)*FQ/2*exp(-2*PI*qr*im);
							 msfdipy+=mom[i](2)*FQ/2*exp(-2*PI*qr*im);
							 msfdipz+=mom[i](3)*FQ/2*exp(-2*PI*qr*im);
                                                         mux=mom[i](1);
                                                         muy=mom[i](2);
                                                         muz=mom[i](3);
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

void neutint(int code,float T,float lambda, float thetamax, float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int & n,Vector * xyz,Vector * mom,float * dwf,float* J,float * gj,float * slr, float * sli,Vector * j0,Vector *  j2,Vector *  j4,Vector *  j6, Vector * Zc, ComplexMatrix ** eigenstates,int & m,Vector *  hkl,float * D,float * theta,float * intmag,float * intmagdip,float * ikern,float * sf,float * lpg,complex <double>*ma,complex <double>*mb,complex <double>*mc,complex <double>*mamb,complex <double>*mamc,complex <double>*mbmc,complex <double>*ma2,complex <double>*mb2,complex <double>*mc2)
{//****************************************************************************
// this routine calculates the intensity of elastic neutrons
// for a given magnetic unit cell (crystal axis orthogonal)
// the magnetic scattering is treated in the dipole approximation
// input:
// code                               governs if a list of hkl given in hkl[] or all hkls should be generated
// lambda                             wavelength[A]
// ovalltemp                          overall temperature factor [A]
// r1(1..3),r2(),r3()                 vectors of primitive unit cell[A]
// n                                  number of atoms per unit cell
// xyz[1...n](1..3)                   atomic positional parameters dr1 dr2 dr3
//                                        '(with respect to primitive lattice)
// dwf[1..n]			      debye waller factors [A]
// slr(1...n),sli(1...n)              nuclear scattering length[10^-12cm]
// mom[1...n](1..3)                   atomic magnetic moment [mb]
//                                       ' (with respect to ortholatice abc)
// gj[1...n]			      Lande factor
// J[1..n]			      total angular momentum (if J=0 dipole approx is used)
// j0[1..n](1..7)                     formfactor j0 for atom 1...n <j0(kr)>-terms A,a,B,b,C,c,D
// j2[1..n](1..7)                     formfactor j2 for atom 1...n <j2(kr)>-terms A,a,B,b,C,c,D
//     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
//     where U(r) is the Radial wave function for the unpaired electrons in the atom
// j4[1..n](1..7)                        formfactor j4 for atom 1...n  (needed to go beyond dipole approx)
// j6[1..n](1..7)                        formfactor j6 for atom 1...n  (needed to go beyond dipole approx)
// T				         temperature [K] (needed to go beyond dipole approx)
// Zc[1...n]			         Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
// eigenstates [1..n](1..2J+1,1..2J+1)   CF+MF eigenstates (needed to go beyond dipole approx)

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
// ma,mb,mc,mamb,mamc,mbmc,ma2mb2mc2[].. fouriertransform of momentunitvectors (for mag xray scattering)

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

          if(getint(hi,ki,li,thetamax,rez1,rez2,rez3,scale,T,lambda,ovalltemp,lorenz,n,xyz,slr,sli,dwf,mom,gj,J,d,j0,j2,j4,j6,Zc,eigenstates,Theta,Imag,Imagdip,inuc,SF,lorentzf,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2))
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
 
                ma[msort] = ma[msort - 1];
                mb[msort] = mb[msort - 1];
                mc[msort] = mc[msort - 1];
                mamb[msort] = mamb[msort - 1];
                mamc[msort] = mamc[msort - 1];
                mbmc[msort] = mbmc[msort - 1];
                ma2[msort] = ma2[msort - 1];
                mb2[msort] = mb2[msort - 1];
                mc2[msort] = mc2[msort - 1];
                  --msort; 
               }
               hkl[msort](1) = hi;hkl[msort](2) = ki; hkl[msort](3) = li; 
               D[msort] = d; theta[msort] = Theta;
               intmag[msort] = Imag;intmagdip[msort] = Imagdip;  ikern[msort] = inuc;
               sf[msort] = SF; lpg[msort] = lorentzf;
               ma[msort]=(double)sqrt(scale)*mqx;
               mb[msort]=(double)sqrt(scale)*mqy;
               mc[msort]=(double)sqrt(scale)*mqz;
               mamb[msort]=(double)sqrt(scale)*mqxy;
               mamc[msort]=(double)sqrt(scale)*mqxz;
               mbmc[msort]=(double)sqrt(scale)*mqyz;
               ma2[msort]=(double)sqrt(scale)*mqx2;
               mb2[msort]=(double)sqrt(scale)*mqy2;
               mc2[msort]=(double)sqrt(scale)*mqz2;
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
          if(!getint((int)hkl[i](1),(int)hkl[i](2),(int)hkl[i](3),thetamax,rez1,rez2,rez3,scale,T,lambda,ovalltemp,lorenz,n,xyz,slr,sli,dwf,mom,gj,J,d,j0,j2,j4,j6,Zc,eigenstates,Theta,Imag,Imagdip,inuc,SF,lorentzf,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2))
                {fprintf(stderr,"ERROR mcdiff: theta for reflection number %i above thetamax=%g\n",i,thetamax);exit(1);}
               D[i] = d; theta[i] = Theta;
               intmag[i] = Imag;intmagdip[i] = Imagdip;  ikern[i] = inuc;
               sf[i] = SF; lpg[i] = lorentzf;
               if(code==1){
               ma[i]=(double)sqrt(scale)*mqx;
               mb[i]=(double)sqrt(scale)*mqy;
               mc[i]=(double)sqrt(scale)*mqz;
               mamb[i]=(double)sqrt(scale)*mqxy;
               mamc[i]=(double)sqrt(scale)*mqxz;
               mbmc[i]=(double)sqrt(scale)*mqyz;
               ma2[i]=(double)sqrt(scale)*mqx2;
               mb2[i]=(double)sqrt(scale)*mqy2;
               mc2[i]=(double)sqrt(scale)*mqz2;
                          }
                  }
 }

return;}


void printeln(int code,const char * filename,const char* infile,char * unitcell,float T, float Ha, float Hb, float Hc,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,Vector * xyz,Vector * mom,float * dwf,float* J,float * gj,float * slr,float * sli,Vector *j0,Vector *j2,Vector *j4,Vector *j6,Vector * Zc, ComplexMatrix ** eigenstates,int m,Vector * hkl,float * ikern,float * intmag,float * intmagdip,float * D,float * theta,float * sf,float * lpg,complex <double>*ma,complex <double>*mb,complex <double>*mc,complex <double>*mamb,complex <double>*mamc,complex <double>*mbmc,complex <double>*ma2,complex <double>*mb2,complex <double>*mc2,float a,float b,float c)
{// ausgabe auf file filename
 FILE * fout;char l[MAXNOFCHARINLINE];
 int i,j,chinr=0;
 double isave=0,isavedip=0;
 time_t curtime;
 struct tm * loctime;
  fout = fopen_errchk (filename, "w");
 fprintf(fout, "#{%s input file: %s ",filename,infile);
 curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
 fprintf(fout,"# unit cell:%s",unitcell);
 fprintf(fout,"#                   / %6.3g A \\     / %6.3g A \\     / %6.3g A \\ \n", r1(1), r2(1), r3(1));
 fprintf(fout,"#                r1=| %6.3g A |  r2=| %6.3g A |  r3=| %6.3g A |\n", r1(2), r2(2), r3(2));
 fprintf(fout,"#                   \\ %6.3g A /     \\ %6.3g A /     \\ %6.3g A /\n", r1(3), r2(3), r3(3));
 fprintf(fout, "# Wavelength=%g A   number of atoms: %i\n",lambda, n);
 fprintf(fout, "# T= %g K Ha= %g T Hb= %g T Hc= %g T\n",T,Ha,Hb,Hc);
 fprintf(fout, "# Overall temperature factor: exp(-2*%g A^2*(sin(theta)/lambda)^2)\n",ovalltemp);

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
 fprintf(fout, "# List of atomic positions dr1 dr2 dr3, moments m scattering lengths sl, Debye Waller (defined as ovalltemp)\n");
 fprintf(fout, "#  and  Lande factors total angular momentum J (=0 if dipole approximation is used) <j0> and <j2> formfactor\n# coefficients\n");
 fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A] gJ     J      <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
 for (i = 1;i<=n;++i)
 {if((double)(i)/50==(double)(i/50))
  {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A] gJ     J      <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");}
  fprintf(fout, "# %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f%+6.3fi %6.3f %6.3f %6.3f ",xyz[i](1),xyz[i](2),xyz[i](3),mom[i](1),mom[i](2),mom[i](3),slr[i],sli[i],dwf[i],gj[i],J[i]);
  for (j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",j0[i](j));}
  for (j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",j2[i](j));}
  fprintf(fout, "\n");
 }
 fprintf(fout, "#}\n");
 hkl[0] = 0; D[0] = 100; theta[0] = 0; ikern[0] = 0; intmag[0] = 0;intmagdip[0] = 0; sf[0] = 0; lpg[0] = 0;ma[0]=0;mb[0]=0;mc[0]=0;ma2[0]=0;mb2[0]=0;mc2[0]=0;mamb[0]=0;mamc[0]=0;mbmc[0]=0;
 double rpvalue=0,chisquared=0,total=0,rpvaluedip=0,chisquareddip=0;
 int imin=1;
 if(code==0)imin=0;
 for(i = imin;i<=m;++i)
 {if(code<2)
   {
    if((double)(i-imin)/50==(double)((i-imin)/50))
    {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t) Itot(2t) |sf|     LF   Imag_dip(2t) F1:max-Isigpi azim Ipisig azim Ipipig azim F2:max-Isigpi azim Ipisig azim Ipipig azim  |^ma_q| |^mb_q| |^mc_q| |^ma^2_q||^mb^2_q||^mc^2_q||(^ma*^mb)_q||(^ma*^mc)_q||(^mb*^mc)_q|}\n");}
    // calculate alpha_i delta_i for reflextion hkl[i](1..3)
    double alpha1,alpha2,alpha3,delta1,delta2,delta3,Q,sqr1,sqr2;
    alpha1=acos(-0.999999*hkl[i](1)*D[i]/a);
    alpha2=acos(-0.999999*hkl[i](2)*D[i]/b);
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
   
      z1=ma[i]*ang(1,1)+mb[i]*ang(1,2)+mc[i]*ang(1,3);
      z2=ma[i]*ang(2,1)+mb[i]*ang(2,2)+mc[i]*ang(2,3);
      z3=ma[i]*ang(3,1)+mb[i]*ang(3,2)+mc[i]*ang(3,3);
   
      z1z2=ma2[i]*ang(1,1)*ang(2,1);
      z1z2+=mb2[i]*ang(1,2)*ang(2,2);
      z1z2+=mc2[i]*ang(1,3)*ang(2,3);
      z1z2+=mamb[i]*(ang(1,1)*ang(2,2)+ang(1,2)*ang(2,1));
      z1z2+=mamc[i]*(ang(1,1)*ang(2,3)+ang(1,3)*ang(2,1));
      z1z2+=mbmc[i]*(ang(1,2)*ang(2,3)+ang(1,3)*ang(2,2));
   
      z2z3=ma2[i]*ang(2,1)*ang(3,1);
      z2z3+=mb2[i]*ang(2,2)*ang(3,2);
      z2z3+=mc2[i]*ang(2,3)*ang(3,3);
      z2z3+=mamb[i]*(ang(2,1)*ang(3,2)+ang(2,2)*ang(3,1));
      z2z3+=mamc[i]*(ang(2,1)*ang(3,3)+ang(2,3)*ang(3,1));
      z2z3+=mbmc[i]*(ang(2,2)*ang(3,3)+ang(2,3)*ang(3,2));
   
      z12=ma2[i]*ang(1,1)*ang(1,1);
      z12+=mb2[i]*ang(1,2)*ang(1,2);
      z12+=mc2[i]*ang(1,3)*ang(1,3);
      z12+=mamb[i]*(ang(1,1)*ang(1,2)+ang(1,2)*ang(1,1));
      z12+=mamc[i]*(ang(1,1)*ang(1,3)+ang(1,3)*ang(1,1));
      z12+=mbmc[i]*(ang(1,2)*ang(1,3)+ang(1,3)*ang(1,2));
   
      z32=ma2[i]*ang(3,1)*ang(3,1);
      z32+=mb2[i]*ang(3,2)*ang(3,2);
      z32+=mc2[i]*ang(3,3)*ang(3,3);
      z32+=mamb[i]*(ang(3,1)*ang(3,2)+ang(3,2)*ang(3,1));
      z32+=mamc[i]*(ang(3,1)*ang(3,3)+ang(3,3)*ang(3,1));
      z32+=mbmc[i]*(ang(3,2)*ang(3,3)+ang(3,3)*ang(3,2));
      st=sin(theta[i]*PI/180);
      ct=cos(theta[i]*PI/180);
      s2t=sin(2.0*theta[i]*PI/180);
   
      f1ps=abs(z1*ct+z3*st); if(f1ps*f1ps>IpsF1){IpsF1=f1ps*f1ps;IpsF1a=azimuth*180/PI;}
      f1sp=abs(z3*st-z1*ct); if(f1sp*f1sp>IspF1){IspF1=f1sp*f1sp;IspF1a=azimuth*180/PI;}
      f1pp=-abs(z2*s2t); if(f1pp*f1pp>IppF1){IppF1=f1pp*f1pp;IppF1a=azimuth*180/PI;}
   
      f2ps=abs(-z1z2*st+z2z3*ct); if(f2ps*f2ps>IpsF2){IpsF2=f2ps*f2ps;IpsF2a=azimuth*180/PI;}
      f2sp=abs(z1z2*st+z2z3*ct); if(f2sp*f2sp>IspF2){IspF2=f2sp*f2sp;IspF2a=azimuth*180/PI;}
      f2pp=abs(-ct*ct*(z12*st*st/ct/ct+z32)); if(f2pp*f2pp>IppF2){IppF2=f2pp*f2pp;IppF2a=azimuth*180/PI;}
      if(code==1)
       {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f        %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
       hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],
       f1ps*f1ps,azimuth*180/PI,
       f1sp*f1sp,azimuth*180/PI,
       f1pp*f1pp,azimuth*180/PI,
       f2ps*f2ps,azimuth*180/PI,
       f2sp*f2sp,azimuth*180/PI,
       f2pp*f2pp,azimuth*180/PI,
       abs(ma[i]),abs(mb[i]),abs(mc[i]),abs(ma2[i]),abs(mb2[i]),abs(mc2[i]),abs(mamb[i]),abs(mamc[i]),abs(mbmc[i]));}
     }

    if(IspF1+IpsF1+IppF1+IspF2+IpsF2+IppF2+ikern[i]+intmag[i]>0.0001)
      {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f        %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
      hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],
      IspF1,IspF1a,IpsF1,IpsF1a,IppF1,IppF1a,IspF2,IspF2a,IpsF2,IpsF2a,IppF2,IppF2a,
       abs(ma[i]),abs(mb[i]),abs(mc[i]),abs(ma2[i]),abs(mb2[i]),abs(mc2[i]),abs(mamb[i]),abs(mamc[i]),abs(mbmc[i]));}
    if(code==1){fprintf(fout,"#\n");}
   }
   if(code==2)//calculate rpvalue and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t) Itot(2t) |sf|     LF   Imag_dip(2t) Iobs\n");}
      fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
      hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],real(ma[i]));
     if(real(ma[i])>=0){
                      rpvalue+=abs(isave+ikern[i]+intmag[i]-abs(ma[i])); total+=abs(ma[i]);
                      rpvaluedip+=abs(isavedip+ikern[i]+intmagdip[i]-abs(ma[i]));
                      isave=0;isavedip=0;
                      }
     else {isave+=ikern[i]+intmag[i];isavedip+=ikern[i]+intmagdip[i];
          }
   }
   if(code==3)//calculate alsorpvalue and chisquared and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t) Itot(2t) |sf|     LF   Imag_dip(2t) Iobs error\n");}
      fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
      hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],ikern[i], intmag[i], ikern[i]+intmag[i],sf[i],lpg[i],intmagdip[i],real(ma[i]),abs(mb[i]));
     if(real(ma[i])>=0){
      rpvalue+=abs(isave+ikern[i]+intmag[i]-abs(ma[i])); total+=abs(ma[i]);
      rpvaluedip+=abs(isavedip+ikern[i]+intmagdip[i]-abs(ma[i]));
      chisquared+=(isave+ikern[i]+intmag[i]-abs(ma[i]))*(isave+ikern[i]+intmag[i]-abs(ma[i]))/abs(mb[i])/abs(mb[i]);
      chisquareddip+=(isavedip+ikern[i]+intmagdip[i]-abs(ma[i]))*(isavedip+ikern[i]+intmagdip[i]-abs(ma[i]))/abs(mb[i])/abs(mb[i]);
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
  float ovalltemp,thetamax,lambda,a=0,b=0,c=0,T=0,Ha=0,Hb=0,Hc=0;
  int i,j,n,lorenz,nat, nofatoms,nna=0,nnb=0,nnc=0,natmagnetic;
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
        while(feof(fin_coq)==false){fgets(instr, MAXNOFCHARINLINE, fin_coq);
                                    // strip /r (dos line feed) from line if necessary
                                    while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}
                                    fprintf(fout,"%s",instr);
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
                                                printf("                 sl=%g%+gi 10^-12cm at %g*r1%+g*r2%+g*r3 DWF=%g\n",sl1r[i],sl1i[i],x1[i],y1[i],z1[i],dwf1[i]);
                                               }
                                    }
              }

// input section 3 *********************************************************
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#'&&nna*nnb*nnc==0) 
  { pos=ftell(fin_coq); 
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   if(a==0)extract(instr, " a", a);
   if(b==0)extract(instr, " b", b);
   if(c==0)extract(instr, " c", c);
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
    extract(instr, "nr1", nna);
    extract(instr, "nr2", nnb);
    extract(instr, "nr3", nnc);
    extract(instr, "nat", natmagnetic);
    extract(instr, "T", T);
    extract(instr, "Ha", Ha);
    extract(instr, "Hb", Hb);
    extract(instr, "Hc", Hc);
  }
  fseek(fin_coq,pos,SEEK_SET); 
if (a == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant a given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
if (b == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant b given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
if (c == 0) {fprintf(stderr,"ERROR mcdiff: no lattice constant c given in section 3 or line does not start with #: \n%s\n",instr);exit (EXIT_FAILURE);}
printf("     section 3 - a=%g A  b=%g A c=%g A\n",a,b,c);
sprintf(unitcellstr," a= %g A  b= %g A c= %g A  alpha=90  beta=90 gamma=90\n",a,b,c);
printf("                    / %5.3ga \\     / %5.3ga \\     / %5.3ga \\ \n", r1(1), r2(1), r3(1));
printf("                 r1=| %5.3gb |  r2=| %5.3gb |  r3=| %5.3gb |\n", r1(2), r2(2), r3(2));
printf("                    \\ %5.3gc /     \\ %5.3gc /     \\ %5.3gc /\n", r1(3), r2(3), r3(3));
r1(1) = a * r1(1) * nna;
r2(1) = a * r2(1) * nnb;
r3(1) = a * r3(1) * nnc;
r1(2) = b * r1(2) * nna;
r2(2) = b * r2(2) * nnb;
r3(2) = b * r3(2) * nnc;
r1(3) = c * r1(3) * nna;
r2(3) = c * r2(3) * nnb;
r3(3) = c * r3(3) * nnc;

// input section 4 *********************************************************

if (nna == 0){fprintf(stderr,"ERROR ELN: nr1 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
if (nnb == 0){fprintf(stderr,"ERROR ELN: nr2 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
if (nnc == 0){fprintf(stderr,"ERROR ELN: nr3 not given or line does not start with # in section 4\n");exit(EXIT_FAILURE);}
printf ("     section 4 - nr1=%i nr2=%i nr3=%i\n",nna,nnb,nnc);
printf ("                 nat=%i magnetic atoms\n",natmagnetic);

n = nna * nnb * nnc * nat + natmagnetic; //atoms in der magnetic unit cell


float slr[n+1], sli[n+1],dwf[n+1], gj[n+1],J[n+1];
Vector * xyz = new Vector [n+1]; for(i=0;i<=n;++i){xyz[i]=Vector (1,3);}
Vector * mom = new Vector[n+1]; for(i=0;i<=n;++i){mom[i]=Vector(1,3);}
Vector * j0 = new Vector[n+1]; for(i=0;i<=n;++i){j0[i]=Vector(1,7);}
Vector * j2 = new Vector[n+1]; for(i=0;i<=n;++i){j2[i]=Vector(1,7);}
Vector * j4 = new Vector[n+1]; for(i=0;i<=n;++i){j4[i]=Vector(1,7);}
Vector * j6 = new Vector[n+1]; for(i=0;i<=n;++i){j6[i]=Vector(1,7);}
Vector * Zc = new Vector[n+1]; for(i=0;i<=n;++i){Zc[i]=Vector(1,7);}
ComplexMatrix ** eigenstates= new ComplexMatrix * [n+1];

printf("                 reading magnetic atoms and moments ...\n");

for(i=1;i<=natmagnetic;++i){ 
                            instr[0]='#';J[i]=-1;
                            while(instr[strspn(instr," \t")]=='#'){pos=ftell(fin_coq);fgets(instr,MAXNOFCHARINLINE,fin_coq);}
			     // get cffilename out of "{filename}   ..."
			    extract(instr,"J",J[i]);
			    if (J[i]>=0) // if J= has been found and J read, overtype it with spaces until a "{" is found
			    {char * token;for(token=instr;*token!='{'&&token<instr+MAXNOFCHARINLINE;++token){*token=' ';};}
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
                             xyz[i](1) = numbers[4] / nna;
			     xyz[i](2) = numbers[5] / nnb;
			     xyz[i](3) = numbers[6] / nnc;
			     dwf[i]=0;slr[i]=0;sli[i]=0;gj[i]=0;

                            fin = fopen_errchk (cffilename, "rb");
                            while(feof(fin)==0)
                                 {instr[0]='#';
                                  while (instr[strspn(instr," \t")]=='#'&&feof(fin)==0) 
                                   {fgets(instr,MAXNOFCHARINLINE,fin);}
                                  extract(instr,"SCATTERINGLENGTHREAL",slr[i]);  
                                  extract(instr,"SCATTERINGLENGTHIMAG",sli[i]);  
                                  extract(instr,"DWF",dwf[i]);  
                                  extract(instr,"GJ",gj[i]);  
                                  extract(instr,"FFj0A",j0[i](1));  
                                  extract(instr,"FFj0a",j0[i](2));  
                                  extract(instr,"FFj0B",j0[i](3));  
                                  extract(instr,"FFj0b",j0[i](4));  
                                  extract(instr,"FFj0C",j0[i](5));  
                                  extract(instr,"FFj0c",j0[i](6));  
                                  extract(instr,"FFj0D",j0[i](7));  

                                  extract(instr,"FFj2A",j2[i](1));  
                                  extract(instr,"FFj2a",j2[i](2));  
                                  extract(instr,"FFj2B",j2[i](3));  
                                  extract(instr,"FFj2b",j2[i](4));  
                                  extract(instr,"FFj2C",j2[i](5));  
                                  extract(instr,"FFj2c",j2[i](6));  
                                  extract(instr,"FFj2D",j2[i](7));  

                                  extract(instr,"FFj4A",j4[i](1));  
                                  extract(instr,"FFj4a",j4[i](2));  
                                  extract(instr,"FFj4B",j4[i](3));  
                                  extract(instr,"FFj4b",j4[i](4));  
                                  extract(instr,"FFj4C",j4[i](5));  
                                  extract(instr,"FFj4c",j4[i](6));  
                                  extract(instr,"FFj4D",j4[i](7));  

                                  extract(instr,"FFj6A",j6[i](1));  
                                  extract(instr,"FFj6a",j6[i](2));  
                                  extract(instr,"FFj6B",j6[i](3));  
                                  extract(instr,"FFj6b",j6[i](4));  
                                  extract(instr,"FFj6C",j6[i](5));  
                                  extract(instr,"FFj6c",j6[i](6));  
                                  extract(instr,"FFj6D",j6[i](7));  
                                   // coefficients of Z(K') according to Lovesey chapter 11.6.1 page 233
                                  extract(instr,"Z1c0",Zc[i](1));  
                                  extract(instr,"Z1c2",Zc[i](2));  
                                  extract(instr,"Z3c2",Zc[i](3));  
                                  extract(instr,"Z3c4",Zc[i](4));  
                                  extract(instr,"Z5c4",Zc[i](5));  
                                  extract(instr,"Z5c6",Zc[i](6));  
                                  extract(instr,"Z7c6",Zc[i](7));  
                                  
                                  }
                             fclose(fin); 
                            if (J[i]>=0) {int dj=(int)(2*J[i]+1);
			                  eigenstates [i]= new ComplexMatrix (0,dj,1,dj);
					  if (myReadComplexMatrix (fin_coq, (*eigenstates[i]))==false){fprintf (stderr, "ERROR mcdiff reading eigenstates for atom %i\n",i);exit(EXIT_FAILURE);}
                                          //myPrintComplexMatrix (stdout,(*eigenstates[i])); 
                                                //calculate partition sum
                                               double z=0;double KBT=T*KB,E0;
    		                               E0=real((*eigenstates[i])(0,1));
                                               for(j=1;j<=dj;++j)
                                               {z+=exp(-((real((*eigenstates[i])(0,j))-E0)/KBT));}
                                               // put boltzmann population into row 0 of eigenstates...
                                               for(j=1;j<=dj;++j)
                                               {(*eigenstates[i])(0,j)=complex<double>(exp(-(real((*eigenstates[i])(0,j))-E0)/KBT)/z,0);}
                             if(Norm(Zc[i])==0){fprintf(stderr,"WARNING mcdiff: Z(K) coefficients not found or zero in file %s\n",cffilename);}
                             if(Norm(j4[i])==0){fprintf(stderr,"WARNING mcdiff: <j4(Q)> coefficients not found or zero in file %s\n",cffilename);}
                             if(Norm(j6[i])==0){fprintf(stderr,"WARNING mcdiff: <j6(Q)> coefficients not found or zero in file %s\n",cffilename);}
 			                 }
                             if(slr[i]==0){fprintf(stderr,"WARNING mcdiff: SCATTERINGLENGTHREAL not found or zero in file %s\n",cffilename);}
                             if(gj[i]==0){fprintf(stderr,"WARNING mcdiff: GJ not found or zero in file %s\n",cffilename);}
                             if(Norm(j0[i])==0){fprintf(stderr,"WARNING mcdiff: <j0(Q)> coefficients not found or zero in file %s\n",cffilename);}
                             if(Norm(j2[i])==0){fprintf(stderr,"WARNING mcdiff: <j2(Q)> coefficients not found or zero in file %s\n",cffilename);}
                              mom[i](1) = numbers[7] * gj[i];
                              mom[i](2) = numbers[8] * gj[i];
                              mom[i](3) = numbers[9] * gj[i];
                           }
  fclose(fin_coq);
printf ("calculating ...\n");  

//now insert also nonmagnetic elements into the unit cell
int ncryst,na,nb,nc;
ncryst = natmagnetic;
for(na = 1;na<=nna;++na){
 for(nb = 1;nb<=nnb;++nb){
  for(nc = 1;nc<=nnc;++nc){
   if(nat!=0){
    for(i=1;i<=nat;++i){
      ++ncryst;
      slr[ncryst]=sl1r[i];
      sli[ncryst]=sl1i[i];
      xyz[ncryst](1) = (na + x1[i] - 1) / nna;
      xyz[ncryst](2) = (nb + y1[i] - 1) / nnb;
      xyz[ncryst](3) = (nc + z1[i] - 1) / nnc;
      mom[ncryst]=0;
      j0[ncryst]=0;
      j2[ncryst]=0;
      j4[ncryst]=0;
      j6[ncryst]=0;
      Zc[ncryst]=0;
      eigenstates[ncryst]=0;
      gj[ncryst]=0;
      J[ncryst]=-1;
      dwf[ncryst]=dwf1[i];
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
complex <double> ma[MAXNOFREFLECTIONS+1];
complex <double> mb[MAXNOFREFLECTIONS+1];
complex <double> mc[MAXNOFREFLECTIONS+1];
complex <double> mamb[MAXNOFREFLECTIONS+1];
complex <double> mamc[MAXNOFREFLECTIONS+1];
complex <double> mbmc[MAXNOFREFLECTIONS+1];
complex <double> ma2[MAXNOFREFLECTIONS+1];
complex <double> mb2[MAXNOFREFLECTIONS+1];
complex <double> mc2[MAXNOFREFLECTIONS+1];

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
                               {hhkkll(1)=nn[1]/a;hhkkll(2)=nn[2]/b;hhkkll(3)=nn[3]/c;++m;
                                code=1;                                
                               // transformieren der millerindizes auf magnetische einheitszelle
                                  hkl[m](1)=hhkkll*r1;
                                  hkl[m](2)=hhkkll*r2;
                                  hkl[m](3)=hhkkll*r3;
                                if(nr>3){ma[m]=complex <double> (nn[4],0);code=2;}// intensities given                                
                                if(nr>4){mb[m]=complex <double> (nn[5],0);code=3;}// errors given                                
                              }}
       fclose(fin);      
           }
neutint(code,T,lambda, thetamax, ovalltemp, lorenz, r1, r2, r3, n, xyz, mom,dwf,J, gj, slr, sli, j0, j2, j4, j6, Zc,eigenstates,m, hkl, D, theta, intmag,intmagdip, ikern, sf, lpg,ma,mb,mc,mamb,mamc,mbmc,ma2,mb2,mc2);



// transformieren der millerindizes auf kristallographische einheitszelle

for(i=1;i<=m;++i){hhkkll=hkl[i];
                  hkl[i]=hhkkll(1)*rez1+hhkkll(2)*rez2+hhkkll(3)*rez3;
                  hkl[i]/=2.0*PI;
                  hkl[i](1)*=a;hkl[i](2)*=b;hkl[i](3)*=c;
                 }


printeln(code,"./results/mcdiff.out","mcdiff.in", unitcellstr,T,Ha,Hb,Hc, lambda, ovalltemp, lorenz, r1, r2, r3, n, xyz, mom,dwf,J, gj , slr,sli,j0,j2,j4,j6,Zc,eigenstates, m, hkl, ikern, intmag,intmagdip, D, theta, sf, lpg,ma,mb,mc,mamb,mamc,mbmc,ma2,mb2,mc2,a,b,c);

fprintf (stderr,"...results written to ./results/mcdiff.out\n");

  delete []j0;
  delete []j2;
  delete []j4;
  delete []j6;
  delete []Zc;
  for (i=1;i<=natmagnetic;++i){if (J[i]>0)delete eigenstates [i];}
  delete []eigenstates;
  delete []mom;
  delete []xyz;
  delete []hkl;
 return 0;
}


