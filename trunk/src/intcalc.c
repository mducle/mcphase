//calculate intensities for given energy
#define PI 3.1415926535
#define KB 0.0862     // Boltzmanns constant in mev/K

double intcalc_approx(ComplexMatrix Tau, int level,double en,inimcdis & ini,par & inputpars,jq & J,Vector & q,Vector & hkl,mdcf & md,int do_verbose)
{//calculates approximate intensity for energylevel i - according to chapter 8.2 mcphas manual

 int i,j,i1,j1,k1,l1,i2,j2,k2,l2,s,ss,stau,sstau,ia,ja;
 double intensity=1.2;
 double QQ,ki,kf;
 complex <double> chileft;

 // determine chi
   ComplexMatrix chi(1,md.nofcomponents*ini.mf.n()*md.nofatoms,1,md.nofcomponents*ini.mf.n()*md.nofatoms);
   
   
    for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
     stau=(ini.mf.nb()*ini.mf.nc()*(i1-1)+ini.mf.nc()*(j1-1)+k1-1)*md.nofatoms;
     s=stau*md.nofcomponents;
     
 for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
     sstau=(ini.mf.nb()*ini.mf.nc()*(i2-1)+ini.mf.nc()*(j2-1)+k2-1)*md.nofatoms;
     ss=sstau*md.nofcomponents;
    for(ia=1;ia<=md.nofatoms;++ia){for(ja=1;ja<=md.nofatoms;++ja){
    for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
  chileft=PI*md.lambda(i1,j1,k1)(md.nofcomponents*ia,md.nofcomponents*ia)*md.U(i1,j1,k1)((ia-1)*md.nofcomponents+i,(ia-1)*md.nofcomponents+md.nofcomponents)*Tau(stau+ia,level);
  chi(s+(ia-1)*md.nofcomponents+i,ss+(ja-1)*md.nofcomponents+j)=
     chileft*conj(Tau(sstau+ja,level))*conj(md.U(i2,j2,k2)((ja-1)*md.nofcomponents+j,(ja-1)*md.nofcomponents+md.nofcomponents))*md.lambda(i2,j2,k2)(md.nofcomponents*ja,md.nofcomponents*ja);
//     chileft*Tau.Conjugate().Transpose()(level,sstau+ja)*md.U(i2,j2,k2).Conjugate().Transpose()((ja-1)*md.nofcomponents+md.nofcomponents,(ja-1)*md.nofcomponents+j)*md.lambda(i2,j2,k2)(md.nofcomponents*ja,md.nofcomponents*ja);
    }}
   }}
  }}}
 }}}


 //  chi'' to  S (bose factor) ... fluctuaion dissipation theorem

   complex<double> im(0,1.0);
   ComplexMatrix S(1,md.nofcomponents*ini.mf.n()*md.nofatoms,1,md.nofcomponents*ini.mf.n()*md.nofatoms);
   double bose;
   if (fabs(en)>SMALL)
   {bose=1.0/(1.0-exp(-en*(1.0/KB/ini.T)));
   }else{//quasielastic needs special treatment
   bose=ini.T*KB/SMALL;
   }
   S=bose*2*chi;
   
 // polarization factor
// neutrons only sense first 3x3 part of S !! - this is taken into account by setting 0 all
// higher components in the polarization factor !!!
 Matrix pol(1,md.nofcomponents,1,md.nofcomponents);
 Vector qxyz(1,3);
 pol=0;
    qxyz(1)=hkl(1)/inputpars.a; // only correct for ortholattices !!!!
    qxyz(2)=hkl(2)/inputpars.b;
    qxyz(3)=hkl(3)/inputpars.c;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qxyz(i)*qxyz(j)/(qxyz*qxyz);
    }}
    QQ=Norm(qxyz)*2*PI;

 //multiply polarization factor, formfactor and debeywallerfactor
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 for(l1=1;l1<=md.nofatoms;++l1){
   s=((((i1-1)*ini.mf.nb()+(j1-1))*ini.mf.nc()+(k1-1))*md.nofatoms+(l1-1))*md.nofcomponents;
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  for(l2=1;l2<=md.nofatoms;++l2){
   ss=((((i2-1)*ini.mf.nb()+(j2-1))*ini.mf.nc()+(k2-1))*md.nofatoms+(l2-1))*md.nofcomponents;
    for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
      S(s+i,ss+j)*=pol(i,j);
      S(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debeywallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
      S(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debeywallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
    }}   
  }
  }}}
 }
 }}}

 // determine dsigma
 //divide by number of crytallographic unit cells  (ini.mf.n()) in magnetic unit cell
intensity=abs(Sum(S))/ini.mf.n()/PI/2.0*3.65/4.0/PI; 

// here should be entered factor  k/k' + absolute scale factor
if (ini.ki==0)
{ki=sqrt(ini.kf*ini.kf+0.4811*en);
intensity*=ini.kf/ki;
}
else
{if (ini.ki*ini.ki-0.4811*en<0)
   {fprintf(stderr,"ERROR mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with ki=const=%g/A at (%g,%g,%g)",en,ini.ki,hkl(1),hkl(2),hkl(3));
                            exit(EXIT_FAILURE);}
kf=sqrt(ini.ki*ini.ki-0.4811*en);
intensity*=kf/ini.ki;
}


return intensity;	
}




double intcalc(double en,inimcdis & ini,par & inputpars,jq & J,Vector & q,Vector & hkl,mdcf & md,int do_verbose,double epsilon)
{int i,j,i1,j1,k1,l1,i2,j2,k2,l2,s,ss;
 double intensity=1.2;
 double QQ,ki,kf;

 complex<double> z(en,epsilon);
 complex<double> eps(epsilon,0);
 // determine chi
   ComplexMatrix chi0c(1,md.nofcomponents*md.nofatoms,1,md.nofcomponents*md.nofatoms);
   ComplexMatrix dd(1,md.nofcomponents*md.nofatoms,1,md.nofcomponents*md.nofatoms);
   ComplexMatrix cc(1,md.nofcomponents*md.nofatoms,1,md.nofcomponents*md.nofatoms);
   ComplexMatrix chi(1,md.nofcomponents*ini.mf.n()*md.nofatoms,1,md.nofcomponents*ini.mf.n()*md.nofatoms);
   ComplexMatrix Ac(1,md.nofcomponents*ini.mf.n()*md.nofatoms,1,md.nofcomponents*ini.mf.n()*md.nofatoms);
   ComplexMatrix Acinv(1,md.nofcomponents*ini.mf.n()*md.nofatoms,1,md.nofcomponents*ini.mf.n()*md.nofatoms);
   ComplexMatrix Bc(1,md.nofcomponents*ini.mf.n()*md.nofatoms,1,md.nofcomponents*ini.mf.n()*md.nofatoms);
   Ac=0;Bc=0;
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
   s=(ini.mf.nb()*ini.mf.nc()*(i1-1)+ini.mf.nc()*(j1-1)+k1-1)*md.nofcomponents*md.nofatoms;
   cc=0; dd=0;
   for(l1=1;l1<=md.nofatoms;++l1){for(i=1;i<=md.nofcomponents;++i){
   Ac(s+md.nofcomponents*(l1-1)+i,s+md.nofcomponents*(l1-1)+i)=1; // set diagonal elements 1 (make Ac a unit matrix)
   if (md.delta(i1,j1,k1)(l1)>SMALL)
    { //normal inelastic intensity
   cc(md.nofcomponents*(l1-1)+i,md.nofcomponents*(l1-1)+i)=1.0/(md.delta(i1,j1,k1)(l1)-z);
   dd(md.nofcomponents*(l1-1)+i,md.nofcomponents*(l1-1)+i)=1.0/(md.delta(i1,j1,k1)(l1)+z);
    }else{ 
     //quasielastic intensity ... artificially we introduce a splitting epsilon !!! compare Jensen 91 p 158
     if (md.delta(i1,j1,k1)(l1)<0) 
           {// this is when transition is between the same states ---> no dd
            cc(md.nofcomponents*(l1-1)+i,md.nofcomponents*(l1-1)+i)=-eps/z;
	    dd(md.nofcomponents*(l1-1)+i,md.nofcomponents*(l1-1)+i)=0.0;
           }else{
            cc(md.nofcomponents*(l1-1)+i,md.nofcomponents*(l1-1)+i)=-eps/z;
            dd(md.nofcomponents*(l1-1)+i,md.nofcomponents*(l1-1)+i)=eps/z;
           }
    }
   }}
   chi0c=md.M(i1,j1,k1)*cc+md.M(i1,j1,k1).Transpose()*dd; 
   for(i=1;i<=md.nofcomponents*md.nofatoms;++i){for(j=1;j<=md.nofcomponents*md.nofatoms;++j){Bc(s+i,s+j)=chi0c(i,j);}}
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
   ss=(ini.mf.nb()*ini.mf.nc()*(i2-1)+ini.mf.nc()*(j2-1)+k2-1)*md.nofcomponents*md.nofatoms;
   cc=chi0c*J.mati(J.in(i1,j1,k1),J.in(i2,j2,k2));
    for(i=1;i<=md.nofcomponents*md.nofatoms;++i){for(j=1;j<=md.nofcomponents*md.nofatoms;++j){
      Ac(s+i,ss+j)-=cc(i,j);
    }}   
  }}}
 }}}


 chi=Ac.Inverse()*Bc;

//myPrintComplexMatrix(stdout,Ac); 

 // determine chi'' and S (bose factor)

   complex<double> im(0,1.0);
   ComplexMatrix S(1,md.nofcomponents*ini.mf.n()*md.nofatoms,1,md.nofcomponents*ini.mf.n()*md.nofatoms);
   complex<double> bose;
   bose=1.0/(1.0-exp(-z*(1.0/KB/ini.T)));
   S=bose/(im)*(chi-chi.Transpose().Conjugate());
   
 // polarization factor
// neutrons only sense first 3x3 part of S !! - this is taken into account by setting 0 all
// higher components in the polarization factor !!!
 Matrix pol(1,md.nofcomponents,1,md.nofcomponents);
 Vector qxyz(1,3);
 pol=0;
    qxyz(1)=hkl(1)/inputpars.a; // only correct for ortholattices !!!!
    qxyz(2)=hkl(2)/inputpars.b;
    qxyz(3)=hkl(3)/inputpars.c;
    for(i=1;i<=3;++i){pol(i,i)=1.0;
    for(j=1;j<=3;++j){pol(i,j)-=qxyz(i)*qxyz(j)/(qxyz*qxyz);
    }}
    QQ=Norm(qxyz)*2*PI;

 //multiply polarization factor, formfactor and debeywallerfactor
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 for(l1=1;l1<=md.nofatoms;++l1){
   s=((((i1-1)*ini.mf.nb()+(j1-1))*ini.mf.nc()+(k1-1))*md.nofatoms+(l1-1))*md.nofcomponents;
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  for(l2=1;l2<=md.nofatoms;++l2){
   ss=((((i2-1)*ini.mf.nb()+(j2-1))*ini.mf.nc()+(k2-1))*md.nofatoms+(l2-1))*md.nofcomponents;
    for(i=1;i<=md.nofcomponents;++i){for(j=1;j<=md.nofcomponents;++j){
      S(s+i,ss+j)*=pol(i,j);
      S(s+i,ss+j)*=(*inputpars.jjj[l1]).gJ/2.0*(*inputpars.jjj[l1]).debeywallerfactor(QQ)*(*inputpars.jjj[l1]).F(QQ); // and formfactor + debey waller factor
      S(s+i,ss+j)*=(*inputpars.jjj[l2]).gJ/2.0*(*inputpars.jjj[l2]).debeywallerfactor(QQ)*(*inputpars.jjj[l2]).F(QQ); // and formfactor + debey waller factor
    }}   
  }
  }}}
 }
 }}}

 // determine dsigma
 //divide by number of crytallographic unit cells  (ini.mf.n()) in magnetic unit cell
intensity=abs(Sum(S))/ini.mf.n()/PI/2.0*3.65/4.0/PI; 

// here should be entered factor  k/k' + absolute scale factor
if (ini.ki==0)
{ki=sqrt(ini.kf*ini.kf+0.4811*en);
intensity*=ini.kf/ki;
}
else
{if (ini.ki*ini.ki-0.4811*en<0)
   {fprintf(stderr,"ERROR mcdisp - calculation of intensity: energy transfer %g meV cannot be reached with ki=const=%g/A at (%g,%g,%g)",en,ini.ki,hkl(1),hkl(2),hkl(3));
                            exit(EXIT_FAILURE);}
kf=sqrt(ini.ki*ini.ki-0.4811*en);
intensity*=kf/ini.ki;
}


return intensity;	
}

