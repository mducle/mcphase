/***********************************************************************
 *
 * mcdisp.c - program to calculate the dispersion of magnetic excitations
 *
 ***********************************************************************/

#define SMALL 1e-4    // deviation from single ion gap delta to take energy into account as not being equal to
                      // delta and therefore being included into output 
#define KB 0.0862     // Boltzmanns constant in mev/K

#include <mcdisp.h>
#include "../version"
#include "myev.c"
#include "intcalc.c"


// procedure to calculate the dispersion
void dispcalc(inimcdis & ini,par & inputpars,int do_Erefine,int do_jqfile, int do_verbose,int maxlevels,double epsilon)
{ int i,j,k,l,ll,s,ss,i1,i2,j1,j2,k1,k2,l1,l2,m,n,nn,norig;
  FILE * fout;
  FILE * fout1;
  FILE * foutds;
  FILE * foutdstot;
  FILE * foutds1;
  FILE * jqfile;
  double E;
  double sta=0;
  Vector hkl(1,3),q(1,3),lambda(1,ini.nofcomponents);
  Vector mf(1,ini.nofcomponents);
  int jmin;
  int noftransitions[inputpars.nofatoms+1]; // vector to remember how many transitions are on each atom
  int offset[inputpars.nofatoms+1]; // vector to remember where higher  transitions are stored 
                                    // (as "separate ions on the same unit cell position")
  mf=0;
   int sort=0;int maxiter=1000000;
  time_t curtime;
  struct tm *loctime;
  float d;
  
  // transition matrix Mij
  ComplexMatrix Mijkl(1,ini.nofcomponents,1,ini.nofcomponents);
  // transformation matrix Uij
  ComplexMatrix Uijkl(1,ini.nofcomponents,1,ini.nofcomponents);

  // enlarge nofatoms in inputpars, so that at each lattice site there  are as many
  // "atoms" as there are single ion transitions and put the corresponding transitionnumber
  // to those atoms
  norig=inputpars.nofatoms;nn=norig;  
  for(l=1;l<=norig;++l)
  {fprintf(stdout,"trying dmcalc for ion %i:\n",l);
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(1,1,1)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s in first 
                                                              //crystallographic unit of magnetic unit cell
   i=(*inputpars.jjj[l]).dmcalc(ini.T,mf,Mijkl,d); 
      // here Mijkl is a nxn matrix n ... numberofcomponents
   noftransitions[l]=1;offset[l]=inputpars.nofatoms;
    while (ini.emin>d||d>ini.emax) //only consider transition if it is in interval emin/emax
     {//first and following  transitions out of energy range ... do not consider them
     ++(*inputpars.jjj[l]).transitionnumber;
     (*inputpars.jjj[l]).dmcalc(ini.T,mf,Mijkl,d);
     if((*inputpars.jjj[l]).transitionnumber>i){fprintf(stderr,"ERROR mcdisp.ini: no transition with energy in range [emin,emax]=[%g,%g] found\n (within first crystallographic unit of magnetic unit cell)\n please increase energy range [emin,emax] in mcdisp.ini\n",ini.emin,ini.emax);
                            exit(EXIT_FAILURE);}
     }
   jmin=(*inputpars.jjj[l]).transitionnumber;
   for(j=jmin+1;j<=i;++j) // for every transition add new "atom" to list
   {
      if (noftransitions[l]==maxlevels){if (do_verbose) fprintf(stdout,"Maximum number of transitions for ion %i reached\n",l);
                       break;} //check if number of transitions  bigger than maximal number 
                                //(given by -max option in command line)
      
        (*inputpars.jjj[l]).transitionnumber=j; // try calculation for transition  j
      (*inputpars.jjj[l]).dmcalc(ini.T,mf,Mijkl,d);
        (*inputpars.jjj[l]).transitionnumber=jmin; // put back transition number for 1st transition
      
      if (ini.emin<d&&d<ini.emax) //only consider transition if it is in interval emin/emax
     {++noftransitions[l]; 
      nn=inputpars.newatom(inputpars.jjj[l]);
      (*inputpars.jjj[nn]).transitionnumber=j;
     }else
     {if (do_verbose)fprintf(stdout," .... this transition is not considered because it is out of interval [emin,emax]=[%g,%g]meV\n",ini.emin,ini.emax);
     }
   }  
  }  

  //calculate single ion properties of every atom in magnetic unit cell
  mdcf md(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),inputpars.nofatoms,ini.nofcomponents);

    // calculate delta(single ion excitation energy), 
    // Malphabeta(transition matrix elements) for every atom in unit cell
    for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
    md.U(i,j,k)=0; // initialize transformation matrix U
    md.M(i,j,k)=0; // initialize matrix M
    md.lambda(i,j,k)=0; // and sqrt(gamma^s) matrix lambda
    ss=norig;
    for(l=1;l<=norig;++l) 
     {
      // do calculation for atom s=(ijkl)
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s
      fprintf(stdout,"transitions of ion %i:\n",l);
      nn=noftransitions[l];
       //enlarge set of atoms by correct number
      for(ll=1;ll<=nn;++ll)
      {s=l;
       if(ll>1) {++ss;s=ss;} //determine s... "ion number" of transition ll in (original) ion l  
     // calculate matrix Ms and gap deltas
     (*inputpars.jjj[s]).dmcalc(ini.T,mf,Mijkl,d);
     md.delta(i,j,k)(s)=d;
     // diagonalizeMs to get unitary transformation matrix Us
     myEigenSystemHermitean (Mijkl,lambda,Uijkl,sort=1,maxiter); 
     if (lambda(ini.nofcomponents)>0&&fabs(lambda(ini.nofcomponents-1))<SMALL) 
                           // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
                           {if (fabs(d)>SMALL)
			    {md.lambda(i,j,k)(ini.nofcomponents*(s-1)+ini.nofcomponents,ini.nofcomponents*(s-1)+ini.nofcomponents)=sqrt(lambda(ini.nofcomponents));// lambda(ini.nofcomponents)=sqr(gamma^s)
                            }else{ //quasielastic line needs lambda=SMALL .... 
			     if (d<0){md.delta(i,j,k)(s)=-SMALL;}else{md.delta(i,j,k)(s)=SMALL;}//set fix delta but keep sign
  			     md.lambda(i,j,k)(ini.nofcomponents*(s-1)+ini.nofcomponents,ini.nofcomponents*(s-1)+ini.nofcomponents)=sqrt(SMALL*lambda(ini.nofcomponents));
			    }
			   }else 
                           {fprintf(stderr,"ERROR eigenvalue of single ion matrix <0: ev1=%g ev2=%g ev3=%g ... evn=%g\n",lambda(1),lambda(2),lambda(3),lambda(ini.nofcomponents));
                            exit(EXIT_FAILURE);}
        for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){
        md.U(i,j,k)(ini.nofcomponents*(s-1)+m,ini.nofcomponents*(s-1)+n)=Uijkl(m,n);
        md.M(i,j,k)(ini.nofcomponents*(s-1)+m,ini.nofcomponents*(s-1)+n)=Mijkl(m,n);
        }}    
      }
     }
    }}}
   
//initialize output files
  errno = 0;
  printf("saving mcdisp.qom\n");
  fout = fopen_errchk ("./results/mcdisp.qom","w");
  fprintf (fout, "#{%s ",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
          fprintf (fout, "#dispersion \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energies[meV] > intensities [barn/sr/f.u.]   f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");

   if (do_Erefine==1){
          errno = 0;
          foutds = fopen_errchk ("./results/mcdisp.dsigma","w");
          printf("saving mcdisp.dsigma\n");
          fprintf (foutds, "#{%s ",MCDISPVERSION);
          curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutds);
          fprintf (foutds, "#Scattering Cross Section \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energy[meV] dsigma/dOmegadE' [barn/mev/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
          foutdstot = fopen_errchk ("./results/mcdisp.dsigma.tot","w");
          printf("saving mcdisp.dsigma.tot\n");
          fprintf (foutdstot, "#{%s ",MCDISPVERSION);
          curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutdstot);
          fprintf (foutdstot, "#Total Scattering Cross Section in energy range [emin;emax]=[%g;%g]\n#Ha[T] Hb[T] Hc[T] T[K] h k l  dsigma/dOmeg [barn/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n",ini.emin,ini.emax);
                     }  
 
// initialize file with jq matrix
if (do_jqfile==1)
{  printf("saving mcdisp.jq\n");
  jqfile = fopen_errchk ("./results/mcdisp.jq","w");
  fprintf (jqfile, "#");
   fputs (asctime(loctime),jqfile);
   fprintf (jqfile, "#q=(hkl)\n #spin s() - spin s'()\n #3x3 matrix jss'(q) real im .... [meV]\n");
} 

//MAIN LOOP - do calculation of excitation energy for every Q vector     
int counter;
for (hkl(1)=ini.qmin(1);hkl(1)<=ini.qmax(1);hkl(1)+=ini.deltaq(1)){
for (hkl(2)=ini.qmin(2);hkl(2)<=ini.qmax(2);hkl(2)+=ini.deltaq(2)){
for (hkl(3)=ini.qmin(3);hkl(3)<=ini.qmax(3);hkl(3)+=ini.deltaq(3)){
 // transform hkl to primitive lattice
 if (ini.hkllist==1){counter=(int)hkl(1);
		     hkl(1)=ini.hkls[counter][1];
		     hkl(2)=ini.hkls[counter][2];
		     hkl(3)=ini.hkls[counter][3];
 		    }
 q=inputpars.r.Transpose()*hkl;
               
fprintf(stdout,"q=(%g,%g,%g)\n",hkl(1),hkl(2),hkl(3));
 if(do_verbose==1){fprintf(stdout,"Setting up J(q) matrix .... \n");}
 // calculate J(q)
 jq J(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),inputpars.nofatoms,ini.nofcomponents);
 jq Jl(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),inputpars.nofatoms,ini.nofcomponents);
 Vector d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
 Vector ij(1,3);
// int signa,signb,signc,sa,sb,sc,
  int sd;
 long int nofneighbours=0;
 ComplexMatrix jsss(1,ini.nofcomponents*inputpars.nofatoms,1,ini.nofcomponents*inputpars.nofatoms);
 complex<double> ipi(0,2*3.1415926535);

 // initialize Js,ss(Q)=0 (see manual for description of this matrix)
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s=J.in(i1,j1,k1);
   for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
   ss=J.in(i2,j2,k2); 
           J.mati(s,ss)= 0;// set (ini.nofcomponents*nofatoms) x (ini.nofcomponents*nofatoms) matrix Js,ss(q)=0
           Jl.mati(s,ss)= 0;// set Js,ss(q)=0 
   }}}
 }}}
 // calculate Js,ss(Q) summing up contributions from the l=1-paranz parameters
   int sl,sll,sublat;
   for(ll=1;ll<=inputpars.nofatoms;++ll){for(l=1;l<=(*inputpars.jjj[ll]).paranz;++l)
   { //sum up l.th neighbour interaction
   // 1. transform dn(l) to primitive lattice
    xyz=(*inputpars.jjj[ll]).dn[l];
    d=inputpars.rez*(const Vector&)xyz;
    for (i=1;i<=3;++i)d_rint(i)=rint(d(i));
   
   //2. sum up parameter
         for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
         s=J.in(i1,j1,k1);
         //next 2 lines special for treatment of assymetry in Rcu2 
	 // i.e. if j1 is 2n then  say: the neighbour is at -d, this ensures that
	 // for half of the atoms the interaction is counted at negative distances
	 // [this makes no error if the lattice is primitive!!!] but for RCu2 it makes
	 // a difference, because atoms at +-0.5c are not equal 
	 sd=-1;if (inputpars.r[2][1]==0.5&&(double)j1/2.0==integer(1.0*j1/2)) sd=1;
  	// if (inputpars.r[2][1]==0.0&&(double)(i1+j1)/2.0==integer(1.0*(i1+j1)/2)) sd=1;
         //calc ss (check in which magnetic cell ss neighbour lies)	 
         i=(int)(i1+sd*d_rint(1)-1); // calculate 
	 j=(int)(j1+sd*d_rint(2)-1);
	 k=(int)(k1+sd*d_rint(3)-1);
	 if (i>=0) ij(1)=integer(1.0*i/ini.mf.na())*ini.mf.na();
	 else      ij(1)=(integer(1.0*(i+1)/ini.mf.na())-1)*ini.mf.na();
	 if (j>=0) ij(2)=integer(1.0*j/ini.mf.nb())*ini.mf.nb();
	 else      ij(2)=(integer(1.0*(j+1)/ini.mf.nb())-1)*ini.mf.nb();
	 if (k>=0) ij(3)=integer(1.0*k/ini.mf.nc())*ini.mf.nc();
	 else      ij(3)=(integer(1.0*(k+1)/ini.mf.nc())-1)*ini.mf.nc();
	 ss=J.in(i-(int)ij(1)+1,j-(int)ij(2)+1,k-(int)ij(3)+1);
	   // sum up 
           jsss=0;
        
	  sl=(*inputpars.jjj[ll]).sublattice[l]; // the whole loop has also to be done 
                                           // for all the other transitions of sublattice sl
          // therefore calculate offset of the set of transitions
          for(sll=1;sll<=noftransitions[sl];++sll)    
          {sublat=sl;if(sll>1){sublat=sll+offset[sl]-1;}   
        
	     for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){ //this should also be ok for nofcomponents > 3 !!! (components 1-3 denote the magnetic moment)
            jsss(ini.nofcomponents*(ll-1)+m,ini.nofcomponents*(sublat-1)+n)=(*inputpars.jjj[ll]).jij[l](m,n);
           }} // but orbitons should be treated correctly by extending 3 to n !!
	  } 
          J.mati(s,ss)+=jsss*exp(ipi*(double)sd*(q*d));// increase Js,ss(q)=0
//          Jl.mati(s,ss)+=jsss*exp(ipi*sd*(q*d));// increase Js,ss(q)=0
          ++nofneighbours; // count neighbours summed up
	 }}}
      
    }}


if (do_jqfile==1){fprintf (jqfile, "#q=(%g, %g, %g) ",hkl(1),hkl(2),hkl(3));
                  fprintf(jqfile,"nofneighbours= %li\n",nofneighbours);}


if(do_verbose==1){fprintf(stdout,"Transform J(q) matrix  with U...\n");}

// transform J(s,ss) (with md.U) and multiply 
// (with eigenvalues sqrt(gamma) [here md.lambda]
// of matrix Malphabeta ) matrix J to J' [here Jl.mati]  ... compare manual
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s=ini.mf.in(i1,j1,k1);
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  ss=ini.mf.in(i2,j2,k2);
//  jsss=Jl.mati(s,ss);
  Jl.mati(s,ss)=md.lambda(i1,j1,k1)*md.U(i1,j1,k1).Conjugate().Transpose()*J.mati(s,ss)*md.U(i2,j2,k2)*md.lambda(i2,j2,k2);
/*  jsss=Jl.mati(ss,s);
  Jl.mati(ss,s)=jsss*md.U(i1,j1,k1)*md.lambda(i1,j1,k1);*/
/*if (do_jqfile==1){fprintf (jqfile, "#spin (%i*r1 %i*r2 %i*r3) - spin (%i*r1 %i*r2 %i*r3)\n",i1,j1,k1,i2,j2,k2);
         myPrintComplexMatrix(jqfile,J.mati(s,ss)); 
         myPrintComplexMatrix(stdout,J.mati(s,ss)); 
	// diagonalize Ac to get eigenvalues (biggest corresponds to Tn) !!!
         Vector Tn(1,3*ini.mf.n());
         myEigenSystemHermitean (J.mati(s,ss),Tn,jsss,sort,maxiter);
        fprintf(jqfile,"#eigenvalues(highest corresponds to Tn, predicted magstructure)\n");
         myPrintVector(jqfile,Tn); 
        fprintf(jqfile,"#eigenvectors(moment direction):\n");
         myPrintComplexMatrix(jqfile,jsss); 
         
	         }*/ //printout here gave errors when unit cell has more than one atom (not diagonalizable atrix
  }}}
 }}}

// calculate Ac
// Ac  is the matrix which is given in manual chapter 8.1 (eq (26) ff) -- diagonalization gives omega_r and Tau

   ComplexMatrix Ac(1,ini.mf.n()*inputpars.nofatoms,1,ini.mf.n()*inputpars.nofatoms);
   ComplexMatrix Tau(1,ini.mf.n()*inputpars.nofatoms,1,ini.mf.n()*inputpars.nofatoms);
   ComplexMatrix J_Q(1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms,1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms);
   Ac=0;J_Q=0;
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
   s=(((i1-1)*ini.mf.nb()+(j1-1))*ini.mf.nc()+(k1-1))*inputpars.nofatoms;
   for(l1=1;l1<=inputpars.nofatoms;++l1){
      Ac(s+l1,s+l1)=md.delta(i1,j1,k1)(l1);}

   for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
   ss=(((i2-1)*ini.mf.nb()+(j2-1))*ini.mf.nc()+(k2-1))*inputpars.nofatoms;

    for(l1=1;l1<=inputpars.nofatoms;++l1){ 
    for(l2=1;l2<=inputpars.nofatoms;++l2){ 
     Ac(s+l1,ss+l2)-=Jl.mati(Jl.in(i1,j1,k1),Jl.in(i2,j2,k2))(ini.nofcomponents*(l1-1)+ini.nofcomponents,ini.nofcomponents*(l2-1)+ini.nofcomponents);
                                             //nofcomponents^th dimension corresponds to 1st in manual 
					     // and it is only necessary to take into 
					     // acount this dimension!!
     }}
    for(i=1;i<=ini.nofcomponents*inputpars.nofatoms;++i){for(j=1;j<=ini.nofcomponents*inputpars.nofatoms;++j){
      J_Q(s+i,ss+j)+=J.mati(J.in(i1,j1,k1),J.in(i2,j2,k2))(i,j);
    }}

  }}}
 }}}

// printout fouriertransform of matrix jq
if (do_jqfile==1){//fprintf (jqfile, "#spin (%i*r1 %i*r2 %i*r3) - spin (%i*r1 %i*r2 %i*r3)\n",i1,j1,k1,i2,j2,k2);
         myPrintComplexMatrix(jqfile,J_Q); 
	// diagonalize JQ to get eigenvalues (biggest corresponds to Tn) !!!
         Vector Tn(1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms);
         ComplexMatrix eigenvectors(1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms,1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms);
         myEigenSystemHermitean (J_Q,Tn,eigenvectors,sort,maxiter);
        fprintf(jqfile,"#eigenvalues(highest corresponds to Tn, predicted magstructure)\n");
         myPrintVector(jqfile,Tn); 
        fprintf(jqfile,"#eigenvectors(moment direction):\n");
         myPrintComplexMatrix(jqfile,eigenvectors); 
	         }


 if(do_verbose==1){fprintf(stdout,"diagonalizing %ix%i matrix...\n",ini.mf.n()*inputpars.nofatoms,ini.mf.n()*inputpars.nofatoms);}
   // diagonalize Ac to get energies  and eigenvectors !!!
   Vector En(1,ini.mf.n()*inputpars.nofatoms);
   Vector ints(1,ini.mf.n()*inputpars.nofatoms);
//   myEigenValuesHermitean (Ac,En,sort,maxiter);
   myEigenSystemHermitean (Ac,En,Tau,sort,maxiter);

 if(do_verbose==1){fprintf(stdout,"saving eigenvalues to mdisp.qom\n");}
   int dim=3;
   if (ini.hkllist==1){dim=(int)ini.hkls[counter][0]-3;}
         Vector dd(1,dim);  dd+=100000.0;
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
   for (i=1;i<=ini.mf.n()*inputpars.nofatoms;++i){
//    i2=1;//check if eigenvalue is equal to some delta(s)
//    for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
//    if(abs(En(i)-md.delta(i1,j1,k1))<SMALL)i2=0;
//    }}} //-if not: print it out and sum up standard deviation
//    if (i2==1){
	       fprintf (fout, " %4.4g ",En(i));
                if(do_verbose==1){fprintf(stdout, " %4.4g ",En(i));}
	       if (En(i)<ini.emin) sta*=1.1-En(i)+ini.emin;
               if (ini.hkllist==1)
	       {double test; // add to sta distance to nearest measured peak squared
	        for (j1=1;j1<=ini.hkls[counter][0]-3;++j1)
	        {if ((test=fabs(En(i)-ini.hkls[counter][j1+3]))<dd(j1))dd(j1)=test;}
	       }
             // }
    }
   sta+=dd*dd;


   // calculate and printout intensities [the energies have already
   // been printed out above, so any refinement of energies during intcalc
   // is not included in the output file]

  if(do_verbose==1){fprintf(stdout,"\ncalculating  intensities approximately ...\n");}
                  fprintf (fout, " > ");
                  for (i=1;i<=ini.mf.n()*inputpars.nofatoms;++i)
		  {//double maxupshift=1e10; if (i<ini.mf.n()) {maxupshift=En(i+1)-En(i);}
		   //double maxdownshift=-En(i); if(i>1){maxdownshift=En(i-1)-En(i);}		 
		   //i2=1;//check if eigenvalue is equal to some delta(s)
                   //for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
                   //for(l1=1;l1<=inputpars.nofatoms;++l1){    
		   //if(abs(En(i)-md.delta(i1,j1,k1))<SMALL)i2=0; 
                   //if(En(i)+maxupshift+SMALL>md.delta(i1,j1,k1)(l1)&&md.delta(i1,j1,k1)(l1)>=En(i)+SMALL){maxupshift=md.delta(i1,j1,k1)(l1)-SMALL-En(i);}
                   //if(En(i)-maxdownshift-SMALL<md.delta(i1,j1,k1)(l1)&&md.delta(i1,j1,k1)(l1)<=En(i)-SMALL){maxdownshift=md.delta(i1,j1,k1)(l1)+SMALL-En(i);}
		   //}
		   //}}} 
		     ints(i)=intcalc_approx(Tau,i,En(i),ini,inputpars,J,q,hkl,md,do_verbose);
                     //printout rectangular function to .mdcisp.qom
	             fprintf (fout, " %4.4g ",ints(i));
                     if(do_verbose==1){fprintf(stdout, " %4.4g ",ints(i));}

		   }

              //initialize output file for display
            errno = 0;
            fout1 = fopen_errchk ("./results/.mcdisp.qom","w");
            fprintf (fout1, "#{%s ",MCDISPVERSION);
            curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout1);
            fprintf (fout1,"\n#Ha[T] Hb[T] Hc[T] T[K] h k l  energies[meV] intensities [barn/meV/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
		     if (do_Erefine==0) epsilon=(Max(En)-Min(En))/100;
                  for (i=1;i<=ini.mf.n()*inputpars.nofatoms;++i)
		    { 
		     if (ints(i)>SMALL)
		      {for (E=0;E<=ints(i)/epsilon;E+=ints(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	                fprintf (fout1, " %4.4g %4.4g \n",En(i)-epsilon+E*epsilon*epsilon/ints(i),E);
		       }
		       for (E=ints(i)/epsilon;E>=0;E-=ints(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	                fprintf (fout1, " %4.4g %4.4g \n",En(i)+epsilon-E*epsilon*epsilon/ints(i),E);
		       }
                       fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	               fprintf (fout1, " %4.4g 0 \n",En(i)+epsilon);
		      }
		    }
	  fclose(fout1);   
	    
		    
                if(do_verbose==1){fprintf(stdout, "\n");}
		    
   // do refinement of energies by output of scattering cross section vs enrgy transfer if required
  if (do_Erefine==1){double totint=0;
                if(do_verbose==1){fprintf(stdout, "refining calculation with exact calculation of energy dependence of scattering cross section\n");}
          errno = 0;
          foutds1 = fopen_errchk ("./results/.mcdisp.dsigma","w");
          fprintf (foutds1, "#{%s ",MCDISPVERSION);
          curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutds1);
          fprintf (foutds1, "#Scattering Cross Section \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energy[meV] dsigma/dOmegadE'[barn/mev/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
		     double intensity;
		     intensity=intcalc(ini.emin,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds1, " %4.4g %4.4g \n",ini.emin,intensity);
		     intensity=intcalc(ini.emax,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds1, " %4.4g %4.4g \n",ini.emax,intensity);
	  fclose(foutds1);
	  for(E=ini.emin;E<=ini.emax;E+=epsilon/2)
	   {
		     intensity=intcalc(E,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
		     totint+=intensity*epsilon/2;

          foutds1 = fopen_errchk ("./results/.mcdisp.dsigma","a");
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds1, " %4.4g %4.4g \n",E,intensity);
          fclose(foutds1);	   
                     fprintf (foutds, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds, " %4.4g %4.4g \n",E,intensity);
	   }
                     fprintf (foutdstot, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutdstot, " %4.4g \n",totint);
                     }  

   
                 
   fprintf (fout, "\n");
   if (ini.hkllist==1){hkl(1)=(double)counter;}
}}}
    fprintf (fout, "#sta= %8.6g ",sta);
    fprintf (stdout, "#sta= %8.6g ",sta);
   
    fclose(fout);if (do_jqfile==1) {fclose(jqfile);}
                 if (do_Erefine==1){fclose(foutdstot);fclose(foutds);}
}

//*************************************************************************************************
// main program
int main (int argc, char **argv)
{int i,do_Erefine=0,do_jqfile=0,do_verbose=0,maxlevels=100;
 char * spinfile="mcdisp.mf"; //default spin-configuration-input file
 double epsilon; //imaginary part of omega to avoid divergence

// check command line and initialize parameters ini
for (i=1;i<=argc-1;++i){
   if(strcmp(argv[i],"-r")==0) {do_Erefine=1; if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -r needs argument epsilon\n");exit(EXIT_FAILURE);}
		                                                epsilon=strtod(argv[i+1],NULL);++i;
							        fprintf(stdout,"epsilon= %g\n",epsilon);
				     }		
         else {if(strcmp(argv[i],"-jq")==0) do_jqfile=1;       
               else {if(strcmp(argv[i],"-v")==0||strcmp(argv[i],"-verbose")==0) do_verbose=1;       
                     else {if(strcmp(argv[i],"-max")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -max needs argument\n");exit(EXIT_FAILURE);}
		                                                maxlevels=(int)strtod(argv[i+1],NULL);++i;
							        fprintf(stdout,"maximum number of single ion excitations taken into account (starting with lowest energy): %i\n",maxlevels);
							       }       
         		     else{spinfile=argv[i];}
		          }
		    }
	      }
	
    }
inimcdis ini("mcdisp.ini",spinfile);

if (argc > 6) {ini.errexit();}
  // as class load  parameters from file
  par inputpars("./mcphas.j");
  

  //calculate dispersion and save on file
  dispcalc(ini,inputpars,do_Erefine,do_jqfile,do_verbose,maxlevels,epsilon);
  exit(0);
 
}





