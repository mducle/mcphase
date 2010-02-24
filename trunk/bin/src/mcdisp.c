/***********************************************************************
 *
 * mcdisp - program to calculate the dispersion of magnetic excitations
 *
 * reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751
 *            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400
 ***********************************************************************/

#define SMALL 1e-6    // deviation from single ion gap delta to take energy into account as not being equal to
                      // delta and therefore being included into output 
		      // transitions of single ions less then SMALL have in Mijkl wn/kT instead of wn-wn'
		      // !!! must match SMALL in jjjpar.cpp !!!!
#define SMALLINT 1e-4 // small intensity treshhold
#define SMALLEDIF 1e-5 // small difference in calculation of transition energy
                       // used to give error if recalculation of mcdisp.trs
		       // energies gives different results than file
#define KB 0.0862     // Boltzmanns constant in mev/K
#define MAXNOFCHARINLINE 1024

#include <mcdisp.h>
#include "../../version"
#include "myev.c"

#ifdef _THREADS
// ----------------------------------------------------------------------------------- //
// Defines to ease interchange between linux and windows thread codes...
// ----------------------------------------------------------------------------------- //
#ifdef __linux__
#include <pthread.h>
#define MUTEX_LOCK    pthread_mutex_lock
#define MUTEX_UNLOCK  pthread_mutex_unlock
#define MUTEX_TYPE    pthread_mutex_t
#define MUTEX_INIT(m) pthread_mutex_init (&m, NULL)
#define EVENT_TYPE    pthread_cond_t
#define EVENT_INIT(e) pthread_cond_init (&e, NULL)    
#define EVENT_SIG(e)  pthread_cond_signal (&e)
#else
#include <windows.h>
#define MUTEX_LOCK    EnterCriticalSection
#define MUTEX_UNLOCK  LeaveCriticalSection
#define MUTEX_TYPE    CRITICAL_SECTION
#define MUTEX_INIT(m) InitializeCriticalSection (&m)
#define EVENT_TYPE    HANDLE
#define EVENT_INIT(e) e = CreateEvent (NULL, TRUE, FALSE, NULL)
#define EVENT_SIG(e)  SetEvent(e)
#endif
#define NUM_THREADS 2

// ----------------------------------------------------------------------------------- //
// Declares a struct to store all the information needed for each disp_calc iteration
// ----------------------------------------------------------------------------------- //
typedef struct{
   ComplexMatrix **chi, **chibey, **S, **Sbey, **Ec;
   Matrix **pol, **polICIC, **polICn, **polnIC;
   mfcf **ev_real, **ev_imag, **eev_real, **eev_imag;
   ComplexMatrix **Tau;
   Vector hkl;
   inimcdis **ini;
   par **inputpars;
   mdcf **md;
   int thread_id;
} intcalcapr_thread_data;
class intcalcapr_input { public:
   Vector hkl;
   int thread_id;
   int dimA, level, do_verbose;
   double En, intensity, intensitybey, QQ;
   intcalcapr_input(int _dimA, int _tid, int _level, int _doverb, double _En, par *pars_in, inimcdis *ini_in) 
   { 
      thread_id = _tid; dimA = _dimA; level = _level; do_verbose = _doverb; En = _En;
   }
};
// ----------------------------------------------------------------------------------- //
// Declares these variables global, so all threads can see them
// ----------------------------------------------------------------------------------- //
intcalcapr_thread_data thrdat;
MUTEX_TYPE mutex_loop;
MUTEX_TYPE mutex_index;
EVENT_TYPE checkfinish;
#endif // if _THREADS

#include "mcdisp_intcalc.c"

int index_s(int i,int j,int k,int l, int t, const mdcf & md, const inimcdis & ini)
{int s=0,i1,j1,k1;
// calculates the index of the Matrix A given 
// ijk ... index of crystallographic unit cell in magnetic unit cell
// l   ... number of atom in crystallographic cell
// t   ... transitionnumber
 for(i1=1;i1<i;++i1){
 for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s+=md.baseindex_max(i1,j1,k1);
 }}}
 i1=i;
 for(j1=1;j1<j;++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s+=md.baseindex_max(i1,j1,k1);
 }}
 j1=j;
 for(k1=1;k1<k;++k1){
 s+=md.baseindex_max(i1,j1,k1);
 }
 s+=md.baseindex(i,j,k,l,t);
 return s;
}

//void s2ijklt(int s, int & i, int & j, int & k, int & l, int & t, mdcf & md, inimcdis & ini)
//{// the reverse of index_s function, returns i j k l t given the index_s
//int ss=0,sss;
// ijk ... index of crystallographic unit cell in magnetic unit cell
// l   ... number of atom in crystallographic cell
// t   ... transitionnumber
 //for(i=1;i<=ini.mf.na();++i){
 //for(j=1;j<=ini.mf.nb();++j){
 //for(k=1;k<=ini.mf.nc();++k){
 //sss=md.baseindex_max(i,j,k);
 //if (ss+sss<	s){ss+=sss;}else{md.baseindex2ltn(s-ss,i,j,k,l,t);return;}
 //}}}
 //fprintf(stderr,"mcdisp internal error: index s too large in function s2ijklt\n");exit(1);
//}

void sortE(Vector & d,ComplexMatrix & z)
{       int i,j,k;
    double p;
    complex <double> p1;
    // lowest and highest column index of the matrix
    int lo = z.Clo();
    int hi = z.Chi();

    for (i = lo; i <= hi; i++) {
	k = i;
	p = d(i);
	for (j = i+1; j <= hi; j++)
	    if (d(j) < p) {
		k = j;
		p = d(j);
	    }
	if (k != i) {
	    d(k) = d(i);
	    d(i) = p;
	    for (j = lo; j <= hi; j++) {
		p1 = z(j,i);
		z(j,i) = z(j,k);
		z(j,k) = p1;
	    }
	}
    }
}

// procedure to calculate the dispersion
void dispcalc(inimcdis & ini,par & inputpars,int do_Erefine,int do_jqfile,int do_createtrs,int do_readtrs, int do_verbose,int maxlevels,double minE,double maxE,double epsilon, const char * filemode)
{ int i,j,k,l,ll,s,ss,i1,i2,j1,j2,k1,k2,l1,l2,t1,t2,b,bb,m,n,tn;
  FILE * fin;
  FILE * fout;
  FILE * foutqei;
  FILE * foutqev;
  FILE * foutqee;
  FILE * fout1;
  FILE * foutds;
  FILE * foutdstot;
  FILE * foutds1;
  FILE * jqfile;
  float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
  int do_gobeyond=1;
  double E;
  double sta=0,sta_int=0;
  double jqsta=-1.0;
  double jq0=0;
  Vector hkl(1,3),q(1,3);
  Vector mf(1,ini.nofcomponents),extmf(1,ini.extended_eigenvector_dimension);
  int jmin,tl,tll;
  IntVector noftransitions(1,inputpars.nofatoms); // vector to remember how many transitions are on each atom
  int offset[inputpars.nofatoms+1]; // vector to remember where higher  transitions are stored 
                                    // (as "separate ions on the same unit cell position")
  mf=0;extmf=0;
   int sort=0;int maxiter=1000000;
  time_t curtime;
  struct tm *loctime;
  float d;
  
  Vector gamma(1,ini.nofcomponents);Vector extgamma(1,ini.extended_eigenvector_dimension);
  complex<double> imaginary(0,1);
  // transition matrix Mij
  ComplexMatrix Mijkl(1,ini.nofcomponents,1,ini.nofcomponents);
 // extended transition matrix Mij
  ComplexMatrix extMijkl(1,ini.extended_eigenvector_dimension,1,ini.extended_eigenvector_dimension);
  // transformation matrix Uij
  ComplexMatrix Uijkl(1,ini.nofcomponents,1,ini.nofcomponents);
// extended transformation matrix Uij
  ComplexMatrix extUijkl(1,ini.extended_eigenvector_dimension,1,ini.extended_eigenvector_dimension);

  //calculate single ion properties of every atom in magnetic unit cell
  mdcf md(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),inputpars.nofatoms,ini.nofcomponents);

  //initialize output of transitions
 if (do_readtrs==0)
 {printf("saving ./results/mcdisp.trs\n");
  fout = fopen_errchk ("./results/mcdisp.trs","w");
   fprintf(fout, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcdisp.trs-->\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(fout,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"#(*)The unpolarized powder average neutron cross section sigma for each transition \n");
   fprintf(fout,"#   is calculated neglectingthe formfactor, the Debye Wallerfactor, factor k'/k.\n");
   fprintf(fout,"# T= %g K\n",ini.T);
   fprintf(fout,"#*********************************************************************\n");
          fprintf (fout, "#i j k ionnr transnr energy |gamma_s|  sigma [barn](*)\n");
 
  for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.nofatoms;++l){
   fprintf(stdout,"trying dmcalc for ion %i in crystallographic unit cell %i %i %i:\n",l,i,j,k);
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s in first 
                                                              //crystallographic unit of magnetic unit cell
   
   md.est_ini(i,j,k,l,(*inputpars.jjj[l]).eigenstates(mf,ini.T)); 
   (*inputpars.jjj[l]).transitionnumber=1;
   fprintf(stdout,"transition number %i: ",(*inputpars.jjj[l]).transitionnumber);
   d=1e10;i1=(*inputpars.jjj[l]).dmcalc(ini.T,mf,Mijkl,d,md.est(i,j,k,l));

      // here Mijkl is a nxn matrix n ... numberofcomponents
   noftransitions(l)=0;
    while ((minE>d||d>maxE)&&(minE>-d||-d>maxE)) //only consider transition if it is in interval minE/maxE
     {//first and following  transitions out of energy range ... do not consider them
     fprintf(stdout," .... transition not stored because out of interval [minE,maxE]=[%g,%g]meV\n",minE,maxE);
     ++(*inputpars.jjj[l]).transitionnumber;
     fprintf(stdout,"transition number %i: ",(*inputpars.jjj[l]).transitionnumber);
     d=1e10;(*inputpars.jjj[l]).dmcalc(ini.T,mf,Mijkl,d,md.est(i,j,k,l));
     if((*inputpars.jjj[l]).transitionnumber>i1){fprintf(stderr,"ERROR mcdisp.par: no transition found within energy in range [minE,maxE]=[%g,%g] found\n (within first crystallographic unit of magnetic unit cell)\n please increase energy range in option -maxE and -minE\n",minE,maxE);
                            exit(EXIT_FAILURE);}
     }
     // calculate powder neutron intensities
     double intensityp=0, intensitym=0;int k1;
     if ((*inputpars.jjj[l]).gJ!=0)
     { for(k1=1;k1<=3;++k1){intensityp+=(*inputpars.jjj[l]).gJ*(*inputpars.jjj[l]).gJ*real(Mijkl(k1,k1));}}
     else
     { for(k1=1;k1<=3;++k1){intensityp+=4*real(Mijkl(2*k1-1,2*k1-1))+real(Mijkl(2*k1,2*k1));}}
     intensityp*=0.61072561;intensitym=intensityp;// prefactor for intensity in barn is 4*PI*2/3*0.54*0.54/4= 0.61072561
     if (d>SMALL){if(d/ini.T/KB<20){intensitym=-intensityp/(1-exp(d/ini.T/KB));intensityp/=(1-exp(-d/ini.T/KB));}}
                                  else{intensityp=intensityp*ini.T*KB;intensitym=intensityp;}
                

   jmin=(*inputpars.jjj[l]).transitionnumber;
//  if (do_verbose==1){fprintf(stdout,"Matrix M(s=%i %i %i)\n",i,j,k);
//                  myPrintComplexMatrix(stdout,Mijkl);fprintf(stdout,"...diagonalising\n");
//                    } 
     // diagonalizeMs to get unitary transformation matrix Us and eigenvalues gamma
     myEigenSystemHermitean (Mijkl,gamma,Uijkl,sort=1,maxiter); 
     if(minE<d&&d<maxE)
    { fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,jmin,d,gamma(ini.nofcomponents),intensityp);
     ++noftransitions(l);}
    if(d>=0&&minE<-d&&-d<maxE) // do not print negative energy transition if d<0 (d<0 means transiton to the same level)
    { fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,jmin,-d,gamma(ini.nofcomponents),intensitym);
    ++noftransitions(l);}

   for(j1=jmin+1;j1<=i1;++j1) 
   // for every transition add new "atom" to list ... changed to "setnumberoftransitions" for
   // ion l in cryst unit cell ijk
   {
      if (noftransitions(l)>=maxlevels){if (do_verbose) fprintf(stdout,"Maximum number of transitions for ion %i reached\n",l);
                       break;} //check if number of transitions  bigger than maximal number 
                                //(given by -max option in command line)
      
        (*inputpars.jjj[l]).transitionnumber=j1; // try calculation for transition  j
      fprintf(stdout,"transition number %i: ",(*inputpars.jjj[l]).transitionnumber);
      d=maxE;(*inputpars.jjj[l]).dmcalc(ini.T,mf,Mijkl,d,md.est(i,j,k,l));
        (*inputpars.jjj[l]).transitionnumber=jmin; // put back transition number for 1st transition
   //printf("noftransitions read by mcdisp: %i",i1);
      
      if ((minE<d&&d<maxE)||(minE<-d&&-d<maxE)) //only consider transition if it is in interval emin/emax
     { 
      // calculate powder neutron intensities
     intensityp=0; intensitym=0;
     if ((*inputpars.jjj[l]).gJ!=0)
     { for(k1=1;k1<=3;++k1){intensityp+=(*inputpars.jjj[l]).gJ*(*inputpars.jjj[l]).gJ*real(Mijkl(k1,k1));}}
     else
     { for(k1=1;k1<=3;++k1){intensityp+=4*real(Mijkl(2*k1-1,2*k1-1))+real(Mijkl(2*k1,2*k1));}}
     intensityp*=0.61072561;intensitym=intensityp;// prefactor for intensity in barn is 4*PI*2/3*0.54*0.54/4= 0.61072561
      if (d>SMALL){if(d/ini.T/KB<20){intensitym=-intensityp/(1-exp(d/ini.T/KB));intensityp/=(1-exp(-d/ini.T/KB));}}
                                  else{intensityp=intensityp*ini.T*KB;intensitym=intensityp;}
               
      // diagonalizeMs to get unitary transformation matrix Us and eigenvalues gamma
     myEigenSystemHermitean (Mijkl,gamma,Uijkl,sort=1,maxiter); 
     if(minE<d&&d<maxE)
     {fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,j1,d,gamma(ini.nofcomponents),intensityp);
      ++noftransitions(l);}
     if(d>=0&&minE<-d&&-d<maxE)// do not print negative energy transition if d<0 (d<0 means transiton to the same level)
     {fprintf(fout,"%i %i %i  %i     %i     %g  %g  %g\n",i,j,k,l,j1,-d,gamma(ini.nofcomponents),intensitym);
      ++noftransitions(l);}
     }else
     {fprintf(stdout," .... transition not stored because  out of interval [minE,maxE]=[%g,%g]meV\n",minE,maxE);
     }


   }  
  }}}}
   fclose(fout);
 }
  if (do_createtrs==1){fprintf(stdout,"single ion transition file ./results/mcdisp.trs created - please comment transitions which should not enter the calculation and restart with option -t\n");exit(0);}

  printf("\n#reading ./results/mcdisp.trs\n\n");
// read transitions to beconsidered from file
 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  fin = fopen_errchk ("./results/mcdisp.trs","rb");
  noftransitions=0;
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3])
    {l=(int)nn[4];++noftransitions(l);}
  }} 
     fclose(fin);
  //just read dimensions of matrices in md
    md.set_noftransitions(i,j,k,noftransitions);
  // for later use:
    for(l=1;l<=inputpars.nofatoms;++l){ //save eigenstates of ions (if possible)
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);} //mf ... mean field vector of atom s in first 
                                                              //crystallographic unit of magnetic unit cell
       if(do_readtrs!=0)md.est_ini(i,j,k,l,(*inputpars.jjj[l]).eigenstates(mf,ini.T)); // initialize ests if not already done above
       }

    md.U(i,j,k)=0; // initialize transformation matrix U
    md.M(i,j,k)=0; // initialize matrix M
    md.sqrt_gamma(i,j,k)=0; // and sqrt(gamma^s) matrix sqrt_gamma
 }}}

// determine the dimension of the matrix Ass' s,s'=1....dimA
int dimA=0;
for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 dimA+=md.baseindex_max(i1,j1,k1);
 }}}

// matrix E^s_alpha' used to store the coefficients for extending the eigenvector (see manual)
ComplexMatrix Ec(1,dimA,1,ini.extended_eigenvector_dimension);Ec=0;
  
 for(i=1;i<=ini.mf.na();++i){for(j=1;j<=ini.mf.nb();++j){for(k=1;k<=ini.mf.nc();++k){
  for(l=1;l<=inputpars.nofatoms;++l){
  fin = fopen_errchk ("./results/mcdisp.trs","rb");
  jmin=0;
  while (feof(fin)==0)
  {if ((i1=inputline(fin,nn))>=5)
   {if(i==(int)nn[1]&&j==(int)nn[2]&&k==(int)nn[3]&&l==(int)nn[4])
    {tn=(int)nn[5];++jmin;  
    // calculate delta(single ion excitation energy), 
    // Malphabeta(transition matrix elements)

      // do calculation for atom s=(ijkl)
      for(ll=1;ll<=ini.nofcomponents;++ll)
       {mf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);    //mf ... mean field vector of atom s
        extmf(ll)=ini.mf.mf(i,j,k)(ini.nofcomponents*(l-1)+ll);}

      fprintf(stdout,"#transition %i of ion %i of cryst. unit cell at pos  %i %i %i in mag unit cell:\n",tn,l,i,j,k);
      if(nn[6]<SMALL){fprintf(stdout,"#-");}else{fprintf(stdout,"#+");}
      
        j1=(*inputpars.jjj[l]).transitionnumber; // try calculation for transition  j
        (*inputpars.jjj[l]).transitionnumber=tn; // try calculation for transition  j
      d=1e10;(*inputpars.jjj[l]).dmcalc(ini.T,mf,Mijkl,d,md.est(i,j,k,l));
       //      myPrintComplexMatrix(stdout,Mijkl);
      
       // here we if required calculate the higher dimension matrix used to do the
       // extension of chi to higher value of (uncoupled) nofcomponents in intcalc_approx ... needed for chargedensityfluctuations, extended eigenvectors ...
             d=1e10;(*inputpars.jjj[l]).dmcalc(ini.T,extmf,extMijkl,d,md.est(i,j,k,l));
      
        (* inputpars.jjj[l]).transitionnumber=j1; // put back transition number for 1st transition

       j1=md.baseindex(i,j,k,l,jmin); 
      
       if(fabs((fabs(d)-fabs(nn[6]))/(fabs(nn[6])+1.0))>SMALLEDIF)
        {fprintf(stderr,"ERROR mcdisp: reading mcdisp.trs with transition energy delta %g meV different from internal calculation %g meV\n",nn[6],d);	 
         exit(EXIT_FAILURE);}
       md.delta(i,j,k)(j1)=nn[6]; // set delta
     // diagonalizeMs to get unitary transformation matrix Us
      myEigenSystemHermitean (Mijkl,gamma,Uijkl,sort=1,maxiter);myEigenSystemHermitean (extMijkl,extgamma,extUijkl,sort=1,maxiter);
      // conjugate:note the eigensystemhgermitean returns eigenvectors as column vectors, but
	// the components need to be complex conjugated 

         // treat correctly case for neutron energy loss
	 if (nn[6]>=0)
	 {
	    for(int ii=Uijkl.Rlo(); ii<Uijkl.Rhi(); ii++)       for(int jj=Uijkl.Clo(); jj<Uijkl.Chi(); jj++)       Uijkl[ii][jj]=conj(Uijkl[ii][jj]);
	    for(int ii=extUijkl.Rlo(); ii<extUijkl.Rhi(); ii++) for(int jj=extUijkl.Clo(); jj<extUijkl.Chi(); jj++) extUijkl[ii][jj]=conj(extUijkl[ii][jj]);
	 }
        // here we should also go for complex conjugate for the vector
         complex <double> extgammas;
     if (gamma(ini.nofcomponents)>=0&&fabs(gamma(ini.nofcomponents-1))<SMALL) 
                           // mind in manual the 1st dimension alpha=1 corresponds
			   // to the nth dimension here, because myEigensystmHermitean
			   // sorts the eigenvalues according to ascending order !!!
                           {if (nn[6]>SMALL)
			    {md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=sqrt(gamma(ini.nofcomponents));// gamma(ini.nofcomponents)=sqr(gamma^s)
                            extgammas=sqrt(extgamma(ini.extended_eigenvector_dimension));}
			    else if (nn[6]<-SMALL)
                            {md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=imaginary*sqrt(gamma(ini.nofcomponents));// gamma(ini.nofcomponents)=sqr(gamma^s)
                            extgammas=imaginary*sqrt(extgamma(ini.extended_eigenvector_dimension));}
 			    else
			    { //quasielastic line needs gamma=SMALL .... because Mijkl and therefore gamma have been set to 
			      // wn/kT instead of wn-wn'=SMALL*wn/kT (in jjjpar.cpp -mdcalc routines)
			      //set fix delta but keep sign
			          if (nn[6]>0){md.delta(i,j,k)(j1)=SMALL;
  			     md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=sqrt(SMALL*gamma(ini.nofcomponents));
                                              extgammas=sqrt(extgamma(ini.extended_eigenvector_dimension));}
				  else        {md.delta(i,j,k)(j1)=-SMALL;
                             md.sqrt_gamma(i,j,k)(ini.nofcomponents*(j1-1)+ini.nofcomponents,ini.nofcomponents*(j1-1)+ini.nofcomponents)=imaginary*sqrt(SMALL*gamma(ini.nofcomponents));
			                      extgammas=imaginary*sqrt(extgamma(ini.extended_eigenvector_dimension));}
			    }
			   }else 
                           {fprintf(stderr,"ERROR eigenvalue of single ion matrix <0: ev1=%g ev2=%g ev3=%g ... evn=%g\n",gamma(1),gamma(2),gamma(3),gamma(ini.nofcomponents));
                            exit(EXIT_FAILURE);}
//printf("extgamma %g %+g i\n",real(extgammas),imag(extgammas));
        for(k1=1;k1<=ini.extended_eigenvector_dimension;++k1){Ec(index_s(i,j,k,l,jmin,md,ini),k1)=extgammas*extUijkl(k1,ini.extended_eigenvector_dimension);}

        for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){
        md.U(i,j,k)(ini.nofcomponents*(j1-1)+m,ini.nofcomponents*(j1-1)+n)=Uijkl(m,n);
        md.M(i,j,k)(ini.nofcomponents*(j1-1)+m,ini.nofcomponents*(j1-1)+n)=Mijkl(m,n);
        }}    
//if (do_verbose==1){
//                  fprintf(stdout,"#Matrix M(s=%i %i %i)\n",i,j,k);
//                  myPrintComplexMatrix(stdout,Mijkl); 
//                  fprintf(stdout,"#Eigenvalues:\n");
//                  myPrintVector(stdout,gamma); 
//                  fprintf(stdout,"#Matrix U(s=%i%i%i)\n",i,j,k);
//                  myPrintComplexMatrix(stdout,Uijkl); 
//                 }

    }}}
    fclose(fin);

  }}}}

//initialize output files
  errno = 0;
if (do_jqfile==0)
{ printf("#saving mcdisp.qom and mcdisp.qei and mcdisp.qev\n");
  fout = fopen_errchk ("./results/mcdisp.qom",filemode);
   fprintf(fout, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcdisp.qom-->\n");
   fprintf(fout,"#*********************************************************************\n");
   fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(fout,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(fout,"#*********************************************************************\n");
          fprintf (fout, "#dispersion \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energies[meV] > intensities [barn/sr/f.u.]   f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
  foutqei = fopen_errchk ("./results/mcdisp.qei",filemode);
   fprintf(foutqei, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutqei);
   fprintf(foutqei,"#!<--mcphas.mcdisp.qei-->\n");
   fprintf(foutqei,"#*********************************************************************\n");
   fprintf(foutqei,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(foutqei,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(foutqei,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(foutqei,"#*********************************************************************\n");
          fprintf (foutqei, "#dispersion displayytext=E(meV)\n#displaylines=false \n#Ha[T] Hb[T] Hc[T] T[K] h k l Q[A^-1] energy[meV] int_dipapprFF) [barn/sr/f.u.] int_beyonddipappr [barn/sr/f.u.]  f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
  foutqev = fopen_errchk ("./results/mcdisp.qev",filemode);
   fprintf(foutqev, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutqev);
   fprintf(foutqev,"#!<--mcphas.mcdisp.qev-->\n");
   fprintf(foutqev,"#*********************************************************************\n");
   fprintf(foutqev,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(foutqev,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(foutqev,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(foutqev,"#*********************************************************************\n");
          fprintf (foutqev, "#!spins_wave_amplitude=1.0\n");
          fprintf (foutqev, "#!spins_show_ellipses=1.0\n");
          fprintf (foutqev, "#!spins_show_static_moment_direction=1.0\n");
          fprintf (foutqev, "#!dispersion displayytext=E(meV)\n#Ha[T] Hb[T] Hc[T] T[K] h k l Q[A^-1] energy[meV] int_dipapprFF) [barn/sr/f.u.] int_beyonddipappr [barn/sr/f.u.]  f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");

  foutqee = fopen_errchk ("./results/mcdisp.qee",filemode);
   fprintf(foutqee, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutqee);
   fprintf(foutqee,"#!<--mcphas.mcdisp.qee-->\n");
   fprintf(foutqee,"#*********************************************************************\n");
   fprintf(foutqee,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(foutqee,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(foutqee,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(foutqee,"#*********************************************************************\n");
          fprintf (foutqee, "#!spins_wave_amplitude=1.0\n");
          fprintf (foutqee, "#!spins_show_ellipses=1.0\n");
          fprintf (foutqev, "#!spins_show_static_moment_direction=1.0\n");
          fprintf (foutqee, "#!extended_eigenvector_dimension=%i\n",ini.extended_eigenvector_dimension); 
          fprintf (foutqee, "#!dispersion displayytext=E(meV)\n#Ha[T] Hb[T] Hc[T] T[K] h k l Q[A^-1] energy[meV] int_dipapprFF) [barn/sr/f.u.] int_beyonddipappr [barn/sr/f.u.]  f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");

  foutdstot = fopen_errchk ("./results/mcdisp.dsigma.tot",filemode);
          printf("#saving mcdisp.dsigma.tot\n");
   fprintf(foutdstot, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutdstot);
   fprintf(foutdstot,"#!<--mcphas.mcdisp.dsigma.tot-->\n");
   fprintf(foutdstot,"#*********************************************************************\n");
   fprintf(foutdstot,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(foutdstot,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(foutdstot,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(foutdstot,"#*********************************************************************\n");
          fprintf (foutdstot, "#!Total Scattering Cross Section in energy range [emin=%g ; emax=%g]\n#Ha[T] Hb[T] Hc[T] T[K] h k l  dsigma/dOmeg [barn/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}",ini.emin,ini.emax);

   if (do_Erefine==1){
          errno = 0;
  foutds = fopen_errchk ("./results/mcdisp.dsigma",filemode);
          printf("#saving mcdisp.dsigma\n");
   fprintf(foutds, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutds);
   fprintf(foutds,"#!<--mcphas.mcdisp.dsigma-->\n");
   fprintf(foutds,"#*********************************************************************\n");
   fprintf(foutds,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(foutds,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(foutds,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(foutds,"#*********************************************************************\n");
          fprintf (foutds, "#Scattering Cross Section \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energy[meV] dsigma/dOmegadE' [barn/mev/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
          fprintf (foutdstot, "for fast algorithm  vs summing dsigma for diff energies");
                     }  
          fprintf (foutdstot, "\n");
}
 
// initialize file with jq matrix
if (do_jqfile==1)
{  printf("#saving mcdisp.jq\n");
 jqfile = fopen_errchk ("./results/mcdisp.jq",filemode);
          printf("#saving mcdisp.jq\n");
   fprintf(jqfile, "#{output file of program %s",MCDISPVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),jqfile);
   fprintf(jqfile,"#!<--mcphas.mcdisp.dsigma.jq-->\n");
   fprintf(jqfile,"#*********************************************************************\n");
   fprintf(jqfile,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
   fprintf(jqfile,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(jqfile,"#            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(jqfile,"#*********************************************************************\n");
  fprintf (jqfile, "#Fourier Transform of 2 Ion Interaction - sta is calculated by comparing the larges eigenvalue\n# to that of the first q vector of the calculation");
   fputs (asctime(loctime),jqfile);
  if (do_verbose==1){   fprintf (jqfile, "#q=(hkl)\n #spin s() - spin s'()\n #3x3 matrix jss'(q) real im .... [meV]\n");}
  else {fprintf(jqfile,"#h  vs  k  vs  l  vs largest eigenvalue of J(hkl) matrix vs components of corresponding eigenvector re im re im re im re im\n");}
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
               
fprintf(stdout,"#q=(%g,%g,%g)\n",hkl(1),hkl(2),hkl(3));
 if(do_verbose==1){fprintf(stdout,"#Setting up J(q) matrix .... \n");}
 // calculate J(q)
 jq J(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),md);
 jq Jl(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),md);

 Vector d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
 Vector ij(1,3);
// int signa,signb,signc,sa,sb,sc,
 int sd;
 long int nofneighbours=0;
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
   { //sum up l.th neighbour interaction of crystallographic atom ll
     // 1. transform dn(l) to primitive lattice and round it to integer value
    sd=(*inputpars.jjj[ll]).sublattice[l];
    xyz=(*inputpars.jjj[ll]).xyz+(*inputpars.jjj[ll]).dn[l]-(*inputpars.jjj[sd]).xyz; // line added 17.6.09 to remove rounding bug in PCSMO calculation
    d=inputpars.rez*(const Vector&)xyz;
    for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); // rint d for loop below to determine crystallographic unit ss ...

    xyz=(*inputpars.jjj[ll]).dn[l];
    d=inputpars.rez*(const Vector&)xyz;// set d to distance for later use to determie phase factor in J(Q) ...

//	  if (do_verbose==1) {printf("#adding neighbor %i (%6.3f %6.3f %6.3f) of atom %i (%6.3f %6.3f %6.3f)- it contributes to J(s,s'):\n",l,xyz(1),xyz(2),xyz(3),ll,(*inputpars.jjj[ll]).xyz[1],(*inputpars.jjj[ll]).xyz[2],(*inputpars.jjj[ll]).xyz[3]);
//                              } 
   //2. in order to sum up we must take into account that the magnetic unit cell is
   //   larger than the crystallographic one - the ll-l neighbor interaction contributes
   //   to many different components of Js,ss(q) ... note s,ss runs over all the atoms
   //   in the magnetic supercell
         for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
         s=J.in(i1,j1,k1); 

         //next 2 lines special for treatment of assymetry in Rcu2 
	 // i.e. if j1 is 2n then  say: the neighbour is at -d, this ensures that
	 // for half of the atoms the interaction is counted at negative distances
	 // [this makes no error if the lattice is primitive!!!] but for RCu2 it makes
	 // a difference, because atoms at +-0.5c are not equal ... deleted 17.6.09 because it might give troubles in general !! 
	 sd=1;//if (inputpars.r[2][1]==0.5&&(double)j1/2.0==integer(1.0*j1/2)) sd=-1;
  	// if (inputpars.r[2][1]==0.0&&(double)(i1+j1)/2.0==integer(1.0*(i1+j1)/2)) sd=1;

         //calc ss (check in which crystallographic unit ss of the magnetic cell the neighbour l-ll lies)	 
         i=(int)(i1+sd*d_rint(1)-1); // calculate 
	 j=(int)(j1+sd*d_rint(2)-1);
	 k=(int)(k1+sd*d_rint(3)-1);
	 if (i>=0) ij(1)=integer(1.0*i/ini.mf.na())*ini.mf.na();
	 else      ij(1)=(integer(1.0*(i+1)/ini.mf.na())-1)*ini.mf.na();
	 if (j>=0) ij(2)=integer(1.0*j/ini.mf.nb())*ini.mf.nb();
	 else      ij(2)=(integer(1.0*(j+1)/ini.mf.nb())-1)*ini.mf.nb();
	 if (k>=0) ij(3)=integer(1.0*k/ini.mf.nc())*ini.mf.nc();
	 else      ij(3)=(integer(1.0*(k+1)/ini.mf.nc())-1)*ini.mf.nc();
//	 if (do_verbose==1) {printf("#ijk=%i %i %i  ij()=%6.3f %6.3f %6.3f ",i,j,k,ij(1),ij(2),ij(3));}
         i=i-(int)ij(1)+1;
	 j=j-(int)ij(2)+1;
	 k=k-(int)ij(3)+1;
	 ss=J.in(i,j,k);
//          if (do_verbose==1) {printf("#s=%i %i %i  s'=%i %i %i\n",i,j,k,i1,j1,k1);}
          // sum up 
           ComplexMatrix jsss(1,ini.nofcomponents*md.baseindex_max(i1,j1,k1),1,ini.nofcomponents*md.baseindex_max(i,j,k));
           jsss=0;
        
	  sl=(*inputpars.jjj[ll]).sublattice[l]; // the whole loop has also to be done 
                                                 // for all the other transitions of sublattice sl

            
          // therefore calculate offset of the set of transitions
          for(tl=1;tl<=md.noft(i1,j1,k1,ll);++tl){
	  for(tll=1;tll<=md.noft(i,j,k,sl);++tll){
	  
       
	     for(m=1;m<=ini.nofcomponents;++m){for(n=1;n<=ini.nofcomponents;++n){ //this should also be ok for nofcomponents > 3 !!! (components 1-3 denote the magnetic moment)
             jsss(ini.nofcomponents*(md.baseindex(i1,j1,k1,ll,tl)-1)+m,ini.nofcomponents*(md.baseindex(i,j,k,sl,tll)-1)+n)=(*inputpars.jjj[ll]).jij[l](m,n);
                                              }                                 } // but orbitons should be treated correctly by extending 3 to n !!
	                                         }} 
// increase Js,ss(q) taking into account the phase factors for the distance l-ll
          // J.mati(s,ss)+=jsss*exp(ipi*(double)sd*(q*d)); // changed
          jsss*=exp(ipi*(double)sd*(q*d)); J.mati(s,ss)+=jsss;
          ++nofneighbours; // count neighbours summed up
	 }}}
      
    }}


if (do_jqfile==1){if (do_verbose==1){fprintf (jqfile, "#q=(%g, %g, %g) ",hkl(1),hkl(2),hkl(3));
                                     fprintf(jqfile,"nofneighbours= %li\n",nofneighbours);
                                    }
                  else
                  {fprintf (jqfile, "%g  %g  %g ",hkl(1),hkl(2),hkl(3));}
                  }


if(do_verbose==1){fprintf(stdout,"#Transform J(q) matrix  with U...\n");}

// transform J(s,ss) (with md.U) and multiply 
// (with eigenvalues sqrt(gamma) [here md.sqrt_gamma]
// of matrix Malphabeta ) matrix J to L [here Jl.mati]  ... compare manual
 for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
 s=ini.mf.in(i1,j1,k1);
  for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
  ss=ini.mf.in(i2,j2,k2);
  Jl.mati(s,ss)=md.sqrt_gamma(i1,j1,k1)*md.U(i1,j1,k1).Conjugate().Transpose()*J.mati(s,ss)*md.U(i2,j2,k2)*md.sqrt_gamma(i2,j2,k2).Conjugate();
//if (do_verbose==1){
//                  fprintf(stdout,"#J(s=%i%i%i,s''=%i%i%i)=\n",i1,j1,k1,i2,j2,k2);
//                  myPrintComplexMatrix(stdout,J.mati(s,ss)); 
//                  fprintf(stdout,"#sqr(gamma_s=%i%i%i)=\n",i1,j1,k1);
//                  myPrintComplexMatrix(stdout,md.sqrt_gamma(i1,j1,k1));
//                  fprintf(stdout,"#U(s=%i%i%i)=\n",i1,j1,k1);
//                  myPrintComplexMatrix(stdout,md.U(i1,j1,k1));
//                  fprintf(stdout,"#sqr(gamma_s=%i%i%i)=\n",i2,j2,k2);
//                  myPrintComplexMatrix(stdout,md.sqrt_gamma(i2,j2,k2));
//                  fprintf(stdout,"#U(s=%i%i%i)=\n",i2,j2,k2);
//                  myPrintComplexMatrix(stdout,md.U(i2,j2,k2));
//                  fprintf(stdout,"#sqr(gamma_s) U(s)T* J(s=%i%i%i,s''=%i%i%i) U(s'') sqr(gamma_s'')*=\n",i1,j1,k1,i2,j2,k2);
//                  myPrintComplexMatrix(stdout,Jl.mati(s,ss)); 
//                 }
  }}}
 }}}


// calculate Ac
if(do_verbose==1){fprintf(stdout,"#calculating matrix A\n");}
// Ac  is the matrix A which is given in manual chapter 9.2.1 (eq (30) ff) 
// -- diagonalization gives omega_r and Tau
   ComplexMatrix Ac(1,dimA,1,dimA);
   ComplexMatrix Lambda(1,dimA,1,dimA);
   ComplexMatrix Tau(1,dimA,1,dimA);
   ComplexMatrix J_Q(1,ini.nofcomponents*dimA,1,ini.nofcomponents*dimA);
   Ac=0;J_Q=0;Lambda=0;
   for(i1=1;i1<=ini.mf.na();++i1){for(j1=1;j1<=ini.mf.nb();++j1){for(k1=1;k1<=ini.mf.nc();++k1){
//   s=(((i1-1)*ini.mf.nb()+(j1-1))*ini.mf.nc()+(k1-1))*inputpars.nofatoms;
//    for(l1=1;l1<=inputpars.nofatoms;++l1){
//      Ac(s+l1,s+l1)=md.delta(i1,j1,k1)(l1);}

    for(l1=1;l1<=inputpars.nofatoms;++l1){
     for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);
      if(md.delta(i1,j1,k1)(b)<0){Lambda(s,s)=-1;}else{Lambda(s,s)=+1;}
      Ac(s,s)=md.delta(i1,j1,k1)(b)*Lambda(s,s);
//      if(do_verbose==1){fprintf(stdout,"#i=%i j=%i k=%i atomnr=%i trans=%i ... s=%i ",i1,j1,k1,l1,t1,s);
//                        fprintf(stdout,"#lambda(%i,%i)xdelta(%i)=%g + i %g\n",s,s,s,real(Ac(s,s)),imag(Ac(s,s)));
//                       }
      }}

   for(i2=1;i2<=ini.mf.na();++i2){for(j2=1;j2<=ini.mf.nb();++j2){for(k2=1;k2<=ini.mf.nc();++k2){
//   ss=(((i2-1)*ini.mf.nb()+(j2-1))*ini.mf.nc()+(k2-1))*inputpars.nofatoms;

    for(l1=1;l1<=inputpars.nofatoms;++l1){ 
     for(t1=1;t1<=md.noft(i1,j1,k1,l1);++t1){
    for(l2=1;l2<=inputpars.nofatoms;++l2){ 
     for(t2=1;t2<=md.noft(i2,j2,k2,l2);++t2){
      s=index_s(i1,j1,k1,l1,t1,md,ini);
      ss=index_s(i2,j2,k2,l2,t2,md,ini);
      b=md.baseindex(i1,j1,k1,l1,t1);
      bb=md.baseindex(i2,j2,k2,l2,t2);
     Ac(s,ss)-=Jl.mati(Jl.in(i1,j1,k1),Jl.in(i2,j2,k2))(ini.nofcomponents*(b-1)+ini.nofcomponents,ini.nofcomponents*(bb-1)+ini.nofcomponents);
                                             //nofcomponents^th dimension corresponds to 1st in manual 
					     // and it is only necessary to take into 
					     // acount this dimension!!
       for(i=1;i<=ini.nofcomponents;++i){for(j=1;j<=ini.nofcomponents;++j){
        J_Q(ini.nofcomponents*(s-1)+i,ini.nofcomponents*(ss-1)+j)+=J.mati(J.in(i1,j1,k1),J.in(i2,j2,k2))(ini.nofcomponents*(b-1)+i,ini.nofcomponents*(bb-1)+j);
       }}

     }}}}

  }}}
 }}}

// printout Fouriertransform of matrix jq
if (do_jqfile==1){
       if (do_verbose==1)
       {//fprintf (jqfile, "#spin (%i*r1 %i*r2 %i*r3) - spin (%i*r1 %i*r2 %i*r3)\n",i1,j1,k1,i2,j2,k2);
         myPrintComplexMatrix(jqfile,J_Q); 
       }

	// diagonalize JQ to get eigenvalues (biggest corresponds to Tn) !!!
         Vector Tn(1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms);
         ComplexMatrix eigenvectors(1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms,1,ini.nofcomponents*ini.mf.n()*inputpars.nofatoms);
         myEigenSystemHermitean (J_Q,Tn,eigenvectors,sort=1,maxiter);

       if(do_verbose==1)
       {fprintf(jqfile,"#eigenvalues(highest corresponds to Tn, predicted magstructure)\n");
         myPrintVector(jqfile,Tn); 
        fprintf(jqfile,"#eigenvectors(moment direction):\n");
         myPrintComplexMatrix(jqfile,eigenvectors); 
       }
       else
       {i2=ini.nofcomponents*ini.mf.n()*inputpars.nofatoms;
        fprintf(jqfile," %g ",Tn(i2));
        for (i1=1;i1<=i2;++i1)
        {fprintf(jqfile," %6.3g ",real(eigenvectors(i1,i2)));
         fprintf(jqfile," %6.3g ",imag(eigenvectors(i1,i2)));
        }
        fprintf(jqfile,"\n");
       }
       if (jqsta<-0.5){jq0=Tn(i2);jqsta=0;}
       else           {if(Tn(i2)>jq0)
                         {jqsta+=(Tn(i2)-jq0)*(Tn(i2)-jq0);}
                      }
 }
 else
 {// no jqfile but excitations to be calculated
 if(do_verbose==1){fprintf(stdout,"#diagonalizing %ix%i matrix A...\n",dimA,dimA);
              //             myPrintComplexMatrix(stdout,Ac); 
                   }
   // diagonalize Ac to get energies  and eigenvectors !!!
   Vector En(1,dimA);
   Vector ints(1,dimA);
   Vector intsbey(1,dimA);
//   myEigenValuesHermitean (Ac,En,sort,maxiter);
  
//   myEigenSystemHermitean (Ac,En,Tau,sort,maxiter);
//    myPrintVector(stdout,En);
   myEigenSystemHermiteanGeneral (Lambda,Ac,En,Tau,sort,maxiter);
   En=1.0/En;
   sortE(En,Tau);
           Tau=Tau.Conjugate();
  	// conjugate inserted 31.10.05, because when calculating simple AF - I noticed
	// that the eigensystemhgermitean returns eigenvectors as column vectors, but
	// the components need to be complex conjugated 

 if(do_verbose==1){// fprintf(stdout,"#eigenvectors (matrix Tau):\n");
                   // myPrintComplexMatrix(stdout,Tau); 
                    fprintf(stdout,"#saving the following eigenvalues (meV) to mcdisp.qom:\n");}
   int dim=3;
   if (ini.hkllist==1){dim=(int)ini.hkls[counter][0]-3;}
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));

   for (i=1;i<=dimA;++i){
	       fprintf (fout, " %4.4g ",En(i));
               if(do_verbose==1){fprintf(stdout, " %4.4g",En(i));}
                         }


   // calculate and printout intensities [the energies have already
   // been printed out above, so any refinement of energies during intcalc
   // is not included in the output file]
#ifndef _THREADS  
  double QQ; mfcf ev_real(ini.mf),ev_imag(ini.mf);
             mfcf eev_real(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ini.extended_eigenvector_dimension);
             mfcf eev_imag(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ini.extended_eigenvector_dimension);
#endif
  double diffint=0,diffintbey=0;
  if(do_verbose==1){fprintf(stdout,"\n#calculating  intensities approximately ...\n");}
                  fprintf (fout, " > ");
diffint=0;diffintbey=0;
                  if(do_gobeyond)do_gobeyond=intcalc_beyond_ini(ini,inputpars,md,do_verbose,hkl);
                  Vector dd(1,dim),dd_int(1,dim);  dd+=100000.0;dd_int+=100000.0;
#ifndef _THREADS  
                     ComplexMatrix chi(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     ComplexMatrix chibey(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     ComplexMatrix S(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     ComplexMatrix Sbey(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     Matrix pol(1,md.nofcomponents,1,md.nofcomponents);
                     Matrix polICIC(1,md.nofcomponents,1,md.nofcomponents);
                     Matrix polICn(1,md.nofcomponents,1,md.nofcomponents);
                     Matrix polnIC(1,md.nofcomponents,1,md.nofcomponents);
#else
                  // Populates the thread data structure
                  intcalcapr_input *tin[NUM_THREADS]; //thrdat.ini = &ini; thrdat.inputpars = &inputpars; thrdat.md = &md;
                  thrdat.ini      = new inimcdis*[NUM_THREADS];      thrdat.inputpars= new par*[NUM_THREADS];
                  thrdat.ev_real  = new mfcf*[NUM_THREADS];          thrdat.ev_imag  = new mfcf*[NUM_THREADS];
                  thrdat.eev_real = new mfcf*[NUM_THREADS];          thrdat.eev_imag = new mfcf*[NUM_THREADS];
                  thrdat.chi      = new ComplexMatrix*[NUM_THREADS]; thrdat.chibey   = new ComplexMatrix*[NUM_THREADS];
                  thrdat.S        = new ComplexMatrix*[NUM_THREADS]; thrdat.Sbey     = new ComplexMatrix*[NUM_THREADS];
                  thrdat.pol      = new Matrix*[NUM_THREADS];        thrdat.polICIC  = new Matrix*[NUM_THREADS];
                  thrdat.polICn   = new Matrix*[NUM_THREADS];        thrdat.polnIC   = new Matrix*[NUM_THREADS];
                  thrdat.Ec       = new ComplexMatrix*[NUM_THREADS]; thrdat.md       = new mdcf*[NUM_THREADS];
                  thrdat.Tau      = new ComplexMatrix*[NUM_THREADS]; thrdat.hkl = hkl; thrdat.thread_id = -1;
                  for (int ithread=0; ithread<NUM_THREADS; ithread++) 
                  {
                     tin[ithread] = new intcalcapr_input(dimA,ithread,1,do_verbose,En,&inputpars,&ini);
                     thrdat.chi[ithread] = new ComplexMatrix(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     thrdat.chibey[ithread] = new ComplexMatrix(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     thrdat.S[ithread] = new ComplexMatrix(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     thrdat.Sbey[ithread] = new ComplexMatrix(1,md.nofcomponents*dimA,1,md.nofcomponents*dimA);
                     thrdat.pol[ithread] = new Matrix(1,md.nofcomponents,1,md.nofcomponents);
                     thrdat.polICIC[ithread] = new Matrix(1,md.nofcomponents,1,md.nofcomponents);
                     thrdat.polICn[ithread] = new Matrix(1,md.nofcomponents,1,md.nofcomponents);
                     thrdat.polnIC[ithread] = new Matrix(1,md.nofcomponents,1,md.nofcomponents);
                     thrdat.ev_real[ithread] = new mfcf(ini.mf); thrdat.ev_imag[ithread] = new mfcf(ini.mf);
                     thrdat.eev_real[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ini.extended_eigenvector_dimension);
                     thrdat.eev_imag[ithread] = new mfcf(ini.mf.na(),ini.mf.nb(),ini.mf.nc(),ini.mf.nofatoms,ini.extended_eigenvector_dimension);
                     thrdat.Ec[ithread] = new ComplexMatrix(1,dimA,1,ini.extended_eigenvector_dimension); *thrdat.Ec[ithread]=Ec;
                     thrdat.Tau[ithread] = new ComplexMatrix(1,dimA,1,dimA); *thrdat.Tau[ithread]=Tau;
                     thrdat.md[ithread] = new mdcf(md); thrdat.ini[ithread] = new inimcdis(ini); thrdat.inputpars[ithread] = new par(inputpars);
                  }
                  // Initialises mutual exclusions and threads
                  MUTEX_INIT(mutex_loop);
                  MUTEX_INIT(mutex_index);
                  EVENT_INIT(checkfinish);
                  #ifdef __linux__
                  pthread_t threads[NUM_THREADS]; int rc; void *status;
                  pthread_attr_t attr;
                  pthread_attr_init(&attr);
                  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                  #else
                  HANDLE threads[NUM_THREADS];
                  DWORD tid[NUM_THREADS], dwError;
                  #endif
                  int ithread=0,oldi=-1,num_threads_started=-1; double QQ; Vector vQQ(1,dimA);
#endif
#ifdef _THREADS  
                  for (i=1;i<=dimA;i+=NUM_THREADS)
#else
                  for (i=1;i<=dimA;++i)
#endif
                  {
#ifdef _THREADS  
                     oldi=i;
                     // Runs threads until all are running - but wait until they are completed in order before printing output.
                     //  This is to ensure that the output is exactly the same as for a single thread (otherwise it would be out of order).
                     for(ithread=0; ithread<NUM_THREADS; ithread++)
                     {
                        i+=ithread; if(i>dimA) break;
                        if(do_gobeyond==0) intsbey(i)=-1.1; else intsbey(i)=+1.1;
                        tin[ithread]->En=En(i); tin[ithread]->intensitybey=intsbey(i); tin[ithread]->level = i;
                        #ifdef __linux__
                        rc = pthread_create(&threads[ithread], &attr, intcalc_approx, (void *) tin[ithread]);
                        if(rc) { printf("Error return code %i from thread %i\n",rc,ithread+1); exit(EXIT_FAILURE); }
                        #else
                        threads[ithread] = CreateThread(NULL, 0, intcalc_approx, (void *) tin[ithread], 0, &tid[ithread]);
                        if(threads[ithread]==NULL) { dwError=GetLastError(); printf("Error code %i from thread %i\n",dwError,ithread+1); exit(EXIT_FAILURE); }
                        #endif
                        num_threads_started = ithread+1;
                     }
                     for(int th=0; th<num_threads_started; th++)
                        #ifdef __linux__
                        rc = pthread_join(threads[th], &status);
                        #else
                        WaitForSingleObject(threads[th],INFINITE);
                        #endif
                     #define ev_real (*thrdat.ev_real[ithread])
                     #define ev_imag (*thrdat.ev_imag[ithread])
                     #define eev_real (*thrdat.eev_real[ithread])
                     #define eev_imag (*thrdat.eev_imag[ithread])
                     for(ithread=0; ithread<num_threads_started; ithread++)
                     {
                        ints(tin[ithread]->level) = tin[ithread]->intensity; vQQ(tin[ithread]->level) = tin[ithread]->QQ;
                        QQ = vQQ(tin[ithread]->level); intsbey(tin[ithread]->level) = tin[ithread]->intensitybey;
                        i = tin[ithread]->level;
#else
                     if(do_gobeyond==0){intsbey(i)=-1.1;}else{intsbey(i)=+1.1;}
                     ints(i)=intcalc_approx(chi,chibey,S,Sbey,pol,polICIC,polICn,polnIC,
                                            intsbey(i),ev_real,ev_imag,eev_real,eev_imag,Ec,dimA,Tau,i,En(i),ini,inputpars,hkl,md,do_verbose,QQ);
#endif
                     if (ini.hkllist==1)
	             {double test; // add to sta distance to nearest measured peak squared
	              for (j1=1;j1<=ini.hkls[counter][0]-3;++j1)
	              {if ((test=fabs(En(i)-ini.hkls[counter][j1+3]))<dd(j1))dd(j1)=test;
                       if ((test=fabs(En(i)-ini.hkls[counter][j1+3]))<dd_int(j1)&&ints(i)>1e-4)dd_int(j1)=test;
                      }
	             }

                     if(intsbey(i)<0)intsbey(i)=-1;
                     //printout rectangular function to .mdcisp.qom
	             fprintf (fout, " %4.4g %4.4g",ints(i),intsbey(i));
                     fprintf (foutqei, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g  %4.4g  %4.4g\n",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3),QQ,En(i),ints(i),intsbey(i));

                     fprintf (foutqev, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g  %4.4g  %4.4g\n",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3),QQ,En(i),ints(i),intsbey(i));
                     fprintf (foutqev, "#eigenvector real part\n");
                     ev_real.print(foutqev); // here we printout the eigenvector of the excitation
                     fprintf (foutqev, "#eigenvector imaginary part\n");
                     ev_imag.print(foutqev); // 
                     fprintf (foutqev, "#\n");

                     fprintf (foutqee, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g  %4.4g  %4.4g\n",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3),QQ,En(i),ints(i),intsbey(i));
                     fprintf (foutqee, "#eigenvector real part\n");
                     eev_real.print(foutqee); // here we printout the eigenvector of the excitation
                     fprintf (foutqee, "#eigenvector imaginary part\n");
                     eev_imag.print(foutqee); // 
                     fprintf (foutqee, "#\n");
                 if(do_verbose==1){fprintf(stdout, "#IdipFF= %4.4g Ibeyonddip=%4.4g\n",ints(i),intsbey(i));}
                     if(En(i)>=ini.emin&&En(i)<=ini.emax){diffint+=ints(i);diffintbey+=intsbey(i);}
#ifdef _THREADS
                     }
                   i=oldi;
#endif
		   }
#ifdef _THREADS
                  #undef ev_real
                  #undef ev_imag
                  #undef eev_real
                  #undef eev_imag
                  for (ithread=0; ithread<NUM_THREADS; ithread++) 
                  {
                     delete thrdat.chi[ithread]; delete thrdat.chibey[ithread]; delete thrdat.S[ithread]; delete thrdat.Sbey[ithread];
                     delete thrdat.pol[ithread]; delete thrdat.polICIC[ithread]; delete thrdat.polICn[ithread]; delete thrdat.polnIC[ithread];
                     delete thrdat.ev_real[ithread]; delete thrdat.ev_imag[ithread]; delete thrdat.eev_real[ithread]; delete thrdat.eev_imag[ithread];
                     delete thrdat.Ec[ithread]; delete thrdat.Tau[ithread]; delete tin[ithread]; 
                     delete thrdat.md[ithread]; delete thrdat.ini[ithread]; delete thrdat.inputpars[ithread];
                  }
                  delete[] thrdat.Ec; delete[] thrdat.Tau; delete[] thrdat.md; delete[] thrdat.ini; delete[] thrdat.inputpars; //delete tin;
                  delete[] thrdat.chi; delete[] thrdat.chibey; delete[] thrdat.S; delete[] thrdat.Sbey;
                  delete[] thrdat.pol; delete[] thrdat.polICIC; delete[] thrdat.polICn; delete[] thrdat.polnIC;
                  delete[] thrdat.ev_real; delete[] thrdat.ev_imag; delete[] thrdat.eev_real; delete[] thrdat.eev_imag; 
#endif
    fprintf (foutdstot, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %4.4g",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3),diffint,diffintbey);
    sta+=dd*dd;sta_int+=dd_int*dd_int;

              //initialize output file for display
            errno = 0;
            fout1 = fopen_errchk ("./results/.mcdisp.qom","w");
            fprintf (fout1, "#{%s ",MCDISPVERSION);
            curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout1);
            fprintf (fout1, "#displayytext=I(barns/meV/sr/f.u.)\n");
            fprintf (fout1, "#displayxtext=E(meV)\n");
            fprintf (fout1, "#displaytitle=(%4.4f %4.4f %4.4f) blue: DMD_Dipapprox red: DMD_exact green: Minv_Dipapprox\n",hkl(1),hkl(2),hkl(3));
            fprintf (fout1,"#Ha[T] Hb[T] Hc[T] T[K] h k l  energies[meV] intensities(dip approx for FF) [barn/meV/sr/f.u.] f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
		     if (do_Erefine==0) epsilon=(Max(En)-Min(En))/100;
		     if (epsilon<=0) epsilon=0.1;
                  for (i=1;i<=dimA;++i)
		    { 
		     if (ints(i)>SMALLINT)  // draw triangles to show calculated intensity
		      {for (E=0;E<=ints(i)/epsilon;E+=ints(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	                fprintf (fout1, " %4.4g %4.4g %4.4g 0\n",En(i)-epsilon+E*epsilon*epsilon/ints(i),E,En(i));
		       }
		       for (E=ints(i)/epsilon;E>=0;E-=ints(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	                fprintf (fout1, " %4.4g %4.4g %4.4g 0\n",En(i)+epsilon-E*epsilon*epsilon/ints(i),E,En(i));
		       }
                       fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	               fprintf (fout1, " %4.4g 0 %4.4g 0\n",En(i)+epsilon,En(i));
		      }
		     if (intsbey(i)>SMALLINT)  // draw triangles to show calculated intensity
		      {for (E=0;E<=intsbey(i)/epsilon;E+=intsbey(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	                fprintf (fout1, " %4.4g 0 %4.4g %4.4g \n",En(i),En(i)-epsilon+E*epsilon*epsilon/intsbey(i),E);
		       }
		       for (E=intsbey(i)/epsilon;E>=0;E-=intsbey(i)/2/epsilon/10)
                       {fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	                fprintf (fout1, " %4.4g 0 %4.4g %4.4g\n",En(i),En(i)+epsilon-E*epsilon*epsilon/intsbey(i),E);
		       }
                       fprintf (fout1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	               fprintf (fout1, " %4.4g 0 %4.4g 0 \n",En(i),En(i)+epsilon);
		      }
		    }
	  fclose(fout1);   
	    
		    
                if(do_verbose==1){fprintf(stdout, "\n");}
		    
   // do refinement of energies by output of scattering cross section vs enrgy transfer if required
  if (do_Erefine==1){double totint=0;
                if(do_verbose==1){fprintf(stdout, "#refining calculation with exact calculation of energy dependence of scattering cross section\n");}
          errno = 0;
          foutds1 = fopen_errchk ("./results/.mcdisp.dsigma","w");
          fprintf (foutds1, "#{%s ",MCDISPVERSION);
          curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),foutds1);
          fprintf (foutds1, "#Scattering Cross Section \n#Ha[T] Hb[T] Hc[T] T[K] h k l  energy[meV] dsigma/dOmegadE'[barn/mev/sr/f.u.] (dipolar approx for FF) f.u.=crystallogrpaphic unit cell (r1xr2xr3)}\n");
		     double intensity;
		     intensity=intcalc(dimA,ini.emin,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds1, " %4.4g %4.4g \n",ini.emin,intensity);
		     intensity=intcalc(dimA,ini.emax,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds1, " %4.4g %4.4g \n",ini.emax,intensity);
	  fclose(foutds1);
	  for(E=ini.emin;E<=ini.emax;E+=epsilon/2)
	   {
		     intensity=intcalc(dimA,E,ini,inputpars,J,q,hkl,md,do_verbose,epsilon);   
		     totint+=intensity*epsilon/2;

          foutds1 = fopen_errchk ("./results/.mcdisp.dsigma","a");
                     fprintf (foutds1, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds1, " %4.4g %4.4g \n",E,intensity);
          fclose(foutds1);	   
                     fprintf (foutds, " %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g ",ini.Ha,ini.Hb,ini.Hc,ini.T,hkl(1),hkl(2),hkl(3));
	             fprintf (foutds, " %4.4g %4.4g \n",E,intensity);
	   }
	             fprintf (foutdstot, " %4.4g ",totint);
                     }  

   fprintf (foutdstot, "\n");              
   fprintf (fout, "\n");
   }
   if (ini.hkllist==1){hkl(1)=(double)counter;}
}}}
    if (do_jqfile==1) 
     {fprintf(jqfile,"#it follows the standard deviation sta defined as:\n");
      fprintf(jqfile,"#the sum of squared differences between the highest eigenvalue\n");
      fprintf(jqfile,"#of a q vector and that of the first q-vector in the list in mcdisp.par.\n");
      fprintf(jqfile,"#only those eigenvalues are taken into account in the sum, which are larger\n");
      fprintf(jqfile,"#than that of the first q-vector in the list in mcdisp.par - this is usefule\n");
      fprintf(jqfile,"#for obtaining an exchange interaction with maximum at the first q-vector\n");
      fprintf(jqfile,"#in the list in mcdisp.par\n");
      fprintf(jqfile,"!sta=%g\n",jqsta);fclose(jqfile);}
    else
     {
    fprintf(fout,"#definitions: sta= sum_i [Eexp(i) - nearestEcalc(i)]^2\n");
    fprintf(fout,"#             sta_int= sum_i [Eexp(i) - nearestEcalc_with_Int>0.1mb/srf.u.(i)]^2\n");
    fprintf (fout, "#!sta= %8.6g \n",sta);
    fprintf (fout, "#!sta_int= %8.6g \n",sta_int);

    fprintf(foutqei,"#definitions: sta= sum_i [Eexp(i) - nearestEcalc(i)]^2\n");
    fprintf(foutqei,"#             sta_int= sum_i [Eexp(i) - nearestEcalc_with_Int>0.1mb/srf.u.(i)]^2\n");
    fprintf (foutqei, "#!sta= %8.6g \n",sta);
    fprintf (foutqei, "#!sta_int= %8.6g \n",sta_int);

    fprintf(stdout,"#definitions: sta= sum_i [Eexp(i) - nearestEcalc(i)]^2\n");
    fprintf(stdout,"#             sta_int= sum_i [Eexp(i) - nearestEcalc_with_Int>0.1mb/srf.u.(i)]^2\n");
    fprintf(stdout, "#!sta= %8.6g \n",sta);
    fprintf(stdout, "#!sta_int= %8.6g \n",sta_int);
   
    fclose(foutqei);
    fclose(foutqev);
    fclose(foutqee);
    fclose(fout);
                 if (do_Erefine==1){fclose(foutds);}
    fclose(foutdstot);
     }
}

//*************************************************************************************************
// main program
int main (int argc, char **argv)
{int i,do_Erefine=0,do_jqfile=0,do_verbose=0,maxlevels=10000000,do_createtrs=0,do_readtrs=0;
 const char * spinfile="mcdisp.mf"; //default spin-configuration-input file
 const char * filemode="w";
 double epsilon; //imaginary part of omega to avoid divergence
 double minE=-100000.0,maxE=+100000.0;
 fprintf(stderr,"***********************************************************************\n");
 fprintf(stderr,"*\n");
 fprintf(stderr,"* mcdisp - program to calculate the dispersion of magnetic excitations\n");
 fprintf(stderr,"*\n");
 fprintf(stderr,"* reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
 fprintf(stderr,"*            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
 fprintf(stderr,"***********************************************************************\n\n");


// check command line and initialize parameters ini
for (i=1;i<=argc-1;++i){
   if(strcmp(argv[i],"-r")==0) {do_Erefine=1; if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -r needs argument epsilon\n");exit(EXIT_FAILURE);}
		                                                epsilon=strtod(argv[i+1],NULL);++i;
							        fprintf(stdout,"#epsilon= %g\n",epsilon);
				     }		
         else {if(strcmp(argv[i],"-jq")==0) {do_jqfile=1;minE=SMALL;maxlevels=1;}       
          else {if(strcmp(argv[i],"-t")==0) do_readtrs=1;       
           else {if(strcmp(argv[i],"-c")==0) do_createtrs=1;       
            else {if(strcmp(argv[i],"-a")==0) filemode="a";       
             else {if(strcmp(argv[i],"-v")==0||strcmp(argv[i],"-verbose")==0) do_verbose=1;       
              else {if(strcmp(argv[i],"-max")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -max needs argument(s)\n");exit(EXIT_FAILURE);}
		                                  maxlevels=(int)strtod(argv[i+1],NULL);++i;
						  fprintf(stdout,"#maximum number of single ion excitations taken into account (starting with lowest energy): %i\n",maxlevels);
					         }       
               else {if(strcmp(argv[i],"-maxE")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -maxE needs argument(s)\n");exit(EXIT_FAILURE);}
		                                  maxE=strtod(argv[i+1],NULL);++i;
						  fprintf(stdout,"#maximum Energy of single ion excitations taken into account: %g\n",maxE);
  					         }       
                else {if(strcmp(argv[i],"-minE")==0) {if(i==argc-1){fprintf(stderr,"Error in command: mcdisp -minE needs argument(s)\n");exit(EXIT_FAILURE);}
 		                                  minE=strtod(argv[i+1],NULL);++i;
						  fprintf(stdout,"#minimum Energy of single ion excitations taken into account: %g\n",minE);
					         }       
           	 else{spinfile=argv[i];}
		   }
		  }          
		 }
		}
	      }
	     }
	    }
           }	
    }
inimcdis ini("mcdisp.par",spinfile);

if (argc > 10) {ini.errexit();}
  // as class load  parameters from file
  par inputpars("./mcphas.j");
  if(ini.nofcomponents!=inputpars.nofcomponents){fprintf(stderr,"Error mcdisp: number of components read from mcdisp.par (%i) and mcphas.j (%i) not equal\n",ini.nofcomponents,inputpars.nofcomponents);exit(1);}
  if(ini.nofatoms!=inputpars.nofatoms){fprintf(stderr,"Error mcdisp: number of atoms in crystal unit cell read from mcdisp.par (%i) and mcphas.j (%i) not equal\n",ini.nofatoms,inputpars.nofatoms);exit(1);}
  inputpars.save("./results/_mcdisp.j");
  inputpars.save_sipfs("./results/_"); 


//calculate dispersion and save to files
dispcalc(ini,inputpars,do_Erefine,do_jqfile,do_createtrs,do_readtrs,do_verbose,maxlevels,minE,maxE,epsilon,filemode);
  
 printf("RESULTS saved in directory ./results/  - files:\n");
   printf("  mcdisp.qei  - T,H,qvector vs energies and neutron intensities\n");
   printf("  mcdisp.qom  - T,H,qvector vs all mode energies in one line (and neutron intensities)\n");
   printf("  mcdisp.qev  - T,H,qvector,E vs eigenvectors\n");
   printf("  mcdisp.qee  - T,H,qvector,E vs extended eigenvectors (more components to plot chrgedens.)\n");
   printf("  mcdisp.dsigma.tot  - T,H,qvector vs total intensity (sum of all modes)\n");
   printf("  mcdisp.dsigma      - (option -r) T,H,qvector,E vs intensity obtained from dyn susz\n");
   printf("  mcdisp.trs  - single ion transitions used\n");
   printf("  _mcdisp.par - input parameters read from mcdisp.par\n");
   printf("  _mcdisp.mf  - input parameters read from mcdisp.mf\n");
   printf("  _mcdisp.j   - input parameters read from mcphas.j\n");
   printf("  ...         - and a copy of the single ion parameter files used.\n\n");
   fprintf(stderr,"************************************************************\n");
   fprintf(stderr,"                    End of Program mcdisp\n");
   fprintf(stderr," reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
   fprintf(stderr,"            M. Rotter J. Comp. Mat. Sci. 38 (2006) 400\n");
   fprintf(stderr,"************************************************************\n");

 exit(0);
 
}





