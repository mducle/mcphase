// *************************************************************************
// ************************ class jjjpar     *******************************
// *************************************************************************
// jjjpar is a class to store all parameters associated
// with a specific ion, for example CEF parameters
// and exchange parameters 
// it is used by many programs in the package
// moreover, it loads also the user defined single ion module functions (linux only)
#include "jjjpar.hpp"

#define MAXNOFNUMBERSINLINE 200
#define MAXNOFCHARINLINE 1024

#define MU_B 0.05788
#define K_B  0.0862
#define SMALL 1e-6   //!!! must match SMALL in mcdisp.c and ionpars.cpp !!!
                     // because it is used to decide wether for small transition
		     // energy the matrix Mijkl contains wn-wn' or wn/kT
#define PI 3.1415926535

#include "jjjpar_intmod_kramer.cpp"   // some functions for intern_mcalc=1
#include "jjjpar_intmod_brillouin.cpp"// some functions for intern_mcalc=3



/****************************************************************************/
// function to calculate magnetisation M from effective field H
// this is the heart of the meanfield algorithm an it is necessary to
// keep this routine as efficient as possible
// at the moment we do only groundstate doublet
/****************************************************************************/
Vector & jjjpar::mcalc (double & T, Vector &  gjmbH, double & lnZ,double & U)
{switch (intern_mcalc)
  {case 1: return kramer(T,gjmbH,lnZ,U);break;
   case 2: return (*iops).cfield(T,gjmbH,lnZ,U);break;
   case 3: return brillouin(T,gjmbH,lnZ,U);break;
   default: (*m)(&mom,&T,&gjmbH,&gJ,&ABC,&cffilename,&lnZ,&U);return mom;
  }
}

/****************************************************************************/
// this function returns n (the number of transitions in the single ion susceptibility)
// the transition matrix mat corresponding to jjjpar.transitionnumber and delta
// for effective field heff and temperature given on input
/****************************************************************************/
int jjjpar::dmcalc(double & T,Vector & gjmbheff,ComplexMatrix & mat,float & delta)
{switch (intern_mcalc)
  {case 1: return kramerdm(transitionnumber,T,gjmbheff,mat,delta);break;
   case 2: return (*iops).cfielddm(transitionnumber,T,gjmbheff,mat,delta);break;
   case 3: return brillouindm(transitionnumber,T,gjmbheff,mat,delta);break;
   default: return (*dm)(&transitionnumber,&T,&gjmbheff,&gJ,&ABC,&cffilename,&mat,&delta);
  }
}


/****************************************************************************/
// returns eigenvalues, boltzmann population and eigenstates matrix parameters of ion
/****************************************************************************/
ComplexMatrix & jjjpar::eigenstates (Vector & gjmbheff,double & T)
{switch (intern_mcalc)
  {case 0:  if(estates!=NULL){(*estates)(&est,&gjmbheff,&gJ,&T,&ABC,&cffilename);}
            return est;break;
   case 2:  est=(*iops).cfeigenstates(gjmbheff,T);return est;break;
   default: fprintf (stderr, "warning class jjjpar: single ion module does not allow to return eigenvalues and eigenstates \n"); 
            est=0;return est;
  }
}


/****************************************************************************/
// returns transition element matrix N(Q) in order to be able to go beyond 

// dipolar approximation in mcdisp - it requires a call to eigenstates first

/****************************************************************************/
int jjjpar::dncalc(Vector & Qvec,double & T, ComplexMatrix & nat,ComplexMatrix & ests)

{double J0,J2,J4,J6;
 double Q,d,s;
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
         double th,ph,Qx,Qy,Qz;							 
	 Qx=Qvec(1);Qy=Qvec(2);Qz=Qvec(3);
	 th=acos(Qz/Q);
	 if(sin(th)>=SMALL){
	                    if(Qx>0){ph=acos(Qx/(Q*sin(th))-SMALL);}
			    else    {ph=acos(Qx/(Q*sin(th))+SMALL);}
			   }
			 else{ph=0;}
	 if (Qy<0){ph=2*PI-ph;} 

 switch (intern_mcalc)
  {static int washere=0;
   
   case 0:if (ddnn!=NULL){return (*ddnn)(&transitionnumber,&th,&ph,&J0,&J2,&J4,&J6,&ests,&T,&nat);break;}
   case 2:return (*iops).cfielddn(transitionnumber,th,ph,J0,J2,J4,J6,Zc,ests,T,nat);break;
   default: if(washere==0){fprintf(stderr,"Warning in scattering operator function dncalc - for ion %s \ngoing beyond dipolar approximation is not implemented\n",cffilename);
                           washere=1;}
            return 0;
  }

}




/****************************************************************************/
// calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
// according to stored eigenstate matrix est
/****************************************************************************/
ComplexVector & jjjpar::MQ(Vector & Qvec)
{double J0,J2,J4,J6;
 double Q,d,s;
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
         double th,ph,Qx,Qy,Qz;							 
	 Qx=Qvec(1);Qy=Qvec(2);Qz=Qvec(3);
	 th=acos(Qz/Q);
	 if(sin(th)>=SMALL){
	                    if(Qx>0){ph=acos(Qx/(Q*sin(th))-SMALL);}
			    else    {ph=acos(Qx/(Q*sin(th))+SMALL);}
			   }
			 else{ph=0;}
	 if (Qy<0){ph=2*PI-ph;} 

 switch (intern_mcalc)
  {case 0: (*mq)(&Mq,&th,&ph,&J0,&J2,&J4,&J6,&est);return Mq;break;
   case 2: return (*iops).MQ(th,ph,J0,J2,J4,J6,Zc,est);break;
   default: fprintf(stderr,"ERROR in scattering operator function M(Q) for ion %s \nM(Q) is currently only implemented for internal module cfield:\n",cffilename);exit(EXIT_FAILURE);
  }
}



//  RETURN TOTAL FORMFACTOR, 
//    however if gJ=0 and Q>0 return spin form factor FS(Q)=<j0(Q)>
//            if gJ=0 and Q<0 return angular  form factor FL(Q)=<j0(Q)>+<j2(Q)>
//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
   double jjjpar::F(double Q)
   {double s,j0,j2;    
    s=Q/4/PI;
    j0=magFFj0(1)*exp(-magFFj0(2)*s*s)+magFFj0(3)*exp(-magFFj0(4)*s*s)+magFFj0(5)*exp(-magFFj0(6)*s*s)+magFFj0(7);
    if(gJ==0&&Q>0){return j0;} // in case of intermediate coupling return spin form factor 
    j2=magFFj2(1)*exp(-magFFj2(2)*s*s)+magFFj2(3)*exp(-magFFj2(4)*s*s)+magFFj2(5)*exp(-magFFj2(6)*s*s)+magFFj2(7);
    j2*=s*s;
    if(gJ==0&&Q<0){return j0+j2;} // in case of intermediate coupling return angular form factor 
   return (j0 + j2 * (2 / gJ - 1)); // formfactor F(Q)

   }

//   debyewallerfactor = exp(-2 * DWF *s*s)
   double jjjpar::debyewallerfactor(double & Q)
   {double s;
    s=Q/4/PI;
    return exp(-2*DWF*s*s);
   }


// returns total angular momentum quantum number J
double jjjpar::J()
{
 switch (intern_mcalc)
  {
   case 2: return (*iops).J;break;
   case 3:  return ABC[1];break;
   default: fprintf (stderr, "error class jjjpar: single ion module does not allow to calculate stevens parameters alpha beta gamma \n"); 
            exit (EXIT_FAILURE);
  }
}

// returns stevens parameters of ion
Vector & jjjpar::tetan ()
{static Vector tt(1,6);
 tt=0;
 switch (intern_mcalc)
  {
   case 2: tt(2)=(*iops).alpha;tt(4)=(*iops).beta;tt(6)=(*iops).gamma;
            return tt;break;
   default: fprintf (stderr, "error class jjjpar: single ion module does not allow to calculate stevens parameters alpha beta gamma \n"); 
            exit (EXIT_FAILURE);
  }
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


//saving parameters to file
void jjjpar::save(FILE * file) 
{ int i,i1,j1;
  saveatom(file);
  fprintf(file,"# x[a]    y[b]      z[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
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
{   fprintf(file,"# da=%4.6g [a] db=%4.6g [b] dc=%4.6g [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g cffilename=%s\n",xyz(1),xyz(2),xyz(3),paranz,diagonalexchange,gJ,cffilename);
}



/*****************************************************************************************/
//constructor with file handle of mcphas.j
jjjpar::jjjpar(FILE * file) 
{ FILE * cf_file;   
  char instr[MAXNOFCHARINLINE];
  cffilename= new char [MAXNOFCHARINLINE];
  int i,j,i1,j1,k1,l;
  double gjcheck;
  float nn[MAXNOFNUMBERSINLINE];
  nn[0]=MAXNOFNUMBERSINLINE;
  xyz=Vector(1,3);
  fgets_errchk (instr, MAXNOFCHARINLINE, file);
  extract(instr,"x",xyz[1]);
  extract(instr,"y",xyz[2]);
  extract(instr,"z",xyz[3]);
  extract(instr,"da",xyz[1]);
  extract(instr,"db",xyz[2]);
  extract(instr,"dc",xyz[3]);
  extract(instr,"nofneighbours",paranz);
  extract(instr,"diagonalexchange",diagonalexchange);
  extract(instr,"gJ",gjcheck);
  extract(instr,"cffilename",cffilename,(size_t)MAXNOFCHARINLINE);
  fgets_errchk (instr, MAXNOFCHARINLINE, file);

// read single ion parameter file and see which type it is (internal module or loadable)
  transitionnumber=1;
  
  //start reading again at the beginning of the file to get formfactors, debye waller factor
  get_parameters_from_sipfile(cffilename);
  if (gJ!=gjcheck){fprintf (stderr, "Error: Lande factor gJ in file mcphas.j and %s are not the same\n",cffilename);
                   exit (EXIT_FAILURE);}
  Mq=ComplexVector(1,3);
//fprintf(stderr,"\nmagnetic formfactors\n");
//fprintf(stderr," FFj0A=%4.4g  FFj0a=%4.4g FFj0B=0%4.4g FFj0b=%4.4g  FFj0C=%4.4g  FFj0c=%4.4g FFj0D=%4.4g\n",magFFj0(1),magFFj0(2),magFFj0(3),magFFj0(4),magFFj0(5),magFFj0(6),magFFj0(7));  
//fprintf(stderr," FFj2A=%4.4g FFj2a=%4.4g FFj2B =%4.4g FFj2b=%4.4g FFj2C=%4.4g  FFj2c=%4.4g FFj2D=%4.4g\n",magFFj2(1), magFFj2(2),magFFj2(3),magFFj2(4),magFFj2(5),magFFj2(6),magFFj2(7));
//fprintf(stderr,"Debey-Waller Factor\n DWF=%g\n\n",DWF);


  nofcomponents=3; // default value for nofcomponents - (important in case nofparameters=0)
// read the exchange parameters from file (exactly paranz parameters!)
  for  (i=1;i<=paranz;++i)
  {while((j=inputline(file, nn))==0&&feof(file)==0){}; // returns 0 if comment line or eof, exits with error, if input string too long
   if(feof(file)!=0){ fprintf (stderr, "Error in jjjpar.cpp: input jjj parameters - \n");
  fprintf(stderr," end of file reached while reading exchange parameter %i(%i)",i,paranz);
      exit (EXIT_FAILURE);
    }
    if(i==1){// determine nofcomponents from number of parameters read in first line of mcphas.j
             if(diagonalexchange==1){nofcomponents=j-3;}else{nofcomponents=(int)sqrt((double)(j-3));}
             if(intern_mcalc==1)
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

// constructor with filename of singleion parameter file
jjjpar::jjjpar(double x,double y,double z, char * sipffile)
{xyz=Vector(1,3);xyz(1)=x;xyz(2)=y;xyz(3)=z;
  mom=Vector(1,9); mom=0; 
  Mq=ComplexVector(1,3);
 get_parameters_from_sipfile(sipffile);

}

void jjjpar::get_parameters_from_sipfile(char * cffilename)
{FILE * cf_file;
 int i,j;
 float nn[MAXNOFNUMBERSINLINE];
 nn[0]=MAXNOFNUMBERSINLINE;
  char modulefilename[MAXNOFCHARINLINE];

 char instr[MAXNOFCHARINLINE];
  cf_file = fopen_errchk (cffilename, "rb");
  fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
  // strip /r (dos line feed) from line if necessary
  char *token;  
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}  
  // determine the file type
  if(strncmp(instr,"#!",2)!=0)
    {fprintf(stderr,"Error: single ion property file %s does not start with '#!'\n",cffilename);
     exit(EXIT_FAILURE);}

  fprintf (stderr,"#parsing single ion property file: %s - loading module %s",cffilename,instr+2);

  if(strncmp(instr,"#!kramer ",9)==0||strncmp(instr,"#!kramer\n",9)==0)
    {intern_mcalc=1;fprintf (stderr,"[internal]\n");
      // input all  lines starting with comments
      while((i=inputparline ("params",cf_file, nn))==0&&feof(cf_file)==false);
      if(i!=3){fprintf(stderr,"Error reading |<+-|Ja|-+>|,|<+-|Jb|-+>|,|<+-|Jc|+->| from file %s\ncorrect file format is:\n");
              fprintf(stderr,"\n#!kramer\n#comment lines ..\n#matrix elements\nparamnames= |<+-|Ja|-+>| |<+-|Jb|-+>| |<+-|Jc|+->|\nparams=2 3 1\n\n",cffilename);exit(EXIT_FAILURE);}
      // now we have the numbers corresponding to the vector ABC() in nn[]
      ABC=Vector(1,i);for(j=1;j<=i;++j){ABC(j)=nn[j];}
      fprintf(stderr," ... kramers doublet with A=<+|Ja|->=%g B=<+-|Jb|+->=+-%g C=<+|Jc|->/i=%g\n",ABC(1),ABC(2),ABC(3));
      est=ComplexMatrix(0,2,1,2); // not used, just initialize to prevent errors
    }
  else
    {if(strncmp(instr,"#!brillouin ",12)==0||strncmp(instr,"#!brillouin\n",12)==0)
     {intern_mcalc=3;fprintf (stderr,"[internal]\n");
      // input all  lines starting with comments
      while((i=inputparline ("params",cf_file, nn))==0&&feof(cf_file)==false);
      if(i!=1){fprintf(stderr,"Error reading spin quantum number J=S from file %s\ncorrect file format is:\n");
              fprintf(stderr,"\n#!brillouin\n#comment lines ..\n#matrix elements\nparamnames= J\nparams=3.5\n\n",cffilename);exit(EXIT_FAILURE);}
      // now we have the numbers corresponding to the vector ABC() in nn[]
      ABC=Vector(1,i);for(j=1;j<=i;++j){ABC(j)=nn[j];}
      fprintf(stderr," ... Brillouin function with J=S=%g\n",ABC(1));
      est=ComplexMatrix(0,2,1,2);// not used, just initialize to prevent errors
      est=0;
     }
     else 
     {if(strncmp(instr,"#!cfield ",9)==0||strncmp(instr,"#!cfield\n",9)==0)
     {intern_mcalc=2;fprintf (stderr,"#[internal]\n");
      iops=new ionpars(cf_file);  
      int dj;dj=(int)(2*J()+1);
      est=ComplexMatrix(0,dj,1,dj);
     // get 1ion parameters - operator matrices
     
     }
     else
     {   
#ifdef __linux__
  instr[1]='=';
  extract(instr,"#",modulefilename,(size_t)MAXNOFCHARINLINE);
  fprintf (stderr,"#[external]\n");
  
       // input all  lines starting with comments
    while((i=inputparline ("params",cf_file, nn))==0&&feof(cf_file)==false);
    // now we have the numbers corresponding to vector ABC() in nn[] - these are the module parameters !
    fprintf(stderr,"#parameters: ");
    if(i>0){
             ABC=Vector(1,i);for(j=1;j<=i;++j){ABC(j)=nn[j];fprintf(stderr,"%g ",nn[j]);} 
            }else{
             ABC=Vector(1,1);
	    } 
    fprintf(stderr,"\n");

  char * error;intern_mcalc=0;
  handle=dlopen (modulefilename,RTLD_NOW | RTLD_GLOBAL);
  if (!handle){fprintf (stderr, "jjjpar::jjjpar - Could not load dynamic library\n");
               if ((error=dlerror())!=NULL) 
	         {fprintf (stderr,"%s\n",error);}
	       exit (EXIT_FAILURE);
	      }
  m=(void(*)(Vector*,double*,Vector*,double*,Vector*,char**,double*,double*))dlsym(handle,"mcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s\n",error);exit (EXIT_FAILURE);}
  dm=(int(*)(int*,double*,Vector*,double*,Vector*,char**,ComplexMatrix*,float*))dlsym(handle,"dmcalc");
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
  fprintf (stderr,"\n Error: non Linux operating system - external loadable modules not supported\n");
  exit(EXIT_FAILURE);
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

  DWF=0;  

  cf_file = fopen_errchk (cffilename, "rb");

  while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
    extract(instr,"SCATTERINGLENGTHREAL",SLR);  
    extract(instr,"SCATTERINGLENGTHIMAG",SLI);  
    extract(instr,"GJ",gJ);  
    extract(instr,"gJ",gJ);  
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
  }    
 }
 
 fclose (cf_file); 
// check gJ
if(intern_mcalc==2&&fabs(gJ-(*iops).gJ)>0.00001)
{fprintf(stderr,"Error internal module cfield : Lande Factor read from %s (gJ=%g) does not conform to internal module value gJ=%g\n",cffilename,gJ,(*iops).gJ);exit(EXIT_FAILURE);}
if (gJ==0){printf("# reading gJ=0 in single ion property file %s -> entering intermediate coupling mode by assigning Ja=Sy Jb=Ly Jc=Sz Jd=Lz Je=Sx Jf=Lx\n",cffilename);
           if (intern_mcalc==1){fprintf(stderr,"Error internal module kramers: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
           if (intern_mcalc==2){fprintf(stderr,"Error internal module cfield : intermediate coupling not supported\n");exit(EXIT_FAILURE);}
           if (intern_mcalc==3){fprintf(stderr,"Error internal module brillouin: intermediate coupling not supported\n");exit(EXIT_FAILURE);}
          }

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

}
//constructor without file
jjjpar::jjjpar(int n,int diag,int nofmom) 
{ cffilename= new char [MAXNOFCHARINLINE];
  diagonalexchange=diag;
  paranz=n;xyz=Vector(1,3);xyz=0;
  int i1;
  intern_mcalc=1;ABC=Vector(1,3);ABC=0;
  transitionnumber=1;
  nofcomponents=nofmom;
  mom=Vector(1,nofcomponents);
  mom=0; 
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
  DWF=0;gJ=0;

}

//kopier-konstruktor
jjjpar::jjjpar (const jjjpar & p)
{ int i;
  xyz=Vector(1,3);
  nofcomponents=p.nofcomponents;
  mom=Vector(1,nofcomponents); 
  xyz=p.xyz;paranz=p.paranz;
  SLR=p.SLR;SLI=p.SLI;

  diagonalexchange=p.diagonalexchange;
  gJ=p.gJ;intern_mcalc=p.intern_mcalc;
  Mq=ComplexVector(1,3);
  Mq=p.Mq;
  
  transitionnumber=p.transitionnumber;
  cffilename= new char [strlen(p.cffilename)+1];
  strcpy(cffilename,p.cffilename);
  if (p.intern_mcalc==1||p.intern_mcalc==0)  ABC=p.ABC;
  if (p.intern_mcalc==2)  {iops=new ionpars((int)(2*(*p.iops).J+1));iops=p.iops;
                           int dj;dj=(int)(2*J()+1);
                           est=ComplexMatrix(0,dj,1,dj);est=p.est;
                           }
//  if (intern_mcalc==2)  iops=new ionpars(4);iops=p.iops;
//  if (intern_mcalc==2)  iops=p.iops;
  
#ifdef __linux__
/*  if (intern_mcalc==0)
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
  if (intern_mcalc==2) delete iops;
#ifdef __linux__
   if (intern_mcalc==0)dlclose(handle);
#endif
}


