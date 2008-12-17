#include "jjjpar.hpp"

#define MAXNOFNUMBERSINLINE 200
#define MAXNOFCHARINLINE 1024

#define MU_B 0.05788
#define K_B  0.0862
#define SMALL 1e-6   //!!! must match SMALL in mcdisp.c and ionpars.cpp !!!
                     // because it is used to decide wether for small transition
		     // energy the matrix Mijkl contains wn-wn' or wn/kT
#define PI 3.1415926535

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



 // *************************************************************************
 // ************************ class parameters *******************************
 // *************************************************************************

//  D = 2 * pi / Q
//  s = 1 / 2 / D: sintheta = lambda * s
   double jjjpar::F(double & Q)
   {double s,j0,j2;    
    s=Q/4/PI;
   j0 = magFF(1) * exp(-magFF(2) * s * s) + magFF(3) * exp(-magFF(4) * s * s);
   j0 = j0 + magFF(5) * exp(-magFF(6) * s * s) + magFF(7);
   j2 = magFF(8) * s * s * exp(-magFF(9) * s * s) + magFF(10) * s * s * exp(-magFF(11) * s * s);
   j2 = j2 + magFF(12) * s * s * exp(-magFF(13) * s * s) + s * s * magFF(14);
   return (j0 + j2 * (2 / gJ - 1)); // formfactor F(Q)

   }

//   debeywallerfactor = exp(-2 * DWF *s*s)
   double jjjpar::debeywallerfactor(double & Q)
   {double s;
    s=Q/4/PI;
    return exp(-2*DWF*s*s);
   }

/****************************************************************************/
// subroutine to calculate magnetisation M from effective field H
// this is the heart of the meanfield algorithm an it is necessary to
// keep this routine as efficient as possible
// at the moment we do only groundstate doublet
Vector & jjjpar::mcalc (double & T, Vector &  gjmbH, double & lnZ,double & U)
{switch (intern_mcalc)
  {case 1: return kramer(T,gjmbH,lnZ,U);break;
   case 2: return (*iops).cfield(T,gjmbH,lnZ,U);break;
   case 3: return brillouin(T,gjmbH,lnZ,U);break;
   default: (*m)(&Jret,&T,&gjmbH,&gJ,&ABC,&lnZ,&U);return Jret;
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

// returns eigenvalues and eigenstates matrix parameters of ion
ComplexMatrix & jjjpar::eigenstates (Vector & gjmbheff)
{switch (intern_mcalc)
  {
   case 2:  return (*iops).cfeigenstates(gjmbheff);break;
   default: fprintf (stderr, "error class jjjpar: single ion module does not allow to return eigenvalues and eigenstates \n"); 
            exit (EXIT_FAILURE);
  }
}

// this function returns n (the number of transitions in the single ion susceptibility)
// the transition matrix mat corresponding to jjjpar.transitionnumber and delta for effective field heff 
int jjjpar::dmcalc(double & T,Vector & gjmbheff,ComplexMatrix & mat,float & delta)
{switch (intern_mcalc)
  {case 1: return kramerdm(transitionnumber,T,gjmbheff,mat,delta);break;
   case 2: return (*iops).cfielddm(transitionnumber,T,gjmbheff,mat,delta);break;
   case 3: return brillouindm(transitionnumber,T,gjmbheff,mat,delta);break;
   default: return (*dm)(&transitionnumber,&T,&gjmbheff,&gJ,&ABC,&mat,&delta);
  }
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



void jjjpar::increase_nofcomponents(int n) // increase nofcomponents by n
{int i,j,k,nold;
  nold=nofcomponents;
  nofcomponents+=n;
  Jret.Resize(1,nofcomponents); 

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

//constructor
jjjpar::jjjpar(FILE * file) 
{ FILE * cf_file;   
  char modulefilename[MAXNOFCHARINLINE];
  char instr[MAXNOFCHARINLINE];
  cffilename= new char [MAXNOFCHARINLINE];
  int i,j,i1,j1,k1,l;
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
  extract(instr,"gJ",gJ);
  extract(instr,"cffilename",cffilename,(size_t)MAXNOFCHARINLINE);
  fgets_errchk (instr, MAXNOFCHARINLINE, file);

// read single ion parameter file and see which type it is (internal module or loadable)
  transitionnumber=1;
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
     }
     else 
     {if(strncmp(instr,"#!cfield ",9)==0||strncmp(instr,"#!cfield\n",9)==0)
     {intern_mcalc=2;fprintf (stderr,"#[internal]\n");
     iops=new ionpars(cf_file);  
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
  m=(void(*)(Vector*,double*,Vector*,double*,Vector*,double*,double*))dlsym(handle,"mcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s\n",error);exit (EXIT_FAILURE);}
  dm=(int(*)(int*,double*,Vector*,double*,Vector*,ComplexMatrix*,float*))dlsym(handle,"dmcalc");
  if ((error=dlerror())!=NULL) {fprintf (stderr,"jjjpar::jjjpar %s -continuing\n",error);dm=NULL;}
#else
  fprintf (stderr,"\n Error: non Linux operating system - external loadable modules not supported\n");
  exit(EXIT_FAILURE);
#endif
    }
   }
  }
  fclose(cf_file);

  magFF=Vector(1,14);
  magFF=0;  magFF[1]=1;
  DWF=0;  
  
  //start reading again at the beginning of the file to get formfactors, debye waller factor
  cf_file = fopen_errchk (cffilename, "rb");

  while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
   // read formfactor if given
    extract(instr,"FFj0A",magFF[1]);
    extract(instr,"FFj0a",magFF[2]);
    extract(instr,"FFj0B",magFF[3]);
    extract(instr,"FFj0b",magFF[4]);
    extract(instr,"FFj0C",magFF[5]);
    extract(instr,"FFj0c",magFF[6]);
    extract(instr,"FFj0D",magFF[7]);
    extract(instr,"FFj2A",magFF[8]);
    extract(instr,"FFj2a",magFF[9]);
    extract(instr,"FFj2B",magFF[10]);
    extract(instr,"FFj2b",magFF[11]);
    extract(instr,"FFj2C",magFF[12]);
    extract(instr,"FFj2c",magFF[13]);
    extract(instr,"FFj2D",magFF[14]);
   // read debeywallerfactor if given    
    extract(instr,"DWF",DWF);
  }    
 }
 
fprintf(stderr,"\nmagnetic formfactors\n");
fprintf(stderr," FFj0A=%4.4g  FFj0a=%4.4g FFj0B=0%4.4g FFj0b=%4.4g  FFj0C=%4.4g  FFj0c=%4.4g FFj0D=%4.4g\n",magFF(1),magFF(2),magFF(3),magFF(4),magFF(5),magFF(6),magFF(7));  
fprintf(stderr," FFj2A=%4.4g FFj2a=%4.4g FFj2B =%4.4g FFj2b=%4.4g FFj2C=%4.4g  FFj2c=%4.4g FFj2D=%4.4g\n",magFF(8), magFF(9),magFF(10),magFF(11),magFF(12),magFF(13),magFF(14));
fprintf(stderr,"Debey-Waller Factor\n DWF=%g\n\n",DWF);
 
 fclose (cf_file); 

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

  Jret=Vector(1,nofcomponents); 

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

//constructor without file
jjjpar::jjjpar(int n,int diag,int nofmom) 
{ cffilename= new char [MAXNOFCHARINLINE];
  diagonalexchange=diag;
  paranz=n;xyz=Vector(1,3);
  int i1;
  intern_mcalc=1;ABC=Vector(1,3);
  transitionnumber=1;
  nofcomponents=nofmom;
  Jret=Vector(1,nofcomponents); 
  dn = new Vector[n+1];for(i1=0;i1<=n;++i1){dn[i1]=Vector(1,3);}
  if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  jij = new Matrix[n+1];for(i1=0;i1<=n;++i1){jij[i1]=Matrix(1,nofcomponents,1,nofcomponents);}
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  magFF=Vector(1,14);
  magFF=0;  magFF[1]=1;
  DWF=0;  

}

//kopier-konstruktor
jjjpar::jjjpar (const jjjpar & p)
{ int i;
  xyz=Vector(1,3);
  nofcomponents=p.nofcomponents;
  Jret=Vector(1,nofcomponents); 
  xyz=p.xyz;paranz=p.paranz;
  diagonalexchange=p.diagonalexchange;
  gJ=p.gJ;intern_mcalc=p.intern_mcalc;
  transitionnumber=p.transitionnumber;
  cffilename= new char [strlen(p.cffilename)+1];
  strcpy(cffilename,p.cffilename);
  if (p.intern_mcalc==1||p.intern_mcalc==0)  ABC=p.ABC;
  if (p.intern_mcalc==2)  iops=new ionpars((int)(2*(*p.iops).J+1));iops=p.iops;
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
/*  }*/
#endif
  magFF=Vector(1,14);
  magFF=p.magFF;
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
{ delete []jij;
  delete []dn;
  delete []sublattice;
  delete []cffilename;
  if (intern_mcalc==2) delete iops;
#ifdef __linux__
   if (intern_mcalc==0)dlclose(handle);
#endif
}

//------------------------------------------------------------------------------------------------
//routine mcalc for kramers doublet
//------------------------------------------------------------------------------------------------
Vector & jjjpar::kramer (double & T, Vector & gjmbH, double & lnZ, double & U)
{ /*on input
    ABC(1...3)  A,M,Ci....saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/
  double alpha, betar, betai, lambdap,lambdap_K_BT, lambdap2, expp, expm, np, nm;
  double nennerp, nennerm, jap, jam, jbp, jbm, jcp, jcm,Z;
  double alpha_lambdap,alphaplambdap,alphaxlambdap;

  static Vector Jret(1,3);
  
  alpha = ABC[2] * gjmbH[2];
  betar = -ABC[1] * gjmbH[1];
  betai = -ABC[3] * gjmbH[3];
  lambdap2 = alpha * alpha + betar * betar + betai * betai;
  lambdap = sqrt (lambdap2);
  lambdap_K_BT=lambdap/K_B/T;
  if (lambdap_K_BT>700){lambdap_K_BT=700;}
  if (lambdap_K_BT<-700){lambdap_K_BT=-700;}
  expm = exp (lambdap_K_BT);
  expp = 1/expm; //=exp (-lambdap_K_BT);
  Z = expp + expm;
  lnZ=log(Z);
  np = expp / Z;
  nm = expm / Z;
  U=lambdap*(np-nm); // energy

//  nennerp = (alpha - lambdap) * (alpha - lambdap) + betar * betar + betai * betai;
//  nennerm = (alpha + lambdap) * (alpha + lambdap) + betar * betar + betai * betai;
    alphaxlambdap=alpha*lambdap;
    alpha_lambdap=alpha-lambdap;
    alphaplambdap=alpha+lambdap;
    nennerp=  2.0*(-alphaxlambdap+lambdap2);    
    nennerm=  2.0*(alphaxlambdap+lambdap2);    

  if (nennerp > SMALL)
    {
      jap = -ABC[1] * 2.0 * betar * (alpha_lambdap) / nennerp;
//      jbp = M * ((alpha_lambdap) * (alpha_lambdap) - (betar * betar + betai * betai)) / nennerp;
      jbp = ABC[2] * (2.0 * alpha*alpha_lambdap) / nennerp;
      jcp = -2.0 * ABC[3] * betai * (alpha_lambdap) / nennerp;
    }
  else
    {
      jap = 0;
      if (alpha * alpha > SMALL)
	{
	  jbp = -copysign (ABC[2], alpha);
	}
      else
	{
	  jbp = ABC[2];
	}
      jcp = 0;
    }

  if (nennerm > SMALL)
    {
      jam = -ABC[1] * 2.0 * betar * (alphaplambdap) / nennerm;
//      jbm = M * ((alpha + lambdap) * (alpha + lambdap) - (betar * betar + betai * betai)) / nennerm;
      jbm = ABC[2] * (2.0 * alpha*alphaplambdap) / nennerm;
      jcm = -2.0 * ABC[3] * betai * (alphaplambdap) / nennerm;
    }
  else
    {
      jam = 0;
      if (alpha * alpha > SMALL)
	{
	  jbm = copysign (ABC[2], alpha);
	}
      else
	{
	  jbm = -ABC[2];
	}
      jcm = 0;
    }

  Jret[1] = np * jap + nm * jam;
  Jret[2] = np * jbp + nm * jbm;
  Jret[3] = np * jcp + nm * jcm;
//  printf ("Ha=%g Hb=%g Hc=%g Ja=%g Jb=%g Jc=%g \n", 
//     gjmbH[1]/MU_B/gjJ, gjmbH[2]/MU_B/gjJ, gjmbH[3]/MU_B/gjJ, J[1], J[2], J[3]);
return Jret;
}

int jjjpar::kramerdm(int & transitionnumber,double & T,Vector & gjmbH,ComplexMatrix & mat,float & delta)
{ 
  /*on input
    transitionnumber ... number of transition to be computed - meaningless for kramers doublet, because there is only 1 transition
    ABC[i]	saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splitting of kramers doublet [meV]
    mat(i,j)	<-|(Ji-<Ji>)|+><+|(Jj-<Jj>|-> tanh(delta/2kT)
*/
  double alpha, betar, betai, lambdap,lambdap_K_BT, lambdap2, expp, expm, np, nm;
  double nennerp, nennerm, nenner;
  complex<double> ja,jb,jc,i(0,1), jap, jam, jbp, jbm, jcp, jcm;
  double alpha_lambdap,alphaplambdap,alphaxlambdap;
  double Z;
  double lnz,u;
  
  static Vector Jret(1,3);
  // clalculate thermal expectation values (needed for quasielastic scattering)
  Jret=kramer(T,gjmbH,lnz,u);
  int pr;
  pr=1;
  if (transitionnumber<0) {pr=0;transitionnumber*=-1;}

  alpha = ABC[2]* gjmbH[2];
  betar = -ABC[1] * gjmbH[1];
  betai = -ABC[3] * gjmbH[3];
  lambdap2 = alpha * alpha + betar * betar + betai * betai;
  lambdap = sqrt (lambdap2);


  
  lambdap_K_BT=lambdap/K_B/T;
  if (lambdap_K_BT>700){lambdap_K_BT=700;}
  if (lambdap_K_BT<-700){lambdap_K_BT=-700;}
  expm = exp (lambdap_K_BT);
  expp = 1/expm; //=exp (-lambdap_K_BT);
  Z = expp + expm;
  np = expp / Z;
  nm = expm / Z;

    alphaxlambdap=alpha*lambdap;
    alpha_lambdap=alpha-lambdap;
    alphaplambdap=alpha+lambdap;
    nennerp=  2.0*(-alphaxlambdap+lambdap2);    
    nennerm=  2.0*(alphaxlambdap+lambdap2);    


if (transitionnumber==2)
{ delta=2*lambdap; //set delta !!!


    nenner=sqrt(nennerp*nennerm);

  if (nenner > SMALL)
    {
      ja = -ABC[1] * 2.0*(alpha * betar+i * betai * lambdap) / nenner;
      jb = -ABC[2] * 2.0 * (betar*betar+betai*betai) / nenner;
      jc = -ABC[3] * 2.0*(alpha*betai -i *betar*lambdap) / nenner;
    }
  else
    {
      if (alpha > SMALL)
	{ja = ABC[1];  // <-| is the ground state
  	 jb = 0;
         jc = -i*ABC[3];
	}
      else
	{ja = ABC[1];  // <+| is the ground state
  	 jb = 0;
         jc = i*ABC[3]; 	
	}
    }
 if (delta>SMALL)
 {// now lets calculate mat
 mat(1,1)=ja*conj(ja)*(nm-np);
 mat(1,2)=ja*conj(jb)*(nm-np);
 mat(1,3)=ja*conj(jc)*(nm-np);
 mat(2,1)=jb*conj(ja)*(nm-np);
 mat(2,2)=jb*conj(jb)*(nm-np);
 mat(2,3)=jb*conj(jc)*(nm-np);
 mat(3,1)=jc*conj(ja)*(nm-np);
 mat(3,2)=jc*conj(jb)*(nm-np);
 mat(3,3)=jc*conj(jc)*(nm-np); 
 }else
 {// quasielastic scattering needs epsilon * nm / KT ....
 mat(1,1)=ja*conj(ja)*nm/K_B/T;
 mat(1,2)=ja*conj(jb)*nm/K_B/T;
 mat(1,3)=ja*conj(jc)*nm/K_B/T;
 mat(2,1)=jb*conj(ja)*nm/K_B/T;
 mat(2,2)=jb*conj(jb)*nm/K_B/T;
 mat(2,3)=jb*conj(jc)*nm/K_B/T;
 mat(3,1)=jc*conj(ja)*nm/K_B/T;
 mat(3,2)=jc*conj(jb)*nm/K_B/T;
 mat(3,3)=jc*conj(jc)*nm/K_B/T; 
 }
}
else
{ delta=-SMALL; // transition within the same level
  if (nennerp > SMALL)
    {
      jap = -ABC[1] * 2.0 * betar * (alpha_lambdap) / nennerp;
//      jbp = M * ((alpha_lambdap) * (alpha_lambdap) - (betar * betar + betai * betai)) / nennerp;
      jbp = ABC[2] * (2.0 * alpha*alpha_lambdap) / nennerp;
      jcp = -2.0 * ABC[3] * betai * (alpha_lambdap) / nennerp;
    }
  else
    {
      jap = 0;
      if (alpha * alpha > SMALL)
	{
	  jbp = -copysign (ABC[2], alpha);
	}
      else
	{
	  jbp = -ABC[2];
	}
      jcp = 0;
    }

  if (nennerm > SMALL)
    {
      jam = -ABC[1] * 2.0 * betar * (alphaplambdap) / nennerm;
//      jbm = M * ((alpha + lambdap) * (alpha + lambdap) - (betar * betar + betai * betai)) / nennerm;
      jbm = ABC[2] * (2.0 * alpha*alphaplambdap) / nennerm;
      jcm = -2.0 * ABC[3] * betai * (alphaplambdap) / nennerm;
    }
  else
    {
      jam = 0;
      if (alpha * alpha > SMALL)
	{
	  jbm = copysign (ABC[2], alpha);
	}
      else
	{
	  jbm = ABC[2];
	}
      jcm = 0;
    }
 if (transitionnumber==1)
 {// now lets calculate mat
 mat(1,1)=(jam-Jret(1))*(jam-Jret(1))*nm/K_B/T;
 mat(1,2)=(jam-Jret(1))*(jbm-Jret(2))*nm/K_B/T;
 mat(1,3)=(jam-Jret(1))*(jcm-Jret(3))*nm/K_B/T;
 mat(2,1)=(jbm-Jret(2))*(jam-Jret(1))*nm/K_B/T;
 mat(2,2)=(jbm-Jret(2))*(jbm-Jret(2))*nm/K_B/T;
 mat(2,3)=(jbm-Jret(2))*(jcm-Jret(3))*nm/K_B/T;
 mat(3,1)=(jcm-Jret(3))*(jam-Jret(1))*nm/K_B/T;
 mat(3,2)=(jcm-Jret(3))*(jbm-Jret(2))*nm/K_B/T;
 mat(3,3)=(jcm-Jret(3))*(jcm-Jret(3))*nm/K_B/T;
 }else{
 // now lets calculate mat
 mat(1,1)=(jap-Jret(1))*(jap-Jret(1))*np/K_B/T;
 mat(1,2)=(jap-Jret(1))*(jbp-Jret(2))*np/K_B/T;
 mat(1,3)=(jap-Jret(1))*(jcp-Jret(3))*np/K_B/T;
 mat(2,1)=(jbp-Jret(2))*(jap-Jret(1))*np/K_B/T;
 mat(2,2)=(jbp-Jret(2))*(jbp-Jret(2))*np/K_B/T;
 mat(2,3)=(jbp-Jret(2))*(jcp-Jret(3))*np/K_B/T;
 mat(3,1)=(jcp-Jret(3))*(jap-Jret(1))*np/K_B/T;
 mat(3,2)=(jcp-Jret(3))*(jbp-Jret(2))*np/K_B/T;
 mat(3,3)=(jcp-Jret(3))*(jcp-Jret(3))*np/K_B/T;
 }
}
if (pr==1) printf ("delta=%4.6g meV\n",delta);

return 3; // kramers doublet has always exactly one transition + 2 levels (quasielastic scattering)!
}
/**************************************************************************/

//------------------------------------------------------------------------------------------------
//routine mcalc for brillouin function
//------------------------------------------------------------------------------------------------
Vector & jjjpar::brillouin (double & T, Vector & gjmbH, double & lnZ, double & U)
{ /*on input
    ABC(1)  J=S....Spin quantum number
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/

static Vector Jret(1,3);

// check dimensions of vector
if(Jret.Hi()!=3||gjmbH.Hi()!=3||ABC.Hi()!=1)
   {fprintf(stderr,"Error loadable module brillouin.so: wrong number of dimensions - check number of columns in file mcphas.j or number of parameters in single ion property file\n");
    exit(EXIT_FAILURE);}
    
double JJ,K_BT,XJ,gmhkt,Jav,gmh,Z,X;

// program brillouin function for S=J=ABC(1)
JJ=ABC[1];
K_BT=T*K_B;
gmh=Norm(gjmbH);
gmhkt=gmh/K_BT;

if(JJ*gmhkt>100||gmhkt>100){Jav=JJ;lnZ=JJ*gmhkt;}
else
{X=exp(gmhkt);
 XJ=exp(JJ*gmhkt);

// printf("1-X=%g gmhkt=%g\n",1-X,gmhkt);

 if (X<=1.000001){Z=2*JJ+1;Jav=0;}
 else
 {Z=(XJ*X-1/XJ)/(X-1.0);
  Jav=JJ*(XJ*X*X-1/XJ)+(JJ+1)*X*(1.0/XJ-XJ);
  Jav/=(X-1);
  Jav/=(XJ*X-1/XJ);
 }
//for (i=-JJ*2;i<=+0.000001;++i)
//{dd=i*gmhkt;
// if (dd<-700){expp=0;}else{expp=exp(dd);}
// Z += expp; //this is not yet Z, a factor exp(J gJ Heff/kT) is missing
//}


//Jav=0;
//for (i=-JJ*2;i<=+0.000001;++i)
//{dd=i*gmhkt;
// if (dd<-700){expp=0;}else{expp=exp(dd);}
// Jav+=(JJ+i)*expp/Z;
//}
//Z*=exp(JJ*gmhkt); //this is now the correct Z
lnZ=log(Z);
}

U=-gmh*Jav;



if (gmh>0)
{ Jret[1] = Jav*gjmbH(1)/gmh;
  Jret[2] = Jav*gjmbH(2)/gmh;
  Jret[3] = Jav*gjmbH(3)/gmh;
 }
 else
 {Jret=0;}
//  printf ("Ha=%g Hb=%g Hc=%g ma=%g mb=%g mc=%g \n", H[1], H[2], H[3], m[1], m[2], m[3]);
return Jret;
}
/**************************************************************************/
// for mcdisp this routine is needed
int jjjpar::brillouindm(int & tn,double & T,Vector & gjmbH,ComplexMatrix & mat,float & delta)
{ 
  /*on input
    tn          transition-number
    ABC(1)      S=J spin quantum number
    g_J		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splittings [meV] 
    mat(i,j)	transition matrix elements ...
*/
// NOT IMPLEMENTED
// printf("Error: external module brillouin.c has not implemented dmcalc function");exit(EXIT_FAILURE);
// return 1; 
static Vector Jret(1,3);
int pr;

  pr=1;
  if (tn<0) {pr=0;tn*=-1;}

  double JJ,K_BT,XJ,gmhkt,gmh,Z,R,X,sinth,hxxyy,jjkt;
  complex <double> i(0,1),bx,by,bz;

// program brillouin function for S=J=ABC(1)
  JJ=ABC[1];
  K_BT=T*K_B;
  gmh=Norm(gjmbH);
  gmhkt=gmh/K_BT;
  X=exp(gmhkt);
  XJ=exp(JJ*gmhkt);
// calculate Z and R
if (X==1.0){Z=2*JJ+1;R=0;}
else
{Z=(XJ*X-1/XJ)/(X-1.0);
 R=JJ*(1/XJ-XJ*X*X)+(JJ+1)*X*(XJ-1.0/XJ);
 R/=0.5*(X-1)*(X-1);
}

// calculate coefficients bx,by,bz
 hxxyy=gjmbH(1)*gjmbH(1)+gjmbH(2)*gjmbH(2);
 if (hxxyy/gjmbH(3)/gjmbH(3)>SMALL*SMALL)
 {sinth=sqrt(hxxyy)/gmh;
  bx=-gjmbH(2)+i*gjmbH(1)*gjmbH(3)/gmh;
  bx/=2*gmh*sinth;
  by=gjmbH(1)+i*gjmbH(2)*gjmbH(3)/gmh;
  by/=2*gmh*sinth;
  }
 else
 {sinth=0;by=0.5;
  if(gjmbH(3)>0)
  {bx=0.5*i;}
  else
  {bx=-0.5*i;}
 }
  bz=-i*sinth*0.5;
// -----------------------------------------

if (tn==2) // transition to finite energy
 {delta=gmh; //set delta !!!

 if (delta>SMALL)
  {// now lets calculate mat
  mat(1,1)=bx*conj(bx)*(-R/Z);
  mat(1,2)=bx*conj(by)*(-R/Z);
  mat(1,3)=bx*conj(bz)*(-R/Z);
  mat(2,1)=by*conj(bx)*(-R/Z);
  mat(2,2)=by*conj(by)*(-R/Z);
  mat(2,3)=by*conj(bz)*(-R/Z);
  mat(3,1)=bz*conj(bx)*(-R/Z);
  mat(3,2)=bz*conj(by)*(-R/Z);
  mat(3,3)=bz*conj(bz)*(-R/Z);
  } else
  {// quasielastic scattering needs epsilon * nm / KT ....
  jjkt=0.6666667*JJ*(JJ+1)/K_BT;
  mat(1,1)=bx*conj(bx)*jjkt;
  mat(1,2)=bx*conj(by)*jjkt;
  mat(1,3)=bx*conj(bz)*jjkt;
  mat(2,1)=by*conj(bx)*jjkt;
  mat(2,2)=by*conj(by)*jjkt;
  mat(2,3)=by*conj(bz)*jjkt;
  mat(3,1)=bz*conj(bx)*jjkt;
  mat(3,2)=bz*conj(by)*jjkt;
  mat(3,3)=bz*conj(bz)*jjkt;
  }
 }
 else
 { delta=-SMALL; // tn=1 ... transition within the same level
   if(X==1.0){jjkt=JJ*(2*JJ*JJ+3*JJ+1)/3/K_BT/(2*JJ+1);}
   else {jjkt=(1-2*JJ-2*JJ*JJ)/XJ;
         jjkt+=JJ*JJ/X/XJ;
	 jjkt+=(JJ*JJ+2*JJ+1)*X/XJ;
	 jjkt-=(JJ+1)*(JJ+1)*XJ;
	 jjkt+=(2*JJ*JJ+2*JJ-1)*XJ*X;
	 jjkt-=JJ*JJ*XJ*X*X;
	 jjkt*=X/(1-X)/(1-X);
	 jjkt/=(1/XJ-X*XJ)*K_BT;}
 // now lets calculate mat
 mat(1,1)=gjmbH(1)*gjmbH(1)*jjkt;
 mat(1,2)=gjmbH(1)*gjmbH(2)*jjkt;
 mat(1,3)=gjmbH(1)*gjmbH(3)*jjkt;
 mat(2,1)=gjmbH(2)*gjmbH(1)*jjkt;
 mat(2,2)=gjmbH(2)*gjmbH(2)*jjkt;
 mat(2,3)=gjmbH(2)*gjmbH(3)*jjkt;
 mat(3,1)=gjmbH(3)*gjmbH(1)*jjkt;
 mat(3,2)=gjmbH(3)*gjmbH(2)*jjkt;
 mat(3,3)=gjmbH(3)*gjmbH(3)*jjkt;
 }
if (pr==1) printf ("delta=%4.6g meV\n",delta);

return 2;
// brillouin function has 2 effective transitions
}
