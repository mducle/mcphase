#include "jjjpar.hpp"

#define MAXNOFNUMBERSINLINE 200
#define MAXNOFCHARINLINE 1024

#define MU_B 0.05788
#define K_B  0.0862
#define SMALL 1e-10
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
   case 2: return cfield(T,gjmbH,lnZ,U);break;
   default: (*m)(&J,&T,&gjmbH,&gJ,&ABC,&lnZ,&U);return J;
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

// this function returns n (the number of transitions in the single ion susceptibility)
// the transition matrix mat corresponding to jjjpar.transitionnumber and delta for effective field heff 
int jjjpar::dmcalc(double & T,Vector & gjmbheff,ComplexMatrix & mat,float & delta)
{switch (intern_mcalc)
  {case 1: return kramerdm(transitionnumber,T,gjmbheff,mat,delta);break;
   case 2: return cfielddm(transitionnumber,T,gjmbheff,mat,delta);break;
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
{   fprintf(file,"# x=%4.6g [a] y=%4.6g [b] z=%4.6g [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g cffilename=%s\n",xyz(1),xyz(2),xyz(3),paranz,diagonalexchange,gJ,cffilename);
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
  jijn = new Matrix[paranz+1](1,nofcomponents,1,nofcomponents);
  dnn = new Vector[paranz+1](1,3);
  
  int i;
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
  dn = new Vector[paranz+1](1,3);
  if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  jij = new Matrix[paranz+1](1,nofcomponents,1,nofcomponents);
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
  extract(instr,"nofneighbours",paranz);
  extract(instr,"diagonalexchange",diagonalexchange);
  extract(instr,"gJ",gJ);
  extract(instr,"cffilename",cffilename,(size_t)MAXNOFCHARINLINE);
  fgets_errchk (instr, MAXNOFCHARINLINE, file);

// read single ion parameter file and see which type it is (internal module or loadable)
  transitionnumber=1;
  cf_file = fopen_errchk (cffilename, "rb");
  fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
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
      // now we have the numbers corresponding to vector ABC() in nn[]
      ABC=Vector(1,i);for(j=1;j<=i;++j){ABC(j)=nn[j];}
      fprintf(stderr," ... kramers doublet with A=<+|Ja|->=%g B=<+-|Jb|+->=+-%g C=<+|Jc|->/i=%g\n",ABC(1),ABC(2),ABC(3));
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
  magFF=Vector(1,14);
  Blm=Vector(1,45);Blm=0;
  magFF=0;  magFF[1]=1;
  DWF=0;  
  
  fclose(cf_file);
  //start reading again at the beginning of the file
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
   // read crystal field parameters 
   if (intern_mcalc==2){
   extract(instr,"B2-2",Blm(1));
   extract(instr,"B2-1",Blm(2));
   extract(instr,"B20",Blm(3));
   extract(instr,"B21",Blm(4));
   extract(instr,"B22",Blm(5));
   
   extract(instr,"B3-3",Blm(6));
   extract(instr,"B3-2",Blm(7));
   extract(instr,"B3-1",Blm(8));
   extract(instr,"B30",Blm(9));
   extract(instr,"B31",Blm(10));
   extract(instr,"B32",Blm(11));
   extract(instr,"B32",Blm(12));

   extract(instr,"B4-4",Blm(13));
   extract(instr,"B4-3",Blm(14));
   extract(instr,"B4-2",Blm(15));
   extract(instr,"B4-1",Blm(16));
   extract(instr,"B40",Blm(17));
   extract(instr,"B41",Blm(18));
   extract(instr,"B42",Blm(19));
   extract(instr,"B43",Blm(20));
   extract(instr,"B44",Blm(21));
  
   extract(instr,"B5-5",Blm(22));
   extract(instr,"B5-4",Blm(23));
   extract(instr,"B5-3",Blm(24));
   extract(instr,"B5-2",Blm(25));
   extract(instr,"B5-1",Blm(26));
   extract(instr,"B50",Blm(27));
   extract(instr,"B51",Blm(28));
   extract(instr,"B52",Blm(29));
   extract(instr,"B53",Blm(30));
   extract(instr,"B54",Blm(31));
   extract(instr,"B55",Blm(32));
 
   extract(instr,"B6-6",Blm(33));
   extract(instr,"B6-5",Blm(34));
   extract(instr,"B6-4",Blm(35));
   extract(instr,"B6-3",Blm(36));
   extract(instr,"B6-2",Blm(37));
   extract(instr,"B6-1",Blm(38));
   extract(instr,"B60",Blm(39));
   extract(instr,"B61",Blm(40));
   extract(instr,"B62",Blm(41));
   extract(instr,"B63",Blm(42));
   extract(instr,"B64",Blm(43));
   extract(instr,"B65",Blm(44));
   extract(instr,"B66",Blm(45));
  }}    
 }
if (intern_mcalc==2){
if((*iops).Hcf==(double)0.0){
// calculation of the cf matrix according 
fprintf(stderr,"crystal field parameters\n");  
for(l=1;l<=45;++l){(*iops).Hcf+=Blm(l)*(*(*iops).Olm[l]);
                   if(Blm(l)!=0){if(l<24){fprintf(stderr,"B%c=%g   ",l+99,Blm(l));}
		                     else{fprintf(stderr,"B(z+%i)=%g   ",l-23,Blm(l));}
		                }
                  }
}}
 
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
             dn = new Vector[paranz+1](1,3);if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
             sublattice = new int[paranz+1];if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
             jij = new Matrix[paranz+1](1,nofcomponents,1,nofcomponents);if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
            }
   //check if correct number of columns has been read	        
    if((diagonalexchange==1&&nofcomponents!=j-3)||(diagonalexchange==0&&nofcomponents!=(int)sqrt((double)(j-3))))
              {fprintf(stderr,"Error reading mcphas.j line %i: check number of columns\n",i);
               exit(EXIT_FAILURE);}

  J=Vector(1,nofcomponents); 

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
  
  intern_mcalc=1;ABC=Vector(1,3);
  transitionnumber=1;
  nofcomponents=nofmom;
  J=Vector(1,nofcomponents); 
  dn = new Vector[n+1](1,3);
  if (dn == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  sublattice = new int[paranz+1];
  if (sublattice == NULL){ fprintf (stderr, "Out of memory\n"); exit (EXIT_FAILURE);}
  jij = new Matrix[n+1](1,nofcomponents,1,nofcomponents);
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
  J=Vector(1,nofcomponents); 
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

//dimension arrays
  jij = new Matrix[paranz+1](1,nofcomponents,1,nofcomponents);
  if (jij == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  dn = new Vector[paranz+1](1,3);
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

  static Vector J(1,3);
  
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
	  jbp = 0;
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
	  jbm = 0;
	}
      jcm = 0;
    }

  J[1] = np * jap + nm * jam;
  J[2] = np * jbp + nm * jbm;
  J[3] = np * jcp + nm * jcm;
//  printf ("Ha=%g Hb=%g Hc=%g Ja=%g Jb=%g Jc=%g \n", 
//     gjmbH[1]/MU_B/gjJ, gjmbH[2]/MU_B/gjJ, gjmbH[3]/MU_B/gjJ, J[1], J[2], J[3]);
return J;
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
  
  static Vector J(1,3);
  // clalculate thermal expectation values (needed for quasielastic scattering)
  J=kramer(T,gjmbH,lnz,u);
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
// now lets calculate mat
mat(1,1)=ja*conj(ja)*(nm-np);
mat(1,2)=ja*conj(jb)*(nm-np);
mat(1,3)=ja*conj(jc)*(nm-np);
mat(2,1)=jb*conj(ja)*(nm-np);
mat(2,2)=jb*conj(jb)*(nm-np);
mat(2,3)=jb*conj(jc)*(nm-np);
mat(3,1)=jc*conj(ja)*(nm-np);
mat(3,2)=jc*conj(jb)*(nm-np);
mat(3,3)=jc*conj(jc)*(nm-np);
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
	  jbp = 0;
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
	  jbm = 0;
	}
      jcm = 0;
    }
 if (transitionnumber==1)
 {// now lets calculate mat
 mat(1,1)=(jam-J(1))*(jam-J(1))*nm/K_B/T;
 mat(1,2)=(jam-J(1))*(jbm-J(2))*nm/K_B/T;
 mat(1,3)=(jam-J(1))*(jcm-J(3))*nm/K_B/T;
 mat(2,1)=(jbm-J(2))*(jam-J(1))*nm/K_B/T;
 mat(2,2)=(jbm-J(2))*(jbm-J(2))*nm/K_B/T;
 mat(2,3)=(jbm-J(2))*(jcm-J(3))*nm/K_B/T;
 mat(3,1)=(jcm-J(3))*(jam-J(1))*nm/K_B/T;
 mat(3,2)=(jcm-J(3))*(jbm-J(2))*nm/K_B/T;
 mat(3,3)=(jcm-J(3))*(jcm-J(3))*nm/K_B/T;
 }else{
 // now lets calculate mat
 mat(1,1)=(jap-J(1))*(jap-J(1))*np/K_B/T;
 mat(1,2)=(jap-J(1))*(jbp-J(2))*np/K_B/T;
 mat(1,3)=(jap-J(1))*(jcp-J(3))*np/K_B/T;
 mat(2,1)=(jbp-J(2))*(jap-J(1))*np/K_B/T;
 mat(2,2)=(jbp-J(2))*(jbp-J(2))*np/K_B/T;
 mat(2,3)=(jbp-J(2))*(jcp-J(3))*np/K_B/T;
 mat(3,1)=(jcp-J(3))*(jap-J(1))*np/K_B/T;
 mat(3,2)=(jcp-J(3))*(jbp-J(2))*np/K_B/T;
 mat(3,3)=(jcp-J(3))*(jcp-J(3))*np/K_B/T;
 }
}
if (pr==1) printf ("delta=%4.6g meV\n",delta);

return 3; // kramers doublet has always exactly one transition + 2 levels (quasielastic scattering)!
}
/**************************************************************************/


//------------------------------------------------------------------------------------------------
// ROUTINE CFIELD for full crystal field + higher order interactions
//------------------------------------------------------------------------------------------------
Vector & jjjpar::cfield(double & T, Vector & gjmbH, double & lnZs, double & U)
//void cf(Vector & J,float * T,Vector & gjmbH, double * gJ,Vector & ABC, double * Z,double * U)
{//ABC not used !!!
    /*on input
    T		temperature[K]
    gJmbH	vector of effective field [meV]
    gJ          Lande factor
    ABC         single ion parameter values (A, B, C corresponding to <+|Ja|->,<-|Jb|->,<+|Jc|->/i
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/

// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error internal module cfield: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)
static Vector J(1,gjmbH.Hi());
   // setup hamiltonian
   int dj,j;
   dj=(*iops).Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
   ComplexMatrix z(1,dj,1,dj);
   ComplexMatrix za(1,dj,1,dj);
   ComplexMatrix zb(1,dj,1,dj);
   ComplexMatrix zc(1,dj,1,dj);
   ComplexMatrix zolm(1,dj,1,dj);    

   Ham=(*iops).Hcf-gjmbH(1)*(*iops).Ja-gjmbH(2)*(*iops).Jb-gjmbH(3)*(*iops).Jc;

   for(j=4;j<=J.Hi();++j){Ham-=gjmbH(j)*(*(*iops).Olm[j-3]);}

/*   int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",(*(*iops).Olm[j])(i1,j1));
    printf ("\n");
    }*/
      
    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=0;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);
     double x,y;int i;
     x=Min(En);
     for (i=1;i<=dj;++i)
     {if ((y=(En(i)-x)/K_B/T)<600) wn[i]=exp(-y); 
      else wn[i]=0.0;
    //  printf("%g\n",En(i));
      }
     double Zs;
     Zs=Sum(wn);wn/=Zs;
     lnZs=log(Zs)-x/K_B/T;
   // calculate U
     U=En*wn;
   // calculate Ja,Jb,Jc
     z=ComplexMatrix(zr,zi);
     
     za=(*iops).Jaa*z;
     zb=(*iops).Jbb*z;
     zc=(*iops).Jcc*z;

    
     J=0;
//    ComplexVector ddd;
    for (i=1;i<=dj;++i)
    {
     J[1]+=wn(i)*real(z.Column(i)*za.Column(i));
     J[2]+=wn(i)*real(z.Column(i)*zb.Column(i));
     J[3]+=wn(i)*real(z.Column(i)*zc.Column(i));
    }
     
   for(j=4;j<=J.Hi();++j)
   {
    zolm=(*(*iops).OOlm[j-3])*z;
    for (i=1;i<=dj;++i) J[j]+=wn(i)*real(z.Column(i)*zolm.Column(i));
   };


return J;
}
/**************************************************************************/

/**************************************************************************/
// for mcdisp this routine is needed
int jjjpar::cfielddm(int & tn,double & T,Vector & gjmbH,ComplexMatrix & mat,float & delta)
{  /*on input
    tn      ... number of transition to be computed 
    ABC[i]	(not used)saturation moment/gJ[MU_B] of groundstate doublet in a.b.c direction
    gJ		lande factor
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta-+	energy of transition [meV]
    mat(i,j)	<-|Ji|+><+|Jj|-> (n+-n-),  n+,n-
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/


// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error loadable module cfield.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)

static Vector J(1,gjmbH.Hi());
double lnz,u;
J=cfield(T,gjmbH,lnz,u);  //expectation values <J>
  int pr;
  pr=1;
  if (tn<0) {pr=0;tn*=-1;}

   // setup hamiltonian
   int dj,j;
   dj=(*iops).Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
    
   Ham=(*iops).Hcf-gjmbH(1)*(*iops).Ja-gjmbH(2)*(*iops).Jb-gjmbH(3)*(*iops).Jc;
 for(j=4;j<=gjmbH.Hi();++j){Ham-=gjmbH(j)*(*(*iops).Olm[j-3]);}

/*   int i1,j1; //printout matrix
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    real(((*(*iops).OOlm[5])-(*iops).Jcc*(*iops).Jcc+(*iops).Jaa*(*iops).Jaa)(i1,j1)));}
//    real(((*iops).Jcc*(*iops).Jaa+(*iops).Jaa*(*iops).Jcc)(i1,j1)));}
//    real((*(*iops).OOlm[1])(i1,j1)));}
    printf ("\n");
    }
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    imag(((*(*iops).OOlm[5])-(*iops).Jcc*(*iops).Jcc+(*iops).Jaa*(*iops).Jaa)(i1,j1)));}
//   imag(((*iops).Jcc*(*iops).Jaa+(*iops).Jaa*(*iops).Jcc)(i1,j1)));}
//   imag((*(*iops).OOlm[1])(i1,j1)));}
    printf ("\n");
    }
exit(0);      
*/    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=1;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
   
   
   
   
   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);double Zs;
     double x,y;int i,k,l,m;
     x=Min(En);
     for (i=1;i<=dj;++i)
     {if ((y=(En(i)-x)/K_B/T)<700) wn[i]=exp(-y); 
      else wn[i]=0.0;
//      printf("%g\n",En(i));
      }
     Zs=Sum(wn);wn/=Zs;  
     Zs*=exp(-x/K_B/T);
   // calculate Ja,Jb,Jc
     ComplexMatrix z(1,dj,1,dj);
     ComplexMatrix * zp[gjmbH.Hi()+1];
     for(l=1;l<=gjmbH.Hi();++l)
      {zp[l]= new ComplexMatrix(1,dj,1,dj);}
     z=ComplexMatrix(zr,zi);
     
     (*zp[1])=(*iops).Jaa*z;
     (*zp[2])=(*iops).Jbb*z;
     (*zp[3])=(*iops).Jcc*z;

     
 for(j=4;j<=gjmbH.Hi();++j)
    {(*zp[j])=(*(*iops).OOlm[j-3])*z;}
     
// calculate mat and delta for transition number tn
// 1. get i and j from tn
k=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=En(j)-En(i);

if (delta<-0.000001){fprintf(stderr,"ERROR module cfield.so - dmcalc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

// 3. set mat
for(l=1;l<=gjmbH.Hi();++l)for(m=1;m<=gjmbH.Hi();++m)
{if(i==j){//take into account thermal expectation values <Jl>
          mat(l,m)=((z.Column(i)*(*zp[l]).Column(j))-J(l))*((z.Column(j)*(*zp[m]).Column(i))-J(m));}
 else    {mat(l,m)=(z.Column(i)*(*zp[l]).Column(j))*(z.Column(j)*(*zp[m]).Column(i));}}



if (delta/K_B/T>0.000001)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Ja|%i>|^2=%4.4g |<%i|Jb|%i>|^2=%4.4g |<%i|Jc|%i>|^2=%4.4g",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
      printf(" n%i-n%i=%4.4g\n",i,j,wn(i)-wn(j));}
    mat*=(wn(i)-wn(j)); // occupation factor    
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Ja-<Ja>|%i>|^2=%4.4g |<%i|Jb-<Jb>|%i>|^2=%4.4g |<%i|Jc-<Jc>|%i>|^2=%4.4g",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
      printf(" n%i=%4.4g\n",i,wn(i));}
    mat*=(wn(i)/K_B/T);
   }

//clean up memory
     for(l=1;l<=gjmbH.Hi();++l)
      {delete zp[l];}
     
// return number of all transitions     
 return (int)(((*iops).J+1)*(2*(*iops).J+1)); 
}
