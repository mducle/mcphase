// methods for class inimcdis 
#include "inimcdis.hpp"
#include <martin.h>
#include <cstring>
#include "../../version"

#define MAXNOFCHARINLINE 1024


 // *************************************************************************
 // ************************ inipar *************************************
 // *************************************************************************
 // class of initial parameters for program mcphas

void inimcdis::errexit() // type info and error exit 
{     printf (" \n			%s \n",MCDISPVERSION);
    printf ("    		    use as: mcdisp\n"); 
    printf ("		 or as: mcdisp [file]\n");
    printf ("                               [file] ... input file with mean field set (default mcdisp.mf)\n");
    printf ("			       \n");
    printf ("	       Note: files which must be in current directory -\n");
    printf ("	             ./mcdisp.par, ./mcphas.j, directory ./results\n");
    printf ("\n");
    printf (" Options: -jq                  ... calculate J(Q) (Fourier transform of exchange)\n");
    printf ("          -max n               ... restrict single ion susceptibility to n lowest\n");
    printf ("		                        lying transitions starting from the ground state\n");
    printf ("	       -minE E              ... an energy range may be given by minE and maxE: only\n");
    printf ("	       -maxE E                  single ion transitions within this energy range will \n");
    printf ("		                        be considered\n");
    printf ("	       -r                   ... refine energies\n");
    printf ("	       -v                   ... verbose\n");
    printf ("	       -a                   ... do not overwrite output files in results - append new results\n");
    printf ("	       -c                   ... only create single ion transition file ./results/mcdisp.trs and exit\n");
    printf ("	       -t                   ... read single ion transition file ./results/mcdisp.trs (do not create it)\n");
      exit (EXIT_FAILURE);
} 
 // *************************************************************************

void inimcdis::save()
{  FILE * fout;int i,j;
  fout=fopen(savfilename,"w");if (fout==NULL) {fprintf(stderr,"ERROR - file %s cannot be opened \n",savfilename);errexit();} 
  fprintf(fout,"# Parameter file  mcdisp.par - read by %s\n",MCDISPVERSION);
  fprintf(fout,"#<!--mcdisp.mcdisp.par>\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
  fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# mcdisp calculates the neutron scattering cross section dsigma/dOmegadE' [barn/sr/meV/f.u.]\n");
  fprintf(fout,"#           f.u.=crystallogrpaphic unit cell (r1xr2xr3) for inelastic and diffuse scattering\n");
  fprintf(fout,"#\n");
  fprintf(fout,"# depending on what is kept constant it follows either kf or ki (1/A)\n");
  if(kf!=0){fprintf(fout,"#!kf=%g\n",kf);}else{fprintf(fout,"#!ki=%g\n",ki);}
  fprintf(fout,"# \n");
  fprintf(fout,"# for full calculation of the dynamical susceptibility (option \"-r\", inversion of the MF-RPA equation \n");
  fprintf(fout,"# for each point in Q-omega space) the minimum and maximum energy has to be given (energy stepwidth is \n");
  fprintf(fout,"# equal to the parameter epsilon given in the command line after \"-r\")\n");

  fprintf(fout,"#\n");
  fprintf(fout,"#!emin=%g\n",emin);

  fprintf(fout,"#!emax=%g\n",emax);

  fprintf(fout,"#\n");
  fprintf(fout,"# optional parameter is extended_eigenvector_dimension\n");
  fprintf(fout,"# which is used to define, how many components of the\n");
  fprintf(fout,"# eigenvector should be in the ouput to file mcdisp.qee\n");
  fprintf(fout,"# (important for charge density movies)\n");
  fprintf(fout,"#!extended_eigenvector_dimension=%i\n",extended_eigenvector_dimension);

  fprintf(fout,"#\n");
  fprintf(fout,"# It follows either \n");
  fprintf(fout,"#\n");
  fprintf(fout,"# (i) a Q vector mesh to be mapped in the calculation\n");
  if(hkllist!=0){fprintf(fout,"#");}else{fprintf(fout,"#!");}
  fprintf(fout,"hmin=%g hmax=%g deltah=%g\n",qmin[1],qmax[1],deltaq[1]);

  if(hkllist!=0){fprintf(fout,"#");}else{fprintf(fout,"#!");}
  fprintf(fout,"kmin=%g kmax=%g deltak=%g\n",qmin[2],qmax[2],deltaq[2]);

  if(hkllist!=0){fprintf(fout,"#");}else{fprintf(fout,"#!");}
  fprintf(fout,"lmin=%g lmax=%g deltal=%g\n",qmin[3],qmax[3],deltaq[3]);

  fprintf(fout,"#\n");
  fprintf(fout,"# or (if no mesh is given the program takes)\n");
  fprintf(fout,"# (ii) a list of Q vectors with (optional) energies of observed excitations to be fitted\n");
  fprintf(fout,"# h k l [E1(meV) E2(meV) ...]\n");
  if(hkllist==0) {fprintf(fout,"0.45 1 1 0.523 0.745\n 0.75 1 1 0.523 0.745\n");}
  else { for (j=1;j<=(int)qmax(1);++j) 
  	      {for(i=1;i<=hkls[j][0];++i)
         	    {fprintf(fout,"%g ",hkls[j][i]);}
               fprintf(fout,"\n");
	      }
       }
  fprintf(fout,"\n");
  fclose (fout);
  fout=fopen("results/_mcdisp.mf","w");
  fprintf(fout,"# Parameter file  mcdisp.mf - read by %s\n",MCDISPVERSION);
  fprintf(fout,"#<!--mcdisp.mcdisp.mf>\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"# mcdisp - program to calculate the dispersion of magnetic excitations\n");
  fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
  fprintf(fout,"#*********************************************************************\n");
  fprintf(fout,"#'T'             temperature T\n");
  fprintf(fout,"#'Ha' 'Hb' 'Hc'  magnetic field\n");
  fprintf(fout,"#'n'             number of atoms in magnetic unit cell\n");
  fprintf(fout,"#'nofatoms'      number of atoms in primitive crystal unit cell\n");
  fprintf(fout,"#'nofcomponents' dimension of moment vector of a magnetic atoms\n");
  fprintf(fout,"T=%g Ha=%g Hb=%g Hc=%g n=%i nofatoms=%i nofcomponents=%i\n",T,Ha,Hb,Hc,mf.n(),nofatoms,nofcomponents);
  mf.print(fout);
  fclose (fout);

}

// *************************************************************************
//constructor ... load initial parameters from file
inimcdis::inimcdis (const char * file,const char * spinfile)
{ float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
  char instr[MAXNOFCHARINLINE];
  int i=0,j=0;
  FILE * fin;
  fin=fopen(spinfile,"rb");if (fin==NULL) {fprintf(stderr,"ERROR - file %s not found \n",spinfile);errexit();}
  instr[0]='#';  
  while(instr[strspn(instr," \t")]=='#'){fgets(instr,MAXNOFCHARINLINE,fin);}
  info= new char [strlen(instr)+1];
  extract(instr,"T",T); 
  extract(instr,"Ha",Ha);
  extract(instr,"Hb",Hb);
  extract(instr,"Hc",Hc);
  strcpy(info,instr);
  printf("#%s \n# reading mean field configuration mf=gj muB heff [meV]\n",instr);
   
  nofatoms=1;nofcomponents=3;
  extract(instr,"nofatoms",nofatoms); 
  extract(instr,"nofcomponents",nofcomponents); 
  mf=mfcf(1,1,1,nofatoms,nofcomponents); 

  if(mf.load(fin)==0)
   {fprintf(stderr,"ERROR loading mean field configuration\n");exit(EXIT_FAILURE);}
  
  savfilename= new char [strlen(file)+11];
  strcpy(savfilename,"results/_");
  strcpy(savfilename+9,file);
  FILE *fin_coq;
  errno = 0;
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  emin=0;emax=10;extended_eigenvector_dimension=nofcomponents;
  printf("reading file %s\n",file);
  fin_coq = fopen(file, "rb"); if (fin_coq==NULL) {fprintf(stderr,"ERROR - file %s not found \n",file);errexit();}   
  // save the parameters read from mcdisp.par into results/mcdisp.par)
  ki=0;kf=0;
  qmin=0;qmax=0;deltaq=0;
  while (fgets(instr,MAXNOFCHARINLINE,fin_coq)!=NULL)
  {++i;
     extract(instr,"hmin",qmin[1]); 
     extract(instr,"kmin",qmin[2]); 
     extract(instr,"lmin",qmin[3]); 
     extract(instr,"hmax",qmax[1]); 
     extract(instr,"kmax",qmax[2]); 
     extract(instr,"lmax",qmax[3]); 
     extract(instr,"deltah",deltaq[1]); 
     extract(instr,"deltak",deltaq[2]); 
     extract(instr,"deltal",deltaq[3]); 
     extract(instr,"emin",emin); 
     extract(instr,"emax",emax); 
     extract(instr,"ki",ki); 
     extract(instr,"kf",kf); 
     extract(instr,"extended_eigenvector_dimension",extended_eigenvector_dimension);
  }
  if (ki==0) {if (kf==0) kf=10;
              fprintf(stdout,"#Calculating intensities for  kf=const=%4.4g/A\n",kf);
	     }
	     else
	     {kf=0;
	      fprintf(stdout,"#Calculating intensities for ki=const=%4.4g/A\n",ki);
	     }
  fclose (fin_coq);
  if (Norm(qmin)+Norm(qmax)+Norm(deltaq)==0)
    {  hkllist=1;
       hkls=new double *[i+10];
       deltaq(1)=1.0;qmin(1)=1.0;
       deltaq(2)=1000.0;deltaq(3)=1000.0;
       fin_coq = fopen(file, "rb");
       while (feof(fin_coq)==0)
          {if ((i=inputline(fin_coq,nn))>=3)
	      {++qmax(1);
	       hkls[(int)qmax(1)]=new double [i+1];
               hkls[(int)qmax(1)][0]=i;    
               for(j=1;j<=i;++j)
         	    {hkls[(int)qmax(1)][j]=nn[j];
	            }
	      }
	  }
     }
   else
     {hkllist=0;}
  save();
}

//kopier-konstruktor 
inimcdis::inimcdis (const inimcdis & p)
{ savfilename= new char [strlen(p.savfilename)+1];
  strcpy(savfilename,p.savfilename);
 info= new char [strlen(p.info)+1];
  strcpy(info,p.info);
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  qmin=p.qmin;
  qmax=p.qmax;
  emin=p.emin;
  emax=p.emax;
  kf=p.kf;
  ki=p.ki;
  extended_eigenvector_dimension=p.extended_eigenvector_dimension;
  deltaq=p.deltaq;  
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  mf=mfcf(1,1,1,nofatoms,nofcomponents);mf=p.mf;T=p.T;
  hkllist=p.hkllist;
  int i,j;
  if (hkllist==1)
   {
      hkls=new double *[(int)qmax(1)+10];
      for (j=1;j<=(int)qmax(1);++j) 
  	      {
	       hkls[j]=new double [(int)p.hkls[j][0]+1];
               for(i=0;i<=p.hkls[j][0];++i)
         	    {hkls[j][i]=p.hkls[j][i];}
	      }
   }
}

//destruktor
inimcdis::~inimcdis ()
{delete []savfilename;delete []info;
 int i;
 if (hkllist==1)
 { for (i=1;i<=(int)qmax(1);++i) 
   { delete []hkls[i];}
   delete []hkls; 
 }
}
