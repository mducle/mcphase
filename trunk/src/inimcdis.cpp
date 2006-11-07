// methods for class inimcdis 
#include "inimcdis.hpp"
#include <martin.h>
#include <cstring>
#include "../version"

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
    printf ("	             ./mcdisp.ini, ./mcphas.j, directory ./results\n");
    printf ("\n");
    printf (" Options: -jq                  ... calculate J(Q) (Fourier transform of exchange)\n");
    printf ("          -max n               ... restrict single ion susceptibility to n lowest\n");
    printf ("		                        lying transitions starting from the ground state\n");
    printf ("	       -minE E              ... an energy range may be given by minE and maxE: only\n");
    printf ("	       -maxE E                  single ion transitions within this energy range will \n");
    printf ("		                        be considered\n");
    printf ("	       -r                   ... refine energies\n");
    printf ("	       -v                   ... verbose\n");
    printf ("	       -c                   ... only create single ion transition file ./results/mcdisp.trs and exit\n");
    printf ("	       -t                   ... read single ion transition file ./results/mcdisp.trs (do not create it)\n");
      exit (EXIT_FAILURE);
} 

//constructor ... load initial parameters from file
inimcdis::inimcdis (char * file,char * spinfile)
{ float nn[MAXNOFCHARINLINE];nn[0]=MAXNOFCHARINLINE;
  char instr[MAXNOFCHARINLINE];
  int i,j;
  FILE * fin;
  fin=fopen(spinfile,"rb");if (fin==NULL) {fprintf(stderr,"ERROR - file %s not found \n",spinfile);errexit();}
  fgets(instr,MAXNOFCHARINLINE,fin);
  info= new char [strlen(instr)+1];
  strcpy(info,instr);
  printf("%s \n reading mean field configuration mf=gj muB heff [meV]\n",instr);
  extract(instr,"T",T); 
  extract(instr,"Ha",Ha);
  extract(instr,"Hb",Hb);
  extract(instr,"Hc",Hc);
   
  nofatoms=1;nofcomponents=3;
  extract(instr,"nofatoms",nofatoms); 
  extract(instr,"nofcomponents",nofcomponents); 
  mf=mfcf(1,1,1,nofatoms,nofcomponents); 

  if(mf.load(fin)==0)
   {fprintf(stderr,"ERROR loading mean field configuration\n");exit(EXIT_FAILURE);}
  mf.print(stdout);
  
  savfilename= new char [strlen(file)+1];
  strcpy(savfilename,file);
  FILE *fin_coq;
  errno = 0;
  qmin=Vector(1,3);qmax=Vector(1,3);deltaq=Vector(1,3);
  emin=0;emax=10;
  printf("reading file %s\n",file);
  fin_coq = fopen(file, "rb"); if (fin_coq==NULL) {fprintf(stderr,"ERROR - file %s not found \n",file);errexit();}
  ki=0;kf=0;
  qmin=0;qmax=0;deltaq=0;
  while (fgets(instr,MAXNOFCHARINLINE,fin_coq)!=NULL)
  {++i;
   if(instr[strspn(instr," \t")]!='#')
   { extract(instr,"hmin",qmin[1]); 
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
   }
  }
  if (ki==0) {if (kf==0) kf=10;
              fprintf(stdout,"Calculating intensities for  kf=const=%4.4g/A\n",kf);
	     }
	     else
	     {kf=0;
	      fprintf(stdout,"Calculating intensities for ki=const=%4.4g/A\n",ki);
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
