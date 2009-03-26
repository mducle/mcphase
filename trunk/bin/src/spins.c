/***********************************************************************
 *
 * spins.c - program to display spinconfiguration at given htpoint
 *
 ***********************************************************************/

#define MAXNOFCHARINLINE 1000
#define MAXNOFATOMS 100
#include "spincf.hpp"
#include "martin.h"
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>

/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ spincf savspins;
 FILE * fin_coq, * fout;
 float delta,dd,ddT,ddHa,ddHb,ddHc,alpha,beta,gamma,T,Ha,Hb,Hc;
 int i,n=0,nofatoms=0,nofcomponents=3;
 long int pos=0,j;
 float numbers[13];numbers[9]=1;numbers[10]=3;
 numbers[0]=13;
 char instr[MAXNOFCHARINLINE];
 char outstr[MAXNOFCHARINLINE];
 float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS],gJ[MAXNOFATOMS];
 char * cffilenames[MAXNOFATOMS];
  Matrix r(1,3,1,3);
  Vector abc(1,3);
// check command line
  if (argc < 5)
    { printf (" program spins - display spincf at HT point\n \
                use as: spins T Ha Hb Hc [file.sps||file.mf]\n");
      exit (1);
    }

 if (argc <6) 
 { fin_coq = fopen_errchk ("./results/mcphas.sps", "rb");}
 else
 { fin_coq = fopen_errchk (argv[5], "rb");}
    
 fout = fopen_errchk ("./results/spins.out", "w");



abc=0;
 // input file header ------------------------------------------------------------------
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
  { pos=ftell(fin_coq); 
   if (pos==-1) 
       {fprintf(stderr,"Error: wrong sps file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq); 
   // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
   for(i=0;i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;} 
   
   if (instr[strspn(instr," \t")]=='#'){fprintf(fout,instr);}
   if(abc[1]==0){extract(instr,"a",abc[1]);extract(instr,"b",abc[2]); extract(instr,"c",abc[3]); 
                 extract(instr,"alpha",alpha);  extract(instr,"beta",beta);extract(instr,"gamma",gamma); 
   }
   extract(instr,"r1x",r[1][1]);extract(instr,"r2x",r[1][2]); extract(instr,"r3x",r[1][3]); 
   extract(instr,"r1y",r[2][1]); extract(instr,"r2y",r[2][2]); extract(instr,"r3y",r[2][3]);
   extract(instr,"r1z",r[3][1]); extract(instr,"r2z",r[3][2]); extract(instr,"r3z",r[3][3]);
   extract(instr,"r1a",r[1][1]);extract(instr,"r2a",r[1][2]); extract(instr,"r3a",r[1][3]); 
   extract(instr,"r1b",r[2][1]); extract(instr,"r2b",r[2][2]); extract(instr,"r3b",r[2][3]);
   extract(instr,"r1c",r[3][1]); extract(instr,"r2c",r[3][2]); extract(instr,"r3c",r[3][3]);
   extract(instr,"nofatoms",nofatoms);    extract(instr,"nofcomponents",nofcomponents); 
   if (nofatoms>0&&(extract(instr,"x",x[n+1])+
                   extract(instr,"y",y[n+1])+
  		       extract(instr,"z",z[n+1])==0)||
		       (extract(instr,"da",x[n+1])+
                   extract(instr,"db",y[n+1])+
		       extract(instr,"dc",z[n+1])==0))
		  {++n;if(n>nofatoms||nofatoms>MAXNOFATOMS)
                    {fprintf(stderr,"ERROR spins.c reading file:maximum number of atoms in unit cell exceeded\n");exit(EXIT_FAILURE);}
                   cffilenames[n]=new char[MAXNOFCHARINLINE];
                   extract(instr,"cffilename",cffilenames[n],(size_t)MAXNOFCHARINLINE);
                   extract(instr,"gJ",gJ[n]);
//		   printf("%s\n",cffilenames[n]);
		  }
  }
  if (alpha!=90||beta!=90||gamma!=90)
  {fprintf(stderr,"ERROR: non orthogonal lattice not supported yet\n");exit(EXIT_FAILURE);}
   Vector gJJ(1,n); for (i=1;i<=n;++i){gJJ(i)=gJ[i];}
  

// load spinsconfigurations and check which one is nearest -------------------------------   

   j=fseek(fin_coq,pos,SEEK_SET); 
    if (j!=0){fprintf(stderr,"Error: wrong sps file format\n");exit (EXIT_FAILURE);}
   
 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      spins.load(fin_coq);
      ddT=strtod(argv[1],NULL)-numbers[3];ddT*=ddT;
      ddHa=strtod(argv[2],NULL)-numbers[5];ddHa*=ddHa;
      ddHb=strtod(argv[3],NULL)-numbers[6];ddHb*=ddHb;
      ddHc=strtod(argv[4],NULL)-numbers[7];ddHc*=ddHc;
      dd=sqrt(ddT+ddHa+ddHb+ddHc+0.000001);
      if (dd<delta)
       {delta=dd;
        sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g n=%g spins nofatoms=%i in primitive basis nofcomponents=%i",numbers[3],numbers[5],numbers[6],numbers[7],numbers[8],(int)numbers[9],(int)numbers[10]);
        savspins=spins;T=numbers[3];Ha=numbers[5];Hb=numbers[6];Hc=numbers[7];
        
       }
    }
  fclose (fin_coq);

  printf("%s - momentum configuration <J(i)>\n",outstr);
  fprintf(fout,"#%s - momentum configuration <J(i)>\n",outstr);
  savspins.printall(fout,abc,r,x,y,z,cffilenames,gJ);
  savspins.print(stdout);
  fclose (fout);
  
// create plot of spinconfiguration -----------------------------------------------------------

    fin_coq = fopen_errchk ("./results/spins.eps", "w");
     savspins.eps(fin_coq,outstr);
    fclose (fin_coq);

// here the 3d file should be created
    fin_coq = fopen_errchk ("./results/spinsab.eps", "w");
     savspins.eps3d(fin_coq,outstr,abc,r,x,y,z,1,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spinsac.eps", "w");
     savspins.eps3d(fin_coq,outstr,abc,r,x,y,z,2,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spinsbc.eps", "w");
     savspins.eps3d(fin_coq,outstr,abc,r,x,y,z,3,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spins3dab.eps", "w");
     savspins.eps3d(fin_coq,outstr,abc,r,x,y,z,4,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spins3dac.eps", "w");
     savspins.eps3d(fin_coq,outstr,abc,r,x,y,z,5,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spins3dbc.eps", "w");
     savspins.eps3d(fin_coq,outstr,abc,r,x,y,z,6,gJJ);
    fclose (fin_coq);

    fin_coq = fopen_errchk ("./results/spins.fst", "w");
     savspins.fst(fin_coq,outstr,abc,r,x,y,z,gJJ);
    fclose (fin_coq);

    
   fin_coq = fopen_errchk ("./results/spins_prim.fst", "w");
     savspins.fstprim(fin_coq,outstr,abc,r,x,y,z,gJJ);
    fclose (fin_coq);

if (argc>=9){// try a spinwave picture
             double h,k,l,E,ddh,ddk,ddl,ddE;
             spincf savev_real(savspins);
             spincf savev_imag(savspins);
             fin_coq = fopen_errchk ("./results/mcdisp.qev", "rb");
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin_coq); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;} 
               // load evs and check which one is nearest -------------------------------   
              }
               j=fseek(fin_coq,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                 spincf ev_real(savspins),ev_imag(savspins);
                 ev_real.load(fin_coq);ev_imag.load(fin_coq);
                 ddT=strtod(argv[1],NULL)-numbers[4];ddT*=ddT;
                 ddHa=strtod(argv[2],NULL)-numbers[1];ddHa*=ddHa;
                 ddHb=strtod(argv[3],NULL)-numbers[2];ddHb*=ddHb;
                 ddHc=strtod(argv[4],NULL)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[6],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[7],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[8],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[9],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;
                  sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  savev_real=ev_real;savev_imag=ev_imag;h=numbers[5];k=numbers[6];l=numbers[7];E=numbers[9];        
                 }
               }
              fclose (fin_coq);
              fprintf(stdout,"#%s - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              savev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              savev_imag.print(stdout);
              spincf spins(savspins);
               sleep(10);                 

              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase=0,amplitude=1.0; // should be controlled
              complex <double> im(0,1);
              while(phase<10000)
              {spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
               fin_coq = fopen_errchk ("./results/spins_prim.fst", "w");
                spins.fstprim(fin_coq,outstr,abc,r,x,y,z,gJJ);
               fclose (fin_coq);
               sleep(1);                 
               phase+=0.1;
              }
          }

  
  for(i=1;i<=nofatoms;++i){  delete cffilenames[i];}
  return 0;
}


