// methods for class parameters 
#include "physprop.hpp"
#include "../../version"
 // *************************************************************************
 // ************************ physproperties *********************************
 // *************************************************************************


void   sort(float * v,int jmin,int jmax,int * jnew) // sorting function
{// input: v[jmin ... jmax] vector to be sorted
 // jnew [jmin .... jmax]   sort index - after sorting the 
 //                         vector jnew contains the indices
 //                         in a sequence such that v[jnew[jmin ... jmax]] is
 //                         sorted in ascending order 
 int gap,i,j,temp,n;
 n=jmax-jmin+1;
 // initialize jnew
 for(i=jmin;i<=jmax;++i){jnew[i]=i;}
 // if nothing is to be done - return
 if (jmax<=jmin) {return;}
 
 for (gap=n/2;gap>0;gap/=2)
   for (i=gap;i<n;++i)
     for (j=i-gap;j>=0 && v[jnew[jmin+j]]>v[jnew[jmin+j+gap]];j-=gap)
      {temp=jnew[jmin+j];jnew[jmin+j]=jnew[jmin+j+gap]; jnew[jmin+j+gap]=temp;}
      
 return;
}
 
//constructor
physproperties::physproperties (int nofspincorrs,int maxnofhkli,int na,int nm)
{washere=0;
 int i;
 nofspincorr=nofspincorrs;
 nofatoms=na;
 nofcomponents=nm;
 
 m=Vector(1,3);
 H=Vector(1,3);
 jj= new Vector [nofspincorrs+1];for(i=0;i<=nofspincorrs;++i){jj[i]=Vector(1,nofcomponents*nofcomponents*nofatoms);} //  ... number of interaction constants (aa bb cc ab ba ac ca bc cb)
   if (jj == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 hkli= new Vector [maxnofhkli+1];for(i=0;i<=maxnofhkli;++i){hkli[i]=Vector(1,7);}
   if (hkli == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 nofhkls=0;
 sps=spincf(1,1,1,nofatoms,nofcomponents);
 mf=mfcf(1,1,1,nofatoms,nofcomponents);
 fe=0;u=0;
}


//kopier-konstruktor
physproperties::physproperties (const physproperties & p)
{int i;
  x=p.x;y=p.y;
  j=p.j;
  T=p.T;
  H=p.H;
  m=p.m;
  nofhkls=p.nofhkls;
  u=p.u;fe=p.fe;
  sps=p.sps;
  mf=p.mf;
  washere=p.washere;
 maxnofhkls=p.maxnofhkls;
 nofspincorr=p.nofspincorr;
 nofatoms=p.nofatoms;
 nofcomponents=p.nofcomponents;
  
 jj= new Vector [nofspincorr+1];for(i=0;i<=nofspincorr;++i){jj[i]=Vector(1,nofcomponents*nofcomponents*nofatoms);} //  ... number of interaction constants (aa bb cc ab ba ac ca bc cb)
   if (jj == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 hkli= new Vector [maxnofhkls+1];for(i=0;i<=maxnofhkls;++i){hkli[i]=Vector(1,7);}
   if (hkli == NULL){fprintf (stderr, "physproperties::physproperties Out of memory\n");exit (EXIT_FAILURE);} 
 for(i=1;i<=nofspincorr;++i)
    {jj[i]=p.jj[i];}
 for(i=1;i<=nofhkls;++i)
    {hkli[i]=p.hkli[i];} 
 
 }


//destruktor
physproperties::~physproperties ()
{delete []jj;delete []hkli;}

void physproperties::update_maxnofhkls(int maxnofhkli)
{delete []hkli;
 maxnofhkls=maxnofhkli;
 int i;
 hkli= new Vector [maxnofhkli+1];for(i=0;i<=maxnofhkli;++i){hkli[i]=Vector(1,7);}
   if (hkli == NULL){fprintf (stderr, "physproperties::update_maxnofhkls - Out of memory\n");exit (EXIT_FAILURE);} 
}


// methode save
double physproperties::save (int verbose, int htfailed, par & inputpars)
{ FILE *fout;
  char filename[50];
  time_t curtime;
  struct tm *loctime;  
  int i,j2,l,i1,j1,nmax;
  Vector null(1,nofcomponents*nofatoms);null=0;
  Vector null1(1,3);null1=0;
  double sta;
  sta=0; 
  float nn[200];nn[0]=199;

  printf("saving properties for T=%g K  Ha= %g Hb= %g Hc= %g T\n", T,H(1),H(2),H(3));

//-----------------------------------------------------------------------------------------  
  errno = 0;
  if (verbose==1) printf("saving mcphas.fum\n");
  if (washere==0)
  {fout = fopen_errchk ("./results/mcphas.fum","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   fprintf (fout, "#1mev/ion=96.48J/mol\n");
   fprintf (fout, "#note: moments and energies are given per ion - not per formula unit !\n");
   fprintf (fout, "#   x    y   T[K] H[T] Ha[T] Hb[T] Hc[T] free energy f[meV/ion] energy u[meV/ion] total moment m     ma mb mc[mb/ion]}\n");
   fclose(fout);
      }
   if (htfailed!=0){fe=0;u=0;m=0;m[1]=0;m[2]=0;m[3]=0;}
   fout = fopen_errchk ("./results/mcphas.fum","a");
   fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g       %8.8g            %8.8g       %4.4g    %4.4g %4.4g %4.4g\n",
            x,y,T,Norm(H),H[1],H[2],H[3],fe,u,Norm(m),m[1],m[2],m[3]);
   fclose(fout);
   fout = fopen_errchk ("./results/.mcphas.fum","a");
   fprintf (fout, "%4.4g %4.4g  %4.4g %4.4g %4.4g %4.4g %4.4g %8.8g %8.8g  %4.4g %4.4g %4.4g %4.4g\n",
            x,y,T,Norm(H),H[1],H[2],H[3],fe,u,Norm(m),m[1],m[2],m[3]);
   fclose(fout);
   if((fout=fopen("./fit/mcphas.fum","rb"))!=NULL)
    {// some measured data should be fitted
     if(washere==0){fprintf(stdout,"#Mcphas- calculating standard deviation to ./fit/mcphas.fum - magnetisation ma mb mc (col 11,12,13)\n");}
     while(feof(fout)==0)
     {
      if (inputline(fout,nn)!=0)
      {if(0.0001>(nn[3]-T)*(nn[3]-T)+(nn[5]-H[1])*(nn[5]-H[1])+(nn[6]-H[2])*(nn[6]-H[2])+(nn[7]-H[3])*(nn[7]-H[3]))
         {sta+=(nn[11]-m[1])*(nn[11]-m[1])+(nn[12]-m[2])*(nn[12]-m[2])+(nn[13]-m[3])*(nn[13]-m[3]);
	 if(verbose==1){fprintf(stdout,"sta_mcphas.fum=%g\n",sta);}
	 }
      }
     }
     fclose(fout);
    }else{errno=0;}
//-----------------------------------------------------------------------------------------  
  errno = 0; Vector totalJ(1,nofcomponents);
  if (verbose==1)printf("saving mcphas.xyt\n");
  if (washere==0)
  {fout = fopen_errchk ("./results/mcphas.xyt","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   fprintf (fout, "#x    y   T[K] H[T] Ha[T] Hb[T] Hc[T] phasnumber-j   period-key ");
           for(i1=1;i1<=nofcomponents;++i1)
	      {fprintf(fout,"<J%c> ",'a'-1+i1);}
	      fprintf(fout,"\n");
   fclose(fout);    
     }
  fout = fopen_errchk ("./results/mcphas.xyt","a");
  totalJ=0; 
  if (htfailed!=0){j=0;}else{totalJ=sps.totalJ();}
   if(j<0){sps.wasstable=j;}
   fprintf (fout, "%4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g       %ip           %ip      ",
            x,y,T,Norm(H),H[1],H[2],H[3],j,sps.wasstable);
           for(i1=1;i1<=nofcomponents;++i1)
	      {fprintf(fout,"%4.4g ",totalJ(i1));}
	      fprintf(fout,"\n");
   fclose(fout);
    if((fout=fopen("./fit/mcphas.xyt","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.xyt not implemented\n");}
     fclose(fout);
    }else{errno=0;}

//-----------------------------------------------------------------------------------------  
 nmax=nofspincorr; // look how many spincorrelationfunction we have indeed calculated - the 
                   // user wanted nofspincorr, but maybe it was fewer ...
 for (l=1;l<=inputpars.nofatoms;++l)
      {if(nmax>(*inputpars.jjj[l]).paranz)
       {nmax=(*inputpars.jjj[l]).paranz;
        fprintf(stderr,"Warning: calculation of nofspincorr=%i correlation functions not possible, \n",nofspincorr);
fprintf(stderr,"         because in mcphas.j for atom %i  only %i neighbours are given.\n",l,(*inputpars.jjj[l]).paranz);
       }
      }
       

  // only output nmax correlation functions ...
 for(i=1;i<=nmax;++i){for(l=1;l<=nofatoms;++l){
  errno = 0;
  if (verbose==1)printf("saving mcphas%i.j%i - spinspin corr for sublattice %i neighbour %i\n",l,i,l,i);
  sprintf(filename,"./results/mcphas%i.j%i",l,i);
  if (washere==0)  //printout file header
  {  fout = fopen_errchk (filename,"w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   fprintf (fout, "# sublattice %i (da=%g a db=%g b dc=%g c)\n",l,(*inputpars.jjj[l]).xyz(1),(*inputpars.jjj[l]).xyz(2),(*inputpars.jjj[l]).xyz(3));
   fprintf (fout, "# correlation fuction <JJ(%g %g %g)>\n",(*inputpars.jjj[l]).dn[i](1),(*inputpars.jjj[l]).dn[i](2),(*inputpars.jjj[l]).dn[i](3));
   fprintf (fout, "#x     y     T[K]  H[T]   Ha[T] Hb[T] Hc[T]  ");
           for(i1=1;i1<=(*inputpars.jjj[l]).nofcomponents;++i1)
	      {fprintf(fout,"<J%cJ%c> ",'a'-1+i1,'a'-1+i1);}
           for(i1=1;i1<=(*inputpars.jjj[l]).nofcomponents-1;++i1)
              {for(j1=i1+1;j1<=(*inputpars.jjj[l]).nofcomponents;++j1)
                        {fprintf(fout,"<J%cJ%c> <J%cJ%c> ",'a'-1+i1,'a'-1+j1,'a'-1+j1,'a'-1+i1);
			}
	      }
                          
   fprintf (fout,"}\n");
   fclose(fout);
      }
  fout = fopen_errchk (filename,"a");
  if (htfailed!=0){jj[i](1)=0;jj[i](2)=0;jj[i](3)=0;}
   fprintf (fout, "%4.4g %4.4g   %4.4g %4.4g   %4.4g %4.4g %4.4g     ",x,y,T,Norm(H),H[1],H[2],H[3]);
        for(j2=1;j2<=nofcomponents*nofcomponents;++j2)               
            {fprintf (fout, "%4.4g ",jj[i](j2+nofcomponents*nofcomponents*(l-1)));
	    }
   fprintf (fout,"\n");
   fclose(fout);
   sprintf(filename,"./fit/mcphas%i.j%i",l,i);
   if((fout=fopen(filename,"rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using %s  not implemented\n",filename);}
     fclose(fout);
    }else{errno=0;}
  }}

//-----------------------------------------------------------------------------------------  
 errno = 0;
  if (verbose==1)printf("saving mcphas*.hkl - neutrons and xrays\n");
  if (washere==0)
  {//neutrons
   fout = fopen_errchk ("./results/mcphas.hkl","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   fprintf (fout, "#Neutron Intensity - Mind: only structure+polarizationfactor+formfactor+debeywallerfactor - no lorentzfactor is  taken into account\n");
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  int       h   k   l   int       h   k   l   int ...}\n");
   fclose(fout);
   //xray a component
   fout = fopen_errchk ("./results/mcphasa.hkl","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - a component\n"); 
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  FT||a       h   k   l   FT||a       h   k   l   FT||a ...}\n");
   fclose(fout);
   //xray b component
   fout = fopen_errchk ("./results/mcphasb.hkl","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - b component\n"); 
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  FT||b       h   k   l   FT||b       h   k   l   FT||b ...}\n");
   fclose(fout);
   //xray a component
   fout = fopen_errchk ("./results/mcphasc.hkl","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   fprintf (fout,"#Absolute Value of the Fourier Transform of the moment configuration - c component\n"); 
   fprintf (fout, "#x   y   T[K]  H[T]  Ha[T] Hb[T] Hc[T]       h   k   l  FT||c       h   k   l   FT||c       h   k   l   FT||c ...}\n");
   fclose(fout);

      }
   int inew[nofhkls];float intensity[nofhkls];
   for (i=1;i<=nofhkls;++i) {intensity[i]=hkli[i](4);}
   sort(intensity,1,nofhkls,inew); // sort according to ascending intensity

  //neutrons
  fout = fopen_errchk ("./results/mcphas.hkl","a");fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",x,y,T,Norm(H),H[1],H[2],H[3]);
   for (i=nofhkls;i>=1;--i)
    {if (htfailed!=0){hkli[inew[i]](1)=0;hkli[inew[i]](2)=0;hkli[inew[i]](3)=0;hkli[inew[i]](4)=0;}
    fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g     ",hkli[inew[i]](1),hkli[inew[i]](2),hkli[inew[i]](3),hkli[inew[i]](4));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray a component
  fout = fopen_errchk ("./results/mcphasa.hkl","a");fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",x,y,T,Norm(H),H[1],H[2],H[3]);
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g     ",hkli[inew[i]](1),hkli[inew[i]](2),hkli[inew[i]](3),hkli[inew[i]](5));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray b component
  fout = fopen_errchk ("./results/mcphasb.hkl","a");fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",x,y,T,Norm(H),H[1],H[2],H[3]);
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g     ",hkli[inew[i]](1),hkli[inew[i]](2),hkli[inew[i]](3),hkli[inew[i]](6));
    } fprintf(fout,"\n");
   fclose(fout);
  //xray c component
  fout = fopen_errchk ("./results/mcphasc.hkl","a");fprintf (fout, " %-4.4g %-4.4g %-4.4g %-4.4g  %-4.4g %-4.4g %-4.4g      ",x,y,T,Norm(H),H[1],H[2],H[3]);
   for (i=nofhkls;i>=1;--i)
    {fprintf (fout, "%4.4g %4.4g %4.4g  %4.4g     ",hkli[inew[i]](1),hkli[inew[i]](2),hkli[inew[i]](3),hkli[inew[i]](7));
    } fprintf(fout,"\n");
   fclose(fout);


    if((fout=fopen("./fit/mcphas.hkl","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.hkl not implemented\n");}
     fclose(fout);
    }else{errno=0;}

//-----------------------------------------------------------------------------------------  
 errno = 0;
  if (verbose==1)printf("saving mcphas.sps - spinconfiguration\n");
  if (washere==0)
  {  fout = fopen_errchk ("./results/mcphas.sps","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   // printout the lattice and atomic positions
   inputpars.savelattice(fout);inputpars.saveatoms(fout);
   fprintf (fout, "#x y T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms(in primitive basis) nofmomentum-components\n");
fprintf (fout, "    #<Ja(1)> <Ja(2)> .... selfconsistent Spinconfiguration  \n");
fprintf (fout, "    #<Jb(1)> <Jb(2)> .... UNITS:  multiply by gJ to get moment [muB]\n");
fprintf (fout, "    #<Jc(1)> <Jc(2)> ....}\n");
    fclose(fout);
   }  
  fout = fopen_errchk ("./results/mcphas.sps","a");
   fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i \n",
            x,y,T,Norm(H),H[1],H[2],H[3],sps.n()*sps.nofatoms,sps.nofatoms,sps.nofcomponents);
   if (htfailed!=0){sps.spinfromq(1,1,1,null1,null,null,null);}
    sps.print(fout);fprintf(fout,"\n");
   fclose(fout);
    if((fout=fopen("./fit/mcphas.sps","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.sps not implemented\n");}
     fclose(fout);
    }else{errno=0;}
 
//-----------------------------------------------------------------------------------------  
 errno = 0;
  if (verbose==1)printf("saving mcphas.mf - mean field configuration\n");
  if (washere==0)
  {  fout = fopen_errchk ("./results/mcphas.mf","w");
   fprintf (fout, "#{%s ",MCPHASVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);//fprintf(fout,"\n");
   // printout the lattice and atomic positions
   inputpars.rems[2]=" ";
   inputpars.savelattice(fout);inputpars.saveatoms(fout);
   fprintf (fout, "#x y T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms(in primitive basis) nofmeanfield-components\n");
   fprintf (fout, "    #mfa(1) mfa(2) .... selfconsistent Mean field configuration \n"); 
   fprintf (fout, "    #mfb(1) mfb(2) .... UNITS: mf(i)=gJ*mu_B*heff(i)[meV] \n"); 
   fprintf (fout, "    #mfc(1) mfc(2) ....         (i.e. divide by gJ and mu_B=0.05788meV/T to get effective field[T]}\n");
    fclose(fout);
   }  
     fout = fopen_errchk ("./results/mcphas.mf","a");
fprintf (fout, " %4.4g %4.4g %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i\n",
            x,y,T,Norm(H),H[1],H[2],H[3],mf.n()*mf.nofatoms,mf.nofatoms,mf.nofcomponents);
   if (htfailed!=0){sps.print(fout);fprintf(fout,"\n");}else
    {mf.print(fout);fprintf(fout,"\n");}
   fclose(fout);
    if((fout=fopen("./fit/mcphas.mf","rb"))!=NULL)
    {// some measured data should be fitted
     if (washere==0){fprintf(stderr,"Warning: Calculation of standard deviation using  ./fit/mcphas.mf not implemented\n");}
     fclose(fout);
    }else{errno=0;}
//-----------------------------------------------------------------------------------------  
 
 washere=1;
return sta;
 }
