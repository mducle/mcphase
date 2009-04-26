// routines for mcphas for calculation of magnetic phases
// htcalc.c
#include "myev.h"

void checkini(testspincf & testspins,qvectors & testqs)
{struct stat filestatus;
 static time_t last_modify_time;
 static int washere=0;
 int loaderr;
  errno = 0;

  if (stat(ini.savfilename,&filestatus)!=0)
    {fprintf (stderr, "Error checking mcphas.ini: Couldn't read status of file %s: %s\n",
              ini.savfilename, strerror (errno));exit (EXIT_FAILURE);
     }

  if(washere==0){washere=1;last_modify_time=filestatus.st_mtime;}
  
   
    if (filestatus.st_mtime!=last_modify_time) //check if file has been modified
    {again:
     last_modify_time=filestatus.st_mtime;
     fprintf(stdout,"mcphas.ini has been modified - reading new mcphas.ini\n");
      sleep(1);
      loaderr=ini.load();
      if(ini.exit_mcphas==1)
        {testspins.save(filemode);  //exit normally
         testqs.save(filemode);
         exit(0);
	}

      while(ini.pause_mcphas==1||loaderr==1) // wait until pause button is released and no loaderror occurs
       {fprintf(stdout,"Pausing ...\n");
        while(filestatus.st_mtime==last_modify_time)  //wait until filestatus changes again
          {sleep(1);
            if (stat(ini.savfilename,&filestatus)!=0)
               {fprintf (stderr, "Error checking file mcphas.ini: Couldn't read status of file %s: %s\n",
                ini.savfilename, strerror (errno));exit (EXIT_FAILURE);
               }
          }
	goto again;  
       }
       

    }
}


int  htcalc (Vector H,double T,par & inputpars,qvectors & testqs,
             testspincf & testspins, physproperties & physprops)
{/* calculates magnetic structure at a given HT- point  
  on input: 
    T	Temperature[K]
    H	Vector of External Magnetic Field [T]
    inputpars	Input parameters (exchange constants etc...)
    testqs	Set of propagation vectors to be tested 
    testspins	Set of Spinconfigurations to be tested
  on return:
    physprops	physical properties at (HT) point (i.e. magnetic structure
		neutron intensities, thermal expansion ...)	
 // returns 0 if successfull
 // returns 1 if too maxnofspinconfigurations is exceeded 
 // returns 2 if no spinconfiguration has been found at ht point
 */
 int i,ii,iii,j,k,tryrandom,nr,rr,ri,is;
 double fe,fered;
 double u,lnz; // free- and magnetic energy per ion [meV]
 Vector momentq0(1,inputpars.nofcomponents*inputpars.nofatoms),phi(1,inputpars.nofcomponents*inputpars.nofatoms),nettom(1,inputpars.nofcomponents*inputpars.nofatoms),q(1,3);
 Vector h1(1,inputpars.nofcomponents),hkl(1,3);
 double femin=10000;char text[100];
 spincf  sps(1,1,1,inputpars.nofatoms,inputpars.nofcomponents),sps1(1,1,1,inputpars.nofatoms,inputpars.nofcomponents);
 spincf  spsmin(1,1,1,inputpars.nofatoms,inputpars.nofcomponents);
 mfcf * mf;
 FILE * felog; // logfile for q dependence of fe
 FILE * fin_coq;

if (T<=0.01){fprintf(stderr," ERROR htcalc - temperature too low - please check mcphas.ini !");exit(EXIT_FAILURE);}

 srand(time(0)); // initialize random number generator
 checkini(testspins,testqs); // check if user pressed a button
 if (ini.logfevsQ==1) {felog=fopen_errchk("./results/mcphas.log","a");
               fprintf(felog,"#Logging of h k l fe[meV] spinconf_nr n1xn2xn3 T=%g Ha=%g Hb=%g Hc=%g\n",T,H(1),H(2),H(3));
               fclose(felog);
	      }
 if (verbose==1)
 { fin_coq= fopen_errchk ("./results/.fe_status.dat","w");
   fprintf(fin_coq,"time(s) log(iteration) log(sta) spinchange stepratio\n");
   fclose(fin_coq);	      
   printf("\n starting T=%g Ha=%g Hb=%g Hc=%g with \n %i spinconfigurations read from mcphas.tst and table \nand\n %i spinconfigurations created from hkl's\n\n",T,H(1),H(2),H(3),testspins.n,testqs.nofqs());
 }

 j=-testqs.nofqs()+(int)rint(rnd(testspins.n+testqs.nofqs())); 
    //begin with j a random number, j<0 means test spinconfigurations 
    //constructed from q vector set testqs, j>0 means test spinconfigurations from
    //set testspins
    //j=0;  //uncomment this for debugging purposes
    
 for (k= -testqs.nofqs();k<=testspins.n;++k)
 {++j; if (j>testspins.n) j=-testqs.nofqs();
  for (tryrandom=0;tryrandom<=ini.nofrndtries&&j!=0;++tryrandom)
   {if (j>0){sps=(*testspins.configurations[j]);// take test-spinconfiguration
	     if (tryrandom==0&&verbose==1) printf ( "conf. no %i (%ix%ix%i spins)"  ,j,sps.na(),sps.nb(),sps.nc());
            }
    else     // take q vector and choose phase and mom dir randomly
            {q=testqs.q(-j);  
	     if (tryrandom==0)
	     {nettom=testqs.nettom(-j);momentq0=testqs.momentq0(-j);phi=testqs.phi(-j);
	     }
	     else
	     {for(i=1;i<=inputpars.nofatoms;++i)
	      {for(ii=1;ii<=inputpars.nofcomponents;++ii)
	        {iii=inputpars.nofcomponents*(i-1)+ii;h1=0;h1(ii)=10*MU_B;
		 nettom(iii)=(*inputpars.jjj[i]).mcalc(T,h1,lnz,u,(*inputpars.jjj[i]).est)(ii)*rnd(1);
	         momentq0(iii)=rnd(1);
	         phi(iii)=rnd(1)*3.1415;
		}
	      }
	     }
	     sps.spinfromq(testqs.na(-j),testqs.nb(-j),testqs.nc(-j),
	                   q,nettom,momentq0,phi);
             hkl=inputpars.rez.Transpose()*q;  
   	     if (tryrandom==0&&verbose==1) printf ( "(hkl)=(%g %g %g)..(%ix%ix%i primitive unit cells) ",hkl(1),hkl(2),hkl(3),sps.na(),sps.nb(),sps.nc());
	    }	 
    if (tryrandom>0){nr=(int)(rint(rnd(1.0)*(sps.n()*inputpars.nofatoms-1)))+1;
	             for (i=1;i<=nr;++i) //MonteCarlo randomize nr spins
                      {rr=(int)rint(rnd(1.0)*(sps.n()-1))+1;
		       ri=inputpars.nofcomponents*(int)rint(rnd(1.0)*(inputpars.nofatoms-1));
	               for(ii=1;ii<=inputpars.nofcomponents;++ii)
		       {sps.mi(rr)(ri+ii)=rnd(1.0) ;}
		       } // randomize spin rr
                    }
     if (H*sps.nettomagmom(inputpars.gJ)<0) //see if nettomoment positiv
        {sps.invert();momentq0=-momentq0;} //if not - invert spinconfiguration

      //!!!calculate free eneregy - this is the heart of this loop !!!!
      mf=new mfcf(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,inputpars.nofcomponents);
      fe=fecalc(H ,T,inputpars,sps,(*mf),u,testspins,testqs);
      delete mf;
     
      // test spinconfiguration  and remember it                                    
      if (fe<femin)
            {               // first - reduce the spinconfiguration if possible
	       sps1=sps;sps1.reduce();
                   mf=new mfcf(sps1.na(),sps1.nb(),sps1.nc(),inputpars.nofatoms,inputpars.nofcomponents);
               if ((fered=fecalc(H ,T,inputpars,sps1,(*mf),u,testspins,testqs))<=fe+0.00001)sps=sps1;
                   delete mf;
	       spsmin=sps;	   
                 // display spinstructure
                if (verbose==1)
                {Vector abc(1,3);
		 abc(1)=inputpars.a;
		 abc(2)=inputpars.b;
		 abc(3)=inputpars.c;
		 float x[inputpars.nofatoms],y[inputpars.nofatoms],z[inputpars.nofatoms];
		 for (is=1;is<=inputpars.nofatoms;++is)
		   {x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     sprintf(text,"fe=%g,fered=%g<femin=%g:T=%gK, |H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT,  %i spins",fe,fered,femin,T,Norm(H),H(1),H(2),H(3),sps.n());
                    fin_coq = fopen_errchk ("./results/.spins3dab.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,4,inputpars.gJ);
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dac.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,5,inputpars.gJ);
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dbc.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,6,inputpars.gJ);
                    fclose (fin_coq);
		   
                    fin_coq = fopen_errchk ("./results/.spins.eps", "w");
                     sps.eps(fin_coq,text);
                    fclose (fin_coq);
		}

	           
                           // see if spinconfiguration is already stored
	     if (0==checkspincf(j,sps,testqs,nettom,momentq0,phi,testspins,physprops))//0 means error in checkspincf/addspincf
	        {fprintf(stderr,"Error htcalc: too many spinconfigurations created");
                 testspins.save(filemode);testqs.save(filemode); 
		 return 1;}
	     femin=fe;	  
            //printout fe
	    if (verbose==1) printf("fe=%gmeV, struc no %i in struct-table (initial values from struct %i)",fe,physprops.j,j);
	     }
            if (tryrandom==ini.nofrndtries&&verbose==1){printf("\n");}
	    
	    
  // log fe if required
   if (ini.logfevsQ==1) {
                 ComplexVector a(1,3*inputpars.nofatoms),b(1,3*inputpars.nofatoms);
                 ComplexVector b1(1,inputpars.nofcomponents*inputpars.nofatoms);
                 float inmax=0;int qh,qk,ql,l;
                 ComplexVector * mq;  
                 mq = new ComplexVector [sps.in(sps.na(),sps.nb(),sps.nc())+2];for(l=0;l<=sps.in(sps.na(),sps.nb(),sps.nc())+1;++l){mq[l]=ComplexVector(1,inputpars.nofcomponents*inputpars.nofatoms);}
                 Vector sq2(1,3*inputpars.nofatoms),qs(1,3),qt(1,3);float in;qs(1)=1000;
                 sps.FT(mq); //Fourier trafo of spincf
		 // get the main propagation vector by looking for the
		 // biggest Fourier component of the magnetic moment arrangement 
                 for(qh=0;qh<=sps.na()/2;++qh){for(qk=0;qk<=sps.nb()/2;++qk){for(ql=0;ql<=sps.nc()/2;++ql)
                  {// get magnetic moment from momentum fouriercomponent into b 
		   b=0;
		   b1 = mq[sps.in(sps.na()-qh,sps.nb()-qk,sps.nc()-ql)];
                   for(l=1;l<=inputpars.nofatoms;++l)
		   {int m1,m1max=3; if ((*inputpars.jjj[l]).gJ==0){m1max=6;}
		    for (m1=1;m1<=m1max;++m1)
		     {if((*inputpars.jjj[l]).gJ==0)
		      {if(m1==2||m1==4||m1==6){b(3*(l-1)+(m1+1)/2)+=b1(inputpars.nofcomponents*(l-1)+m1);}
		       else                   {b(3*(l-1)+(m1+1)/2)+=2.0*b1(inputpars.nofcomponents*(l-1)+m1);}
		      }
		      else
		      {b(3*(l-1)+m1)=b1(inputpars.nofcomponents*(l-1)+m1)*(*inputpars.jjj[l]).gJ;
		      }
		     }    
		    }
		   a = b.Conjugate();
		   b1 = mq[sps.in(qh,qk,ql)];
		   b=0;
                   for(l=1;l<=inputpars.nofatoms;++l)
		   {int m1,m1max=3; if ((*inputpars.jjj[l]).gJ==0){m1max=6;}
		    for (m1=1;m1<=m1max;++m1)
		     {if((*inputpars.jjj[l]).gJ==0)
		      {if(m1==2||m1==4||m1==6){b(3*(l-1)+(m1+1)/2)+=b1(inputpars.nofcomponents*(l-1)+m1);}
		       else                   {b(3*(l-1)+(m1+1)/2)+=2.0*b1(inputpars.nofcomponents*(l-1)+m1);}
		      }
		      else
		      {b(3*(l-1)+m1)=b1(inputpars.nofcomponents*(l-1)+m1)*(*inputpars.jjj[l]).gJ;
		      }
		     }    
		    }                   
		   // inner product
                   sq2=Abs(b+a)/(double)sps.n()/(double)inputpars.nofatoms;
                   Vector q(1,3);
		   q(1)=1.0*qh/sps.na();
	           q(2)=1.0*qk/sps.nb();
                   q(3)=1.0*ql/sps.nc();
                   qt=inputpars.rez.Transpose()*q;
		   in=Norm(sq2)*Norm(sq2);
	           if ((in>inmax-0.01&&Norm(qt)<Norm(qs))
		      ||in>inmax)
                    {inmax=in;qs=qt;}
                   }}}
                   felog=fopen_errchk("./results/mcphas.log","a");
                   if (verbose==1||fe>10000){fprintf(felog,"#");}
                   fprintf(felog,"%10.6g %10.6g %10.6g %10.6g %3i %3i %3i %3i \n",qs(1),qs(2),qs(3),fe,j,sps.na(),sps.nb(),sps.nc());
                   if (verbose==1&&fe<20000){sps.print(felog);}
	           fclose(felog);
                  delete []mq;
                 }

      }
  }

if (femin>=10000) // did we find a stable structure ??
 {fprintf(stderr,"Warning propcalc: femin positive ... no stable structure found at  T= %g K / H= %g T\n",
                 physprops.T,Norm(physprops.H));return 2;}
else // if yes ... then
 {// calculate physical properties
 if (physprops.j>0){ // take spinconfiguration ----
                     sps=(*testspins.configurations[physprops.j]);
                       if (sps.wasstable==0)
                       {// go through qvectors and spinfconfigurations and see if periodicity matches
                        for (i=1;i<=testqs.nofqs();++i)
                         {if (testqs.na(i)==sps.na()&&testqs.nb(i)==sps.nb()&&testqs.nc(i)==sps.nc())
                             {sps.wasstable=-i;break;}
                         }
		        if (sps.wasstable==0)
                         {for (i=1;i<=testspins.n;++i)
                          {if ((*testspins.configurations[i]).na()==sps.na()&&
			       (*testspins.configurations[i]).nb()==sps.nb()&&
			       (*testspins.configurations[i]).nc()==sps.nc())
                             {sps.wasstable=i;break;}
                          }
                         }
			if (sps.wasstable==0){fprintf(stderr,"internal ERROR htcalc - calculating periodicity not possible");exit(EXIT_FAILURE);}
			//---mark it as stable with periodicity key---
			(*testspins.configurations[physprops.j]).wasstable=sps.wasstable;    
                       }
	      }
    else     // ---- or take q vector
            { sps.spinfromq(testqs.na(-physprops.j),testqs.nb(-physprops.j),
	              testqs.nc(-physprops.j),testqs.q(-physprops.j),
		      testqs.nettom(-physprops.j),testqs.momentq0(-physprops.j),
		      testqs.phi(-physprops.j));
	      }

     sps=spsmin;//take spinconfiguration which gave minimum free energy as starting value
     if (H*sps.nettomagmom(inputpars.gJ)<0)   //see if nettomoment positiv
        {sps.invert();} //if not - invert spinconfiguration
  // now really calculate the physical properties
      mf=new mfcf(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,inputpars.nofcomponents);
      physprops.fe=fecalc(H ,T,inputpars,sps,(*mf),physprops.u,testspins,testqs); 
             // display spinstructure
                if (verbose==1)
                {Vector abc(1,3);
		 abc(1)=inputpars.a;
		 abc(2)=inputpars.b;
		 abc(3)=inputpars.c;
		 float x[inputpars.nofatoms],y[inputpars.nofatoms],z[inputpars.nofatoms];
		 for (is=1;is<=inputpars.nofatoms;++is)
		   {x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
                     sprintf(text,"recalculated: fe=%g,femin=%g:T=%gK,|H|=%gT,Ha=%gT, Hb=%gT, Hc=%gT, %i spins",physprops.fe,femin,T,Norm(H),H(1),H(2),H(3),sps.n());
                    fin_coq = fopen_errchk ("./results/.spins3dab.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,4,inputpars.gJ);
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dac.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,5,inputpars.gJ);
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dbc.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,6,inputpars.gJ);
                    fclose (fin_coq);
		 fin_coq = fopen_errchk ("./results/.spins.eps", "w");
                     sps.eps(fin_coq,text);
                    fclose (fin_coq);
		}
 //check if fecalculation gives again correct result
   if (physprops.fe>femin+(0.00001*fabs(femin))){fprintf(stderr,"Warning htcalc.c: at T=%g K /  H= %g Tfemin=%4.9g was calc.(conf no %i),\n but recalculation  gives fe= %4.9gmeV -> no structure saved\n",
                            T,Norm(H),femin,physprops.j,physprops.fe);delete mf;return 2;}
 physpropclc(H,T,sps,(*mf),physprops,inputpars);
      delete mf;
 }

return 0; // ok we are done with this (HT) point- return ok
}




/************************************************************************/

void physpropclc(Vector H,double T,spincf & sps,mfcf & mf,physproperties & physprops,par & inputpars)
{ int i,j,k,l,n;div_t result;float mmax; char text[100];FILE * fin_coq;
 //save fe and u
 // calculate nettomoment from spinstructure
     physprops.m=sps.nettomagmom(inputpars.gJ)/(double)sps.n()/(double)sps.nofatoms;
 

// thermal expansion - magnetostricton correlation-functions
// according to each neighbour given in mcphas.j a correlation
// function is calculated - up to ini.nofspincorrs neighbours
 int nmax,i1,j1,k1,l1,i2,j2,k2,is;Vector xyz(1,3),d(1,3),d_rint(1,3);
 nmax=ini.nofspincorrs;
 for (l=1;l<=inputpars.nofatoms;++l){if(nmax>(*inputpars.jjj[l]).paranz){nmax=(*inputpars.jjj[l]).paranz;}}

 for(n=1;n<=nmax;++n){physprops.jj[n]=0;for(l=1;l<=inputpars.nofatoms;++l){
    //calculate spincorrelation function of neighbour n of sublattice l
    // 1. transform dn(n) to primitive lattice
     xyz=(*inputpars.jjj[l]).dn[n]-(*inputpars.jjj[l]).xyz; 
     d=inputpars.rez*(const Vector&)xyz;
     for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); //round relative position to integer numbers (to do 
                                             // something sensible if not integer, i.e. if sublattice
					     // of neighbour has not been identified by par.cpp)

     // go through magnetic unit cell and sum up the contribution of every atom
      for(i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k){
         // now we are at atom ijk of sublattice l: the spin is  sps.m(i,j,k)(1..3+3*(l-1))
	 // which i1 j1 k1 l1 is the neighbour in distance d ?
	l1=(*inputpars.jjj[l]).sublattice[n];  // ... yes, this is the sublattice of the neighbour

	i1=i+(int)(d_rint(1));
	j1=j+(int)(d_rint(2));
	k1=k+(int)(d_rint(3));
        while (i1<=0) i1+=sps.na();result=div(i1,sps.na());i1=result.rem;
        while (j1<=0) j1+=sps.nb();result=div(j1,sps.nb());j1=result.rem;
        while (k1<=0) k1+=sps.nc();result=div(k1,sps.nc());k1=result.rem;
        result=div(i1,sps.na());i1=result.rem; if(i1==0)i1=sps.na();
        result=div(j1,sps.nb());j1=result.rem; if(j1==0)j1=sps.nb();
        result=div(k1,sps.nc());k1=result.rem; if(k1==0)k1=sps.nc();

	// sum up correlation function
           for(i2=1;i2<=inputpars.nofcomponents;++i2)
               {physprops.jj[n](i2+inputpars.nofcomponents*inputpars.nofcomponents*(l-1))+=
	        (sps.m(i,j,k)(i2+inputpars.nofcomponents*(l-1)))*
		(sps.m(i1,j1,k1)(i2+inputpars.nofcomponents*(l1-1))); //<JaJa>,<JbJb> ...
               }
               k2=inputpars.nofcomponents;       
           for(i2=1;i2<=inputpars.nofcomponents-1;++i2)
              {for(j2=i2+1;j2<=inputpars.nofcomponents;++j2)
                        {++k2;physprops.jj[n](k2+inputpars.nofcomponents*inputpars.nofcomponents*(l-1))+=
			 (sps.m(i,j,k)(i2+inputpars.nofcomponents*(l-1)))*
			 (sps.m(i1,j1,k1)(j2+inputpars.nofcomponents*(l1-1))); //<JaJb>
			++k2;physprops.jj[n](k2+inputpars.nofcomponents*inputpars.nofcomponents*(l-1))+=
			 (sps.m(i,j,k)(j2+inputpars.nofcomponents*(l-1)))*
			 (sps.m(i1,j1,k1)(i2+inputpars.nofcomponents*(l1-1))); //<JbJa>
			}
	      }
      }}}
 }physprops.jj[n]/=sps.n(); // divide by number of basis sets in magnetic unit cell
 }
     
// neutron intensities (structure factor)
      // check if maxnofhklis was modified by user 
  if (ini.maxnofhkls!=physprops.maxnofhkls){physprops.update_maxnofhkls(ini.maxnofhkls);}
  mmax=0; int qh,qk,ql;j=0;
  ComplexVector a(1,3),qeuklid(1,3);
    complex<double> piq(0,2*3.1415926535),g;
  ComplexVector * mq;  
  Vector ri(1,3);
  double QQ;
 mq = new ComplexVector [sps.in(sps.na(),sps.nb(),sps.nc())+2];for(i=0;i<=sps.in(sps.na(),sps.nb(),sps.nc())+1;++i){mq[i]=ComplexVector(1,sps.nofcomponents*sps.nofatoms);}
  float in;
  sps.FT(mq); //Fourier trafo of momentum configuration
     // mq[n]=mq[0];

  for(qh=0;qh<=sps.na();++qh){for(qk=0;qk<=sps.nb();++qk){for(ql=0;ql<=sps.nc();++ql)
   { 
      Vector q(1,3),hkl(1,3),Q(1,3);
      // this is q- Vector
      q(1)=(double)qh/sps.na();q(2)=(double)qk/sps.nb();q(3)=(double)ql/sps.nc();

      // try different Q vectors corresponding to q !!
     int i1,j1,k1,i2,k2,j2,i1max,j1max,k1max;
     ri(1)=inputpars.a*inputpars.r(1,1);
     ri(2)=inputpars.b*inputpars.r(2,1);
     ri(3)=inputpars.c*inputpars.r(3,1);
     i1max=(int)(ini.maxQ*abs(ri)/2/3.1415926535)+1;
     ri(1)=inputpars.a*inputpars.r(1,2);
     ri(2)=inputpars.b*inputpars.r(2,2);
     ri(3)=inputpars.c*inputpars.r(3,2);
     j1max=(int)(ini.maxQ*abs(ri)/2/3.1415926535)+1;
     ri(1)=inputpars.a*inputpars.r(1,3);
     ri(2)=inputpars.b*inputpars.r(2,3);
     ri(3)=inputpars.c*inputpars.r(3,3);
     k1max=(int)(ini.maxQ*abs(ri)/2/3.1415926535)+1;



              for (i1=0;i1<=i1max;++i1){for(i2=-1;i2<=1&&i2-i1!=1;i2+=2){
              for (j1=0;j1<=j1max;++j1){for(j2=-1;j2<=1&&j2-j1!=1;j2+=2){
              for (k1=0;k1<=k1max;++k1){for(k2=-1;k2<=1&&k2-k1!=1;k2+=2){
       Q(1)=q(1)+i1*i2;Q(2)=q(2)+j1*j2;Q(3)=q(3)+k1*k2;
        //project back to big lattice  
       hkl=inputpars.rez.Transpose()*Q;
      qeuklid(1)=hkl(1)/inputpars.a;
      qeuklid(2)=hkl(2)/inputpars.b;
      qeuklid(3)=hkl(3)/inputpars.c;
      QQ=Norm(qeuklid)*2*3.1415;

     a=0;
     for(l=1;l<=inputpars.nofatoms;++l)
      {//multiply mq by lattice positions exp(iqr_i) and sum into a
       ri=inputpars.rez*(const Vector&)(*inputpars.jjj[l]).xyz; // ri ... atom position with respect to primitive lattice
       g=exp(piq*(Q*ri))*(*inputpars.jjj[l]).debyewallerfactor(QQ)*(*inputpars.jjj[l]).F(QQ)/2.0; // and formfactor + debey waller factor
       if((*inputpars.jjj[l]).gJ==0)
        {for(l1=1;l1<=6&&l1<=inputpars.nofcomponents;++l1){
         if(l1==2||l1==4||l1==6){a(l1/2)+=g*mq[sps.in(qh,qk,ql)](inputpars.nofcomponents*(l-1)+l1);}
         else                   {a((l1+1)/2)+=2.0*g*mq[sps.in(qh,qk,ql)](inputpars.nofcomponents*(l-1)+l1);}
                                                           }
        }
       else
        {for(l1=1;l1<=3&&l1<=inputpars.nofcomponents;++l1){
         a(l1)+=g*mq[sps.in(qh,qk,ql)](inputpars.nofcomponents*(l-1)+l1)*(*inputpars.jjj[l]).gJ;
                                                           }
        }
      } 
      if (QQ<ini.maxQ||(i1==0&&j1==0&&k1==0&&i2==-1&&j2==-1&&k2==-1&&ini.maxQ==0))
          {
	   //calculate intensity using polarization factor (a^*.a-|a.q/|q||^2)
           
	   if (Q==0.0) // tests showed that a*a is square of norm(a) !!! (therefore no conjugate !!!) changed 12.4.03 !!
	    {in=abs((complex<double>)(a*a));}    
           else       
	    {in=abs((complex<double>)(a*a)-(a.Conjugate()*qeuklid)*(a*qeuklid)/(qeuklid*qeuklid));}

           in/=sps.n()*sps.n()*sps.nofatoms*sps.nofatoms;
           //fprintf(stdout,"%i %i %i %g %g %g\n",qh,qk,ql,in,real(a(1)),imag(a(1)));
           //not implemented: lorentzfactor

           if (in>0.0001&&ini.maxnofhkls>0)
	    {if(j<ini.maxnofhkls)
        	{++j; physprops.hkli[j](1)=hkl(1);
	         physprops.hkli[j](2)=hkl(2);
                 physprops.hkli[j](3)=hkl(3);
	         physprops.hkli[j](4)=in;
		 physprops.hkli[j](5)=abs((complex<double>)a(1)/(double)sps.n()/(double)sps.nofatoms);
		 physprops.hkli[j](6)=abs((complex<double>)a(2)/(double)sps.n()/(double)sps.nofatoms);
		 physprops.hkli[j](7)=abs((complex<double>)a(3)/(double)sps.n()/(double)sps.nofatoms);

// if (abs((qeuklid*qeuklid)*2.0*2.0*3.1415*3.1415-1.92247*1.92247)<0.01) {fprintf(stderr,"a1=%g+i%g  a2=%g+i%g  a3=%g+i%g\n",
//	      real(a(1)),imag(a(1)),real(a(2)),imag(a(2)),real(a(3)),imag(a(3)));
//            fprintf(stderr,"in=%g a*a=%g+i%g\n",in,real((complex<double>)(a.Conjugate()*a)),imag((complex<double>)(a.Conjugate()*a)));
//            fprintf(stderr,"in=%g aa=%g+i%g\n",in,real((complex<double>)(a*a)),imag((complex<double>)(a*a)));
//            fprintf(stderr,"q^2=%g+i%g\n",real((qeuklid*qeuklid)),imag((qeuklid*qeuklid)));
//            fprintf(stderr,"(a*q)=%g+i%g (aq)=%g+i%g\n",real((a.Conjugate()*qeuklid)),imag((a.Conjugate()*qeuklid)),
//              real((a*qeuklid)),imag((a*qeuklid)));
//           }
		 }
	    else
	        {int jmin,jj; //determine minimal intensity in list
		 jmin=1;for(jj=1;jj<=j;++jj){if(physprops.hkli[jj](4)<physprops.hkli[jmin](4)){jmin=jj;}}
		 if (in>physprops.hkli[jmin](4)) //if calc int is bigger than smallest value in list
		    { physprops.hkli[jmin](1)=hkl(1);
	              physprops.hkli[jmin](2)=hkl(2);
                      physprops.hkli[jmin](3)=hkl(3);
	              physprops.hkli[jmin](4)=in;                                 
	              physprops.hkli[jmin](5)=abs((complex<double>)a(1)*(1.0/(double)sps.n()/sps.nofatoms));                                 
	              physprops.hkli[jmin](6)=abs((complex<double>)a(2)*(1.0/(double)sps.n()/sps.nofatoms));                                 
	              physprops.hkli[jmin](7)=abs((complex<double>)a(3)*(1.0/(double)sps.n()/sps.nofatoms));}                                 
		}	 
            }
           }
      }}}}}}

   }}}physprops.nofhkls=j;
//save spinarrangement 
      physprops.sps=sps;

//save mfarrangement 
      physprops.mf=mf;

// display spinstructure
           Vector abc(1,3);
		 abc(1)=inputpars.a;
		 abc(2)=inputpars.b;
		 abc(3)=inputpars.c;
		 float x[inputpars.nofatoms],y[inputpars.nofatoms],z[inputpars.nofatoms];
		 for (is=1;is<=inputpars.nofatoms;++is)
		   {x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
    sprintf(text,"physpropclc:T=%gK, |H|=%gT, Ha=%gT, Hb=%gT, Hc=%gT  %i spins",T,Norm(H),H(1),H(2),H(3),sps.n());
                    fin_coq = fopen_errchk ("./results/.spins3dab.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,4,inputpars.gJ);
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dac.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,5,inputpars.gJ);
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dbc.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,6,inputpars.gJ);
                    fclose (fin_coq);
   fin_coq = fopen_errchk ("./results/.spins.eps", "w");
    sps.eps(fin_coq,text);
   fclose (fin_coq);

  //sps.display(text);
delete []mq; 
}




/*****************************************************************************/
// this sub checks if a spinconfiguration has already been added to
// table testspins and adds it if necessary
int checkspincf(int j,spincf & sps,qvectors & testqs,Vector & nettom,
		     Vector & momentq0, Vector & phi, 
                     testspincf & testspins,physproperties & physprops)
{ int i;
  spincf spq(1,1,1,sps.nofatoms,sps.nofcomponents);
 
// compare spinconfigurations stabilized by 
// index j with existing spinconfigurations in testspins

      

   // 1. check initial config: take just used nettom,momentq0,phi for comparison
   if (j<0)
   {spq.spinfromq(testqs.na(-j),testqs.nb(-j),testqs.nc(-j),testqs.q(-j),
                  nettom,momentq0,phi);
    if (spq==sps) {physprops.j=j;testqs.nettom(-j)=nettom;
                  testqs.momentq0(-j)=momentq0;testqs.phi(-j)=phi;return 1;} //ok
   } 


// compare new configuration to all stored configurations 
//check all spinconfigurations

 for (i= -testqs.nofqs();i<=testspins.n;++i)
 {
  if (i>0) 
   {if (sps==(*testspins.configurations[i])) 
	 {
	 physprops.j=i;return 1;} //ok
   }
  if (i<0)
  { spq.spinfromq(testqs.na(-i),testqs.nb(-i),testqs.nc(-i),testqs.q(-i),
                   testqs.nettom(-i),testqs.momentq0(-i),testqs.phi(-i));
    if (spq==sps){
    physprops.j=i;return 1;} //ok
   
  }
 } 

// if it gets here, the spins sps configuration has not been found
// -. add configuration to testspins
return (physprops.j=testspins.addspincf(sps));  //ok=1
}




/*****************************************************************************/
// here the free energy is calculated for a given (initial) spinconfiguration
// using the meanfield algorithm 
double fecalc(Vector  Hex,double T,par & inputpars,
             spincf & sps,mfcf & mf,double & u,testspincf & testspins, qvectors & testqs)
{/*on input:
    T		Temperature[K]
    Hex		Vector of external magnetic field [T]
    inputpars	exchange and other parameters
    sps		initial spinconfiguration
    testspins	all other testspinconfigurations
  on return:
    returns free energy[meV]
    sps		selfconsistently stabilized spinconfiguration (may be different
		from initial spinconfiguration)
    u		mangetic energy[meV]
        
 */
 double fe; // free energy 
 Vector diff(1,inputpars.nofcomponents*inputpars.nofatoms),d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
 Vector meanfield(1,inputpars.nofcomponents),moment(1,inputpars.nofcomponents),d1(1,inputpars.nofcomponents);
 char text[60]; // some text variable
 int i,j,k,i1,j1,k1,di,dj,dk,l,r,s,sdim,m,n,m1;
 div_t result; // some modulo variable
 float    sta=1000000; // initial value of standard deviation
 float staold=2000000; 
 float bigstep;
 float smallstep;
 int slowct=10;
 float stepratio=1.0;

 float spinchange=0; // initial value of spinchange
 sdim=sps.in(sps.na(),sps.nb(),sps.nc()); // dimension of spinconfigurations
 Vector  * lnzi; lnzi=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){lnzi[i]=Vector(1,inputpars.nofatoms);} // partition sum for every atom
 Vector  * ui; ui=new Vector [sdim+2];for(i=0;i<=sdim+1;++i){ui[i]=Vector(1,inputpars.nofatoms);} // magnetic energy for every atom
 ComplexMatrix ** ests;ests=new ComplexMatrix*[inputpars.nofatoms*sdim+2];
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {for (l=1;l<=inputpars.nofatoms;++l){
  ests[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1]=new ComplexMatrix((*inputpars.jjj[l]).est);
  }}}}
 int diagonalexchange=1;
 FILE * fin_coq;
 time_t time_of_last_output=0; 
 
 spincf  spsold(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,inputpars.nofcomponents); // spinconf variable to store old sps
 mfcf  mfold(mf.na(),mf.nb(),mf.nc(),inputpars.nofatoms,inputpars.nofcomponents); // spinconf variable to store old mf
 spsold=sps;
 
// coupling coefficients jj[](a-c) berechnen
// for (r=0;r<=sdim;++r)

 Matrix * jj; //if (inputpars.diagonalexchange()==0){i=9;}else{i=3;}
  jj= new Matrix [(sdim+1)+1];for(i=0;i<=sdim+1;++i){jj[i]=Matrix(1,inputpars.nofcomponents*inputpars.nofatoms,1,inputpars.nofcomponents*inputpars.nofatoms);} // coupling coeff.variable
   if (jj == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

   // initialize mfold with zeros 
   for(s=0;s<=mfold.in(mfold.na(),mfold.nb(),mfold.nc());++s){mfold.mi(s)=0;} 
   for(s=0;s<=sdim;++s){jj[s]=0;} //clear jj(j,...)
   
   for(m=1;m<=inputpars.nofatoms;++m)
   {if ((*inputpars.jjj[m]).diagonalexchange==0){diagonalexchange=0;} // if any ion has anisotropic exchange - calculate anisotropic
    for(l=1;l<=(*inputpars.jjj[m]).paranz;++l)
    {//sum up l.th neighbour interaction of atom m
                                             // atom m = sublattice m
	n=(*inputpars.jjj[m]).sublattice[l]; // n set to sublattice of neighbor l

    // determine s (index of difference between crystal unit cells in the magnetic supercell)
    // start with calculating the difference vector of origins of crystal unit cells
                   // bugfix GdVO3: sign of 2nd term changed and last term added 12.12.07
     xyz=(*inputpars.jjj[m]).dn[l]+(*inputpars.jjj[m]).xyz-(*inputpars.jjj[n]).xyz; 
 
    // transform distance vector to primitive lattice
     d=inputpars.rez*(const Vector&)xyz;
     
     for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); //round relative position to integer numbers (to do 
                                             // something sensible if not integer, i.e. if sublattice
					     // of neighbour has not been identified by par.cpp)
        
        i=(int)(d_rint(1));
	j=(int)(d_rint(2));
	k=(int)(d_rint(3));
        // here we have the difference between crystal unitc cells ijk in the magnetic
        // supercell given by the indices i j k: if they point out of the magnetic supercell
        // they are folded back into it in the next 3 lines: this is allowed  because it is 
        // irrelevant for the mean field summation
        // where the neighbor actually sits, but only on which sublattice it sits...
        while (i<=0) i+=sps.na();result=div(i,sps.na());i=result.rem; // only distance is important ...
        while (j<=0) j+=sps.nb();result=div(j,sps.nb());j=result.rem;
        while (k<=0) k+=sps.nc();result=div(k,sps.nc());k=result.rem;
      // s is determined from a vector ijk connecting the different crystal unit cells
	s=sps.in(i,j,k); //ijk range here from 0 to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!

        //     myPrintMatrix(stdout,(*inputpars.jjj[m]).jij[l]);

	// sum up the contribution of the interaction parameter to the interaction matrix jj[s] to be
        // used in the meanfield calculation below
	for(i=1;i<=inputpars.nofcomponents;++i){for(j=1;j<=inputpars.nofcomponents;++j){
	  jj[s](inputpars.nofcomponents*(m-1)+i,inputpars.nofcomponents*(n-1)+j)+=(*inputpars.jjj[m]).jij[l](i,j); 

	//remark: function par:jij(l) returns exchange constants (*inputpars.jjj[1]).jij[l](1-9)	   	
        }}

    } 
   }
  
   
   
if (ini.displayall==1)   // display spincf if button is pressed
 {
     fin_coq = fopen_errchk ("./results/.spins.eps", "w");
     sprintf(text,"fecalc:%i spins, iteration 0",sps.n());
     sps.eps(fin_coq,text);
     fclose (fin_coq);
     sleep(2);

// sps.display(text);
 }


// loop for selfconsistency
for (r=1;sta>ini.maxstamf;++r)
{if (r>ini.maxnofmfloops)
    {delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
      delete ests[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []ests;

     if (verbose==1) fprintf(stderr,"feDIV!MAXlooP");
     return 20000;}
 if (spinchange>ini.maxspinchange)
    {delete []jj;delete []lnzi;delete []ui;
          for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
      delete ests[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []ests;
     if (verbose==1) fprintf(stderr,"feDIV!MAXspinchangE");
     return 20001;}

 //1. calculate mf from sps (and calculate sta)
 sta=0;
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {mf.mf(i,j,k)=0;
  for (l=1;l<=inputpars.nofatoms;++l){
    if(inputpars.gJ(l)==0)              {
     for (i1=1;i1<=6&&i1<=inputpars.nofcomponents;++i1){
            if(i1==2||i1==4||i1==6){ mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+i1]=Hex(i1/2)*MU_B;}
	    else                   { mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+i1]=2*Hex((i1+1)/2)*MU_B;}
  				                       }
            				}
    else                                {
     for (i1=1;i1<=3&&i1<=inputpars.nofcomponents;++i1){
               mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+i1]=Hex(i1)*inputpars.gJ(l)*MU_B;
  				                       }
					}
				     }
  for (i1=1;i1<=sps.na();++i1){if (i<i1){di=i-i1+sps.na();}else{di=i-i1;} 
                               for (j1=1;j1<=sps.nb();++j1){if (j<j1){dj=j-j1+sps.nb();}else{dj=j-j1;}
			                                    for (k1=1;k1<=sps.nc();++k1){if (k<k1){dk=k-k1+sps.nc();}else{dk=k-k1;}
    l=sps.in(di,dj,dk);//di dj dk range from 0 to to sps.na()-1,sps.nb()-1,sps.nc()-1 !!!!
                       // and index a difference between crystal unit cell positions in the
                       // magnetic supercell
     
     // here the contribution of the crystal unit cell i1 j1 k1 (i1,j1,k1 indicate the
     // position of the crystal unit cell in the magnetic supercell) to the mean field
     // of the crystal unit cell i j k is calculated by one matrix multiplication
     if (diagonalexchange==0||inputpars.nofatoms>1)
     {mf.mf(i,j,k)+=jj[l]*(const Vector&)sps.m(i1,j1,k1);
     }else
     {//do the diagonal elements separately to accellerate the sum
      for(m1=1;m1<=inputpars.nofatoms*inputpars.nofcomponents;++m1)
         {mf.mf(i,j,k)(m1)+=sps.m(i1,j1,k1)(m1)*jj[l](m1,m1);}
     }
    }}}
  diff=mf.mf(i,j,k)-mfold.mf(i,j,k);sta+=diff*diff;
  diff*=stepratio;mf.mf(i,j,k)=mfold.mf(i,j,k)+diff;//step gently ... i.e. scale change of MF with stepratio
  }}}
  mfold=mf;      
  sta=sqrt(sta/sps.n()/inputpars.nofatoms);
  
  bigstep=fmodf(ini.bigstep-0.0001,1.0);
  if (ini.bigstep>1.0){smallstep=bigstep/(ini.bigstep-bigstep);}else{smallstep=bigstep/5;}
  if (r==1) {stepratio=smallstep;} //in first loop initialize stepratio to smallstep
  if (staold<sta&&stepratio==bigstep){stepratio=smallstep;slowct=10;}//if sta increases then set stepratio to bigstep
  if (staold>sta&&stepratio<bigstep){--slowct;if (slowct<=0)stepratio=bigstep;} // at least for 10 cycles
  staold=sta;

//2. calculate sps from mf
 for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
 {diff=sps.m(i,j,k);s=sps.in(i,j,k);
  for(l=1;l<=inputpars.nofatoms;++l)
  {int lm1m3;
   lm1m3=inputpars.nofcomponents*(l-1);
   for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {d1[m1]=mf.mf(i,j,k)[lm1m3+m1];}
   moment=(*inputpars.jjj[l]).mcalc(T,d1,lnzi[s][l],ui[s][l],(*ests[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1]));
   for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {sps.m(i,j,k)(lm1m3+m1)=moment[m1];}
  }
  diff-=sps.m(i,j,k);
  spinchange+=sqrt(diff*diff)/sps.n();  
  }}}

  //treat program interrupts 
  checkini(testspins,testqs); 
if (ini.displayall==1)  // if all should be displayed - write sps picture to file .spins.eps
 {
      fin_coq = fopen_errchk ("./results/.spins.eps", "w");
     sprintf(text,"fecalc:%i spins, iteration %i sta=%g",sps.n(),r,sta);
     sps.eps(fin_coq,text);
     fclose (fin_coq);
     sleep(2);
 }
   //for verbose mode do some outputs
 if (verbose==1)
 {if (time(0)-time_of_last_output>2)
  {time_of_last_output=time(0);
   fin_coq= fopen_errchk ("./results/.fe_status.dat","a");
   fprintf(fin_coq,"%i %g %g %g %g\n",time(0),log((double)r)/log(10.0),log(sta)/log(10.0),spinchange,stepratio);
   fclose(fin_coq);
  }
 }

}

// calculate free energy fe and energy u
fe=0;u=0;
for (i=1;i<=sps.na();++i){for (j=1;j<=sps.nb();++j){for (k=1;k<=sps.nc();++k)
{s=sps.in(i,j,k);
 for(l=1;l<=inputpars.nofatoms;++l)
 {fe-=K_B*T*lnzi[s][l];
  u+=ui[s][l]; 
// correction term
  for(m1=1;m1<=inputpars.nofcomponents;++m1)
   {d1[m1]=sps.m(i,j,k)[inputpars.nofcomponents*(l-1)+m1];
  meanfield[m1]=mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+m1];}
  // subtract external field (only necessary for magnetic field, not for quadrupolar fields,
  // because the Cf parameters are treated separately in mcalc and not as part of the quadrupolar
  // field)
  if(inputpars.gJ(l)==0)
  {for(m1=1;m1<=6&&m1<=inputpars.nofcomponents;++m1)
    {if(m1==2||m1==4||m1==6) {meanfield[m1]-=Hex[m1/2]*MU_B;}
     else                    {meanfield[m1]-=2*Hex[(m1+1)/2]*MU_B;}
    }
  }
  else
  {
  for(m1=1;m1<=3&&m1<=inputpars.nofcomponents;++m1){meanfield[m1]-=Hex[m1]*inputpars.gJ(l)*MU_B;}
  }
  
  // add correction term
  fe+=0.5*(meanfield*d1);
  u+=0.5*(meanfield*d1);
//  printf ("Ha=%g Hb=%g Hc=%g ma=%g mb=%g mc=%g \n", H[1], H[2], H[3], m[1], m[2], m[3]);
 }
}}}
fe/=(double)sps.n()*sps.nofatoms; //norm to formula unit
u/=(double)sps.n()*sps.nofatoms;

if (ini.displayall==1)
 {
      fin_coq = fopen_errchk ("./results/.spins.eps", "w");
       sprintf(text,"fecalc:%i spins, iteration %i, fe=%gmeV",sps.n(),i,fe);
       sps.eps(fin_coq,text);
      fclose (fin_coq);
      sleep(2);
// sps.display(text);
  }

 delete []jj;delete []lnzi;delete []ui;
     for (i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k)
     {for (l=1;l<=inputpars.nofatoms;++l){
 delete ests[inputpars.nofatoms*sps.in(i-1,j-1,k-1)+l-1];
     }}}} delete []ests;

return fe;     
}

