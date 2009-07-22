//  class chargedensity for the storage of a chargedensity surface
// used in program mcphase

#include"chargedensity.hpp"
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cerrno>
#include<martin.h>
#include<vector.h>

/* private:
  double dtheta;
  double dfi;
  int nofpoints;  
  Vector **rfiteta;
  
  */
  Vector & chargedensity::rtf(int p) // pointer to rft (i=1 ... nofpts)
{if (p>(*this).nofpts){fprintf(stderr,"Error class chargedensity: point index %i exceed maximum nofpts %i \n",p,nofpts);exit (EXIT_FAILURE);}
 return (*rtetafi[p]);
}

void chargedensity::calc_cd_surface(Vector & moments,jjjpar & ionpar,  double ccc)
{ // double ccc = surface value of chargedensity
// here the set of points for this ion should be created corresponding to its
// charge density: 
//REM dieses programm berechnet aus alm die ladungsdichteverteilung
//REM und zeichnet diese auf
int anzahl = 0;
int imax = 3,iii,iv;
Vector rp(1,imax),rrttff(1,3);
double rmax = 0;
double rstp = .1;
double max = .01 * ccc;  //end of intervalschachtelung to find r(ro=ccc)
double theta,fi,R,rin,ro,deltaa,rstpp,delta1,dd;
int nt,tt,ff;
for(tt=0;tt<=3.1415/dtheta;++tt){for(ff=0;ff<=2*3.1415/dfi;++ff){
   rp=0;theta=(double)tt*dtheta;fi=(double)ff*dfi;
   nt=0;
   for(rin=0.1;rin<=2.2;rin+=0.2){
   R=rin;ro=ionpar.rocalc(theta,fi,R,moments);
   deltaa=fabs(ro-ccc);rstpp=rstp;delta1=1e4;
   for(iii=1;iii<=100&delta1>=max;++iii)
     {R+=rstpp;ro=ionpar.rocalc(theta,fi,R,moments);
      delta1=fabs(ro-ccc);
      if(delta1>=max){
                      if(delta1<deltaa){deltaa=delta1;}
		      else{R-=rstpp;rstpp*=-0.5;}
                     }
     }
     if(delta1<max)
     {for(iv=1;iv<=nt&fabs(rp(iv)-R)>=0.1;++iv);
      if(fabs(rp(iv)-R)>=0.1)
       {++nt;rp(nt)=R;
        if (R>rmax)rmax=R;
       }
     }
   }//next rin
   if(nt==0) {rp(1)=0.05;rp(2)=0.06;}
   if(nt==1) {rp(2)=rp(1)+0.0001;}
   // select most outsinde value of rp
   R=0;
   for(iv=1;iv<=imax;++iv){if(rp(iv)>R)R=rp(iv);}
	++anzahl;
	rrttff(1)=R;rrttff(2)=theta;rrttff(3)=fi;
   rtf(anzahl)=rrttff;
   if(tt==0){ff=(int)(2*3.1415/dfi+1);}
  }}

        
}

 int chargedensity::nofpoints() // returns number of points on chargedensity surface
 {return nofpts;}
 
//constructor
chargedensity::chargedensity(double dt,double df)
  {nofpts=(int)(2*3.1415/df+1)*(int)(3.1415/dt)+1;
   rtetafi = new Vector * [1+nofpts];
   int i;
   for (i=1;i<=nofpts;++i){rtetafi[i]=new Vector(1,3);(*rtetafi[i])=0;}     
   dtheta=dt;
   dfi=df;
 }
  
chargedensity::chargedensity(const chargedensity & p) //kopier-konstruktor
  {dtheta=p.dtheta;dfi=p.dfi;
   nofpts=p.nofpts;
   rtetafi = new Vector * [1+nofpts];
   int i;
   for (i=1;i<=nofpts;++i){rtetafi[i]=new Vector(1,3);(*rtetafi[i])=(*p.rtetafi[i]);}     
  }

chargedensity::~chargedensity () //destruktor
  {int i;
   for (i=1;i<=nofpts;++i){delete rtetafi[i];}     
  delete []rtetafi;
  }

