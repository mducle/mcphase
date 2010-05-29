//  class density for the storage of a density surface
// used in program mcphase

#include"density.hpp"
#include<martin.h>

/* private:
  double dtheta;
  double dfi;
  int nofpoints;  
  Vector **rfiteta;
  
  */
  Vector & density::rtf(int p) // pointer to rft (i=1 ... nofpts)
{if (p>(*this).nofpts){fprintf(stderr,"Error class density: point index %i exceed maximum nofpts %i \n",p,nofpts);exit (EXIT_FAILURE);}
 return (*rtetafi[p]);
}

void density::calc_cd_surface(Vector & moments,jjjpar & ionpar,  double ccc)
{ // double ccc = surface value of density
// here the set of points for this ion should be created corresponding to its
// charge density: 
//REM dieses programm berechnet aus alm die ladungsdichteverteilung
//REM und zeichnet diese auf
int anzahl = 0;
Vector momentsspin(1,49);
Vector momentsorb(1,49);
if(type==3){int i;for(i=1;i<=49;++i){momentsspin(i)=moments(i);momentsorb(i)=moments(i+49);}}
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
   R=rin; 
    switch (type) 
     {case 1:      ro=ionpar.spindensity_calc(theta,fi,R,moments);break;
      case 2:      ro=ionpar.orbmomdensity_calc(theta,fi,R,moments);break;
      case 3:      ro=ionpar.spindensity_calc(theta,fi,R,momentsspin)+
                 ionpar.orbmomdensity_calc(theta,fi,R,momentsorb);break;
      default:ro=ionpar.rocalc(theta,fi,R,moments);
     }
   deltaa=fabs(ro-ccc);rstpp=rstp;delta1=1e4;
   for(iii=1;iii<=100&delta1>=max;++iii)
     {R+=rstpp;
    switch (type) 
     {case 1:      ro=ionpar.spindensity_calc(theta,fi,R,moments);break;
      case 2:      ro=ionpar.orbmomdensity_calc(theta,fi,R,moments);break;
      case 3:      ro=ionpar.spindensity_calc(theta,fi,R,momentsspin)+
                 ionpar.orbmomdensity_calc(theta,fi,R,momentsorb);break;
      default:ro=ionpar.rocalc(theta,fi,R,moments);
     }
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
   // select most outsinde value of rp and set rrttff accordingly
   R=0;
   for(iv=1;iv<=imax;++iv){if(rp(iv)>R)R=rp(iv);}
	++anzahl;
	rrttff(1)=R;rrttff(2)=theta;rrttff(3)=fi;
   rtf(anzahl)=rrttff;
   if(tt==0){ff=(int)(2*3.1415/dfi+1);}// this is to avoid more than one point for tt=0
  }}

        
}

 int density::nofpoints() // returns number of points on density surface
 {return nofpts;}
 
//constructor
density::density(char * title,double dt,double df)
  {type=0;
   if(strncmp(title,"spindensity",10)==0){type=1;}
   if(strncmp(title,"orbmomdensity",10)==0){type=2;}
   if(strncmp(title,"momdensity",10)==0){type=3;}
   nofpts=(int)(2*3.1415/df+1)*(int)(3.1415/dt)+1;
   rtetafi = new Vector * [1+nofpts];
   int i;
   for (i=1;i<=nofpts;++i){rtetafi[i]=new Vector(1,3);(*rtetafi[i])=0;}     
   dtheta=dt;
   dfi=df;
 }
  
density::density(const density & p) //kopier-konstruktor
  {dtheta=p.dtheta;dfi=p.dfi;
   nofpts=p.nofpts;
   type=p.type;
   rtetafi = new Vector * [1+nofpts];
   int i;
   for (i=1;i<=nofpts;++i){rtetafi[i]=new Vector(1,3);(*rtetafi[i])=(*p.rtetafi[i]);}     
  }

density::~density () //destruktor
  {int i;
   for (i=1;i<=nofpts;++i){delete rtetafi[i];}     
  delete []rtetafi;
  }

