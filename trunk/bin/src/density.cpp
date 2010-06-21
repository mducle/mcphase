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
Vector momentsx(1,49);
Vector momentsy(1,49);
Vector momentsz(1,49);
Vector momentlx(1,49);
Vector momently(1,49);
Vector momentlz(1,49);
int i;Vector mom(1,3);

if(moments.Hi()>=3*49&&type<=3){type+=4;} // now the type of density is clear

    switch (type)
   {case 1: break; // component of spindensity
    case 2: break; // component of orbmomdensity
    case 3: for(i=1;i<=49;++i){momentsx(i)=moments(i);momentlx(i)=moments(i+49);}
            break; // component of total momdensity
    case 8:
    case 4: for(i=1;i<=49;++i){momentlx(i)=moments(i);momently(i)=moments(i+49);momentlz(i)=moments(i+2*49);}
            break; // absvalue of current density
    case 5: for(i=1;i<=49;++i){momentsx(i)=moments(i);momentsy(i)=moments(i+49);momentsz(i)=moments(i+2*49);}
            break; // absvalue of spindensity
    case 6: for(i=1;i<=49;++i){momentlx(i)=moments(i);momently(i)=moments(i+49);momentlz(i)=moments(i+2*49);}
            break; // absvalue of orbmomdensity
    case 7: for(i=1;i<=49;++i){momentsx(i)=moments(i);momentsy(i)=moments(i+49);momentsz(i)=moments(i+2*49);
                                     momentlx(i)=moments(i+3*49);momently(i)=moments(i+4*49);momentlz(i)=moments(i+5*49);}
            break; // absvalue of tot momdensity
    default: ;     // chargedensity
   }


int imax = 3,iii,iv;
Vector rp(1,imax),rrttff(1,6);
double rmax = 0;
double rstp = .1;
double max = .01 * ccc;  //end of intervalschachtelung to find r(ro=ccc)
double theta,fi,R,rin,ro,deltaa,rstpp,delta1;
int nt,tt,ff;
for(tt=0;tt<=3.1415/dtheta;++tt){for(ff=0;ff<=2*3.1415/dfi;++ff){
   rp=0;theta=(double)tt*dtheta;fi=(double)ff*dfi;
   nt=0;
   for(rin=0.1;rin<=3.0;rin+=0.2){
   R=rin; 
    switch (type) 
     {case 1:      ro=ionpar.spindensity_calc(theta,fi,R,moments);break;
      case 2:      ro=ionpar.orbmomdensity_calc(theta,fi,R,moments);break;
      case 3:      ro=ionpar.spindensity_calc(theta,fi,R,momentsx)+
                      ionpar.orbmomdensity_calc(theta,fi,R,momentlx);break;
      case 4:      ro=Norm(ionpar.currdensity_calc(theta,fi,R,momentlx,momently,momentlz));break;
      case 5:      ro=Norm(ionpar.spindensity_calc(theta,fi,R,momentsx,momentsy,momentsz));break;
      case 6:      ro=Norm(ionpar.orbmomdensity_calc(theta,fi,R,momentlx,momently,momentlz));break;
      case 7:      ro=Norm(ionpar.spindensity_calc(theta,fi,R,momentsx,momentsy,momentsz)+
                           ionpar.orbmomdensity_calc(theta,fi,R,momentlx,momently,momentlz));break;
      case 8:      ro=pr*ionpar.currdensity_calc(theta,fi,R,momentlx,momently,momentlz);break;
      default:ro=ionpar.rocalc(theta,fi,R,moments);
     }

   deltaa=fabs(ro-ccc);rstpp=rstp;delta1=1e4;
   for(iii=1;(iii<=100)&(delta1>=max);++iii)
     {R+=rstpp;
    switch (type) 
     {case 1:      ro=ionpar.spindensity_calc(theta,fi,R,moments);break;
      case 2:      ro=ionpar.orbmomdensity_calc(theta,fi,R,moments);break;
      case 3:      ro=ionpar.spindensity_calc(theta,fi,R,momentsx)+
                      ionpar.orbmomdensity_calc(theta,fi,R,momentlx);break;
      case 4:      ro=Norm(ionpar.currdensity_calc(theta,fi,R,momentlx,momently,momentlz));break;
      case 5:      ro=Norm(ionpar.spindensity_calc(theta,fi,R,momentsx,momentsy,momentsz));break;
      case 6:      ro=Norm(ionpar.orbmomdensity_calc(theta,fi,R,momentlx,momently,momentlz));break;
      case 7:      ro=Norm(ionpar.spindensity_calc(theta,fi,R,momentsx,momentsy,momentsz)+
                           ionpar.orbmomdensity_calc(theta,fi,R,momentlx,momently,momentlz));break;
      case 8:      ro=pr*ionpar.currdensity_calc(theta,fi,R,momentlx,momently,momentlz);break;
      default:ro=ionpar.rocalc(theta,fi,R,moments);
     }
     delta1=fabs(ro-ccc);
      if(delta1>=max){
                      if(delta1<deltaa){deltaa=delta1;}
		      else{R-=rstpp;rstpp*=-0.5;}
                     }
     }
     if(delta1<max)
     {for(iv=1;(iv<=nt)&(fabs(rp(iv)-R)>=0.1);++iv);// check if this minimum was already stored
      if(fabs(rp(iv)-R)>=0.1) // if not store it (up to nt=imax minima will be stored)
       {++nt;rp(nt)=R;
        if (R>rmax)rmax=R;
       }
     }
   }//next rin
   if(nt==0) {rp(1)=0.05;rp(2)=0.06;} // no minimum found
   if(nt==1) {rp(2)=rp(1)+0.0001;} // one minimum found
   // select most outside value of rp and set rrttff accordingly
   R=0;
   for(iv=1;iv<=imax;++iv){if(rp(iv)>R)R=rp(iv);}// get the maximal radius (outside surface)
	++anzahl;
	rrttff(1)=R;rrttff(2)=theta;rrttff(3)=fi;
        if(type>=4){ switch (type)
     {case 8:
      case 4:      mom=ionpar.currdensity_calc(theta,fi,R,momentlx,momently,momentlz);break;
      case 5:      mom=ionpar.spindensity_calc(theta,fi,R,momentsx,momentsy,momentsz);break;
      case 6:      mom=ionpar.orbmomdensity_calc(theta,fi,R,momentlx,momently,momentlz);break;
      case 7:      mom=ionpar.spindensity_calc(theta,fi,R,momentsx,momentsy,momentsz)+
                           ionpar.orbmomdensity_calc(theta,fi,R,momentlx,momently,momentlz);break;
         }
   if (nt==0)mom=0;// if no minimum has been found, set moment vector to zero.
                    rrttff(4)=mom(1);rrttff(5)=mom(2);rrttff(6)=mom(3);}
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
   if(strncmp(title,"currdensityabsvalue",15)==0){type=4;}
   if(strncmp(title,"currdensityprojection",15)==0){type=8;}
   pr=Vector(1,3);extract(title,"i",pr(1));extract(title,"j",pr(2));extract(title,"k",pr(3));
   nofpts=(int)(2*3.1415/df+1)*(int)(3.1415/dt)+1;
   rtetafi = new Vector * [1+nofpts];
   int i;
   for (i=1;i<=nofpts;++i){rtetafi[i]=new Vector(1,6);(*rtetafi[i])=0;}
   dtheta=dt;
   dfi=df;
 }
  
density::density(const density & p) //kopier-konstruktor
  {dtheta=p.dtheta;dfi=p.dfi;
   nofpts=p.nofpts;
   type=p.type;
   pr=Vector(1,3);pr=p.pr;
   rtetafi = new Vector * [1+nofpts];
   int i;
   for (i=1;i<=nofpts;++i){rtetafi[i]=new Vector(1,6);(*rtetafi[i])=(*p.rtetafi[i]);}
  }

density::~density () //destruktor
  {int i;
   for (i=1;i<=nofpts;++i){delete rtetafi[i];}     
  delete []rtetafi;
  }

