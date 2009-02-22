//  class chargedensity for the storage of a chargedensity surface
// used in program mcphase

#include"chargedensity.hpp"
#include<cstdio>
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

