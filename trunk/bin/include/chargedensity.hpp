//  class chargedensity for the storage of a chargedensity surface
// used in program mcphase

#ifndef CHARGEDENSITY
#define CHARGEDENSITY

#include<cstdio>
#include<cstring>
#include<cerrno>
#include<martin.h>
#include<vector.h>

class chargedensity
{ private:
  double dtheta;
  double dfi;
  int nofpts;  
  Vector **rtetafi;
  
  public:
  Vector & rtf(int i); // pointer to rft (i=1 ... nofpoints)
  int nofpoints(); // returns nofpoints in chargedensitysurface  
//constructor
  chargedensity(double dtheta,double dfi);
  chargedensity(const chargedensity & p);//kopier-konstruktor
 ~chargedensity ();//destruktor
};

#endif

