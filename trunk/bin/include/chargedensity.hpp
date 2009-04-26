//  class chargedensity for the storage of a chargedensity surface
// used in program mcphase

#ifndef CHARGEDENSITY
#define CHARGEDENSITY

#include<cstdio>
#include<cstring>
#include<cerrno>
#include<martin.h>
#include<vector.h>
#include<ionpars.hpp>

class chargedensity
{ private:
  int nofpts;  
  Vector **rtetafi;
  
  public:
  double dtheta;
  double dfi;
  Vector & rtf(int i); // pointer to rft (i=1 ... nofpoints)
  int nofpoints(); // returns nofpoints in chargedensitysurface  
  
  // sub for calculation of charge density given a radiu R and polar angles teta, 
  // fi and expansion coeff. alm
  double rocalc (double & teta,double & fi,double & R, Vector & moments,ionpars & ionpar);

  void calc_cd_surface(Vector & moments,ionpars & ionpar,double ccc);
  //ccc ... surface value of chargedensity


  //constructor
  chargedensity(double dtheta,double dfi);
  chargedensity(const chargedensity & p);//kopier-konstruktor
 ~chargedensity ();//destruktor
};

#endif

