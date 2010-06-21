//  class density for the storage of a charge/spin/orbmom...density surface
// used in program mcphase

#ifndef CHARGEDENSITY
#define CHARGEDENSITY

#include<martin.h>
#include<jjjpar.hpp>

class density
{ private:
  int nofpts;  
  Vector **rtetafi;
  int type;
   Vector pr;

  public:
  double dtheta;
  double dfi;
  Vector & rtf(int i); // pointer to rft (i=1 ... nofpoints)
  int nofpoints(); // returns nofpoints in densitysurface  
  
  // sub for calculation of charge density given a radiu R and polar angles teta, 
  // fi and expansion coeff. alm
  double rocalc (double & teta,double & fi,double & R, Vector & moments,ionpars & ionpar);

  void calc_cd_surface(Vector & moments,jjjpar & ionpar,double ccc);
  //ccc ... surface value of density


  //constructor
  density(char*title,double dtheta,double dfi);
  density(const density & p);//kopier-konstruktor
 ~density ();//destruktor
};

#endif

