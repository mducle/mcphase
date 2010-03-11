//class for the generation and storage of a spinconfiguration
// used in mcphase


#include<cstdlib>
#include<cerrno>
#include<vector.h>
#include<complex>
#include<cstdio>

#ifndef SPINCF_H
#define SPINCF_H

class spincf 
{
  private:
 // number of spins  
   int nofa,nofb,nofc;
 // this subtracts n2 if n1>n2
   int mod(int n1,int n2);
   int mxa,mxb,mxc;
 // frame of display
   void epsarrow(FILE * fout,Vector a,Vector b);
   Vector * mom; // momentums <J>
   int iv[4];
   Vector xy(Vector xyz,int orientation,Vector min,Vector max,float bbwidth,float bbheight);
   int spequal(Vector a,Vector b);// routine to compare spins
     
   void calc_minmax(Vector & min,Vector & max,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc);
   void calc_minmax_scale(Vector & min,Vector & max,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc,double scale_view_1,double scale_view_2,double scale_view_3);
   Vector magmom(int i,int j,int k,int l,Vector & gJ); // returns magnetic moment (1,3)
   Vector moment(int i,int j,int k,int l); // returns moment of atom l (1,nofcomponents)
   double nndist(float * x, float * y, float * z,Vector & abc,Matrix & p,Vector &dd);

 public:
 // array of spins 
   int in(int i,int j, int k); 
    int wasstable; // index to remember if it was stable: if a sinconfiguration is set stable, its periodicity key is stored in wasstable
   
    int nofatoms;
    int nofcomponents;
    Vector & m(int i,int j,int k); // returns pointer to spin (ijk) 
    Vector & mi(int in); // returns pointer to spin i
    void  FT(ComplexVector * mq); // returns Fourier transform mq of spins (for use see htcalc.c)
    int * ijk(int in);  // returns spin indizes (ijk)(in): in=0,...,n(=na*nb*nc)
    
    int n(); // returns total number of spins
    int na(); // returns number of spinsl
    int nb(); // returns number of spins
    int nc(); // returns number of spins
 
    void calc_prim_mag_unitcell(Matrix & p,Vector & abc, Matrix & r);
    Vector nettomagmom ( Vector & gJ); // returns nettomagneticmoment [muB]
    Vector totalJ (); // returns nettomoment <J>
    Vector pos(int i, int j, int k, int l,Vector & abc,Matrix & r,float * x,float *y,float*z); 
                      //returns position of atom l at lattice site (i j k) (Angstrom)
    void invert();// inverts all spins (AND higher order moments)
    void reduce();// reduces spinconfiguration
    void spinfromq (int n1,int n2, int n3,Vector & qvector, Vector & nettom,Vector & momentq0, Vector & phi);

    void print(FILE * fout);

//print list of atoms + positions + moments
    void printall(FILE * fout,Vector & abc,Matrix & r,float * x,float *y,float*z, char ** cffilenames,float * gJ);
    void eps(FILE * fout);
    void eps(FILE * fout,const char * text);
    void eps3d(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z,int orientation, Vector & gJ);
    void fst(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ);
    void jvx(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ,
                 double show_abc_unitcell,double show_primitive_crystal_unitcell,double show_magnetic_unitcell,double show_atoms,double scale_view_1,double scale_view_2,double scale_view_3,
                 int showprim,double phase,spincf & savev_real,spincf & savev_imag,double amplitude,Vector & hkl,
                 double spins_show_ellipses,double spins_scale_moment,double spins_show_static_moment_direction);

              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase

    void jvx_cd(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ,
                 double show_abc_unitcell,double show_primitive_crystal_unitcell,double show_magnetic_unitcell,double show_atoms,double scale_view_1,double scale_view_2,double scale_view_3,
                 int showprim,double phase,spincf & savev_real,spincf & savev_imag,double amplitude,Vector & hkl,
                 double spins_show_ellipses,double spins_scale_moment,char ** cffilenames,double show_chargedensity,double spins_show_static_moment_direction, double threshhold);

     void cd(FILE * fout,Vector & abc,Matrix & r,float * x,float *y,float*z,char ** cffilenames,int showprim,
                int nofpointsi, int nofpointsj, int nofpointsk,double scale_view_1,double scale_view_2,double scale_view_3,
                spincf & savev_real,spincf & savev_imag,double phase,double amplitude,Vector & hkl);

    void fstprim(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ);
    int  load(FILE * fin_coq);	
     
    spincf & operator + (const spincf & op2); // addition    
    spincf & operator * (const double factor); // multiplication with constant
    spincf & operator= (const spincf & op2); // zuweisung
    int operator== (spincf & op2); // vergleich

   
spincf (int n1=1,int n2=1,int n3=1,int nofatoms=1,int nofcomponents=3);
                                   //konstruktor mit initialisierung (wenn noetig)
   spincf (const spincf & spins);	// kopier-konstruktor
   
~spincf ();		//destruktor

};

#endif
