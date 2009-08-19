/* ic1ion_module.cpp
 *
 * Functions:
 *   void myPrintMatrix(FILE * file,sMat<double> & M,int d)                        // Prints out full matrix
 *   bool checkmat(ComplexMatrix &cmat, complexdouble *fmat,int r, int c)          // Compares Matpack and fortran matrices
 *   void mcalc(Vector &J, double *T, Vector &gjmbH, double *gJ, Vector &ABC,      // Calculates the meanfield moment
 *                 char **sipffile, double *lnZ, double *U, ComplexMatrix &est)    //
 *   int dmcalc(int &tn, double &T, Vector &gjmbH, double &g_J, Vector &ABC,       // Calculates the transition
 *                 char **sipffilename, ComplexMatrix &mat, float &delta,          //   matrix elements
 *                 ComplexMatrix &est)                                             //
 *   void estates(ComplexMatrix *est, Vector &gjmbheff, double *g_J, double &T,    // Calculates the energy and wavefunctions
 *                 Vector &ABC, char **sipffilename)                               //   "estates" matrix
 *   bool get_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,         // Calculates/loads the Q_q operators
 *                 std::vector<double> &Jvec)                                      //   for beyond dipole calculations
 *   void save_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l,        // Saves the Q_q operators to a file
 *                 std::vector<double> &Jvec)                                      //
 *   void mq(ComplexVector &Mq, double &th, double &ph, double &J0, double &J2,    // Calculates the thermal expectation
 *                 double &J4, double &J6, ComplexMatrix &est)                     //   of the magnetisation density
 *   int dncalc(int &tn, double &th, double &ph, double &J0, double &J2,           // Calculates the transition matrix
 *                 double &J4, double &J6, ComplexMatrix &est, double &T,          //   elements beyond the dipole
 *                 ComplexMatrix &mat)                                             //   approximation.
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include "vector.h"          // MatPack vector class
#include <fstream>
#include <ctime>
#define SMALL 1e-6   // must match SMALL in mcdisp.c and ionpars.cpp because it is used to decide wether for small
		     // transition, energy the matrix Mijkl contains wn-wn' or wn/kT

void myPrintMatrix(FILE * file,sMat<double> & M,int d)
{int i1,j1;
    fprintf (file,"Matrix\n");
   for (i1=0;i1<=d;++i1){
    for (j1=0;j1<=d;++j1) fprintf (file,"%6.3f ",M(i1,j1));
    fprintf (file,"\n");
    }
}    

// --------------------------------------------------------------------------------------------------------------- //
// Checks whether the Matpack matrix is the same as the c-array
// --------------------------------------------------------------------------------------------------------------- //
bool checkmat(ComplexMatrix &cmat, complexdouble *fmat,int r, int c)
{
   int i,j;
   for(i=0; i<(cmat.Rows()-r); i++)
      for(j=0; j<(cmat.Cols()-c); j++)
      { 
       //std::cout << "cmat["<<i+r<<"]["<<j+c<<"]=" << cmat[i+1][j+1] << "\tfmat["<<j<<"*" << cmat.Cols()-c << "+"<<i<<"]=";
       //   std::cout << fmat[j*(cmat.Cols()-c)+i].r << "+" << fmat[j*(cmat.Cols()-c)+i].i << "i\n";
       //std::cout << "cmat["<<j+c<<"]["<<i+r<<"]=" << cmat[j+c][i+r] << "\tfmat["<<j<<"*" << cmat.Cols()-c << "+"<<i<<"]=";
       //   std::cout << fmat[j*(cmat.Cols()-c)+i].r << "+" << fmat[j*(cmat.Cols()-c)+i].i << "i\n";
    //   if(real(cmat[i+r][j+c])!=fmat[j*(cmat.Cols()-c)+i].r || imag(cmat[i+r][j+c])!=fmat[j*(cmat.Cols()-c)+i].i) return false; }
         if(real(cmat[j+c][i+r])!=fmat[j*(cmat.Cols()-c)+i].r || imag(cmat[j+c][i+r])!=fmat[j*(cmat.Cols()-c)+i].i) return false; }
   return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the magnetic moment at a particular temperature and field
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           void mcalc(Vector &J,          // Output single ion momentum vector <Ja>,<Jb>,<Jc>, etc.
                      double *T,          // Input scalar temperature
                      Vector &gjmbH,      // Input vector of mean fields (meV) 
 /* Not Used */       double *gJ,         // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      double *lnZ,        // Output scalar logarithm of partition function
                      double *U,          // Output scalar internal energy 
                      ComplexMatrix &est) // Input/output eigenstate matrix (initialized in estates)                                          
{
   // Parses the input file for parameters
   icpars pars; 
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // Calculates the IC Hamiltonian matrix
   int i,k,q,Hsz=getdim(pars.n,pars.l);
   complexdouble *H=0,*Jm=0;
   bool Hicnotcalc = false;
   std::vector<double> parval; parval.reserve(35);
   if(pars.n==1) { parval.push_back(pars.xi); for(k=2; k<=(2*pars.l); k+=2) for(q=-k; q<=k; q++) parval.push_back(pars.B(k,q)); }
   else {
      for(i=0; i<4; i++) parval.push_back(pars.F[i]); parval.push_back(pars.xi); for(i=0; i<3; i++) parval.push_back(pars.alpha[i]);
      for(k=2; k<=(2*pars.l); k+=2) for(q=-k; q<=k; q++) parval.push_back(pars.B(k,q)); }
   if(parval.size()%2==1) parval.push_back(0.);

   if(est.Cols()!=(Hsz+1) || est.Rows()!=(Hsz+1)) Hicnotcalc = true;
   else if(real(est[0][0])==-0.1 && imag(est[0][0])==-0.1)  // Hic previously calculated
   {
      for(i=0; i<(int)(parval.size()/2); i++) if(real(est[0][i+1])!=parval[2*i] && imag(est[0][i+1])!=parval[2*i+1]) { Hicnotcalc = true; break; }
   }
   else Hicnotcalc = true;
   if(Hicnotcalc)
   {
      for(i=0; i<(int)(parval.size()/2); i++) est[0][i+1] = complex<double> (parval[2*i],parval[2*i+1]);
      if(pars.spectrelevels==-1)
      {
         sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); Hic/=MEV2CM; iHic/=MEV2CM; H = zmat2f(Hic,iHic);
//       est = ComplexMatrix(0,Hsz,0,Hsz); I comment this out - you should not reinitialize est !!!
         if(est.Rhi()!=Hsz||est.Chi()!=Hsz)
         {
            std::cerr << "ERROR module ic1ion - mcalc: Hsz recalculation does not agree with eigenstates matrix dimension\n"; exit(EXIT_FAILURE);
         }
         est[0][0] = complex<double> (-0.1,-0.1);
         if ((int)(parval.size()/2)>Hsz)
         {
            std::cerr << "ERROR module ic1ion - mcalc: storing Hamiltonian in est matrix not possible - parval size too big\n"; exit(EXIT_FAILURE);
         } 
         for(i=1; i<=Hsz; i++) memcpy(&est[i][1],&H[(i-1)*Hsz],Hsz*sizeof(complexdouble));
         free(H);
      }
      else  // Calculates using the Spectre method...                                       // lvl = 10;     % Number of |LSJ> levels to keep
      {
         complexdouble *Vcso, *Hrot, *zmt, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0; double *Ecso;
         sMat<double> Hcso = ic_Hcso(pars); Hcso/=MEV2CM; int Hcso_sz = Hcso.nc();          // Hf = full(ic_freehmltn(conf,racah_FtoF_k(F),xi,alpha));
         Ecso = new double[Hcso_sz]; Vcso = new complexdouble[Hcso_sz*Hcso_sz];
         sMat<double> zeros(Hcso_sz,Hcso_sz); int j = ic_diag(Hcso,zeros,Vcso,Ecso);        // [Vcso,Ecso] = eig(Hf);
         if(j!=0) { std::cerr << "error: specre\n"; return; }
         for(i=0; i<Hcso_sz; i++) for(j=0; j<Hcso_sz; j++) 
            if(fabs(Vcso[i+j*Hcso_sz].r)<DBL_EPSILON*50) Vcso[i+j*Hcso_sz].r=0.;
         char notranspose = 'N', transpose = 'C'; int cb = pars.spectrelevels; 
         complexdouble *zHcso = zmat2f(Hcso,zeros);
         zmt = new complexdouble[Hsz*cb]; Hrot = new complexdouble[cb*cb];
         F77NAME(zgemm)(&notranspose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,zHcso,      // Hrot = Vcso(:,1:lvl)'*Hcso*Vcso(:,1:lvl);
                           &Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
         F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,
                           &Hcso_sz,&zbeta,Hrot,&cb); free(zHcso); delete[]zmt;
         std::vector<int> J2r(pars.spectrelevels,0), Jz2r, id; 
         int imax,icv=0,mj; fconf conf(pars.n,0,pars.l); int incx = 1;
         std::vector< std::vector<int> > cvSO2CF; std::vector<int> cvEl;
         for(i=0; i<pars.spectrelevels; i++)                                                // for i=1:lvl; 
         {
            // IZAMAX returns fortran-style index (starts at 1)
            imax = F77NAME(izamax)(&Hcso_sz,&Vcso[i*Hcso_sz],&incx); 
            J2r[i] = conf.states[imax-1].J2;                                                //    J(i) = st_SO{find(abs(Vcso(:,i))==max(abs(Vcso(:,i))))}{5}; 
            cvEl.push_back(icv); icv += J2r[i]+1; cvEl.push_back(icv-1);                    //    cvSO2CF(i,1) = icv+1; icv=icv+(2*J(i)+1);
            cvSO2CF.push_back(cvEl); cvEl.clear();                                          //    cvSO2CF(i,2) = icv;
            for(mj=-J2r[i]; mj<=J2r[i]; mj+=2) { Jz2r.push_back(mj); id.push_back(i); }     //    for mj = 1:(2*J(i)+1); id(mj+icv) = i; Jz(mj+icv) = mj-J(i)-1; end
         }                                                                                  // end
         complexdouble *Hrmj = new complexdouble[icv*icv]; 
         convH2H(Hrot,Hrmj,cb,icv,cvSO2CF); delete[]Hrot;                                   // Hrot = full(convH2H(Hf,icv,cvSO2CF));
         // Calculates the CF rotated matrix, and expands to |LSJ,mJ>
         int k,iq,qp,qm; sMat<double> UJ; complexdouble *fUJ, *UJr; 
         double em,ep,icfact[]={0,0,-sqrt(15./7.)/2.,0,sqrt(11./14.),0,-sqrt(429./7.)/10.};
         for(k=2; k<=6; k+=2)
         {
            icfact[k] *= MEV2CM;
            UJ = racah_uJ(pars.n,k,pars.l); fUJ = zmat2f(UJ,zeros);                         // [U4,H,HZ,st_SO] = nshell_rme84(conf);
            UJr = new complexdouble[cb*cb]; zmt = new complexdouble[Hsz*cb];
            F77NAME(zgemm)(&notranspose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,fUJ,     // for l = 1:3; UCOEFF{l} = (RVEC(:,NKPT)'*U4{l}*RVEC(:,NKPT)) ./ icfact(l); end
                             &Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
            F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,
                             &Hcso_sz,&zbeta,UJr,&cb); free(fUJ); delete[]zmt;              // for k = 2:2:6
            if(pars.B(k,0)!=0)                                                              //    for iq = find(B{k/2})
            {                                                                               //        Uj = zeros(icv); q = iq-1-k;
               for(i=0; i<icv; i++) for(j=0; j<icv; j++)                                    //        for i = 1:icv; for j = 1:icv
               {
                  Hrmj[i+j*icv].r += pow(-1.,(J2r[id[i]]-Jz2r[i])/2.)                       //           Uj(i,j) = (-1)^(J(id(i))-Jz(i)-q)  
                                     * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],0,Jz2r[j]) //                       * threej([J(id(i)) k J(id(j)); -Jz(i) q Jz(j)])
                                     * UJr[id[i]+id[j]*cb].r * pars.B(k,0)/icfact[k];       //                       * UCOEFF{k/2}(id(i),id(j));
               }                                                                            //        end; end
            }                                                                               //        Hrmj = Hrmj + Uj.*B{k/2}(iq);      % Note have to give -q and +q pars!
            for(iq=0; iq<k; iq++) if(pars.B(k,iq-k)!=0 || pars.B(k,-(iq-k))!=0)             //    end
            {                                                                               // end
               qm = iq-k; qp = -qm;
               for(i=0; i<icv; i++) for(j=0; j<icv; j++)
               {
                  ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qp,Jz2r[j]) * UJr[id[i]+id[j]*cb].r;
                  em = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qm,Jz2r[j]) * UJr[id[i]+id[j]*cb].r;
                  if(pars.B(k,qp)!=0) { if(pars.n>(2*pars.l+1)) Hrmj[i+j*icv].r -= (ep + em*pow(-1.,qp)) * (pars.B(k,qp)/icfact[k]);
                                        else                    Hrmj[i+j*icv].r += (ep + em*pow(-1.,qp)) * (pars.B(k,qp)/icfact[k]); }
                  if(pars.B(k,qm)!=0) { if(pars.n>(2*pars.l+1)) Hrmj[i+j*icv].i -= (ep - em*pow(-1.,qm)) * (pars.B(k,qm)/icfact[k]);
                                        else                    Hrmj[i+j*icv].i += (ep - em*pow(-1.,qm)) * (pars.B(k,qm)/icfact[k]); }
               }
            }
         }
         est[0][0] = complex<double> (-0.1,-0.1); est[0][parval.size()/2] = complex<double> (icv,0.);
         for(i=1; i<=icv; i++) memcpy(&est[i][1],&Hrmj[(i-1)*icv],icv*sizeof(complexdouble)); delete[]Ecso; delete[]UJr; //delete[]Hrmj;
         // Calculates the matrix elements of the magnetic operators L and S
         complexdouble *Lrm,*Srm,*Lrj,*Srj; Lrm = new complexdouble[Hcso_sz*Hcso_sz]; Srm = new complexdouble[Hcso_sz*Hcso_sz]; int S2,L2,S2p,L2p,J2,J2p;
         for(i=0; i<Hcso_sz; i++) for(j=0; j<Hcso_sz; j++)
         {
            S2 = conf.states[i].S2; S2p = conf.states[j].S2; L2 = abs(conf.states[i].L)*2; L2p = abs(conf.states[j].L)*2; J2 = conf.states[i].J2; J2p = conf.states[j].J2;
            Lrm[i+j*Hcso_sz].i = 0.; Srm[i+j*Hcso_sz].i = 0.; Lrm[i+j*Hcso_sz].r = 0.; Srm[i+j*Hcso_sz].r = 0.;
            if(abs(J2-J2p)>2 || S2!=S2p || L2!=L2p || conf.states[i].v!=conf.states[j].v || conf.states[i].U!=conf.states[j].U) continue;
            Lrm[i+j*Hcso_sz].r = pow(-1.,(S2p+L2p+J2+2.)/2.)  * sixj(L2,J2,S2p,J2p,L2p,2) * sqrt((L2p+1.)*(J2+1.)*(J2p+1.)*((L2p/2.)*(L2p/2.+1)));
            Srm[i+j*Hcso_sz].r = pow(-1.,(S2p+L2p+J2p+2.)/2.) * sixj(S2,J2,L2p,J2p,S2p,2) * sqrt((S2p+1.)*(J2+1.)*(J2p+1.)*((S2p/2.)*(S2p/2.+1)));
         }
         zmt = new complexdouble[Hsz*cb]; Lrj = new complexdouble[cb*cb]; Srj = new complexdouble[cb*cb];
         F77NAME(zgemm)(&transpose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,Lrm,&Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
         F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,&Hcso_sz,&zbeta,Lrj,&cb); delete[]Lrm;
         F77NAME(zgemm)(&transpose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,Srm,&Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
         F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,&Hcso_sz,&zbeta,Srj,&cb); delete[]Srm; delete[]zmt;

         double sqrt2 = sqrt(2);
         char fSx[] = "results/ic1ion.m1"; std::fstream FSx; FSx.open(fSx, std::fstream::out); FSx.precision(24);
         char fLx[] = "results/ic1ion.m2"; std::fstream FLx; FLx.open(fLx, std::fstream::out); FLx.precision(24);
         char fSy[] = "results/ic1ion.m3"; std::fstream FSy; FSy.open(fSy, std::fstream::out); FSy.precision(24);
         char fLy[] = "results/ic1ion.m4"; std::fstream FLy; FLy.open(fLy, std::fstream::out); FLy.precision(24);
         char fSz[] = "results/ic1ion.m5"; std::fstream FSz; FSz.open(fSz, std::fstream::out); FSz.precision(24);
         char fLz[] = "results/ic1ion.m6"; std::fstream FLz; FLz.open(fLz, std::fstream::out); FLz.precision(24);
         for(i=0; i<icv; i++) for(j=0; j<icv; j++)
         {
            ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2,J2r[id[j]],-Jz2r[i], 2,Jz2r[j])/sqrt2;
            em = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2,J2r[id[j]],-Jz2r[i],-2,Jz2r[j])/sqrt2;
            FLx << (em-ep)*Lrj[id[i]+id[j]*cb].r << "\n"; FLy << (em+ep)*Lrj[id[i]+id[j]*cb].r << "\n";
            FSx << (em-ep)*Srj[id[i]+id[j]*cb].r << "\n"; FSy << (em+ep)*Srj[id[i]+id[j]*cb].r << "\n";
            ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2,J2r[id[j]],-Jz2r[i], 0,Jz2r[j]);
            FLz << ep*Lrj[id[i]+id[j]*cb].r << "\n"; FSz << ep*Srj[id[i]+id[j]*cb].r << "\n";
         }
         FSx.close(); FLx.close(); FSy.close(); FLy.close(); FSz.close(); FLz.close(); delete[]Lrj; delete[]Srj; delete[]Vcso;
      }
   }
   if(pars.spectrelevels==-1)
   {
      // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
      icmfmat mfmat(pars.n,pars.l,J.Hi()-J.Lo()+1);
      std::vector<double> vgjmbH((J.Hi()-J.Lo()+1),0.); for(i=J.Lo(); i<=J.Hi(); i++) vgjmbH[i-J.Lo()] = -gjmbH[i];
      sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH); 
      complex<double> a(1.,0.); int incx = 1;
      Jm = zmat2f(Jmat,iJmat); for(i=1; i<=Hsz; i++) F77NAME(zaxpy)(&Hsz,(complexdouble*)&a,(complexdouble*)&est[i][1],&incx,&Jm[(i-1)*Hsz],&incx);

      // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
      iceig VE; if(pars.partial) VE.lcalc(pars,Jm); else if(pars.arnoldi) VE.acalc(pars,Jm); else VE.calc(Hsz,Jm); free(Jm);

      // Calculates the expectation values sum_n{ <n|Ja|n> exp(-En/kT) }
      std::vector< std::vector<double> > matel; std::vector<double> vJ = mfmat.expJ(VE,*T,matel);
      for(i=J.Lo(); i<=J.Hi(); i++) J[i] = vJ[i-J.Lo()]; 
      *lnZ = vJ[i-1]; *U = vJ[i];
   }
   else
   {
      // Reloads the magnetic operator matrix in the rotated basis
      int icv = (int)real(est[0][parval.size()/2]); char fn[] = "results/ic1ion.mq"; Jm = new complexdouble[icv*icv];
      int Jhi = (J.Hi()>6) ? 6 : J.Hi(); memset(Jm,0,icv*icv*sizeof(complexdouble)); std::fstream FILEIN; 
      complexdouble *LS[6]; int iq,j; for(iq=0; iq<6; iq++) { LS[iq] = new complexdouble[icv*icv]; memset(LS[iq],0,icv*icv*sizeof(complexdouble)); }
      for(iq=J.Lo(); iq<=Jhi; iq++)
      {
         fn[16]=iq+48; FILEIN.open(fn,std::fstream::in); FILEIN.precision(24);  // In ascii, 48=="0"
         if(iq==3||iq==4) for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> LS[iq-1][i+j*icv].i; Jm[i+j*icv].i += LS[iq-1][i+j*icv].i*(-gjmbH[iq]); }
         else             for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> LS[iq-1][i+j*icv].r; Jm[i+j*icv].r += LS[iq-1][i+j*icv].r*(-gjmbH[iq]); }
         FILEIN.close();
      }
      complex<double> a(1.,0.); int incx = 1;
      for(i=1; i<=icv; i++) F77NAME(zaxpy)(&icv,(complexdouble*)&a,(complexdouble*)&est[i][1],&incx,&Jm[(i-1)*icv],&incx);
      // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
      iceig VE; VE.calc(icv,Jm); delete[]Jm; 
      // Calculates the moments, energy and partition function.
      std::vector<double> E(icv,0.),eb; int Esz; *U = 0.; double Z = 0.; char uplo = 'U';
      for(Esz=0; Esz<icv; Esz++) { E[Esz] = VE.E(Esz)-VE.E(0); eb.push_back(exp(-E[Esz]/(KB**T))); if(eb[Esz]<DBL_EPSILON) break; }
      complexdouble *zt = new complexdouble[icv]; complexdouble zalpha,zbeta,zme; zalpha.r=1.; zalpha.i=0.; zbeta.r=0.; zbeta.i=0.;
      for(iq=J.Lo(); iq<=Jhi; iq++) 
      { 
         J[iq] = 0.;
         for(i=0; i<Esz; i++)
         {
            F77NAME(zhemv)(&uplo, &icv, &zalpha, LS[iq-1], &icv, VE.zV(i), &incx, &zbeta, zt, &incx); 
            zme = F77NAME(zdotc)(&icv, VE.zV(i), &incx, zt, &incx);
            J[iq]+=zme.r*eb[i]; if(iq==J.Lo()) { Z+=eb[i]; *U+=(E[i]+VE.E(0))*eb[i]; }
         } 
         J[iq]/=Z;
      }
      for(i=0; i<6; i++) delete[]LS[i]; delete[]zt;
      *lnZ = log(Z)-VE.E(0)/(KB**T); *U/=Z;
   }

}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate transition matrix elements
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           int dmcalc(int &tn,            // Input transition number; if tn>0, print debug info
                      double &T,          // Input temperature
                      Vector &gjmbH,      // Input vector of mean fields (meV) 
 /* Not Used */       double &g_J,        // Input Lande g-factor
 /* Not Used */       Vector &ABC,        // Input vector of parameters from single ion property file
                      char **sipffilename,// Single ion properties filename
                      ComplexMatrix &mat, // Output M or Q matrix
                      float &delta,       // Output transition energy
                      ComplexMatrix &est) // Input eigenstate matrix (stored in estates)
                                          // Returns total number of transitions
{ 
   int i,j,k;

   // check if printout should be done and make tn positive
   int pr=1; if (tn<0) { pr=0; tn*=-1; }

   // Parses the input file for parameters
   icpars pars; const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
   int num_op = gjmbH.Hi()-gjmbH.Lo()+1; icmfmat mfmat(pars.n,pars.l,(num_op>6?num_op:6));

   // Copies the already calculated energy levels / wavefunctions from *est
   if(est.Rows()!=est.Cols()) { std::cerr << "dmcalc(): Input rows and columns of eigenstates matrix don't match.\n"; return 0; }
   int Hsz = est.Rows()-1; double *en = new double[Hsz]; for(i=0; i<Hsz; i++) en[i] = est[0][i+1].real();
   iceig VE(Hsz,en,(complexdouble*)&est[1][0],1);
 
   // Calculates the transition matrix elements:
   //    M_ab = <i|Ja|j><j|Jb|i> * (exp(-Ei/kT)-exp(-Ej/kT)) / Z    if delta > small
   //    M_ab = <i|Ja|j><j|Jb|i> * (exp(-Ei/kT)) / kTZ              if delta < small (quasielastic scattering)
   //    See file icpars.cpp, function mfmat::Mab() to see the actual code to calculate this.
   j=0;k=0; for(i=0; i<Hsz; ++i){ for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   sMat<double> Mab, iMab; mfmat.Mab(Mab,iMab,VE,T,i,j,pr,delta);

   for(i=1; i<=(num_op>6?num_op:6); i++)
      for(j=1; j<=(num_op>6?num_op:6); j++)
         mat(i,j) = complex<double> (Mab(i,j), iMab(i,j));

   return Hsz*(Hsz-1)/2;
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the eigenstates of Hic+effective_mean_fields
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
          void estates(ComplexMatrix *est,// Output Eigenstates matrix (row 0: real==Eigenvalues;imag==population)
                      Vector &gjmbheff,   // Input  Effective mean fields (meV)
 /* Not Used */       double *g_J,        // Input  Lande g-factor
                      double &T,          // Input  temperature
 /* Not Used */       Vector &ABC,        // Input  Vector of parameters from single ion property file
                      char **sipffilename)// Input  Single ion properties filename
{
   clock_t start,end; start = clock();

   // Parses the input file for parameters
   icpars pars;
   const char *filename = sipffilename[0];
   ic_parseinput(filename,pars);

   // Calculates the IC Hamiltonian matrix
   sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars); int Hsz = Hic.nr();
 
   // Calculates the mean field matrices <Sx>, <Lx>, etc. and the matrix sum_a(gjmbH_a*Ja)
   int num_op = gjmbheff.Hi()-gjmbheff.Lo()+1; icmfmat mfmat(pars.n,pars.l,(num_op>6?num_op:6));
   int i,j; std::vector<double> vgjmbH(6,0.); for(i=0; i<6; i++) vgjmbH[i] = gjmbheff[i+1];
   sMat<double> Jmat,iJmat; mfmat.Jmat(Jmat,iJmat,vgjmbH); 

   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
   Hic/=MEV2CM; Hic+=Jmat; if(!iHic.isempty()) iHic/=MEV2CM; if(!iJmat.isempty()) iHic+=iJmat; 
   iceig VE; if(iHic.isempty()) VE.calc(Hic); else VE.calc(Hic,iHic);

   // Initialises the output matrix
   (*est) = ComplexMatrix(0,Hsz,0,Hsz);

   // Stores the number of electrons and the orbital number in element (0,0)
   (*est)(0,0) = complex<double> (pars.n, pars.l);

   // Puts eigenvectors/values into the est matrix
   for(i=0; i<Hsz; i++) (*est)(0,i+1) = complex<double> (VE.E(i), exp(-(VE.E(i)-VE.E(0))/(KB*T)));   // Row 0

   if(VE.iscomplex())
      for(i=1; i<=Hsz; i++) memcpy(&(*est)[i][1],VE.zV(i-1),Hsz*sizeof(complexdouble)); // std::cout << "\nest==VE.zV = " << checkmat((*est),VE.zV(0),1,1) << endl; }
   else
      for(i=0; i<Hsz; i++){//printf("\n");
         for(j=0; j<Hsz; j++) 
            {(*est)(i+1,j+1) = complex<double> (VE.V(i)[j], 0.);
//                         printf("%6.3f %+6.3f i  ",VE.V(i)[j],0.0);
             if(VE.V(i)[j]!=VE.V(j,i)){fprintf(stderr,"compiler problem: bad memory mapping of vectors\n");exit(EXIT_FAILURE);}
             if(VE.V(i)[j]!=(*est)[i+1][j+1]){fprintf(stderr,"compiler problem: bad memory mapping of vectors\n");exit(EXIT_FAILURE);}
            }}

   end = clock(); std::cerr << "Time to do estates() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";

}

// --------------------------------------------------------------------------------------------------------------- //
// Loads a Q_q matrix from file if the file exists and has the same parameters n,l,Jvec
// --------------------------------------------------------------------------------------------------------------- //
bool get_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   int i, j, mn, r, c, sz, ml;
   std::vector<double> mJv(6,0.); 
   char filename[] = "results/mcphas.Qq"; filename[16]=q+120;               // 120==x, 121==y, 122==z
   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) return false;
   FILEIN >> mn >> ml; for(i=0; i<6; i++) FILEIN >> mJv[i]; FILEIN >> r >> c;
   if(mn!=n || ml!=(int)l) return false; for(i=0; i<6; i++) if(fabs(mJv[i]-Jvec[i])>1e-4) return false;
   Qq.clear(); sMat<double> emptymat(r,c); double Qt;
   for(i=0; i<6; i++) 
   {
      Qq.push_back(emptymat); FILEIN >> sz; for(j=0; j<sz; j++) {  FILEIN >> r >> c >> Qt; Qq[i](r,c) = Qt; }
   }
   return true;
}

// --------------------------------------------------------------------------------------------------------------- //
// Saves a Q_q matrix to a temporary file in the results/ directory
// --------------------------------------------------------------------------------------------------------------- //
void save_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   int i,j,sz;
   std::vector< std::vector<int> > nz;
   char filename[] = "results/mcphas.Qq"; filename[16]=q+120;               // 120==x, 121==y, 122==z
   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out);
   FILEOUT << n << " " << (int)l << " "; for(i=0; i<6; i++) FILEOUT << Jvec[i] << " "; FILEOUT << "\n";
   FILEOUT << Qq[0].nr() << " " << Qq[0].nc() << " ";
   for(i=0; i<6; i++)
   {
      nz = Qq[i].find(); sz = nz.size(); FILEOUT << sz << "\n";
      for(j=0; j<sz; j++) FILEOUT << nz[j][0] << " " << nz[j][1] << " " << Qq[i](nz[j][0],nz[j][1]) << "\n";
   }
   FILEOUT.close();
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the thermal expectation value of the FT of the magnetisation density -2Q in Bohr magnetons
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
          void mq(ComplexVector &Mq,      // Output expectation values -2[<Q>_{x} <Q>_y <Q>_z]
                  double &th, double &ph, // Input polar and azimuth angles theta and phi
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est)     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
{
   int i,q,n=1,Hsz=est.Cols()-1; orbital l;
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); l = (i==2) ? D : F;
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp, Qm; 
   std::vector< std::vector< sMat<double> > > Qmat; for(i=0; i<3; i++) Qmat.push_back(Qp);
   complexdouble *zQmat, *zt, zme, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
   double zMqr,zMqi,Z=0.;
   char trans = 'U'; int incx=1;

   Mq = ComplexVector(1,3);

   if(!get_Qq(Qmat[0],0,n,l,Jvec) || !get_Qq(Qmat[1],1,n,l,Jvec))            // Qmat[0]==Qx, Qmat[1]==Qy, Qmat[2]==Qz
   {
      lovesey_Qq(Qm,-1,n,l,Jvec); lovesey_Qq(Qp,1,n,l,Jvec);
      for(i=0; i<6; i++)  
      {
         Qmat[0].push_back( (Qp[i]-Qm[i]) * (-1/sqrt(2)) );                  // Qx = -1/sqrt(2) * (Q_{+1} - Q_{-1})
         if(i%2==0) Qmat[1].push_back( (Qp[i+1]+Qm[i+1]) * (-1/sqrt(2)) );   // real(Qy) = i^2/sqrt(2) * imag(Q_{+1}+Q_{-1})
         else       Qmat[1].push_back( (Qp[i-1]+Qm[i-1]) *  (1/sqrt(2)) );   // imag(Qy) = i/sqrt(2) * real(Q_{+1}+Q_{-1})
      }
      save_Qq(Qmat[0],0,n,l,Jvec); save_Qq(Qmat[1],1,n,l,Jvec); 
   }

   if(!get_Qq(Qmat[2],2,n,l,Jvec))                                           // Loads the Q_q matrix if prev. saved
   {
      lovesey_Qq(Qmat[2],0,n,l,Jvec);                                        // Calcs. scattering operator matrix Q_q
      save_Qq(Qmat[2],2,n,l,Jvec);
  // myPrintMatrix(stdout,Qmat[2][0],Hsz-1);
   }
   for(q=0; q<3; q++)
   {
      zQmat = zmat2f(Qmat[q][0],Qmat[q][1]);
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      zMqr = 0.; zMqi = 0.;
      for(i=1; i<=Hsz; i++)
      {
         F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
         zme = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx);
//         printf ("%i zme=%g %+g i  Ei=%6.3f ni=%6.3f \n",i,zme.r,zme.i,est[0][i].real(),est[0][i].imag());
         zMqr += (-2.)*zme.r*est[0][i].imag(); zMqi += (-2.)*zme.i*est[0][i].imag(); if(q==0) Z += est[0][i].imag();
      }
      free(zQmat); free(zt); Mq[q+1] = complex<double> (zMqr, zMqi)/Z;
   }
//   printf("MQ=(%g %+g i, %g %+g i,%g %+g i)\n",real(Mq(1)),imag(Mq(1)),real(Mq(2)),imag(Mq(2)),real(Mq(3)),imag(Mq(3)));
}

// --------------------------------------------------------------------------------------------------------------- //
// Routine to calculate the transition matrix using the scattering operator of Balcar and Lovesey
// --------------------------------------------------------------------------------------------------------------- //
extern "C"
#ifdef _WINDOWS
__declspec(dllexport)
#endif
           int                            // Returns total number of transitions
               dncalc(int &tn,            // Input transition number |tn|. If tn<0 omit printout. If tn>0 print info.
                  double &th,             // Input zenith angle (with the z==b axis) in radians.
                  double &ph,             // Input azimuth angle (with the x==a axis, to projection in x-y plane).
                  double &J0, double &J2, // Input radial parameters <j_0>, <j_2>
                  double &J4, double &J6, // Input radial parameters <j_4>, <j_6>
                  ComplexMatrix &est,     // Input eigenvalues/vectors of the system Hamiltonian, H_SI+H_mf 
                  double &T,              // Input temperature (K)
                  ComplexMatrix &mat)     // Output transition matrix, N(alpha,beta) = <-|Qa|+><+|Qb|-> (n- - n+)
/* 
     Note on Qalpha (Qa or Qb)
       if gJ>0:
        Kartesian components of the scattering operator Qalpha, alpha=1,2,3=a,b,c
        according to Lovesey Neutron Scattering equation 6.87b 
        scattering operator is given in  spherical coordinates Q-1,Q0,Q+1 (introduced
        as described above on input of th and ph) these are related to Qa,Qb,Qc by
        Q1=Qbx(axb)=Qy= i/sqrt(2)(Q+1 + Q-1) 
        Q2=Qb      =Qz= Q0                   
        Q3=Qaxb    =Qx=-1/sqrt(2)(Q+1 - Q-1)
                   
       if gJ=0 
        the orbital and spin contributions have to be given as separate components of Qalpha=1,2,3,4,5,6
        according to Lovesey Neutron Scattering equations 11.55 and 11.71 (the spin part 11.71 has to be
        divided by 2), i.e.
        <-|Q1,3,5|+>=<-|QSa,b,c|+>=
          =<-|sum_i exp(i k ri) s_(a,b,c)|+> /2                   as defined by 11.71 / 2
				   
        <-|Q2,4,6|+>=<-|QLa,b,c|+>=
          =<-|sum_i exp(i k ri) (-(k x grad_i)_(a,b,c)/|k|)|+>     as defined by 11.54 /(-|k|)
*/
{
   // check if printout should be done and make tn positive
   int pr=1; if (tn<0) { pr=0; tn*=-1; }

   int i,iJ,q,n=1,Hsz=est.Cols()-1; orbital l;//=D; find_nl_from_dim(Hsz,*&n,*&l,(complexdouble*)&est[1][1]);
   n = (int)est[0][0].real(); i = (int)est[0][0].imag(); l = (i==2) ? D : F;
   std::vector<double> E,Jvec(6,0.); Jvec[0]=th; Jvec[1]=ph; Jvec[2]=J0; Jvec[3]=J2; Jvec[4]=J4; Jvec[5]=J6;
   std::vector< sMat<double> > Qp1, Qm1; 
   std::vector< std::vector< sMat<double> > > Qq; for(i=0; i<3; i++) Qq.push_back(Qp1);
   complexdouble *zQmat, *zt, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
   std::vector<complexdouble> zij(7,zbeta), zji(7,zbeta);
   double Z=0., therm;
   char trans = 'U'; int incx=1;

   // Calculates the scattering operator, Q.
   if(!get_Qq(Qq[0],0,n,l,Jvec) || !get_Qq(Qq[1],1,n,l,Jvec))                // Qq[0]==Qx, Qq[1]==Qy, Qq[2]==Qz
   {
      lovesey_Qq(Qm1,-1,n,l,Jvec); lovesey_Qq(Qp1,1,n,l,Jvec);
      for(i=0; i<6; i++)  
      {
         Qq[0].push_back( (Qp1[i]-Qm1[i]) * (-1.0/sqrt(2.0)) );                  // Qx = -1/sqrt(2) * (Q_{+1} - Q_{-1})
         if(i%2==0) Qq[1].push_back( (Qp1[i+1]+Qm1[i+1]) * (-1.0/sqrt(2.0)) );   // real(Qy) = i^2/sqrt(2) * imag(Q_{+1}+Q_{-1})
         else       Qq[1].push_back( (Qp1[i-1]+Qm1[i-1]) *  (1.0/sqrt(2.0)) );   // imag(Qy) = i/sqrt(2) * real(Q_{+1}+Q_{-1})
      }
      save_Qq(Qq[0],0,n,l,Jvec); save_Qq(Qq[1],1,n,l,Jvec); 
     // myPrintMatrix(stdout,Qq[0][2],Hsz-1);
     // myPrintMatrix(stdout,Qq[0][3],Hsz-1);
     // myPrintMatrix(stdout,Qq[0][4],Hsz-1);
    //  myPrintMatrix(stdout,Qq[0][5],Hsz-1);
   }
   if(!get_Qq(Qq[2],2,n,l,Jvec))                                             // Loads the Q_q matrix if prev. saved
   {
      lovesey_Qq(Qq[2],0,n,l,Jvec);                                          // Calcs. scattering operator matrix Q_q
      save_Qq(Qq[2],2,n,l,Jvec);
   }

   for(i=0; i<Hsz; i++) { 
      therm = exp(-(est[0][i+1].real()-est[0][1].real())/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   int a,b,j=0,k=0; for(i=0; i<Hsz; ++i) {for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   ++i;++j; // because in est i and j start from 1...Hsz
 
   for(q=0; q<3; q++)
   {  zQmat = zmat2f(Qq[q][2],Qq[q][3]);     // Spin part
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[j][1], &incx, &zbeta, zt, &incx);
      zij[2*q+1] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx) ;
//         int k;for(k=0;k<Hsz;++k)printf("%6.3f %+6.3f i  ",real(est(j,1+k)),imag(est(j,1+k)));
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
      zji[2*q+1] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx) ;

      if(i==j)                               //subtract thermal expectation value from zij=zii
      {
         complexdouble expQ; double thexp=0;
         for(iJ=1;iJ<=Hsz;++iJ)    
         {
            therm = exp(-(est[0][iJ].real()-est[0][1].real())/(KB*T)); if(therm<DBL_EPSILON) break; 
            F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[iJ][1], &incx, &zbeta, zt, &incx);
            expQ = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
            thexp += expQ.r * therm / Z;
         }
         zij[2*q+1].r-=thexp;zji[2*q+1].r-=thexp;
      }
      free(zQmat); free(zt);

      zQmat = zmat2f(Qq[q][4],Qq[q][5]);     // orbital part
      zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[j][1], &incx, &zbeta, zt, &incx);
      zij[2*q+2] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[i][1], &incx, zt, &incx);
      F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[i][1], &incx, &zbeta, zt, &incx);
      zji[2*q+2] = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[j][1], &incx, zt, &incx);
      if(i==j)                               //subtract thermal expectation value from zij=zii
      {
         complexdouble expQ;double thexp=0;
         for(iJ=1;iJ<=Hsz;++iJ)    
         {
            therm = exp(-(est[0][iJ].real()-est[0][1].real())/(KB*T)); if(therm<DBL_EPSILON) break; 
            F77NAME(zhemv)(&trans, &Hsz, &zalpha, zQmat, &Hsz, (complexdouble*)&est[iJ][1], &incx, &zbeta, zt, &incx);
            expQ = F77NAME(zdotc)(&Hsz, (complexdouble*)&est[iJ][1], &incx, zt, &incx);
            thexp += expQ.r * therm / Z;
         }
         zij[2*q+2].r-=thexp;zji[2*q+2].r-=thexp;
      }
      free(zQmat); free(zt);
   }

   // check if zij are complex conjugate
   for(iJ=1;iJ<=6;++iJ)
      if(fabs(zij[iJ].i+zji[iJ].i)>SMALL) { std::cerr << "ERROR module ic1ion - dncalc: <i|Qalpha|j>not hermitian\n"; exit(EXIT_FAILURE); }

                
   complex<double> im(0,1);
   ComplexVector iQalphaj(1,6);
   
   for(a=1; a<=6; a++){iQalphaj(a) = complex<double> (zij[a].r,zij[a].i);if(a%2==1){iQalphaj(a)*=0.5;}} 
                                                                         // divide spin part by 2
   mat = ComplexMatrix(1,6,1,6);
   for(a=1; a<=6; a++)
      for(b=1; b<=6; b++)
         mat(a,b) = iQalphaj(a)*conj(iQalphaj(b));

   double delta;
   delta = est[0][j].real()-est[0][i].real();
   if(delta<-0.000001) { std::cerr << "ERROR module ic1ion - dncalc: energy gain delta gets negative\n"; exit(EXIT_FAILURE); }

   if(j==i) delta = -SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(est[0][i].real()-est[0][1].real())/(KB*T)) - exp(-(est[0][j].real()-est[0][1].real())/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV",i,j,delta);
         printf(" |<%i|Qa|%i>|^2=%6.3f |<%i|Qb|%i>|^2=%6.3f |<%i|Qc|%i>|^2=%6.3f",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
         printf(" |<%i|Qd|%i>|^2=%6.3f |<%i|Qe|%i>|^2=%6.3f |<%i|Qf|%i>|^2=%6.3f",i,j,real(mat(4,4)),i,j,real(mat(5,5)),i,j,real(mat(6,6)));
         printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
      }
   }
   else
   {
      therm = exp(-(est[0][i].real()-est[0][1].real())/(KB*T))/(KB*T);
          // quasielastic scattering has not wi-wj but wj*epsilon/kT
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV",i,j,delta);
         printf(" |<%i|Qa-<Qa>|%i>|^2=%6.3f |<%i|Qb-<Qb>|%i>|^2=%6.3f |<%i|Qc-<Qc>|%i>|^2=%6.3f",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
         printf(" |<%i|Qd-<Qd>|%i>|^2=%6.3f |<%i|Qe-<Qe>|%i>|^2=%6.3f |<%i|Qf-<Qf>|%i>|^2=%6.3f",i,j,real(mat(4,4)),i,j,real(mat(5,5)),i,j,real(mat(6,6)));
         printf(" n%i=%6.3f\n",i,therm/Z);
      }
   }
         mat *= therm / Z;

   return Hsz*(Hsz-1)/2;
}