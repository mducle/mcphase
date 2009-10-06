/* spectre.cpp
 *
 * Routines to calculate the IC Hamiltonian using the method in the Spectra/XTAL programs of Hannah Crosswhite.
 * References: https://chmwls.chm.anl.gov/downloads/spectra/about/matrix.html
 *             Carnall, Crosswhite, Crosswhite, Conway, J. Chem. Phys., v64, p3582, 1976
 *
 * Functions:
 * void spectre_hmltn(pars, est, parvalsize)                         // Calculates the IC Hamiltonian
 * vector<double> spectre_expJ(pars,est,parvalsize,gjmbH,Jhi,Jlo,T)  // Calculates the expectation values
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2009 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include "vector.h"          // MatPack vector class
#include <fstream>

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the IC Hamiltonian so that's basis states are diagonal in the free-ion part (Hc+Hso)
// --------------------------------------------------------------------------------------------------------------- //
void spectre_hmltn(icpars &pars, ComplexMatrix &est, int parvalsize)
{
   int i,Hsz=getdim(pars.n,pars.l);

   complexdouble *Vcso, *Hrot, *zmt, zalpha, zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0; double *Ecso;
   sMat<double> Hcso = ic_Hcso(pars); Hcso/=MEV2CM; int Hcso_sz = Hcso.nc();          // Hf = full(ic_freehmltn(conf,racah_FtoF_k(F),xi,alpha));
   if (pars.spectrelevels > Hcso.nr()) pars.spectrelevels=Hcso.nr(); 
   Ecso = new double[Hcso_sz]; Vcso = new complexdouble[Hcso_sz*Hcso_sz];
   sMat<double> zeros(Hcso_sz,Hcso_sz); int j = ic_diag(Hcso,zeros,Vcso,Ecso);        // [Vcso,Ecso] = eig(Hf);
   if(j!=0) { std::cerr << "error: spectre\n"; return; }
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
   double em,ep,epe,eme,p = 1./( pow(-1.,(double)abs(pars.l))*(2.*pars.l+1.) );
   double icfact[] = {0,0,p/threej(2*pars.l,4,2*pars.l,0,0,0),0,p/threej(2*pars.l,8,2*pars.l,0,0,0),0,p/threej(2*pars.l,12,2*pars.l,0,0,0)};
   for(k=2; k<=6; k+=2)
   {
      icfact[k] *= MEV2CM;
      UJ = racah_uJ(pars.n,k,pars.l); fUJ = zmat2f(UJ,zeros);                         // [U4,H,HZ,st_SO] = nshell_rme84(conf);
      UJr = new complexdouble[cb*cb]; zmt = new complexdouble[Hsz*cb];
      F77NAME(zgemm)(&notranspose,&notranspose,&Hcso_sz,&cb,&Hcso_sz,&zalpha,fUJ,     // for l = 1:3; UCOEFF{l} = (RVEC(:,NKPT)'*U4{l}*RVEC(:,NKPT)) ./ icfact(l); end
                       &Hcso_sz,Vcso,&Hcso_sz,&zbeta,zmt,&Hcso_sz);
      F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&Hcso_sz,&zalpha,Vcso,&Hcso_sz,zmt,
                       &Hcso_sz,&zbeta,UJr,&cb); free(fUJ); delete[]zmt;              // for k = 2:2:6
/*    if(pars.B(k,0)!=0)                                                              //    for iq = find(B{k/2})
      {                                                                               //        Uj = zeros(icv); q = iq-1-k;
         for(i=0; i<icv; i++) for(j=0; j<icv; j++)                                    //        for i = 1:icv; for j = 1:icv
         {
            Hrmj[i+j*icv].r += pow(-1.,(J2r[id[i]]-Jz2r[i])/2.)                       //           Uj(i,j) = (-1)^(J(id(i))-Jz(i)-q)  
                               * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],0,Jz2r[j]) //                       * threej([J(id(i)) k J(id(j)); -Jz(i) q Jz(j)])
                               * UJr[id[i]+id[j]*cb].r * pars.B(k,0)/icfact[k];       //                       * UCOEFF{k/2}(id(i),id(j));
         }                                                                            //        end; end
      }                                                                               //        Hrmj = Hrmj + Uj.*B{k/2}(iq);      % Note have to give -q and +q pars!
*/    for(iq=0; iq<=k; iq++) // if(pars.B(k,iq-k)!=0 || pars.B(k,-(iq-k))!=0)         //    end
      {                                                                               // end
         qm = iq-k; qp = -qm;
         if(qm==0)
         {
            char fh[] = "results/ic1ion.uX0"; fh[16]=k+48; std::fstream FH; FH.open(fh, std::fstream::out); FH.precision(24);
            for(i=0; i<icv; i++) for(j=0; j<icv; j++)
            {
               em = ( pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qm,Jz2r[j]) * UJr[id[i]+id[j]*cb].r ) /icfact[k];
               Hrmj[i+j*icv].r += em * pars.B(k,qp); FH << em << "\n";
            }
            FH.close();
         }
         else
         {
            char fhm[] = "results/ic1ion.uXmX"; fhm[16]=k+48; fhm[18]=qp+48; std::fstream FM; FM.open(fhm, std::fstream::out); FM.precision(24);
            char fhp[] = "results/ic1ion.uXX";  fhp[16]=k+48; fhp[17]=qp+48; std::fstream FP; FP.open(fhp, std::fstream::out); FP.precision(24);
            for(i=0; i<icv; i++) for(j=0; j<icv; j++)
            {
               ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qp,Jz2r[j]) * UJr[id[i]+id[j]*cb].r;
               em = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2*k,J2r[id[j]],-Jz2r[i],2*qm,Jz2r[j]) * UJr[id[i]+id[j]*cb].r;
               epe = (ep + em*pow(-1.,qp)) / icfact[k]; eme = (ep - em*pow(-1.,qm)) / icfact[k];
               if(pars.n>(2*pars.l+1))
               {
                  if(fabs(eme)>DBL_EPSILON) { FM << -eme << "\n"; if(fabs(pars.B(k,qm))>DBL_EPSILON) Hrmj[i+j*icv].i -= eme * pars.B(k,qm); } else FM << 0 << "\n";
                  if(fabs(eme)>DBL_EPSILON) { FP << -epe << "\n"; if(fabs(pars.B(k,qp))>DBL_EPSILON) Hrmj[i+j*icv].r -= epe * pars.B(k,qp); } else FP << 0 << "\n";
               }
               else
               {
                  if(fabs(eme)>DBL_EPSILON) { FM <<  eme << "\n"; if(fabs(pars.B(k,qm))>DBL_EPSILON) Hrmj[i+j*icv].i += eme * pars.B(k,qm); } else FM << 0 << "\n";
                  if(fabs(eme)>DBL_EPSILON) { FP <<  epe << "\n"; if(fabs(pars.B(k,qp))>DBL_EPSILON) Hrmj[i+j*icv].r += epe * pars.B(k,qp); } else FP << 0 << "\n";
               }
            }
            FM.close(); FP.close();
         }
      }
   }
   est[0][0] = complex<double> (-0.1,-0.1); est[0][parvalsize/2+1] = complex<double> (icv,0.);
   for(i=1; i<=icv; i++) memcpy(&est[i][1],&Hrmj[(i-1)*icv],icv*sizeof(complexdouble)); delete[]Ecso; delete[]UJr; delete[]Hrmj;
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

   double sqrt2 = sqrt(2.),elem;
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
      elem = (em-ep)*Lrj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FLx << elem << "\n"; else FLx << 0 << "\n";
      elem = (em+ep)*Lrj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FLy << elem << "\n"; else FLy << 0 << "\n";
      elem = (em-ep)*Srj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FSx << elem << "\n"; else FSx << 0 << "\n";
      elem = (em+ep)*Srj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FSy << elem << "\n"; else FSy << 0 << "\n";
      ep = pow(-1.,(J2r[id[i]]-Jz2r[i])/2.) * threej(J2r[id[i]],2,J2r[id[j]],-Jz2r[i], 0,Jz2r[j]);
      elem = ep*Lrj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FLz << elem << "\n"; else FLz << 0 << "\n";
      elem = ep*Srj[id[i]+id[j]*cb].r; if(fabs(elem)>DBL_EPSILON) FSz << elem << "\n"; else FSz << 0 << "\n";
   }
   FSx.close(); FLx.close(); FSy.close(); FLy.close(); FSz.close(); FLz.close(); delete[]Lrj; delete[]Srj; delete[]Vcso;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the expectation values of the multipolar operators for the spectre-like rotated Hamiltonian
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> spectre_expJ(icpars &pars, ComplexMatrix &est, int parvalsize, Vector &gjmbH, int Jlo, int Jhi, double T)
{
   std::vector<double> J(Jhi+2,0.);
   // Reloads the magnetic operator matrix in the rotated basis
   int i,icv = (int)real(est[0][parvalsize/2+1]); char fn[] = "results/ic1ion.mq"; complexdouble *Jm = new complexdouble[icv*icv]; double elem;
 /*int Jhi = (J.Hi()>6) ? 6 : J.Hi();*/ memset(Jm,0,icv*icv*sizeof(complexdouble)); std::fstream FILEIN; 
   int k[] = {0,1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   complexdouble **LS = new complexdouble*[6]; 
   int iq,j; for(iq=0; iq<6; iq++) { LS[iq] = new complexdouble[icv*icv]; memset(LS[iq],0,icv*icv*sizeof(complexdouble)); }
   for(iq=Jlo; iq<=Jhi; iq++)
   {
      if(iq<=6)
      {
         fn[16]=iq+48; FILEIN.open(fn,std::fstream::in); FILEIN.precision(24);  // In ascii, 48=="0"
         if(iq==3||iq==4) for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> LS[iq-1][i+j*icv].i; Jm[i+j*icv].i += LS[iq-1][i+j*icv].i*(-gjmbH[iq]); }
         else             for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> LS[iq-1][i+j*icv].r; Jm[i+j*icv].r += LS[iq-1][i+j*icv].r*(-gjmbH[iq]); }
         FILEIN.close();
      }
      else
      {
         if(fabs(gjmbH[iq])>(DBL_EPSILON*1000))
         {
            char fh[] = "results/ic1ion.uXmX"; fh[16]=k[iq]+48; 
            if(q[iq]<0) 
            {
               fh[18]=abs(q[iq])+48; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24);
               for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].i += elem*(-gjmbH[iq]); } FILEIN.close();
            }
            else 
            { 
               fh[17]=q[iq]+48; fh[18]=0; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24);
               for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].r += elem*(-gjmbH[iq]); } FILEIN.close();
            }
         }
      }
   }
   complex<double> a(1.,0.); int incx = 1;
   for(i=1; i<=icv; i++) F77NAME(zaxpy)(&icv,(complexdouble*)&a,(complexdouble*)&est[i][1],&incx,&Jm[(i-1)*icv],&incx);
   // Diagonalises the Hamiltonian H = Hic + sum_a(gjmbH_a*Ja)
   iceig VE; VE.calc(icv,Jm);
   // Calculates the moments, energy and partition function.
   std::vector<double> E(icv,0.),eb; int Esz; double U = 0.; double Z = 0.; char uplo = 'U';
   for(Esz=0; Esz<icv; Esz++) { E[Esz] = VE.E(Esz)-VE.E(0); eb.push_back(exp(-E[Esz]/(KB*T))); if(eb[Esz]<DBL_EPSILON) break; }
   complexdouble *zt = new complexdouble[icv]; complexdouble zalpha,zbeta,zme; zalpha.r=1.; zalpha.i=0.; zbeta.r=0.; zbeta.i=0.;
   for(iq=Jlo; iq<=Jhi; iq++) 
   { 
      if((iq>6 && k[iq]%2==1) || (k[iq]>4 && pars.l==D)) continue;
      J[iq] = 0.;
      if(iq>6)
      { 
         memset(Jm,0,icv*icv*sizeof(complexdouble)); char fh[] = "results/ic1ion.uXmX"; fh[16]=k[iq]+48; 
         if(q[iq]<0) 
         {
            fh[18]=abs(q[iq])+48; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24); 
            for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].i = elem; } FILEIN.close();
         }
         else
         {
            fh[17]=q[iq]+48; fh[18]=0; FILEIN.open(fh, std::fstream::out); FILEIN.precision(24);
            for(i=0; i<icv; i++) for(j=0; j<icv; j++) { FILEIN >> elem; Jm[i+j*icv].r = elem; } FILEIN.close();
         }
      }
      for(i=0; i<Esz; i++)
      {
         if(iq<=6)
            F77NAME(zhemv)(&uplo, &icv, &zalpha, LS[iq-1], &icv, VE.zV(i), &incx, &zbeta, zt, &incx); 
         else
            F77NAME(zhemv)(&uplo, &icv, &zalpha, Jm, &icv, VE.zV(i), &incx, &zbeta, zt, &incx); 
#ifdef _G77
         F77NAME(zdotc)(&zme, &icv, VE.zV(i), &incx, zt, &incx);
#else
         zme = F77NAME(zdotc)(&icv, VE.zV(i), &incx, zt, &incx);
#endif
         J[iq]+=zme.r*eb[i]; if(iq==Jlo) { Z+=eb[i]; U+=(E[i]+VE.E(0))*eb[i]; }
      }
      J[iq]/=Z; if(fabs(J[iq])<DBL_EPSILON) J[iq]=0.;
   }
   for(i=0; i<6; i++) delete[]LS[i]; delete[]LS;
   delete[]zt; delete[]Jm; 
   J[iq] = log(Z)-VE.E(0)/(KB*T); J[iq+1] = U/Z;
   return J;
}
