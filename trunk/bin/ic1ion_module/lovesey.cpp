/* lovesey.cpp
 *
 * Function:  bool lovesey_akk(aKK,K,Kp,n,l)                    //  Calculates the a(K,K') matrix
 *            bool lovesey_ckk(cKK,K,Kp,n,l)                    //  Calculates the c(K,K') matrix
 *            void lovesey_Qq(Qq,iQq,q,n,l,Jvec)                //  Calculates the scattering operator Q_q
 *            complexdouble spherical_harmonics(k,q,theta,phi)  //  Calculates Yq^k(theta,phi)
 *
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2009 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 * Ref: E. Balcar and S.W. Lovesey, Theory of Magnetic Neutron and Photon Scattering, OUP, 1989
 * The formulae for the matrices are reproduced below, but may also be obtained from equations 9.8-9.16 of
 * Lovesey and Rimmer, Rep. Prog. Phys. v32, p333-394, 1969
 *
 */

#include "ic1ion.hpp"
#include <fstream>

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the a(K,K') matrix - the routine returns a true if the matrix is real, and false if it is imaginary
//    The matrix is returned in the argument aKK
// --------------------------------------------------------------------------------------------------------------- //
bool lovesey_aKK(sMat<double> &aKK, int K, int Kp, int n, orbital l)
{
   // Loads a previously save matrix if it exists
   char nstr[6]; char filename[255]; char basename[255], *mcphasedir = getenv("MCPHASE_DIR");
   if(mcphasedir==NULL) strcpy(basename,"mms/"); else { strcpy(basename,mcphasedir); strcat(basename,"/bin/ic1ion_module/mms/"); }
   nstr[0] = (l==F?102:100); if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 65;   // 65 is ASCII for "A", 100=="d" and 102=="f"
   nstr[1] = K+48; nstr[2] = Kp+48; nstr[3] = 0; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
   aKK = mm_gin(filename); if(!aKK.isempty()) return 1;

   bool df;
   if(l==D) df = true;
   else if(l==F) df = false;
   else {  std::cerr << "racah_mumat(): Only d- and f- configurations are implemented.\n"; return false; }
   int np = (n==1)?n:(n-1); 
   fconf confp(np,l);
   fconf conf(n,l);
   int num_states = (int)conf.states.size(), ns=0;
   int i,j,k,kk,v,vp,isz,jsz,iJ,jJ;
   int j2min,j2max,j2pmin,j2pmax;
   int L2,L2p,S2,S2p,J2,J2p,L2b,S2b;
   std::vector<int> indexJstart,indexJstop;
   std::vector<cfpls> cfpsi,cfpsj;
   double rmLS,sumcfp;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2) ns++;
      indexJstop.push_back(ns);
   }

   aKK.zero(ns,ns);

   // Selection rules on K.
   if(K<0  || K>(2*abs(l))   || (K%2)!=0)      return 1;        // Eqn 4.2.2
   if(Kp<1 || Kp>=(2*abs(l)) || ((Kp+1)%2)!=0) return 1;        // Eqn 4.2.3
   if(K!=(Kp+1) && K!=(Kp-1))                  return 1;        // Eqn 4.2.4

   // Calculates the matrix elements:
   // Eqn. 3.6.8 (d==delta)
   //
   //                                  K'+1           J'+S+L+L'+l  1/2    2                  1/2
   // A(K,K') = ( <j    > + <j    > ) i      d    (-1)            2    [l]  [L,L',K',K',K,J']
   //               K'-1      K'+1            SS'
   //
   //              \/   ( 1 K K' )  { 1  1  1  }  { J' K' J  } 
   //              /\   ( 0 0 0  )  { K' K  K' }  { L  S  L' }  A(K',K',l)
   //                                            _
   //                     ---     _   _          L                                _
   //              \/   n >   (t{|t) (t|}t') (-1)  { L' K' L }                Lb==L
   //              /\     ---                      { l  Lb l }
   //                      t

   // Where: (Eqn 3.4.9)
   //                                    1/2
   //                 l+1 [ (2l+3)(l+1) ]    ( l K' l+1 )  { 1  K'  K }
   // A(K,K',l) = (-1)    [ ----------- ]    ( 0 0  0   )  { l l+1  l }
   //                     [   (2l+1)    ]
   //
   
   // Calculate the coefficient A(K,K',l)
   double AKKl = pow(-1.,l+1.) * sqrt( (2.*l+3.)*(l+1.)/(2.*l+1.) ) * threej(2*l,2*Kp,2*(l+1),0,0,0) * sixj(2,2*Kp,2*Kp,2*l,2*(l+1),2*l);

   // Calculate the none state dependent part of the matrix elements
   double redmat = sqrt(2) * (2*l+1)*(2*l+1) * (2*Kp+1.) * sqrt(2*K+1.) * threej(2,2*K,2*Kp,0,0,0) * sixj(2,2,2,2*Kp,2*K,2*Kp) * AKKl;

   // Determines the factor i^(K'-1) - because K' is always odd, matrix is always real
   if((Kp-1)%4==0) redmat = -redmat;

   // Calculates the matrix elements at particular |J,M>, |J',M'>
   for(i=0; i<num_states; i++)
   {
      L2 = abs(conf.states[i].L)*2; S2 = conf.states[i].S2; v = conf.states[i].v;
      if(n!=1) 
      {
         if(df) cfpsi = racah_parents(n,v,S2,conf.states[i].L); 
         else   cfpsi = racah_parents(n,v,conf.states[i].U,S2,conf.states[i].L);
      }

      for(j=0; j<num_states; j++)
      {
         L2p = abs(conf.states[j].L)*2; S2p = conf.states[j].S2; vp = conf.states[j].v;

         if(S2!=S2p) continue;                                  // delta_SS'
         if(abs(L2-L2p)>(2*Kp)) continue;                       // Triangular condition on 6j symbol in 3.6.9

	 if(n==1) 
            rmLS = -redmat;
         else
         {
            if(df) cfpsj = racah_parents(n,vp,S2p,conf.states[j].L); 
            else   cfpsj = racah_parents(n,vp,conf.states[j].U,S2p,conf.states[j].L);

            sumcfp = 0.; isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(k=0; k<isz; k++)
               for(kk=0; kk<jsz; kk++)
                  if(cfpsi[k].ind==cfpsj[kk].ind) {
                     L2b = 2*abs(confp.states[cfpsi[k].ind].L); S2b = confp.states[cfpsi[k].ind].S2;
                     sumcfp += cfpsi[k].cfp*cfpsj[kk].cfp * pow(-1.,L2b/2) * sixj(L2p,2*Kp,L2,2*l,L2b,2*l); }
            rmLS = sqrt( (L2+1.)*(L2p+1.) ) * n * sumcfp * redmat;
         }

         // Caculate the J-dependent reduced matrix elements 
         j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2p); j2pmax = L2p+S2p;
         iJ=indexJstart[i]-1; jJ=indexJstart[j]-1;
         for(J2=j2min; J2<=j2max; J2+=2)
         {
            for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
            {
               if(abs(J2-J2p)<=(2*Kp))                          // Triangular condition on 6j symbol in 3.6.9
                  aKK(iJ,jJ) = pow(-1.,(J2p+S2+L2+L2p)/2.+l) * sqrt(J2p+1.) * sixj(J2p,2*Kp,J2,L2,S2,L2p) * rmLS;
               jJ++;
            }
            iJ++; jJ=indexJstart[j]-1;
         }
      }
   }

   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out); FILEOUT.close();
   rmzeros(aKK); mm_gout(aKK,filename); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the c(K,K') matrix, returning true if matrix is real, false if imaginary
// --------------------------------------------------------------------------------------------------------------- //
bool lovesey_cKK(sMat<double> &cKK, int K, int Kp, int n, orbital l)
{
   // Loads a previously save matrix if it exists
   char nstr[6]; char filename[255]; char basename[255], *mcphasedir = getenv("MCPHASE_DIR");
   if(mcphasedir==NULL) strcpy(basename,"mms/"); else { strcpy(basename,mcphasedir); strcat(basename,"/bin/ic1ion_module/mms/"); }
   nstr[0] = (l==F?102:100); if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 67;   // 67 is ASCII for "C", 100=="d" and 102=="f"
   nstr[1] = K+48; nstr[2] = Kp+48; nstr[3] = 0; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
   cKK = mm_gin(filename); if(!cKK.isempty()) return 1;

   bool df;
   if(l==D) df = true;
   else if(l==F) df = false;
   else { std::cerr << "racah_mumat(): Only d- and f- configurations are implemented.\n"; return false; }
   int np = (n==1)?n:(n-1);
   fconf confp(np,l);
   fconf conf(n,l);
   int num_states = (int)conf.states.size(), ns=0;
   int i,j,k,kk,j2min,j2max,j2pmin,j2pmax;
   int L2,L2p,S2,S2p,J2,J2p;
   std::vector<cfpls> cfpsi,cfpsj;
   std::vector<int> indexJstart,indexJstop,irm;
   sMat<double> mJmat_i(1,1);
   sMat<double> rmJ;
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   double rmLS,sumcfp;
   int v,vp,isz,jsz,L2b,S2b,iJ,jJ; 

   // Eqn. 3.6.11
   //
   //                  K      1/2                            1/2      1/2+S'+L'
   // C(K,K') = <j >  i  (1/2)     [l,l,S,S',L,L',J',K,K',K']     (-1)
   //             K                                                       _ _
   //                              { 1  K  K' }    ---     _   _          S+L
   //              \/   ( l K l )  { S' L' J' }  n >   (t{|t) (t|}t') (-1)    { S  1  S' }  { L  K  L' }
   //              /\   ( 0 0 0 )  { S  L  J  }    ---                        { s  Sb s  }  { l  Lb l  }
   //                                               t

   // Calculates part of the prefactor which only depends on K,K',l
   double redmat = sqrt(1./2) * (2*l+1.) * sqrt(2*K+1.) * (2*Kp+1.) * threej(2*l,2*K,2*l,0,0,0); 
   if(fabs(redmat)<DBL_EPSILON) return 1;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2) { ns++; irm.push_back(i); }
      indexJstop.push_back(ns);
   }

   cKK.zero(ns,ns);

   // Determines the factor i^K - because K is always even, matrix is always real
   if(K%4==2) redmat = -redmat;

   rmJ.zero(num_states,num_states);
   // Calculates the matrix elements at particular |J,M>, |J',M'>
   for(i=0; i<num_states; i++)
   {
      L2 = abs(conf.states[i].L)*2; S2 = conf.states[i].S2; v = conf.states[i].v;
      if(n!=1) 
      {
         if(df) cfpsi = racah_parents(n,v,S2,conf.states[i].L); 
         else   cfpsi = racah_parents(n,v,conf.states[i].U,S2,conf.states[i].L); 
      }

      for(j=0; j<num_states; j++)
      {
         L2p = abs(conf.states[j].L)*2; S2p = conf.states[j].S2; vp = conf.states[j].v;

         if(abs(S2-S2p)>2)     continue;                        // Triangular condition on 6j symbol in 3.6.11
         if(abs(L2-L2p)>(2*K)) continue;                        // Triangular condition on 6j symbol in 3.6.11

	 if(n==1)
	    rmLS = redmat;
	 else
	 {
            if(df) cfpsj = racah_parents(n,vp,S2p,conf.states[j].L); 
            else   cfpsj = racah_parents(n,vp,conf.states[j].U,S2p,conf.states[j].L);

            sumcfp = 0.; isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(k=0; k<isz; k++)
               for(kk=0; kk<jsz; kk++)
                  if(cfpsi[k].ind==cfpsj[kk].ind) {
                     L2b = 2*abs(confp.states[cfpsi[k].ind].L); S2b = confp.states[cfpsi[k].ind].S2;
                     sumcfp += cfpsi[k].cfp*cfpsj[kk].cfp * pow(-1.,(1+S2p+L2p+L2b+S2b)/2.) 
                                  * sixj(S2,2,S2p,1,S2b,1)*sixj(L2,2*K,L2p,2*l,L2b,2*l); }

            rmLS = sqrt( (S2+1.)*(S2p+1.)*(L2+1.)*(L2p+1.) ) * n * sumcfp * redmat;
         }

         // Caculate the J-dependent reduced matrix elements 
         j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2p); j2pmax = L2p+S2p;
         iJ=indexJstart[i]-1; jJ=indexJstart[j]-1;
         for(J2=j2min; J2<=j2max; J2+=2)
         {
            for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
            {
               if(abs(J2-J2p)<=(2*Kp))                          // Triangular condition on 9j symbol in 3.6.11
                  cKK(iJ,jJ) = sqrt(J2p+1.) * ninej(2,2*K,2*Kp,S2p,L2p,J2p,S2,L2,J2) * rmLS;
               jJ++;
            }
            iJ++; jJ=indexJstart[j]-1;
         }
      }
   }

   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out); FILEOUT.close();
   rmzeros(cKK); mm_gout(cKK,filename); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the spherical harmonic functions Y_q^k(theta,phi)
// --------------------------------------------------------------------------------------------------------------- //
complexdouble spherical_harmonics(int k, int q, double theta, double phi)
{
   complexdouble Ykq; Ykq.r=0; Ykq.i=0;
   double c=0.; double ct = cos(theta), st = sin(theta);
   if(k<0) k=-k;
   if(q<-k || q>k) return Ykq;
   if(k%2!=0) return Ykq;
   
#define SQ switch(q)
#define CQ(M,V) case M: c=V; break
#define SA switch(abs(q))
#define TQ(M,T) case M: c*=T; break
   switch(k)
   {
      case 0: Ykq.r = (1/2.)*sqrt(1./PI); break;
      case 2: SQ { CQ(-2,(1/4.)*sqrt(15/(2*PI))); CQ(-1,(1/2.)*sqrt(15/(2*PI))); CQ(0,(1/4.)*sqrt(5/PI));
                   CQ(1,(-1/2.)*sqrt(15/(2*PI))); CQ(2,(1/4.)*sqrt(15/(2*PI))); }
              SA { TQ(2,pow(st,2.)); TQ(1,st*ct); TQ(0,3*pow(ct,2.)-1.); } break;
      case 4: SQ { CQ(-4,(3/16.)*sqrt(35/(2*PI))); CQ(-3,(3/8.)*sqrt(35/PI)); CQ(-2,(3/8.)*sqrt(5/(2*PI)));
                   CQ(-1,(3/8.)*sqrt(5/PI)); CQ(0,(3/16.)*sqrt(1/PI)); CQ(1,(-3/8.)*sqrt(5/PI));
                   CQ(2,(3/8.)*sqrt(5/(2*PI))); CQ(3,(-3/8.)*sqrt(35/PI)); CQ(4,(3/16.)*sqrt(35/(2*PI))); }
              SA { TQ(4,pow(st,4.)); TQ(3,pow(st,3.)*ct); TQ(2,pow(st,2.)*(7*pow(ct,2.)-1.));
                   TQ(1,st*(7*pow(ct,3.)-3*ct)); TQ(0,35*pow(ct,4.)-30*pow(ct,2.)+3); } break;
      case 6: SQ { CQ(-6,(1/64.)*sqrt(3003/PI)); CQ(-5,(3/32.)*sqrt(1001/PI)); CQ(-4,(3/32.)*sqrt(91/(2*PI)));
                   CQ(-3,(1/32.)*sqrt(1365/PI)); CQ(-2,(1/64.)*sqrt(1365/PI)); CQ(-1,(1/16.)*sqrt(273/(2*PI)));
                   CQ(0,(1/32.)*sqrt(13/PI)); CQ(1,(-1/16.)*sqrt(273/(2*PI))); CQ(2,(1/64.)*sqrt(1365/PI));
                   CQ(3,(-1/32.)*sqrt(1365/PI)); CQ(4,(3/32.)*sqrt(91/(2*PI))); CQ(5,(-3/32.)*sqrt(1001/PI));
                   CQ(6,(1/64.)*sqrt(3003/PI)); }
              SA { TQ(6,pow(st,6.)); TQ(5,pow(st,5.)*ct); TQ(4,pow(st,4.)*(11*pow(ct,2.)-1.));
                   TQ(3,pow(st,3.)*(11*pow(ct,3.)-3*ct)); TQ(2,pow(st,2.)*(33*pow(ct,4.)-18*pow(ct,2.)+1));
                   TQ(1,st*(33*pow(ct,5.)-30*pow(ct,3.)+5*ct)); TQ(0,231*pow(ct,6.)-315*pow(ct,4.)+105*pow(ct,2.)-5); }
              break;
   }

   if(k>0) { Ykq.r = c*cos(q*phi); Ykq.i = c*sin(q*phi); }
   return Ykq;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the transition operator Q_q
// --------------------------------------------------------------------------------------------------------------- //
void lovesey_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   int i,j,ns,nn=(n>(4*abs(l)+2))?(4*abs(l)+2-n):n; 
   ns=1; for(i=(4*abs(l)+2-nn+1); i<=(4*abs(l)+2); i++) ns*=i; 
   j=1; for(i=n; i>1; i--) j*=i; ns/=j;                         // Number of states = ^{4l+2}C_{n}

   std::string errormsg("lovesey_Qq(): Unable to calculate the A(K,K') or B(K,K') matrix\n");
   double theta = Jvec[0], phi = Jvec[1], J[]={Jvec[2],0,Jvec[3],0,Jvec[4],0,Jvec[5]};
   int K,Kp,Q,Qp;
   sMat<double> Akk,Bkk,ckk,Qmat,QLmat,QSmat;
   Qmat.zero(ns,ns); Qq.clear(); for(i=0; i<6; i++) Qq.push_back(Qmat);

   bool iA,ic;
   complexdouble Ykq; double Tj,Tj2;

   int k,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int imJ_i,imJ_j,J2,J2p,Jz2,Jz2p;
   std::vector< std::vector<int> > nz;

   sMat<double> mJmat_i(1,1);
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   std::vector<int> iJst,iJen;
   fconf conf(n,0,l); minJ2=99; maxJ2=0; k=0;
   for(i=0; i<(int)conf.states.size(); i++) {
      j = conf.states[i].J2; if(j<minJ2) minJ2 = j; if(j>maxJ2) maxJ2 = j; iJst.push_back(k+1); k+=j+1; iJen.push_back(k); }

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // The scattering operator is given by: (eqn. 3.7.9 and 3.7.11 for the orbital and spin parts respectively)
   //
   //              ____ ---  {  K'-1 ^    2K'+1                                                  K'                      }
   // <i|Q |j> = \/4*pi >    { Y    (k) [ ------ ] [ A(K'-1,K') + B(K'-1,K') ] ( K'-1 K' 1  ) + Y  B(K',K') ( K' K' 1  ) }
   //     q             ---  {  Q    -     K'+1                                ( Q    Q' -q )    Q          ( Q  Q' -q ) }
   //                  K'QQ'
   //                              ( K' J' J )
   //                           *  ( Q' M' M )                   with |i>==|avUSLJM> and |j>==|a'v'U'S'L'J'M'> 
   //
   // The matrices A(K,K') and B(K,K') are:                                                              1/2 
   //                                                                K+2 [                        ( K+1 )            ]
   //  A(K,K') = ( <j    > + <j    > ) a(K,K')           B(K,K+1) = ---- [ <j > c(K,K+1) + <j   > ( --- ) c(K+2,K+1) ]
   //                K'-1      K'+1                                 2K+3 [   K               K+2  ( K+2 )            ]
   //
   //                                                                                                    1/2
   //                                                                K-1 [                        (  K  )            ]
   //  B(K,K) = <j > c(K,K)                              B(K,K-1) = ---- [ <j > c(K,K-1) + <j   > ( --- ) c(K-2,K-1) ]
   //             K                                                 2K-1 [   K               K-2  ( K-1 )            ]
   //
   //
   //  The 3j symbols with the matrices a(K,K') and c(K,K') restrict K to even integers up to 2l, and K' to odd integers
   //  less than 2l. I.e. for f-electrons, K=0,2,4,6 and K'=1,3,5. The summation over K' in the first case however is over
   //  all integers up to 2l, so that in the case of even K', the first term (with the A+B) is zero and in the case of odd
   //  K', the second term (with just the B term) is zero.
   //
   //  The matrices a(K,K') and c(K,K') are solely functions of |i> and |j> and do not depend on the bessel functions
   //  <j_k>, or spherical harmonics Y^K_Q, and are calculated in the functions lovesey_aKK() and loveset_cKK() above.

   for(Kp=0; Kp<=(2*l+1); Kp++)                                 // K may take even values between [0,2l]; K' odd values [1,(2l+1)]
   {
      if(Kp%2==1)                                               // Calculates the matrices for the first term with A+B
      {
         iA = lovesey_aKK(Akk,Kp-1,Kp,n,l); if(iA) Akk*= (J[Kp-1]+J[Kp+1])*((2*Kp+1)/(Kp+1.)); else { std::cerr << errormsg; return; }
         ic = lovesey_cKK(ckk,Kp-1,Kp,n,l); if(ic) Bkk = ckk*J[Kp-1];                          else { std::cerr << errormsg; return; }
         ic = lovesey_cKK(ckk,Kp+1,Kp,n,l); if(ic) Bkk+= ckk*(J[Kp+1]*sqrt(Kp/(Kp+1.)));       else { std::cerr << errormsg; return; }
         Bkk *= ((Kp+1.)/(2*Kp+1.)) * ( (2*Kp+1)/(Kp+1.) ); K = Kp-1; ckk = Bkk+Akk;
      }
      else                                                      // Calculates the matrices for the second term with just B
      {
         ic = lovesey_cKK(ckk,Kp,Kp,n,l); if(ic) ckk *= J[Kp]; else { std::cerr << errormsg; return; } K = Kp;
      }

      for(Q=-K; Q<=K; Q++)
      {
         Qp = -(Q-q); if(Qp<-Kp || Qp>Kp) continue;             // 3j symbol requires Q+Q'-q = 0
 
         Ykq = spherical_harmonics(K,Q,theta,phi); 
         Tj = pow(-1.,K-Kp+q) * sqrt(3) * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q);
         if((fabs(Ykq.r)<DBL_EPSILON && fabs(Ykq.i)<DBL_EPSILON) || fabs(Tj)<DBL_EPSILON) continue;

         Qmat.zero(ns,ns); QLmat.zero(ns,ns); QSmat.zero(ns,ns);

         nz = ckk.find();
         for(i=0; i<(int)nz.size(); i++)                        // The matrices above already contain the dependence on |vULSJ>
         {                                                      //    we now use the W-E theorem to add the Jz, Q, dependence
            J2 = conf.states[nz[i][0]].J2;
            J2p = conf.states[nz[i][1]].J2;
            valJ_i = (J2-minJ2)/2; valJ_j = (J2p-minJ2)/2;
            if(mJmat[valJ_i][valJ_j].isempty())
            {
               mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
               for(imJ_i=0; imJ_i<=J2; imJ_i++)
               {
                  Jz2 = imJ_i*2-J2;
                  for(imJ_j=0; imJ_j<=J2p; imJ_j++)
                  {
                     Jz2p = imJ_j*2-J2p; Tj2 = threej(2*Kp,J2p,J2,2*Qp,Jz2p,-Jz2);
                     if(fabs(Tj2)>DBL_EPSILON) mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,Kp+(-J2p+Jz2)/2.) * sqrt(J2+1.) * Tj2;
                  }
               }
            }
            Qmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*ckk(nz[i][0],nz[i][1]));
            if(Kp%2==1)
            {
               QLmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Akk(nz[i][0],nz[i][1]));
               QSmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Bkk(nz[i][0],nz[i][1]));
            }
         }
         for(i=0; i<=valJ; i++) for(j=0; j<=valJ; j++) mJmat[i][j].clear();

         Qq[0] += Qmat*(Ykq.r*Tj); Qq[1] += Qmat*(Ykq.i*Tj);
         if(Kp%2==1) { Qq[2] += QSmat*(Ykq.r*Tj); Qq[3] += QSmat*(Ykq.i*Tj);  Qq[4] += QLmat*(Ykq.r*Tj); Qq[5] += QLmat*(Ykq.i*Tj); }
         else        { Qq[2] += Qmat*(Ykq.r*Tj); Qq[3] += Qmat*(Ykq.i*Tj); }
      }
   }

   for(i=0; i<6; i++) { rmzeros(Qq[i]); Qq[i] *= sqrt(4*PI); }
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile:
//    g++ -g states.cpp cfp.cpp njsyms.cpp so_cf.cpp maths.cpp mmio.cpp lovesey.cpp -llapack -lblas; ./a.out
/* --------------------------------------------------------------------------------------------------------------- //
int main(int argc, char *argv[])
{
   int n,K,Kp;
   orbital l=F;
   if(argc>1) n = atoi(argv[1]); else n = 2;
   if(argc>2) K = atoi(argv[2]); else K = 2;
   if(argc>3) Kp = atoi(argv[3]); else Kp = 1;

   clock_t start,end; start = clock();

   sMat<double> aKK,cKK;

   bool ia = lovesey_aKK(aKK,K,Kp,n,l);
   end = clock(); std::cerr << "Time to calculate a(K,K') = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   std::cerr << "a(" << K << "," << Kp << ") is "; if(ia) std::cerr << "real\n"; else std::cerr << "imaginary\n";
   std::cout << aKK << "\naKK=x; clear x;\n";

   start = clock();
   bool ic = lovesey_cKK(cKK,K,Kp,n,l);
   end = clock(); std::cerr << "Time to calculate c(K,K') = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   std::cerr << "c(" << K << "," << Kp << ") is "; if(ia) std::cerr << "real\n"; else std::cerr << "imaginary\n";
   std::cout << cKK << "\ncKK=x; clear x;\n";

   start = clock();
   std::vector<double> Jvec(6,1.);
   lovesey_Qq(aKK,cKK,0,n,l,Jvec);
   end = clock(); std::cerr << "Time to calculate Q = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   std::cout << aKK << "Qq=x; clear x;\n" << cKK << "iQq=x; clear x;\n";

   return 0;
}*/
