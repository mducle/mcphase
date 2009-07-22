/* maths.cpp
 *
 * This file implements mathematical and matrix classes with similar syntax to Matlab, and may be used to
 * translate programs from Matlab to C++
 *
 * Class:					 		//
 * 								//
 * Function:							//
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 */

#include "maths.hpp"
#include<ctime>
#include<cstdlib>                      // For rand()

#ifdef USE_LAPACK
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the eigenvalues and eigenvectors of a symmetric matrix using the LAPACK routine dsyev
// --------------------------------------------------------------------------------------------------------------- //
eigVE<double> eig(const sMat<double> & M)
{
   eigVE<double> retval;
   if(!M.issymm()) return retval;     // Checks that matrix is symmetric first otherwise returns empty

   // We must pass all parameters to Fortran routines by reference. Hence we must first set up variables for each
   /* parameters. The following are parameters to work out a good block size for the algorithm using ILAENV
   int ispec = 1;          // Calculates optimum block size
   char *name = "SSYTRD";  // Calculates using the tri-diagonal reduction of a real symmetric matrix
   int unused = -1;        // For unused parameters N3, N4. We use opts==jobz, N1==lda, N2==n (see below) */
   
   // Paramaters for DSYEV itself
   int lda = M.nc();
   int n = M.nr();
   char jobz = 'V';        // Job type - 'N' is eigenvalues only; 'V' calculates eigenvectors as well
   char uplo = 'U';        // States whether to use the upper ('U') or lower triangle ('L')
   int info = 0;           // Success flag on return. 0 is success. <0 is wrong input. >0 means not converged
   double *work;           // A C-array for the workspace
   int lwork;              // The length of the workspace

   // Declaration for the C-array holding the matrix on input and the eigenvectors on output
   double *m = M.f_array();
   double *eigval = (double*)malloc(n*sizeof(double));

   /* Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   lwork = ( F77NAME(ilaenv)(&ispec, name, &jobz, &lda, &n, &unused, &unused) +2 ) * n; */
   lwork = 4*n;
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dsyev)(&jobz, &uplo, &n, m, &lda, eigval, work, &lwork, &info);

   if(info==0)
   {
      retval.E = f2vec(eigval,n);
      retval.V = f2mat(m,lda,n);
   }
   else
   {
      std::cerr << "eig(): DSYEV return error code " << info << "\n";
   }

   free(eigval); free(m); free(work);
   return retval;
}
#else   // USE_LAPACK
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the eigenvalues and eigenvectors of a matrix using tri-diagonalisation and QL-reduction
// --------------------------------------------------------------------------------------------------------------- //
eigVE<double> eig(const sMat<double> & M)
{
   eigVE<double> retval;
   if(!M.issymm()) return retval;     // Checks that matrix is symmetric first otherwise returns empty

   // Calculates the Householder tri-diagonal reduction of the matrix M, using the tred2 algorithm of Wilkinson,
   //    Numerische Mathematik vol 4, pp. 354-361 (1962), and converted from the Algol routine printed in Reinsch
   //    and Wilkinson (eds.), Handbook for Automatic Computation, vol 2: Linear Algebra (Springer-Verlag 1971).

   int i,j,k,l;
   int n = M.nr();
   double f,g,h,hh;
   double tol = DBL_EPSILON;		// z holds the eigenvectors after execution of imtql2. Before, its upper part is
   sMat<double> z = M;			//    the original matrix M, and lower tri. has information on the transformation.
   std::vector<double> d(n,0.),e(n,0.);	// d holds the diagonal and e the off-diagonal elements of the tri-diagonal matrix

   for(i=n-1; i>0; i--)
   {
      l = i-2; f = z(i,i-1); g = 0.;
      for(k=0; k <= l; k++)
         g += z(i,k)*z(i,k);
      h = g + f*f;
      // If g is too small for orthogonality to be guaranteed, the transformation is skipped
      if(g<=tol)
      {
         e[i] = f; h = 0.;
      }
      else
      {
         l++; g = (f>=0.) ? -sqrt(h) : sqrt(h); e[i] = g;
	 h -= f*g; z(i,i-1) = f-g; f = 0.;
	 for(j=0; j<=l; j++)
	 {
	    z(j,i) = z(i,j)/h; g = 0.;
	    // Form elements of A x u
            for(k=0; k<=j; k++)
               g += z(j,k)*z(i,k);
            for(k=j+1; k<=l; k++)
               g += z(k,j)*z(i,k);
            // Form elements of _p_
            e[j] = g/h; f += g*z(j,i);
         }
         // Form K
         hh = f/(h+h);
         // Form reduced A
         for(j=0; j<=l; j++)
	 {
            f = z(i,j); g = e[j] - hh*f; e[j] = g;
	    for(k=0; k<=j; k++)
               z(j,k) = z(j,k) - f*e[k] - g*z(i,k);
         }
      }
      d[i] = h;
   }

   d[0] = 0.; e[0] = 0.;
   // Accumulation of transformation matrices
   for(i=0; i<n; i++)
   {
      l = i-1;
      if(d[i]!=0)
         for(j=0; j<=l; j++)
         {
            g = 0.;
            for(k=0; k<=l; k++)
               g += z(i,k)*z(k,j);
            for(k=0; k<=l; k++)
               z(k,j) = z(k,j) - g*z(k,i);
         }
      d[i] = z(i,i); z(i,i) = 1.;
      for(j=0; j<=l; j++)
      {
         z.del(i,j); z.del(j,i);
      }
   }
 
   //std::cout << "Tridiagonal matrix, diagonal is: " << dispvect(d) << "; subdiagonal is " << dispvect(e) << "\n";

   // Calculates the eigenvalues and eigenvectors of the tridiagonal matrix represented in e[] and d[] by implicit QL
   //    factorisation, using imqtl2 algorithm of Dubrulle, Martin and Wilkinson, Num. Math. 12, 377-383 (1968), 
   //    which is a modification of the implicit QR factorisation algorithm of Francis, Comput. J. 4, 265-271 and
   //    332-345 (1961-1962), as converted from the Algol routine in the Handbook of Automatic Computation (1971).

   int ia,m,its;
   double c,p,q,s,t,u;
   bool doiter = true;
   bool failflag = false;

   for(i=1; i<n; i++)
      e[i-1] = e[i];
   e[n-1] = 0.; k = n-2;
   for(j=0; j<n; j++)
   {
      its = 0;
      while(doiter)
      {
         // Look for single small sub-diagonal element
         for(m=j; m<=k; m++)
            if( fabs(e[m]) <= (DBL_EPSILON*(fabs(d[m])+fabs(d[m+1]))) ) break;
       //m = n-1;        // If the loop finishes without condition satisfied, m=k+1 = n-1;
         if(m==k) break;
         u = d[j];
         if(m!=j)
         {
            if(its==100) { failflag = true; break; }
            its++;
            // Form shift
            q = (d[j+1]-u) / (2*e[j]); t = sqrt(1.+q*q);
            q = d[m] - u + e[j] / ( (q<0.) ? q-t : q+t );
            u = 0.; s = 1.; c = 1.;
            for(i=m-1; i>=j; i--)
            {
               p = s*e[i]; h = c*e[i];
               if(fabs(p)>=fabs(q))
               {
                  c = q/p; t = sqrt(c*c+1.);
                  e[i+1] = p*t; s = 1./t; c *= s;
               }
               else
               {
                  s = p/q; t = sqrt(s*s+1.);
                  e[i+1] = q*t; c = 1./t; s *= c;
               }
               q = d[i+1] - u; t = (d[i]-q)*s + 2.*c*h;
               u = s*t; d[i+1] = q+u; q = c*t-h;
               // Form vector
               for(ia=0; ia<n; ia++)
               {
                  p = z(ia,i+1);
                  z(ia,i+1) = s*z(ia,i) + c*p;
                  z(ia,i) = c*z(ia,i) - s*p;
               }
            }
            d[j] = d[j] - u; e[j] = q; e[m] = 0.;
         }
	 else { doiter = false; }
      }
      doiter = true;
      if(failflag) break;
   }

   if(!failflag)
   {  // Order eigenvalues and eigenvectors
      for(i=0; i<n; i++)
      {
         k = i; p = d[i];
         for(j=i+1; j<n; j++)
            if(d[j]<p)
            {
               k = j; p = d[j];
            }
         if(k!=i)
         {
	       d[k] = d[i]; d[i] = p;
            for(j=0; j<n; j++)
            {
               p = z(j,i); z(j,i) = z(j,k); z(j,k) = p;
            }
         }
      }

      // Puts results into structure
      retval.E = d;
      retval.V = z;

   }

   return retval;
}
#endif  // USE_LAPACK

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the eigenvalues and eigenvectors of a tridiagonal matrix using QL-reduction
// --------------------------------------------------------------------------------------------------------------- //
eigVE<double> eigT(const sMat<double> & M)
{
   eigVE<double> retval;
   if(!M.issymm()) return retval;           // Checks that matrix is symmetric first otherwise returns empty

   // Calculates the eigenvalues and eigenvectors of the tridiagonal matrix represented in e[] and d[] by implicit QL
   //    factorisation, using imqtl2 algorithm of Dubrulle, Martin and Wilkinson, Num. Math. 12, 377-383 (1968), 
   //    which is a modification of the implicit QR factorisation algorithm of Francis, Comput. J. 4, 265-271 and
   //    332-345 (1961-1962), as converted from the Algol routine in the Handbook of Automatic Computation (1971).

   int n = M.nr();
   int i,j,k,ia,m,its;
   double c,p,q,s,t,u,h;
   bool doiter = true;
   bool failflag = false;

   std::vector<double> d = diag(M);         // Diagonal elements
   std::vector<double> e = diag(M,1);       // Off-diagonal elements
   e.push_back(0.);                         // Algorithm expects n elements - otherwise get a segfault!
   sMat<double> z(n,n); for(i=0; i<n; i++) z(i,i) = 1.;

   k = n-2;
   for(j=0; j<n; j++)
   {
      its = 0;
      while(doiter)
      {
         // Look for single small sub-diagonal element
         for(m=j; m<=k; m++)
            if( fabs(e[m]) <= (DBL_EPSILON*(fabs(d[m])+fabs(d[m+1]))) ) break;
       //m = n-1;        // If the loop finishes without condition satisfied, m=k+1 = n-1;
         if(m==k) break;
         u = d[j];
         if(m!=j)
         {
            if(its==100) { failflag = true; break; }//{ return retval; }
            its++;
            // Form shift
            q = (d[j+1]-u) / (2*e[j]); t = sqrt(1.+q*q); q = d[m] - u + e[j] / ( (q<0.) ? q-t : q+t );
            u = 0.; s = 1.; c = 1.;
            for(i=m-1; i>=j; i--)
            {
               p = s*e[i]; h = c*e[i];
               if(fabs(p)>=fabs(q))
               { c = q/p; t = sqrt(c*c+1.); e[i+1] = p*t; s = 1./t; c *= s; }
               else
               { s = p/q; t = sqrt(s*s+1.); e[i+1] = q*t; c = 1./t; s *= c; }
               q = d[i+1] - u; t = (d[i]-q)*s + 2.*c*h; u = s*t; d[i+1] = q+u; q = c*t-h;
               // Form vector
               for(ia=0; ia<n; ia++)
               { p = z(ia,i+1); z(ia,i+1) = s*z(ia,i) + c*p; z(ia,i) = c*z(ia,i) - s*p; }
            }
            d[j] = d[j] - u; e[j] = q; e[m] = 0.;
         }
	 else { doiter = false; }
      }
      doiter = true;
      if(failflag) break;
   }
   if(!failflag)
   {  // Order eigenvalues and eigenvectors
      for(i=0; i<n; i++)
      {
         k = i; p = d[i];
         for(j=i+1; j<n; j++) if(d[j]<p) { k = j; p = d[j]; }
         if(k!=i)
         {
            d[k] = d[i]; d[i] = p;
            for(j=0; j<n; j++) { p = z(j,i); z(j,i) = z(j,k); z(j,k) = p; }
         }
      }
      // Puts results into structure
      retval.E = d; retval.V = z;
   }
   return retval;
}

#ifdef USE_LAPACK
std::vector<double> svd(const sMat<double> & M)
{
   std::vector<double> retval;

   // The singular values q[] of the matrix A is given by:		       T
   // 			 T     T				A = U diag(q) V
   // 	 with:		U U = V V = I

   // Paramaters for DGESVD
   char jobu = 'S';        // Job type - 'A' returns all columns of U or V^T [jobu/jobvt]
   char jobvt = 'N';       //          - 'S' first min(m,n) columns of U or V^T ("economy" mode)
                           //          - 'O' fisrt min(m,n) columns of U or V^T overwritten on A
                           //          - 'N' don't compute U or V^T
   int m = M.nr();         // Number of rows of input matrix A
   int n = M.nc();         // Number of columns of input matrix A
   double *a = M.f_array();// Input matrix A (C-array version of M)
   int lda = m;            // Leading dimension of A >=max(1,m)
   double *s;              // On output the singular values of A, dimension min(m,n)
   double *u;              // On output the left singular vectors, U
   int ldu = m;            // Leading dimension of U
 //double *vt;             // On output the right singular vectors V^T
   int ldvt = n;           // Leading dimension of V^T
   double *work;           // A C-array for the workspace
   int lwork;              // The length of the workspace, LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
   int info = 0;           // Success flag on return. 0 is success. <0 is wrong input. >0 means not converged

   // Allocates memory to output arrays
   s  = (double*)malloc( (n<m?n:m) * sizeof(double));    // dimension min(m,n)
   u  = (double*)malloc( m*m * sizeof(double));          // dimension m*m

   // Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   lwork = 3*(m<n?m:n)+(m>n?m:n); lwork = (lwork>(5*(m>n?m:n))) ? lwork : (5*(m>n?m:n));
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, u, &ldvt, work, &lwork, &info);

   if(info==0)
      retval = f2vec(s,n<m?n:m);
   else
      std::cerr << "svd(): DGESVD return error code " << info << "\n";

   // Frees memory of output arrays
   free(s); free(u); free(work);
   return retval;
}
std::vector<double> svd(const sMat<double> & M, sMat<double> & V)
{
   std::vector<double> retval;

   // Calculates the right singular vectors as well this time

   // Paramaters for DGESVD
   char jobu = 'N',jobvt = 'A';
   int m = M.nr(), n = M.nc(), lda = m, ldu = m, ldvt = n, info = 0;
   double *a = M.f_array();
   double *s,*vt,*work;
   // Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   int lwork = 3*(m<n?m:n)+(m>n?m:n); lwork = (lwork>(5*(m>n?m:n))) ? lwork : (5*(m>n?m:n));

   // Allocates memory to arrays
   s  = (double*)malloc( (n<m?n:m) * sizeof(double));    // dimension min(m,n)
   vt = (double*)malloc( n*n * sizeof(double));          // dimension n*n
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, vt, &ldu, vt, &ldvt, work, &lwork, &info);

   if(info==0)
   {  retval = f2vec(s,n<m?n:m); V = f2mat(vt,n,n); }
   else
      std::cerr << "svd(): DGESVD return error code " << info << "\n";

   // Frees memory of output arrays
   free(s); free(vt); free(work);
   return retval;
}
std::vector<double> svd(const sMat<double> & M, sMat<double> & U, sMat<double> & V)
{
   std::vector<double> retval;

   // Calculates both the left and right singular vectors this time

   // Paramaters for DGESVD
   char jobu = 'A',jobvt = 'A';
   int m = M.nr(), n = M.nc(), lda = m, ldu = m, ldvt = n, info = 0;
   double *a = M.f_array();
   double *s,*u,*vt,*work;
   // Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   int lwork = 3*(m<n?m:n)+(m>n?m:n); lwork = (lwork>(5*(m>n?m:n))) ? lwork : (5*(m>n?m:n));

   // Allocates memory to arrays
   s  = (double*)malloc( (n<m?n:m) * sizeof(double));    // dimension min(m,n)
   u  = (double*)malloc( m*m * sizeof(double));          // dimension m*m
   vt = (double*)malloc( n*n * sizeof(double));          // dimension n*n
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

   if(info==0)
   {  retval = f2vec(s,n<m?n:m); U = f2mat(u,m,m); V = f2mat(vt,n,n); }
   else
      std::cerr << "svd(): DGESVD return error code " << info << "\n";

   // Frees memory of output arrays
   free(s); free(u); free(vt); free(work);
   return retval;
}
#else  // USE_LAPACK
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the single value decomposition of a real matrix
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> svd(const sMat<double> & M)
{
   // Calculates the singular value decomposition of a real rectangular matrix, M, by reducing M to bidiagonal 
   //    form using Householder rotations, then using the QR algorithm to find the singular values of this matrix.
   //    The proceedure used is that due to G.H. Golub and C. Reinsch, Num. Math. 14, 403-420 (1970), as converted
   //    from the Algol routine published in the Handbook of Automatic Computation (1971).
   
   // The singular values q[] of the matrix A is given by:		       T
   // 			 T     T				A = U diag(q) V
   // 	 with:		U U = V V = I

   int m = M.nr();
   int n = M.nc();
   //double tol = DBL_MIN/DBL_EPSILON;
   double tol = DBL_EPSILON;
   sMat<double> a = M;
   if(m<n) { a.transpose(); m = n; n = M.nr(); }
   std::vector<double> q(n,0.); //d.reserve(n+1);  // To store the eigenvalues
   sMat<double> u = a;
   sMat<double> v(n,n);

   int i,j,k,l,l1;
   double c,f,g,h,s,x,y,z;
   std::vector<double> e(n,0.);

   // Householder's reduction to bidiagonal form
   g = 0.; x = 0.;
   for(i=0; i<n; i++)
   {
      e[i] = g; s = 0.; l = i+1;
      for(j=i; j<m; j++)
         s += u(j,i)*u(j,i);
      if(s<tol) g = 0.;
      else
      {
         f = u(i,i); g = (f<0) ? sqrt(s) : -sqrt(s);
	 h = f*g - s; u(i,i) = f-g;
	 for(j=l; j<n; j++)
	 {
	    s = 0.;
	    for(k=i; k<m; k++)
	       s += u(k,i)*u(k,j);
	    f = s/h;
	    for(k=i; k<m; k++)
	       u(k,j) = u(k,j) + f*u(k,i);
	 }
      }
      q[i] = g; s = 0.;
      for(j=l; j<n; j++)
         s += u(i,j)*u(i,j);
      if(s<tol) g = 0.;
      else
      {
         f = u(i,i+1); g = (f<0) ? sqrt(s) : -sqrt(s);
	 h = f*g - s; u(i,i+1) = f-g;
	 for(j=l; j<n; j++)
	    e[j] = u(i,j)/h;
	 for(j=l; j<m; j++)
	 {
	    s = 0.;
	    for(k=l; k<n; k++)
	       s += u(j,k)*u(i,k);
	    for(k=l; k<n; k++)
	       u(j,k) = u(j,k) + s*e[k];
	 }
      }
      y = fabs(q[i]) + fabs(e[i]); if(y>x) x=y;
   }

   // Diagonalisation of the bidiagonal form
   tol *= x;
   for(k=n-1; k>=0; k--)
   {
      // test f splitting:
      while(1)
      {
         for(l=k; l>=0; l--)
         {
            if(fabs(e[l])<=tol)   // goto test f convergence
               goto test_f_convergence;
            if(fabs(q[l-1])<=tol) // goto cancellation
               break;
         }
   
         // Cancellation of e[l] if l>1
         // cancellation:
         c = 0.; s = 1.; l1 = l-1;
         for(i=l; i<=k; i++)
         {
            f = s*e[i]; e[i] = c*e[i];
            if(fabs(f)<=tol)      // goto test f convergence
	       break;
            g = q[i]; h = sqrt(f*f + g*g); q[i] = h; c = g/h; s = -f/h;
         }
         test_f_convergence:
         z = q[k];
         //if(l==k) // goto convergence
         if(l!=k)
         {
            // Shift from bottom 2x2 minor
            x = q[l]; y = q[k-1]; g = e[k-1]; h = e[k];
            f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y); g = sqrt(f*f+1.);
            f = ((x-z)*(x+z) + h*(y / ( (f<0)?(f-g):(f+g) ) - h)) / x;
   
            // Next QR transformation
            c = 1.; s = 1.;
            for(i=l+1; i<=k; i++)
            {
               g = e[i]; y = q[i]; h = s*g; g *= c;
               z = sqrt(f*f + h*h); e[i-1] = z; c = f/z; s = h/z;
               f = x*c + g*s; g = -x*s + g*c; h = y*s; y *= c;
               z = sqrt(f*f + h*h); q[i-1] = z; c = f/z; s = h/z;
               f = g*c + y*s; x = -g*s + y*c;
            }
            e[l] = 0.; e[k] = f; q[k] = x; // goto test f splitting
            continue;
         }
   
         // convergence:
         if(z<0)
            q[k] = -z;   // q[k] is made non-negative
         break;
      }
   }

   return q;
}
std::vector<double> svd(const sMat<double> & M, sMat<double> & V)
{
   // This version of the function also returns the U and V matrices.
   // The singular values q[] of the matrix A is given by:		       T
   // 			 T     T				A = U diag(q) V
   // 	 with:		U U = V V = I

   int m = M.nr();
   int n = M.nc();
   //double tol = DBL_MIN/DBL_EPSILON;
   double tol = DBL_EPSILON;
   sMat<double> a = M; a.transpose();
   std::vector<double> q(n,0.);                          // To store the eigenvalues
   sMat<double> u = M;
   sMat<double> v(n,n);

   int i,j,k,l=0,l1;//,pass;
   double c,f,g,h,s,x,y,z;
   std::vector<double> e(n,0.);

   // Householder's reduction to bidiagonal form
   g = 0.; x = 0.;
   for(i=0; i<n; i++)
   {
      e[i] = g; s = 0.; l = i+1;
      for(j=i; j<m; j++)
         s += u(j,i)*u(j,i);
      if(s<tol) g = 0.;
      else
      {
         f = u(i,i); g = (f<0) ? sqrt(s) : -sqrt(s);
	 h = f*g - s; u(i,i) = f-g;
	 for(j=l; j<n; j++)
	 {
	    s = 0.;
	    for(k=i; k<m; k++)
	       s += u(k,i)*u(k,j);
	    f = s/h;
	    for(k=i; k<m; k++)
	       u(k,j) = u(k,j) + f*u(k,i);
	 }
      }
      q[i] = g; s = 0.;
      for(j=l; j<n; j++)
         s += u(i,j)*u(i,j);
      if(s<tol) g = 0.;
      else
      {
         f = u(i,i+1); g = (f<0) ? sqrt(s) : -sqrt(s);
	 h = f*g - s; u(i,i+1) = f-g;
	 for(j=l; j<n; j++)
	    e[j] = u(i,j)/h;
	 for(j=l; j<m; j++)
	 {
	    s = 0.;
	    for(k=l; k<n; k++)
	       s += u(j,k)*u(i,k);
	    for(k=l; k<n; k++)
	       u(j,k) = u(j,k) + s*e[k];
	 }
      }
      y = fabs(q[i]) + fabs(e[i]); if(y>x) x=y;
   }

   // Accumulation of right-hand transformations
   for(i=n-1; i>=0; i--)   // if withv
   {
      if(g!=0)
      {
         h = u(i,i+1)*g;
	 for(j=l; j<n; j++)
	    v(j,i) = u(i,j)/h;
	 for(j=l; j<n; j++)
	 {
	    s = 0.;
	    for(k=l; k<n; k++)
	       s += u(i,k)*v(k,j);
	    for(k=l; k<n; k++)
	       v(k,j) = v(k,j) + s*v(k,i);
	 }
      }
      for(j=l; j<n; j++)
      {
         //v.del(i,j); v.del(j,i);
	 v(i,j) = 0.; v(j,i) = 0.;
      }
      v(i,i) = 1.; g = e[i]; l = i;
   }

   // Accumulation of left-hand transformations
   //for(i=n-1; i>=0; i--)   // if withu
   //{
   //   l = i + 1; g = q[i];
   //   for(j=l; j<n; j++)
   //      u(i,j) = 0.; //u.del(i,j);
   //   if(g!=0)
   //   {
   //      h = u(i,i)*g;
   //      for(j=l; j<n; j++)
   //      {
   //         s = 0.;
   //         for(k=l; k<m; k++)
   //            s += u(k,i)*u(k,j);
   //         f = s/h;
   //         for(k=i; k<m; k++)
   //            u(k,j) = u(k,j) + f*u(k,i);
   //      }
   //      for(j=i; j<m; j++)
   //         u(j,i) = 0.; //u.del(j,i);
   //      }
   //   else for(j=i; j<m; j++)
   //      u(j,i) = 0.; //u.del(j,i);
   //   u(i,i) = u(i,i) + 1;
   //}

   // Diagonalisation of the bidiagonal form
   tol *= x;
   for(k=n-1; k>=0; k--)
   {
      // test f splitting:
      while(1)
      {
         for(l=k; l>=0; l--)
         {
            if(fabs(e[l])<=tol)   // goto test f convergence
               goto test_f_convergence;
            if(fabs(q[l-1])<=tol) // goto cancellation
               break;
         }
   
         // Cancellation of e[l] if l>1
         // cancellation:
         c = 0.; s = 1.; l1 = l-1;
         for(i=l; i<=k; i++)
         {
            f = s*e[i]; e[i] = c*e[i];
            if(fabs(f)<=tol)      // goto test f convergence
	       break;
            g = q[i]; h = sqrt(f*f + g*g); q[i] = h; c = g/h; s = -f/h;
            //for(j=0; j<m; j++) // if withu
            //{
            //   y = u(j,l1); z = u(j,i);
            //   u(j,l1) = y*c + z*s; u(j,i) = -y*s + z*c;
            //}
         }
         test_f_convergence:
         z = q[k];
         //if(l==k) // goto convergence
         if(l!=k)
         {
            // Shift from bottom 2x2 minor
            x = q[l]; y = q[k-1]; g = e[k-1]; h = e[k];
            f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y); g = sqrt(f*f+1.);
            f = ((x-z)*(x+z) + h*(y / ( (f<0)?(f-g):(f+g) ) - h)) / x;
   
            // Next QR transformation
            c = 1.; s = 1.;
            for(i=l+1; i<=k; i++)
            {
               g = e[i]; y = q[i]; h = s*g; g *= c;
               z = sqrt(f*f + h*h); e[i-1] = z; c = f/z; s = h/z;
               f = x*c + g*s; g = -x*s + g*c; h = y*s; y *= c;
               for(j=0; j<n; j++)  // if withv
               {
                  x = v(j,i-1); z = v(j,i);
                  v(j,i-1) = x*c + z*s; v(j,i) = -x*s + z*c;
               }
               z = sqrt(f*f + h*h); q[i-1] = z; c = f/z; s = h/z;
               f = g*c + y*s; x = -g*s + y*c;
               //for(j=0; j<m; j++)  // if withu
               //{
               //   y = u(j,i-1); z = u(j,i);
               //   u(j,i-1) = y*c + z*s; u(j,i) = -y*s + z*c;
               //}
            }
            e[l] = 0.; e[k] = f; q[k] = x; // goto test f splitting
            continue;
         }
   
         // convergence:
         if(z<0)
         {   // q[k] is made non-negative
            q[k] = -z;
            for(j=0; j<n; j++)  // if withv
               v(j,k) = -v(j,k);
         }
         break;
      }
   }

   V = v; 

   return q;
}
std::vector<double> svd(const sMat<double> & M, sMat<double> & U, sMat<double> & V)
{
   int m = M.nr();
   int n = M.nc();
   std::vector<double> d;
   sMat<double> Mt = M; Mt.transpose();

   if(m>n)
   {
      d = svd(M,V);
      d = svd(Mt,U);
   }
   else if (m<n)
   {
      d = svd(Mt,V);
      d = svd(M,U);
   }   
   else
   {
      d = svd(M,V);
      U = V;
   }
   
   return d;
}
#endif // USE_LAPACK

#ifdef USE_LAPACK
sMat<double> qr(const sMat<double> & M, sMat<double> & Q)
{
   sMat<double> retval;    // The matrix R to be returned.

   // Parameters for DGEQRF
   int m = M.nr();         // Number of rows of input matrix A
   int n = M.nc();         // Number of columns of A
   double *a = M.f_array();// Input matrix A as an array
   int lda = m;            // Leading dimension of A
   double *tau;            // On output, the scalar factors of elementary reflectors
   double *work;           // Workspace
   int lwork = -1;         // Length of workspace. If lwork=-1, calculates optimal workspace
   int info = 0;           // Success flag on return. 0 is success. <0 is wrong input. >0 means not converged
   // Extra parameters for DORGQR
   int k = m<n?m:n;        // The number of elementary reflectors

   // Allocates memory to array initially
   tau = (double*)malloc( (m<n?m:n) * sizeof(double) );
   work = (double*)malloc( 4*n * sizeof(double) );      // Initial value only. We use the workspace query to find optimum
   // Calculates optimum workspace size
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
   if(info==0) lwork = work[0]; else { std::cerr << "qr(): DGEQRF could not determine optimum workspace size\n"; return retval; }
   work = (double*)realloc( work, lwork*sizeof(double) );
   if(work==NULL) { std::cerr << "qr(): error reallocating workspace memory\n"; return retval; }

   // Calculates the QR factorisation
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);

   if(info==0)
   {
      retval = triu(f2mat(a,m,n));
      // Calculates the orthogonal matrix Q from the output of DGEQRF
      F77NAME(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
      if(info==0) Q = f2mat(a,m,n); else std::cerr << "qr(): DGEQRF return error code " << info << "\n";
   }
   else
      std::cerr << "qr(): DGEQRF return error code " << info << "\n";

   free(tau); free(work);
   return retval;
}
sMat<double> qr(const sMat<double> & M, sMat<double> & Q, int zero)
{
   sMat<double> retval;    // The matrix R to be returned.

   // Parameters for DGEQRF and DORGQR
   int m = M.nr(), n = M.nc(), lda = m, lwork = -1, info = 0, k = m<n?m:n;
   double *a, *tau, *work;
   if(m>n) { m = n; lda = n; retval = M; retval.resize(1,n,1,n); a = retval.f_array(); } else a = M.f_array();
   // Allocates memory to array initially
   tau = (double*)malloc( (m<n?m:n) * sizeof(double) );
   work = (double*)malloc( 4*n * sizeof(double) );      // Initial value only. We use the workspace query to find optimum
   // Calculates optimum workspace size
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
   if(info==0) lwork = work[0]; else { std::cerr << "qr(): DGEQRF could not determine optimum workspace size\n"; return retval; }
   work = (double*)realloc( work, lwork*sizeof(double) );
   if(work==NULL) { std::cerr << "qr(): error reallocating workspace memory\n"; return retval; }

   // Calculates the QR factorisation
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);

   if(info==0)
   {
      retval = triu(f2mat(a,m,n));
      // Calculates the orthogonal matrix Q from the output of DGEQRF
      F77NAME(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
      if(info==0) Q = f2mat(a,m,n); else std::cerr << "qr(): DGEQRF return error code " << info << "\n";
   }
   else
      std::cerr << "qr(): DGEQRF return error code " << info << "\n";

   return retval;
}
sMat<double> qr(sMat<double> & M, int zero)
{
   sMat<double> retval;    // The matrix R to be returned.

   // Parameters for DGEQRF
   int m = M.nr(), n = M.nc(), lda = m, lwork = -1, info = 0;
   double *a, *tau, *work;
   if(m>n) { m = n; lda = n; retval = M; retval.resize(1,n,1,n); a = retval.f_array(); } else a = M.f_array();
   // Allocates memory to array initially
   tau = (double*)malloc( (m<n?m:n) * sizeof(double) );
   work = (double*)malloc( 4*n * sizeof(double) );      // Initial value only. We use the workspace query to find optimum
   // Calculates optimum workspace size
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
   if(info==0) lwork = work[0]; else { std::cerr << "qr(): DGEQRF could not determine optimum workspace size\n"; return retval; }
   work = (double*)realloc( work, lwork*sizeof(double) );
   if(work==NULL) { std::cerr << "qr(): error reallocating workspace memory\n"; return retval; }

   // Calculates the QR factorisation
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);

   if(info==0)
      retval = triu(f2mat(a,m,n));
   else
      std::cerr << "qr(): DGEQRF return error code " << info << "\n";

   free(tau); free(work);
   return retval;
}
#else  // USE_LAPACK
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the QR decomposition of a rectangular matrix
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> qr(const sMat<double> & M, sMat<double> & Q)
{
   // Calculates the QR decomposition using Householder reflections with column pivoting, after Businger and Golub.
   //    This function is a translation from the Algol proceedure "decompose", which is part of the proceedure
   //    "least squares solution", by P. Businger and G.H. Golub, Num. Math. 7, 269-276 (1965), as reproduced in the 
   //    Handbook of Automatic Computation (1971).

 //double eta = DBL_EPSILON;
   int m = M.nr();
   int n = M.nc();
   sMat<double> qr = M;
   int i,j,jbar,k;
   double sigma,alphak,qrkk;
   std::vector<double> y(n,0.),sum(n,0.),beta(n,0.);
   std::vector<int> pivot(n,0);
   double ip,ipp;                             // Temporary variable for the inner_product macro
   sMat<double> alpha(m,n);                   // The upper triangular matrix R
   Q.zero(m,m);                               // The orthonormal matrix Q
   sMat<double> P(m,m);
   sMat<double> eye(m,m); for(i=0; i<m; i++) eye(i,i) = 1.;

#define inner_product(I,M,N,A,B,C,D) ip=C; for(I=M; I<=N; I++) ip += A*B; D = ip

   for(j=0; j<n; j++)
   {  // j-th column sum
      /* sum[j] = */ inner_product(i,0,m-1,qr(i,j),qr(i,j),0.,sum[j]); 
      pivot[j] = j;
   }
   for(k=0; k<n; k++)
   {  // k-th Householder transformation
      sigma = sum[k]; jbar = k;
      for(j=k+1; j<n; j++)
         if(sigma<sum[j])
         {
            sigma = sum[j]; jbar = j;
         }
    /*if(jbar!=k)
      {  // column interchange
         i = pivot[k]; pivot[k] = pivot[jbar]; pivot[jbar] = i;
         sum[jbar] = sum[k]; sum[k] = sigma;
         for(i=0; i<m; i++)
         {
            sigma = qr(i,k); qr(i,k) = qr(i,jbar); 
            qr(i,jbar) = sigma;
         }
      }*/
      /* sigma = */ inner_product(i,k,m-1,qr(i,k),qr(i,k),0.,sigma);
      if(sigma==0) // go to singular
      {
	 break;
      }
      qrkk = qr(k,k);
      alphak = (qrkk<0) ? sqrt(sigma) : -sqrt(sigma); alpha(k,k) = alphak;
      beta[k] = 1/(sigma - qrkk*alphak);
      qr(k,k) = qrkk-alphak;
      for(j=k+1; j<n; j++)
      {
         inner_product(i,k,m-1,qr(i,k),qr(i,j),0.,ipp); y[j] = beta[k]*ipp;
      }
      for(j=k+1; j<n; j++)
      {
         for(i=k; i<m; i++)
            qr(i,j) = qr(i,j) - qr(i,k)*y[j];
         sum[j] = sum[j] - qr(k,j)*qr(k,j);
      }
   }
   // End of decompose

   alpha = alpha + triu(qr,1);

   // Works out the Q matrix, from the factors in the lower half of qr()
   //     T  ___                       (k) (k)T
   //    Q = | | P   where P  = I  -  v   v     * beta_k  are the Householder transformation matrices
   //         k   k         k    m
   Q = eye;
   for(k=m-1; k>=0; k--)
   {
      P = eye;
      for(i=k; i<m; i++)
      {
         for(j=k; j<=i; j++)
         {
            P(i,j) = P(i,j) - qr(i,k)*qr(j,k)*beta[k];
            P(j,i) = P(i,j);
         }
      }
      Q *= P;
   }
   Q.transpose();

   // This algorithm actually calculates the decomposition: AP = QR, where P is a permutation matrix
   //    P = [e_p1, e_p2, ..., e_pn] where p1,...,pn == pivot. Thus we need to back convert it to get
   //    A = QR
   //P.zero(n,n);
   //for(k=0; k<n; k++)
   //   P(k,pivot[k]) = 1.;

   return alpha;
}
sMat<double> qr(const sMat<double> & M, sMat<double> & Q, int zero)
{
   // Calculates only the first n columns of Q

 //double eta = DBL_EPSILON;
   int m = M.nr();
   int n = M.nc();
   sMat<double> qr = M;
   int i,j,jbar,k;
   double sigma,alphak,qrkk;
   std::vector<double> y(n,0.),sum(n,0.),beta(n,0.);
   std::vector<int> pivot(n,0);
   double ip,ipp;                             // Temporary variable for the inner_product macro
   sMat<double> alpha(m,n);                   // The upper triangular matrix R
   Q.zero(m,m);                               // The orthonormal matrix Q
   sMat<double> P(m,m);
   sMat<double> eye(m,m); for(i=0; i<m; i++) eye(i,i) = 1.;

#define inner_product(I,M,N,A,B,C,D) ip=C; for(I=M; I<=N; I++) ip += A*B; D = ip

   for(j=0; j<n; j++)
   {  // j-th column sum
      /* sum[j] = */ inner_product(i,0,m-1,qr(i,j),qr(i,j),0.,sum[j]); 
      pivot[j] = j;
   }
   for(k=0; k<n; k++)
   {  // k-th Householder transformation
      sigma = sum[k]; jbar = k;
      for(j=k+1; j<n; j++)
         if(sigma<sum[j])
         {
            sigma = sum[j]; jbar = j;
         }
      /* sigma = */ inner_product(i,k,m-1,qr(i,k),qr(i,k),0.,sigma);
      if(sigma==0) // go to singular
      {
	 break;
      }
      qrkk = qr(k,k);
      alphak = (qrkk<0) ? sqrt(sigma) : -sqrt(sigma); alpha(k,k) = alphak;
      beta[k] = 1/(sigma - qrkk*alphak);
      qr(k,k) = qrkk-alphak;
      for(j=k+1; j<n; j++)
      {
         inner_product(i,k,m-1,qr(i,k),qr(i,j),0.,ipp); y[j] = beta[k]*ipp;
      }
      for(j=k+1; j<n; j++)
      {
         for(i=k; i<m; i++)
            qr(i,j) = qr(i,j) - qr(i,k)*y[j];
         sum[j] = sum[j] - qr(k,j)*qr(k,j);
      }
   }
   // End of decompose

   alpha = alpha + triu(qr,1);

   Q = eye; 
   for(k=((m>n)?n:m)-1; k>=0; k--)
   {
      P = eye;
      for(i=k; i<m; i++)
      {
         for(j=k; j<=i; j++)
         {
            P(i,j) = P(i,j) - qr(i,k)*qr(j,k)*beta[k];
            P(j,i) = P(i,j);
         }
      }
      Q *= P;
   }
   Q.transpose();

   return alpha;
}

sMat<double> qr(sMat<double> & M, int zero)
{
   int m = M.nr();
   int n = M.nc();
   sMat<double> qr = M;
   int i,j,jbar,k;
   double sigma,alphak,qrkk;
   std::vector<double> y(n,0.),sum(n,0.),beta(n,0.);
   std::vector<int> pivot(n,0);
   double ip,ipp;                             // Temporary variable for the inner_product macro
   sMat<double> alpha(m,n);                   // The upper triangular matrix R
   sMat<double> P(m,m);
   sMat<double> eye(m,m); for(i=0; i<m; i++) eye(i,i) = 1.;

#define inner_product(I,M,N,A,B,C,D) ip=C; for(I=M; I<=N; I++) ip += A*B; D = ip

   for(j=0; j<n; j++)
   {  // j-th column sum
      /* sum[j] = */ inner_product(i,0,m-1,qr(i,j),qr(i,j),0.,sum[j]); 
      pivot[j] = j;
   }
   for(k=0; k<n; k++)
   {  // k-th Householder transformation
      sigma = sum[k]; jbar = k;
      for(j=k+1; j<n; j++)
         if(sigma<sum[j])
         {
            sigma = sum[j]; jbar = j;
         }
      /* sigma = */ inner_product(i,k,m-1,qr(i,k),qr(i,k),0.,sigma);
      if(sigma==0) // go to singular
      {
	 break;
      }
      qrkk = qr(k,k);
      alphak = (qrkk<0) ? sqrt(sigma) : -sqrt(sigma); alpha(k,k) = alphak;
      beta[k] = 1/(sigma - qrkk*alphak);
      qr(k,k) = qrkk-alphak;
      for(j=k+1; j<n; j++)
      {
         inner_product(i,k,m-1,qr(i,k),qr(i,j),0.,ipp); y[j] = beta[k]*ipp;
      }
      for(j=k+1; j<n; j++)
      {
         for(i=k; i<m; i++)
            qr(i,j) = qr(i,j) - qr(i,k)*y[j];
         sum[j] = sum[j] - qr(k,j)*qr(k,j);
      }
   }
   // End of decompose

   alpha = alpha + triu(qr,1);

   return alpha;
}
#endif // USE_LAPACK

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the absolute values of each non zero element in a matrix
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> mabs(const sMat<double> & M)
{
   std::vector< std::vector<int> > nz;
   int i,sz;
   sMat<double> r = M;

   nz = M.find(); sz = (int)nz.size();
   for (i=0; i<sz; i++)
      r(nz[i][0],nz[i][1]) = fabs(M(nz[i][0],nz[i][1]));

   return r;
}
std::vector<double> vabs(const std::vector<double> & V)
{
   std::vector<double> r(V.size(),0.);
   int i;

   for (i=0; i<(int)V.size(); i++)
      r[i] = fabs(V[i]);

   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the norm of a matrix
// --------------------------------------------------------------------------------------------------------------- //
double norm(const sMat<double> & M)
{
   return vmax<double>(svd(M));
}
double norm(const sMat<double> & M, normp p)
{
   sMat<double> Mt = M; Mt.transpose();
   switch(p)                                                       // To be compatible with the matlab function
   {
      case 1: return vmax<double>(msum<double>(mabs(M))); break;   // The 1-norm
      case 2: return norm(M); break;                               // Largest singular value
      case 3: return vmax<double>(msum<double>(mabs(Mt)));         // The infinity norm (normp inf==3, see maths.hpp)
              break;
      case 4: return sqrt(vsum<double>(diag<double>(Mt*M)));       // The frobenius norm (normp fro==4, see maths.hpp)
              break;
      default: return 0.;
   }
   return 0.;
}

// --------------------------------------------------------------------------------------------------------------- //
// Returns a (mxn) matrix of normally distributed random numbers
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> randn(int m, int n)
{
   sMat<double> M(m,n);
   int i,j,k=0;
   srand(time(NULL));
   double x1=0,x2=0,y=0;
   double maxrand = RAND_MAX/2.;

   // This function uses the C++ standard function rand() to generate uniformly distributed random numbers
   //    between [0,1], and then the polar Box-Muller transformation to generate pairs of normally distributed
   //    random numbers (corresponding to sin(angle) and cos(angle) respectively). Hence we need some trickery to
   //    use both the pair, and not waste half the computing time! [ mean=0, s.d.=1 ]

   for(i=0; i<m; i++)
      for(j=0; j<n; j++)
      {
         if(k%2 == 0)		// First of pair, generate new pair.
	 {
	    do
	    {
               x1 = rand()/maxrand - 1.;
               x2 = rand()/maxrand - 1.;
               y  = x1*x1 + x2*x2;
            } while (y>=1.);
	    y = sqrt( (-2. * log(y)) / y);
	    M(i,j) = x1*y;
	 }
	 else			// Second of pair.
	    M(i,j) = x2*y;
	 k++;
      }

   return M;
}

// --------------------------------------------------------------------------------------------------------------- //
// Returns an orthonormal basis for a matrix. Translation of the MATLAB routine using svd.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> orth(const sMat<double> & M)
{
   sMat<double> U,Q,V;
   int i,r=0;
   sMat<double> Mt = M; Mt.transpose();
   int m = M.nr(); int n = M.nc();                       // [m,n] = size(A);
   std::vector<double> s = svd(Mt,U); U.resize(1,m,1,n); // [U,S,V] = svd(A,0);
   double tol;

   if(m>1) {}                                            // if m > 1, s = diag(S);
      else if(m==1) { s.erase(s.begin()+1,s.end()); }    //    elseif m == 1, s = S(1);
      else { s.clear(); s.push_back(0.); }               //    else s = 0;
                                                         // end
   tol = (double)vmax<int>(M.size()) * vmax<double>(s)   // tol = max(m,n) * max(s) * eps(class(A)); 
           * DBL_EPSILON;
   for(i=0; i<(int)s.size(); i++)
      if(s[i]>tol) r++;                                  // r = sum(s > tol);

   Q = U.setp(-1,-1,1,r);                                // Q = U(:,1:r);

   return Q;
}

// --------------------------------------------------------------------------------------------------------------- //
// Estimates the number of linearly independent rows or columns of a matrix. Translation of the MATLAB function.
// --------------------------------------------------------------------------------------------------------------- //
int rank(const sMat<double> & M, double tol)
{
   std::vector<double> s = svd(M);
   int i,r=0;
   for(i=0; i<(int)s.size(); i++)
      if(s[i]>tol) r++;

   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Removes entries in matrix less than machine epsilon
// --------------------------------------------------------------------------------------------------------------- //
void rmzeros(sMat<double> & M)
{
   std::vector< std::vector<int> > nz;
   int i,sz;
   nz = M.find(); sz = (int)nz.size();
   for (i=0; i<sz; i++)
      if(fabs(M(nz[i][0],nz[i][1]))<DBL_EPSILON*10) M.del(nz[i][0],nz[i][1]);
}

// --------------------------------------------------------------------------------------------------------------- //
// Converts a 1D C-array into a std::vector
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> f2vec(double *v, int n)
{
   int i;
   std::vector<double> retval(n,0.);
   for (i=0; i<n; i++) retval[i] = v[i];
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Converts a 2D C-array into an sMat
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> f2mat(double *M, int m, int n)
{
   int i,j;
   sMat<double> retval(m,n);
   for(j=0; j<n; j++)
      for(i=0; i<m; i++)
         if(fabs(M[m*j+i])>DBL_EPSILON) retval(i,j) = M[m*j+i];
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Converts two double precision sMat to a complex C-array
// --------------------------------------------------------------------------------------------------------------- //
complexdouble* zmat2f(sMat<double> &r, sMat<double> &i)
{
   sMat<double> tmp = r+i;
   std::vector< std::vector<int> > u = tmp.find/*upper*/();
   int j, rows=tmp.nr();
   // Allocates an _r*_c array and initiallises all elements to zero.
   complexdouble x, *retval = (complexdouble*)calloc(rows*tmp.nc(),sizeof(complexdouble));
   for (j=0; j<(int)u.size(); j++)
   {
      x.r = r(u[j][0],u[j][1]); x.i = i(u[j][0],u[j][1]);
      retval[rows*u[j][1]+u[j][0]] = x;
   }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: g++ maths.cpp && ./a.out
// --------------------------------------------------------------------------------------------------------------- //
/*int main(int argc, char *argv[])
{
   sMat<double> l,m,n,a,b,c;//(3,3);
   std::vector< std::vector<int> > nz;
   std::vector<int> size;
   std::vector<double> v,v2;
   int i,sz;
   eigVE<double> dg;
  
   std::cout << "Empty matrix is " << m; // << "\n";

   // Tests mset("") member function
   m.mset("1 0 0 ; 0 1 0 ; 0 0 1");
   nz = m.find();
   std::cout << "m.mset(\"1 0 0 ; 0 1 0 ; 0 0 1\"); m = " << m;
   std::cout << "Setting (2,3) to 2:\n"; m(1,2) = 2; std::cout << m;

   // Test find() member function.
   std::cout << "Testing find(): nz = m.find(); gives:";
   sz = (int)nz.size(); //std::cout << sz << "\n";
   for (i=0; i<sz; i++)
   {
      std::cout << "(" << nz[i][0] << "," << nz[i][1] << ")\t" << "\n";
   }
   m.mset("1/3 sqrt(2) 2 ; 2 1/2 2 ; sqrt(2/4) 2 1");
   nz = m.find(); sz = (int)nz.size();
   for (i=0; i<sz; i++)
   {
      std::cout << "(" << nz[i][0] << "," << nz[i][1] << ")\t" << m(nz[i][0],nz[i][1]) << "\n";
   }

   // Tests the iostream output friend function.
   n.mset("1*3 sqrt(2) 2/6 ; 2 1/2 2 ; sqrt(2/4) 2+1 1-6");
   std::cout << "Testing iostream: cout << n gives:\n" << n;
   //n >> l;

   // Tests size() and issquare()
   size = m.size(); std::cout << "Matrix size is " << size[0] << "x" << size[1] << "\n";
   std::cout << "Matrix issquare is " << (m.issquare() ? "true" : "false" ) << "\n";
  
   // Tests overloaded operators
   m+=n;
   std::cout << m; 
   m-=m;
   std::cout << m; 
   m.mset("1 0 2 ; 0-1 3 1"); n.mset("3 1 ; 2 1 ; 1 0");
   std::cout << m; std::cout << " x \n"; std::cout << n;
   std::cout << "null value is " << m(0,1) << "\n";
   m*=n;
   std::cout << m;
   //m.mset("1 0 2 ; 0-1 3 1"); n.mset("3 ; 2 ; 1");
   //std::cout << n;
   //m*=n; std::cout << m;
   m.mset("1 0 2 ; 0-1 3 1"); std::cout << "Full display gives:\n" << m.display_full();
   m(0,1) = 11; std::cout << "m(0,1)=11 gives " << m(0,1) << " gives:\n" << m.display_full();
   v.push_back(3); v.push_back(2); v.push_back(1);
   v = m*v; std::cout << "V size is " << v.size() << ". V = [" << v[0] << "," << v[1] << "]\n";
   m*=5; std::cout << "m=\n" << m.display_full(); 
   std::cout << "n=\n" << n.display_full();
   n.transpose(); a = m % n; std::cout << "a = m % n' =\n" << a.display_full();
   a = m ^ n; std::cout << "a = m ^ n' =\n" << a.display_full();
   
   // Tests numdenom
   //numdenom nd,nd3;
   //std::cout << "before: " << nd << " after: ";
   //nd = 1; std::cout << "nd = 1 is " << nd << "\n";
   //std::cout << "(double)nd gives "<< (double)nd << "\n";
   //numdenom nd2(1,2,2,1); std::cout << "nd = sqrt(2)/2 is " << nd2 << "\n";
   //std::cout << "double(nd) gives "<< double(nd2) << "\n";
   //nd = 2; std::cout << "nd = 2 is " << nd << "\n";
   //nd3 = sqrt(nd);
   //std::cout << "nd = sqrt(2) gives "<< nd3 << " = " << double(nd3) << "\n";

   std::cout << "#rows in m is: " << m.nr() << "\n";

   // Tests pset
   a.mset("6 0-2 5 7 ; 8 0-7 4 2 ; 1 0-3 1 5 ; 2 5 3 1"); b.mset("6 0-3 ; 2 9");
   std::cout << "a =\n" << a.display_full();
   std::cout << "b =\n" << b.display_full();
   a.pset(2,3,2,3,b); std::cout << "a.pset(2,3,2,3,b) gives:\n" << a.display_full();

   // Tests diagonalisation routines
   c = mabs(a); std::cout << "mabs(a) =\n" << c.display_full();
   v = mmax(a); std::cout << "mmax(a) = " << dispvect(v) << "\t vmax(mmax(a)) = " << vmax(v) << "\n";
   v = msum(a); std::cout << "msum(a) = " << dispvect(v) << "\t vsum(msum(a)) = " << vsum(v) << "\n";
   c = a; c.transpose(); std::cout << "a.transpose() =\n" << c.display_full();
   std::cout << "a.issymm() = " << a.issymm() << "\n";
   b = triu(a); c = tril(c); c = b+c; std::cout << "c = triu(a)+ tril(a.transpose) =\n" << c.display_full();
   std::cout << "c.issymm() = " << c.issymm() << "\n";
   // Matrix from the Handbook of Automatic Computation, pp.223
   a.mset("10 1 2 3 4 ; 1 9 0-1 2 0-3 ; 2 0-1 7 3 0-5 ; 3 2 3 12 0-1 ; 4 0-3 0-5 0-1 15");
   b = triu(a,1); std::cout << "triu(a,1) =\n" << b.display_full();
   b.mset("5 1 0-2 0 0-2 5 ; 1 6 0-3 2 0 6 ; 0-2 0-3 8 0-5 0-6 0 ; 0 2 0-5 5 1 0-2 ; 0-2 0 0-6 1 6 0-3 ; 5 6 0 0-2 0-3 8");
   std::cout << "a =\n" << a.display_full();
   dg = eig(a); std::cout << "dg=eig(a); dg.E = " << dispvect(dg.E) << "\n" << "dg.V =\n" << dg.V.display_full();
   std::cout << "b =\n" << b.display_full();
   dg = eig(b); std::cout << "dg=eig(b); dg.E = " << dispvect(dg.E) << "\n" << "dg.V =\n" << dg.V.display_full();
   v = svd(a); std::cout << "svd(a) = " << dispvect(v) << "\n";
   v = svd(b); std::cout << "svd(b) = " << dispvect(v) << "\n";
   a.mset("1 1 9 ; 5 8 5"); std::cout << "a is now\n" << a.display_full();
   v2 = svd(a); std::cout << "svd(a) = " << dispvect(v) << "\n";
   a.mset("5 2 ; 4 4 ; 5 6"); std::cout << "a is now\n" << a.display_full();
   v = svd(a,m,n); std::cout << "svd(a,u,v) = " << dispvect(v) << "\nU =\n" << m.display_full() << "V =\n" << n.display_full();
   a.transpose(); std::cout << "a is now\n" << a.display_full();
   v = svd(a,m); std::cout << "svd(a,u) = " << dispvect(v2) << "\nU =\n" << m.display_full();
   std::cout << "norm(a) = " << norm(a) << "\n";
   std::cout << "rank(a) = " << rank(a,(double)vmax<int>(a.size())*vmax<double>(v)*DBL_EPSILON) << "\n";
   m = orth(a); std::cout << "orth(a) =\n" << m.display_full();
   a.transpose(); std::cout << "a is now\n" << a.display_full();
   m = orth(a); std::cout << "orth(a) =\n" << m.display_full();
   b.mset("7 ; 8 ; 9"); a.pset(1,a.nr(),a.nc()+1,a.nc()+1,b); std::cout << "a is now\n" << a.display_full();
   a = randn(2,3); std::cout << "randn(2,3) =\n" << a.display_full(); 
   //dg = jacobi(a); std::cout << "E = " << dispvect(dg.E) << "\n" << "V =\n" << dg.V.display_full();
   v.resize(5,0.); v[0]=-1; v[1]=0; v[2]=2; v[3]=4; v[4]=6; v2.resize(4,0.); v2[0]=-1; v2[1]=0; v2[2]=1; v2[3]=3;
   std::cout << "setunion(" << dispvect(v) << "," << dispvect(v2) << ") = " << dispvect(setunion(v,v2)) << "\n";
   std::cout << "setdiff(" << dispvect(v) << "," << dispvect(v2) << ") = " << dispvect(setdiff(v,v2)) << "\n";
   std::cout << "setxor(" << dispvect(v) << "," << dispvect(v2) << ") = " << dispvect(setxor(v,v2)) << "\n";
   a((int)(rand()%2),(int)(rand()%3)) = 0; std::cout << "a is now\n" << a.display_full();
   v = all(a); std::cout << "all(a) = " << dispvect(v) << "\n";
   v = all(a,2); std::cout << "all(a,2) = " << dispvect(v) << "\n";
   b = qr(a,c); std::cout << "qr(a,c) =\n" << b.display_full() << "c =\n" << c.display_full();
   //a.mset("5 1 0-2 0 0-2 5 ; 1 6 0-3 2 0 6 ; 0-2 0-3 8 0-5 0-6 1 ; 0 2 0-5 5 1 0-2 ; 0-2 0 0-6 1 6 0-3 ; 5 6 0 0-2 0-3 8");
   a.mset("1 2 3 ; 1 5 6 ; 1 8 9 ; 1 11 12");
   std::cout << "a is now\n" << a.display_full();
   b = qr(a,c); std::cout << "qr(a,c) =\n" << b.display_full() << "c =\n" << c.display_full();
   a.mset("1 2 0 0 ; 2 -3 5 0 ; 0 5 4 -1 ; 0 0 -1 2");
   std::cout << "a is now\n" << a.display_full();
   dg = eigT(a); std::cout << "dg=eigT(a); dg.E = " << dispvect(dg.E) << "\n dg.V =\n" <<dg.V.display_full();
}*/
