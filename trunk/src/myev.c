#define VERYSMALL 1e-04

// subs to be able to check and directly diagonalize hermitean
// matrizes, inverse a nearly singular matrix
void myPrintComplexMatrix(FILE * file,ComplexMatrix & M)
{int i1,j1;
 fprintf (file,"#Real Part\n");
   
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
    for (j1=M.Clo();j1<=M.Chi();++j1) fprintf (file,"%6.3g ",real(M(i1,j1)));
    fprintf (file,"\n");
    }
    fprintf (file,"#Imaginary Part\n");
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
   for (j1=M.Clo();j1<=M.Chi();++j1)fprintf (file,"%6.3g ",imag(M(i1,j1)));
    fprintf (file,"\n");
    }
}    

void myPrintMatrix(FILE * file,Matrix & M)
{int i1,j1;
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
    for (j1=M.Clo();j1<=M.Chi();++j1) fprintf (file,"%6.3g ",M(i1,j1));
    fprintf (file,"\n");
    }
}    

void myPrintVector(FILE * file,Vector & M)
{int j1;
 fprintf (file,"#Components:\n");
   
    for (j1=M.Lo();j1<=M.Hi();++j1) fprintf (file,"%6.3g ",M(j1));
    fprintf (file,"\n");    
}    

void myEigenValuesHermitean (ComplexMatrix & M,Vector & lambda,int & sort,int & maxiter)
{ // this sub diagonalizes M and puts eigenvalues to lambda
  Matrix mat1(M.Rlo(),M.Rhi(),M.Clo(),M.Chi());
  int i1,j1;

  //check if M it is hermitean
   if (NormFro(M-M.Conjugate().Transpose())>VERYSMALL)
   {fprintf(stderr,"myEigenSystemHermitean: ERROR-matrix not hermitian\n");
    //printout matrix
    myPrintComplexMatrix(stderr,M);
    getchar();
   }


  // put matrix to format needed for library diagonalize function
   for(i1=M.Rlo();i1<=M.Rhi();++i1){for(j1=M.Clo();j1<=M.Chi();++j1){
    mat1(j1,i1)=imag(M(i1,j1)); 
    mat1(i1,j1)=real(M(i1,j1));
   }}
   EigenValuesHermitean (mat1,lambda,sort,maxiter);
   return;

}

void myEigenSystemHermitean (ComplexMatrix & M,Vector & lambda,ComplexMatrix & l,int & sort,int & maxiter)
{ // this sub diagonalizes M and puts eigenvalues to lambda, the eigenvectors to l
  Matrix mat1(M.Rlo(),M.Rhi(),M.Clo(),M.Chi());
  Matrix lr(l.Rlo(),l.Rhi(),l.Clo(),l.Chi());
  Matrix li(l.Rlo(),l.Rhi(),l.Clo(),l.Chi());
  complex<double> ii(0,1);  
  int i1,j1;

  //check if M it is hermitean
   if (NormFro(M-M.Conjugate().Transpose())>VERYSMALL)
   {fprintf(stderr,"myEigenSystemHermitean: ERROR-matrix not hermitian\n");
    //printout matrix
    myPrintComplexMatrix(stderr,M);
    getchar();
   }


  // put matrix to format needed for library diagonalize function
   for(i1=M.Rlo();i1<=M.Rhi();++i1){for(j1=M.Clo();j1<=M.Chi();++j1){
    mat1(j1,i1)=imag(M(i1,j1)); 
    mat1(i1,j1)=real(M(i1,j1));
   }}
 // setup matrix and diagonalize
 //  Driver routine to compute the  eigenvalues and normalized eigenvectors 
 //  of a complex Hermitian matrix z.The real parts of the elements must be
 //  stored in the lower triangle of z,the imaginary parts (of the elements
 //  corresponding to the lower triangle) in the positions
 //  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
 //  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
 //  True, otherwise  not ordered for sort = False. The real  and imaginary
 //  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
 //  The storage requirement is 3*n*n + 4*n complex numbers. 
 //  All matrices and vectors have to be allocated and removed by the user.
 //  They are checked for conformance !
 // void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
 // 			   int sort, int maxiter)
   EigenSystemHermitean (mat1,lambda,lr,li,sort,maxiter);
  // l=li;l*=ii;l+=lr;
  l=ComplexMatrix(li,lr);
   return;
}

