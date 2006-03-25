#include "ionpars.hpp"
#include "martin.h"
#include <cstring>
// ionpars: class to load and store matrices for internal module cfield
#include "ionpars.h"

#define NOF_OLM_MATRICES 45
#define MAXNOFCHARINLINE 1024

 ionpars::ionpars (const ionpars & p) //copy constructor
 {J=p.J;
  Ja=p.Ja; Jb=p.Jb; Jc=p.Jc;Hcf=p.Hcf;
  Jaa=p.Jaa; Jbb=p.Jbb; Jcc=p.Jcc;
  alpha=p.alpha;beta=p.beta;gamma=p.gamma;
  
 int i;
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 

 for(i=1;i<=NOF_OLM_MATRICES;++i)
 { Olm [i]= new Matrix(1,(*p.Olm[i]).Rhi(),1,(*p.Olm[i]).Chi()); 
 // define first matrix 
   OOlm [i] = new ComplexMatrix(1,(*p.OOlm[i]).Rhi(),1,(*p.OOlm[i]).Chi()); 
   (*Olm[i])=(*p.Olm[i]);
   (*OOlm[i])=(*p.OOlm[i]);
 }  

}
 
ionpars::ionpars (int dimj) // constructor from dimj
 {
  J=((double)dimj-1)/2;
  Ja=Matrix(1,dimj,1,dimj);
  Jb=Matrix(1,dimj,1,dimj);
  Jc=Matrix(1,dimj,1,dimj);
  Hcf=Matrix(1,dimj,1,dimj);
  Jaa=ComplexMatrix(1,dimj,1,dimj);
  Jbb=ComplexMatrix(1,dimj,1,dimj);
  Jcc=ComplexMatrix(1,dimj,1,dimj);

   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 

 int i;   
 for(i=1;i<=NOF_OLM_MATRICES;++i)
 { Olm [i]= new Matrix(1,dimj,1,dimj); 
 // define first matrix 
   OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
 }  

}

ionpars::~ionpars(){
 int i;
 for (i=1;i<=NOF_OLM_MATRICES;++i)
 {delete Olm[i];delete OOlm[i];}
  delete []Olm;
  delete []OOlm;
 } //destructor

ionpars::ionpars(FILE * cf_file) 
//constructor with commands from file handle (filename of cf parameters etc)
{ 
static int pr=1;
  FILE * tryfile;
  int dimj;complex<double> im(0,1);
  int i,j,dj=30; //30 ... maximum number of 2j+1
  char instr[MAXNOFCHARINLINE];
  char iontype[MAXNOFCHARINLINE];
//  char * cf_filename;

  // read in lines until  IONTYPE=
  while(fgets(instr,MAXNOFCHARINLINE,cf_file)!=NULL&&
        (i=extract(instr,"IONTYPE",iontype,(size_t)MAXNOFCHARINLINE))==1);
  
  
 // instr[strspn(instr," \t")]=='#');
  // take parameters from standard input file if no filename is given
 // if(instr[strspn(instr," \t")]=='#'||strlen(instr)-strspn(instr," \t")<=1)
 if(i==1){fprintf(stderr,"Error: no line in single ion property file contains IONTYPE field, e.g. IONTYPE=Nd3+\n");exit(EXIT_FAILURE);}
// get filename of parameter file out of first uncommented line in FILE * cf_file   
//  cf_filename=strtok(instr," \t\n");
  
  double ** hcfr,**hcfi,**Jxr,**Jxi,**Jyr,**Jyi,**Jzr,**Jzi;

  double ** mo22sr,**mo22si;
  double ** mo21sr,**mo21si;
  double ** mo20cr,**mo20ci;
  double ** mo21cr,**mo21ci;
  double ** mo22cr,**mo22ci;
  
  double ** mo33sr,**mo33si;
  double ** mo32sr,**mo32si;
  double ** mo31sr,**mo31si;
  double ** mo30cr,**mo30ci;
  double ** mo31cr,**mo31ci;
  double ** mo32cr,**mo32ci;
  double ** mo33cr,**mo33ci;

  double ** mo44sr,**mo44si;
  double ** mo43sr,**mo43si;
  double ** mo42sr,**mo42si;
  double ** mo41sr,**mo41si;
  double ** mo40cr,**mo40ci;
  double ** mo41cr,**mo41ci;
  double ** mo42cr,**mo42ci;
  double ** mo43cr,**mo43ci;
  double ** mo44cr,**mo44ci;
  
  double ** mo55sr,**mo55si;
  double ** mo54sr,**mo54si;
  double ** mo53sr,**mo53si;
  double ** mo52sr,**mo52si;
  double ** mo51sr,**mo51si;
  double ** mo50cr,**mo50ci;
  double ** mo51cr,**mo51ci;
  double ** mo52cr,**mo52ci;
  double ** mo53cr,**mo53ci;
  double ** mo54cr,**mo54ci;
  double ** mo55cr,**mo55ci;

  double ** mo66sr,**mo66si;
  double ** mo65sr,**mo65si;
  double ** mo64sr,**mo64si;
  double ** mo63sr,**mo63si;
  double ** mo62sr,**mo62si;
  double ** mo61sr,**mo61si;
  double ** mo60cr,**mo60ci;
  double ** mo61cr,**mo61ci;
  double ** mo62cr,**mo62ci;
  double ** mo63cr,**mo63ci;
  double ** mo64cr,**mo64ci;
  double ** mo65cr,**mo65ci;
  double ** mo66cr,**mo66ci;


    Jxr=new double*[dj+1];Jxi=new double*[dj+1];
    Jyr=new double*[dj+1];Jyi=new double*[dj+1];
    Jzr=new double*[dj+1];Jzi=new double*[dj+1];
    hcfr=new double*[dj+1];hcfi=new double*[dj+1];

  mo22sr=new double*[dj+1];mo22si=new double*[dj+1];
  mo21sr=new double*[dj+1];mo21si=new double*[dj+1];
  mo20cr=new double*[dj+1];mo20ci=new double*[dj+1];
  mo21cr=new double*[dj+1];mo21ci=new double*[dj+1];
  mo22cr=new double*[dj+1];mo22ci=new double*[dj+1];

  mo33sr=new double*[dj+1];mo33si=new double*[dj+1];
  mo32sr=new double*[dj+1];mo32si=new double*[dj+1];
  mo31sr=new double*[dj+1];mo31si=new double*[dj+1];
  mo30cr=new double*[dj+1];mo30ci=new double*[dj+1];
  mo31cr=new double*[dj+1];mo31ci=new double*[dj+1];
  mo32cr=new double*[dj+1];mo32ci=new double*[dj+1];
  mo33cr=new double*[dj+1];mo33ci=new double*[dj+1];

  mo44sr=new double*[dj+1];mo44si=new double*[dj+1];
  mo43sr=new double*[dj+1];mo43si=new double*[dj+1];
  mo42sr=new double*[dj+1];mo42si=new double*[dj+1];
  mo41sr=new double*[dj+1];mo41si=new double*[dj+1];
  mo40cr=new double*[dj+1];mo40ci=new double*[dj+1];
  mo41cr=new double*[dj+1];mo41ci=new double*[dj+1];
  mo42cr=new double*[dj+1];mo42ci=new double*[dj+1];
  mo43cr=new double*[dj+1];mo43ci=new double*[dj+1];
  mo44cr=new double*[dj+1];mo44ci=new double*[dj+1];

  mo55sr=new double*[dj+1];mo55si=new double*[dj+1];
  mo54sr=new double*[dj+1];mo54si=new double*[dj+1];
  mo53sr=new double*[dj+1];mo53si=new double*[dj+1];
  mo52sr=new double*[dj+1];mo52si=new double*[dj+1];
  mo51sr=new double*[dj+1];mo51si=new double*[dj+1];
  mo50cr=new double*[dj+1];mo50ci=new double*[dj+1];
  mo51cr=new double*[dj+1];mo51ci=new double*[dj+1];
  mo52cr=new double*[dj+1];mo52ci=new double*[dj+1];
  mo53cr=new double*[dj+1];mo53ci=new double*[dj+1];
  mo54cr=new double*[dj+1];mo54ci=new double*[dj+1];
  mo55cr=new double*[dj+1];mo55ci=new double*[dj+1];

  mo66sr=new double*[dj+1];mo66si=new double*[dj+1];
  mo65sr=new double*[dj+1];mo65si=new double*[dj+1];
  mo64sr=new double*[dj+1];mo64si=new double*[dj+1];
  mo63sr=new double*[dj+1];mo63si=new double*[dj+1];
  mo62sr=new double*[dj+1];mo62si=new double*[dj+1];
  mo61sr=new double*[dj+1];mo61si=new double*[dj+1];
  mo60cr=new double*[dj+1];mo60ci=new double*[dj+1];
  mo61cr=new double*[dj+1];mo61ci=new double*[dj+1];
  mo62cr=new double*[dj+1];mo62ci=new double*[dj+1];
  mo63cr=new double*[dj+1];mo63ci=new double*[dj+1];
  mo64cr=new double*[dj+1];mo64ci=new double*[dj+1];
  mo65cr=new double*[dj+1];mo65ci=new double*[dj+1];
  mo66cr=new double*[dj+1];mo66ci=new double*[dj+1];

    for (i=1;i<=dj;++i)
     {Jxr[i]=new double [dj];Jxi[i]=new double [dj];
      Jyr[i]=new double [dj];Jyi[i]=new double [dj];
      Jzr[i]=new double [dj];Jzi[i]=new double [dj];
      hcfr[i]=new double [dj];hcfi[i]=new double [dj];

    mo22sr[i]=new double [dj];mo22si[i]=new double [dj];
    mo21sr[i]=new double [dj];mo21si[i]=new double [dj];
    mo20cr[i]=new double [dj];mo20ci[i]=new double [dj];
    mo21cr[i]=new double [dj];mo21ci[i]=new double [dj];
    mo22cr[i]=new double [dj];mo22ci[i]=new double [dj];

    mo33sr[i]=new double [dj];mo33si[i]=new double [dj];
    mo32sr[i]=new double [dj];mo32si[i]=new double [dj];
    mo31sr[i]=new double [dj];mo31si[i]=new double [dj];
    mo30cr[i]=new double [dj];mo30ci[i]=new double [dj];
    mo31cr[i]=new double [dj];mo31ci[i]=new double [dj];
    mo32cr[i]=new double [dj];mo32ci[i]=new double [dj];
    mo33cr[i]=new double [dj];mo33ci[i]=new double [dj];

    mo44sr[i]=new double [dj];mo44si[i]=new double [dj];
    mo43sr[i]=new double [dj];mo43si[i]=new double [dj];
    mo42sr[i]=new double [dj];mo42si[i]=new double [dj];
    mo41sr[i]=new double [dj];mo41si[i]=new double [dj];
    mo40cr[i]=new double [dj];mo40ci[i]=new double [dj];
    mo41cr[i]=new double [dj];mo41ci[i]=new double [dj];
    mo42cr[i]=new double [dj];mo42ci[i]=new double [dj];
    mo43cr[i]=new double [dj];mo43ci[i]=new double [dj];
    mo44cr[i]=new double [dj];mo44ci[i]=new double [dj];

    mo55sr[i]=new double [dj];mo55si[i]=new double [dj];
    mo54sr[i]=new double [dj];mo54si[i]=new double [dj];
    mo53sr[i]=new double [dj];mo53si[i]=new double [dj];
    mo52sr[i]=new double [dj];mo52si[i]=new double [dj];
    mo51sr[i]=new double [dj];mo51si[i]=new double [dj];
    mo50cr[i]=new double [dj];mo50ci[i]=new double [dj];
    mo51cr[i]=new double [dj];mo51ci[i]=new double [dj];
    mo52cr[i]=new double [dj];mo52ci[i]=new double [dj];
    mo53cr[i]=new double [dj];mo53ci[i]=new double [dj];
    mo54cr[i]=new double [dj];mo54ci[i]=new double [dj];
    mo55cr[i]=new double [dj];mo55ci[i]=new double [dj];

    mo66sr[i]=new double [dj];mo66si[i]=new double [dj];
    mo65sr[i]=new double [dj];mo65si[i]=new double [dj];
    mo64sr[i]=new double [dj];mo64si[i]=new double [dj];
    mo63sr[i]=new double [dj];mo63si[i]=new double [dj];
    mo62sr[i]=new double [dj];mo62si[i]=new double [dj];
    mo61sr[i]=new double [dj];mo61si[i]=new double [dj];
    mo60cr[i]=new double [dj];mo60ci[i]=new double [dj];
    mo61cr[i]=new double [dj];mo61ci[i]=new double [dj];
    mo62cr[i]=new double [dj];mo62ci[i]=new double [dj];
    mo63cr[i]=new double [dj];mo63ci[i]=new double [dj];
    mo64cr[i]=new double [dj];mo64ci[i]=new double [dj];
    mo65cr[i]=new double [dj];mo65ci[i]=new double [dj];
    mo66cr[i]=new double [dj];mo66ci[i]=new double [dj];
    }
     

if (pr==1) printf("#using cfield ...\n");

  
  fprintf(stderr,"# module cfield ... for ion %s\n",iontype);
  cfield_mcphasnew(iontype,Jxr,Jxi,  Jyr, Jyi, Jzr, Jzi,
  mo22sr,mo22si,
  mo21sr,mo21si,
  mo20cr,mo20ci,
  mo21cr,mo21ci,
  mo22cr,mo22ci,

  mo33sr,mo33si,
  mo32sr,mo32si,
  mo31sr,mo31si,
  mo30cr,mo30ci,
  mo31cr,mo31ci,
  mo32cr,mo32ci,
  mo33cr,mo33ci,

  mo44sr,mo44si,
  mo43sr,mo43si,
  mo42sr,mo42si,
  mo41sr,mo41si,
  mo40cr,mo40ci,
  mo41cr,mo41ci,
  mo42cr,mo42ci,
  mo43cr,mo43ci,
  mo44cr,mo44ci,

  mo55sr,mo55si,
  mo54sr,mo54si,
  mo53sr,mo53si,
  mo52sr,mo52si,
  mo51sr,mo51si,
  mo50cr,mo50ci,
  mo51cr,mo51ci,
  mo52cr,mo52ci,
  mo53cr,mo53ci,
  mo54cr,mo54ci,
  mo55cr,mo55ci,

  mo66sr,mo66si,
  mo65sr,mo65si,
  mo64sr,mo64si,
  mo63sr,mo63si,
  mo62sr,mo62si,
  mo61sr,mo61si,
  mo60cr,mo60ci,
  mo61cr,mo61ci,
  mo62cr,mo62ci,
  mo63cr,mo63ci,
  mo64cr,mo64ci,
  mo65cr,mo65ci,
  mo66cr,mo66ci,

  &dimj,&alpha,&beta,&gamma);

  Hcf= Matrix(1,dimj,1,dimj); 
  Hcf=0;


if (pr==1) printf("#end using cfield\n");

   J=((double)dimj-1)/2; //momentum quantum number

if (pr==1) printf("#J=%g\n",J);

   Ja = Matrix(1,dimj,1,dimj); 
   Jaa = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jaa(i,j)=im*(Jxi[i])[j]+(Jxr[i])[j];
    if(i<j){Ja(i,j)=(Jxi[j])[i];}else{Ja(i,j)=(Jxr[i])[j];}
   }

   Jb = Matrix(1,dimj,1,dimj); 
   Jbb = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jbb(i,j)=im*(Jyi[i])[j]+(Jyr[i])[j];
    if(i<j){Jb(i,j)=(Jyi[j])[i];}else{Jb(i,j)=(Jyr[i])[j];}
   }

   Jc = Matrix(1,dimj,1,dimj); 
   Jcc = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jcc(i,j)=im*(Jzi[i])[j]+(Jzr[i])[j];
    if(i<j){Jc(i,j)=(Jzi[j])[i];}else{Jc(i,j)=(Jzr[i])[j];}
   }

//---------------------------------------------------------------------------

   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
   
   for(i=1;i<=NOF_OLM_MATRICES;++i)   
    {   Olm[i]= new Matrix(1,dimj,1,dimj); 
 // define memory for all matrices 
        OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
    }   
    
    
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {(*OOlm[1])(i,j)=im*(mo22si[i])[j]+(mo22sr[i])[j];
if(i<j){(*Olm[1])(i,j)=(mo22si[j])[i];}else{(*Olm[1])(i,j)=(mo22sr[i])[j];}
    (*OOlm[2])(i,j)=im*(mo21si[i])[j]+(mo21sr[i])[j];
if(i<j){(*Olm[2])(i,j)=(mo21si[j])[i];}else{(*Olm[2])(i,j)=(mo21sr[i])[j];}
    (*OOlm[3])(i,j)=im*(mo20ci[i])[j]+(mo20cr[i])[j];
if(i<j){(*Olm[3])(i,j)=(mo20ci[j])[i];}else{(*Olm[3])(i,j)=(mo20cr[i])[j];}
    (*OOlm[4])(i,j)=im*(mo21ci[i])[j]+(mo21cr[i])[j];
if(i<j){(*Olm[4])(i,j)=(mo21ci[j])[i];}else{(*Olm[4])(i,j)=(mo21cr[i])[j];}
    (*OOlm[5])(i,j)=im*(mo22ci[i])[j]+(mo22cr[i])[j];
if(i<j){(*Olm[5])(i,j)=(mo22ci[j])[i];}else{(*Olm[5])(i,j)=(mo22cr[i])[j];}
    
    (*OOlm[6])(i,j)=im*(mo33si[i])[j]+(mo33sr[i])[j];
if(i<j){(*Olm[6])(i,j)=(mo33si[j])[i];}else{(*Olm[6])(i,j)=(mo33sr[i])[j];}
    (*OOlm[7])(i,j)=im*(mo32si[i])[j]+(mo32sr[i])[j];
if(i<j){(*Olm[7])(i,j)=(mo32si[j])[i];}else{(*Olm[7])(i,j)=(mo32sr[i])[j];}
    (*OOlm[8])(i,j)=im*(mo31si[i])[j]+(mo31sr[i])[j];
if(i<j){(*Olm[8])(i,j)=(mo31si[j])[i];}else{(*Olm[8])(i,j)=(mo31sr[i])[j];}
    (*OOlm[9])(i,j)=im*(mo30ci[i])[j]+(mo30cr[i])[j];
if(i<j){(*Olm[9])(i,j)=(mo30ci[j])[i];}else{(*Olm[9])(i,j)=(mo30cr[i])[j];}
    (*OOlm[10])(i,j)=im*(mo31ci[i])[j]+(mo31cr[i])[j];
if(i<j){(*Olm[10])(i,j)=(mo31ci[j])[i];}else{(*Olm[10])(i,j)=(mo31cr[i])[j];}
    (*OOlm[11])(i,j)=im*(mo32ci[i])[j]+(mo32cr[i])[j];
if(i<j){(*Olm[11])(i,j)=(mo32ci[j])[i];}else{(*Olm[11])(i,j)=(mo32cr[i])[j];}
    (*OOlm[12])(i,j)=im*(mo33ci[i])[j]+(mo33cr[i])[j];
if(i<j){(*Olm[12])(i,j)=(mo33ci[j])[i];}else{(*Olm[12])(i,j)=(mo33cr[i])[j];}
    
    (*OOlm[13])(i,j)=im*(mo44si[i])[j]+(mo44sr[i])[j];
if(i<j){(*Olm[13])(i,j)=(mo44si[j])[i];}else{(*Olm[13])(i,j)=(mo44sr[i])[j];}
    (*OOlm[14])(i,j)=im*(mo43si[i])[j]+(mo43sr[i])[j];
if(i<j){(*Olm[14])(i,j)=(mo43si[j])[i];}else{(*Olm[14])(i,j)=(mo43sr[i])[j];}
    (*OOlm[15])(i,j)=im*(mo42si[i])[j]+(mo42sr[i])[j];
if(i<j){(*Olm[15])(i,j)=(mo42si[j])[i];}else{(*Olm[15])(i,j)=(mo42sr[i])[j];}
    (*OOlm[16])(i,j)=im*(mo41si[i])[j]+(mo41sr[i])[j];
if(i<j){(*Olm[16])(i,j)=(mo41si[j])[i];}else{(*Olm[16])(i,j)=(mo41sr[i])[j];}
    (*OOlm[17])(i,j)=im*(mo40ci[i])[j]+(mo40cr[i])[j];
if(i<j){(*Olm[17])(i,j)=(mo40ci[j])[i];}else{(*Olm[17])(i,j)=(mo40cr[i])[j];}
    (*OOlm[18])(i,j)=im*(mo41ci[i])[j]+(mo41cr[i])[j];
if(i<j){(*Olm[18])(i,j)=(mo41ci[j])[i];}else{(*Olm[18])(i,j)=(mo41cr[i])[j];}
    (*OOlm[19])(i,j)=im*(mo42ci[i])[j]+(mo42cr[i])[j];
if(i<j){(*Olm[19])(i,j)=(mo42ci[j])[i];}else{(*Olm[19])(i,j)=(mo42cr[i])[j];}
    (*OOlm[20])(i,j)=im*(mo43ci[i])[j]+(mo43cr[i])[j];
if(i<j){(*Olm[20])(i,j)=(mo43ci[j])[i];}else{(*Olm[20])(i,j)=(mo43cr[i])[j];}
    (*OOlm[21])(i,j)=im*(mo44ci[i])[j]+(mo44cr[i])[j];
if(i<j){(*Olm[21])(i,j)=(mo44ci[j])[i];}else{(*Olm[21])(i,j)=(mo44cr[i])[j];}
    
    (*OOlm[22])(i,j)=im*(mo55si[i])[j]+(mo55sr[i])[j];
if(i<j){(*Olm[22])(i,j)=(mo55si[j])[i];}else{(*Olm[22])(i,j)=(mo55sr[i])[j];}
    (*OOlm[23])(i,j)=im*(mo54si[i])[j]+(mo54sr[i])[j];
if(i<j){(*Olm[23])(i,j)=(mo54si[j])[i];}else{(*Olm[23])(i,j)=(mo54sr[i])[j];}
    (*OOlm[24])(i,j)=im*(mo53si[i])[j]+(mo53sr[i])[j];
if(i<j){(*Olm[24])(i,j)=(mo53si[j])[i];}else{(*Olm[24])(i,j)=(mo53sr[i])[j];}
    (*OOlm[25])(i,j)=im*(mo52si[i])[j]+(mo52sr[i])[j];
if(i<j){(*Olm[25])(i,j)=(mo52si[j])[i];}else{(*Olm[25])(i,j)=(mo52sr[i])[j];}
    (*OOlm[26])(i,j)=im*(mo51si[i])[j]+(mo51sr[i])[j];
if(i<j){(*Olm[26])(i,j)=(mo51si[j])[i];}else{(*Olm[26])(i,j)=(mo51sr[i])[j];}
    (*OOlm[27])(i,j)=im*(mo50ci[i])[j]+(mo50cr[i])[j];
if(i<j){(*Olm[27])(i,j)=(mo50ci[j])[i];}else{(*Olm[27])(i,j)=(mo50cr[i])[j];}
    (*OOlm[28])(i,j)=im*(mo51ci[i])[j]+(mo51cr[i])[j];
if(i<j){(*Olm[28])(i,j)=(mo51ci[j])[i];}else{(*Olm[28])(i,j)=(mo51cr[i])[j];}
    (*OOlm[29])(i,j)=im*(mo52ci[i])[j]+(mo52cr[i])[j];
if(i<j){(*Olm[29])(i,j)=(mo52ci[j])[i];}else{(*Olm[29])(i,j)=(mo52cr[i])[j];}
    (*OOlm[30])(i,j)=im*(mo53ci[i])[j]+(mo53cr[i])[j];
if(i<j){(*Olm[30])(i,j)=(mo53ci[j])[i];}else{(*Olm[30])(i,j)=(mo53cr[i])[j];}
    (*OOlm[31])(i,j)=im*(mo54ci[i])[j]+(mo54cr[i])[j];
if(i<j){(*Olm[31])(i,j)=(mo54ci[j])[i];}else{(*Olm[31])(i,j)=(mo54cr[i])[j];}
    (*OOlm[32])(i,j)=im*(mo55ci[i])[j]+(mo55cr[i])[j];
if(i<j){(*Olm[32])(i,j)=(mo55ci[j])[i];}else{(*Olm[32])(i,j)=(mo55cr[i])[j];}
    
    (*OOlm[33])(i,j)=im*(mo66si[i])[j]+(mo66sr[i])[j];
if(i<j){(*Olm[33])(i,j)=(mo66si[j])[i];}else{(*Olm[33])(i,j)=(mo66sr[i])[j];}
    (*OOlm[34])(i,j)=im*(mo65si[i])[j]+(mo65sr[i])[j];
if(i<j){(*Olm[34])(i,j)=(mo65si[j])[i];}else{(*Olm[34])(i,j)=(mo65sr[i])[j];}
    (*OOlm[35])(i,j)=im*(mo64si[i])[j]+(mo64sr[i])[j];
if(i<j){(*Olm[35])(i,j)=(mo64si[j])[i];}else{(*Olm[35])(i,j)=(mo64sr[i])[j];}
    (*OOlm[36])(i,j)=im*(mo63si[i])[j]+(mo63sr[i])[j];
if(i<j){(*Olm[36])(i,j)=(mo63si[j])[i];}else{(*Olm[36])(i,j)=(mo63sr[i])[j];}
    (*OOlm[37])(i,j)=im*(mo62si[i])[j]+(mo62sr[i])[j];
if(i<j){(*Olm[37])(i,j)=(mo62si[j])[i];}else{(*Olm[37])(i,j)=(mo62sr[i])[j];}
    (*OOlm[38])(i,j)=im*(mo61si[i])[j]+(mo61sr[i])[j];
if(i<j){(*Olm[38])(i,j)=(mo61si[j])[i];}else{(*Olm[38])(i,j)=(mo61sr[i])[j];}
    (*OOlm[39])(i,j)=im*(mo60ci[i])[j]+(mo60cr[i])[j];
if(i<j){(*Olm[39])(i,j)=(mo60ci[j])[i];}else{(*Olm[39])(i,j)=(mo60cr[i])[j];}
    (*OOlm[40])(i,j)=im*(mo61ci[i])[j]+(mo61cr[i])[j];
if(i<j){(*Olm[40])(i,j)=(mo61ci[j])[i];}else{(*Olm[40])(i,j)=(mo61cr[i])[j];}
    (*OOlm[41])(i,j)=im*(mo62ci[i])[j]+(mo62cr[i])[j];
if(i<j){(*Olm[41])(i,j)=(mo62ci[j])[i];}else{(*Olm[41])(i,j)=(mo62cr[i])[j];}
    (*OOlm[42])(i,j)=im*(mo63ci[i])[j]+(mo63cr[i])[j];
if(i<j){(*Olm[42])(i,j)=(mo63ci[j])[i];}else{(*Olm[42])(i,j)=(mo63cr[i])[j];}
    (*OOlm[43])(i,j)=im*(mo64ci[i])[j]+(mo64cr[i])[j];
if(i<j){(*Olm[43])(i,j)=(mo64ci[j])[i];}else{(*Olm[43])(i,j)=(mo64cr[i])[j];}
    (*OOlm[44])(i,j)=im*(mo65ci[i])[j]+(mo65cr[i])[j];
if(i<j){(*Olm[44])(i,j)=(mo65ci[j])[i];}else{(*Olm[44])(i,j)=(mo65cr[i])[j];}
    (*OOlm[45])(i,j)=im*(mo66ci[i])[j]+(mo66cr[i])[j];
if(i<j){(*Olm[45])(i,j)=(mo66ci[j])[i];}else{(*Olm[45])(i,j)=(mo66cr[i])[j];}
    
   }}

// ------------------------------------------------------------
   


   
    for (i=1;i<=dj;++i)
     {delete[]Jxr[i];delete[]Jxi[i];
      delete[]Jyr[i];delete[]Jyi[i];
      delete[]Jzr[i];delete[]Jzi[i];
      delete[]hcfr[i];delete[]hcfi[i];

   delete[]mo22sr[i];delete[]mo22si[i];
   delete[]mo21sr[i];delete[]mo21si[i];
   delete[]mo20cr[i];delete[]mo20ci[i];
   delete[]mo21cr[i];delete[]mo21ci[i];
   delete[]mo22cr[i];delete[]mo22ci[i];

   delete[]mo33sr[i];delete[]mo33si[i];
   delete[]mo32sr[i];delete[]mo32si[i];
   delete[]mo31sr[i];delete[]mo31si[i];
   delete[]mo30cr[i];delete[]mo30ci[i];
   delete[]mo31cr[i];delete[]mo31ci[i];
   delete[]mo32cr[i];delete[]mo32ci[i];
   delete[]mo33cr[i];delete[]mo33ci[i];

   delete[]mo44sr[i];delete[]mo44si[i];
   delete[]mo43sr[i];delete[]mo43si[i];
   delete[]mo42sr[i];delete[]mo42si[i];
   delete[]mo41sr[i];delete[]mo41si[i];
   delete[]mo40cr[i];delete[]mo40ci[i];
   delete[]mo41cr[i];delete[]mo41ci[i];
   delete[]mo42cr[i];delete[]mo42ci[i];
   delete[]mo43cr[i];delete[]mo43ci[i];
   delete[]mo44cr[i];delete[]mo44ci[i];

   delete[]mo55sr[i];delete[]mo55si[i];
   delete[]mo54sr[i];delete[]mo54si[i];
   delete[]mo53sr[i];delete[]mo53si[i];
   delete[]mo52sr[i];delete[]mo52si[i];
   delete[]mo51sr[i];delete[]mo51si[i];
   delete[]mo50cr[i];delete[]mo50ci[i];
   delete[]mo51cr[i];delete[]mo51ci[i];
   delete[]mo52cr[i];delete[]mo52ci[i];
   delete[]mo53cr[i];delete[]mo53ci[i];
   delete[]mo54cr[i];delete[]mo54ci[i];
   delete[]mo55cr[i];delete[]mo55ci[i];

   delete[]mo66sr[i];delete[]mo66si[i];
   delete[]mo65sr[i];delete[]mo65si[i];
   delete[]mo64sr[i];delete[]mo64si[i];
   delete[]mo63sr[i];delete[]mo63si[i];
   delete[]mo62sr[i];delete[]mo62si[i];
   delete[]mo61sr[i];delete[]mo61si[i];
   delete[]mo60cr[i];delete[]mo60ci[i];
   delete[]mo61cr[i];delete[]mo61ci[i];
   delete[]mo62cr[i];delete[]mo62ci[i];
   delete[]mo63cr[i];delete[]mo63ci[i];
   delete[]mo64cr[i];delete[]mo64ci[i];
   delete[]mo65cr[i];delete[]mo65ci[i];
   delete[]mo66cr[i];delete[]mo66ci[i];
   }

     delete[]Jxr;delete[]Jxi;
     delete[]Jyr;delete[]Jyi;
     delete[]Jzr;delete[]Jzi;
     delete[]hcfr;delete[]hcfi;

   delete []mo22sr;delete []mo22si;
   delete []mo21sr;delete []mo21si;
   delete []mo20cr;delete []mo20ci;
   delete []mo21cr;delete []mo21ci;
   delete []mo22cr;delete []mo22ci;

   delete []mo33sr;delete []mo33si;
   delete []mo32sr;delete []mo32si;
   delete []mo31sr;delete []mo31si;
   delete []mo30cr;delete []mo30ci;
   delete []mo31cr;delete []mo31ci;
   delete []mo32cr;delete []mo32ci;
   delete []mo33cr;delete []mo33ci;

   delete []mo44sr;delete []mo44si;
   delete []mo43sr;delete []mo43si;
   delete []mo42sr;delete []mo42si;
   delete []mo41sr;delete []mo41si;
   delete []mo40cr;delete []mo40ci;
   delete []mo41cr;delete []mo41ci;
   delete []mo42cr;delete []mo42ci;
   delete []mo43cr;delete []mo43ci;
   delete []mo44cr;delete []mo44ci;

   delete []mo55sr;delete []mo55si;
   delete []mo54sr;delete []mo54si;
   delete []mo53sr;delete []mo53si;
   delete []mo52sr;delete []mo52si;
   delete []mo51sr;delete []mo51si;
   delete []mo50cr;delete []mo50ci;
   delete []mo51cr;delete []mo51ci;
   delete []mo52cr;delete []mo52ci;
   delete []mo53cr;delete []mo53ci;
   delete []mo54cr;delete []mo54ci;
   delete []mo55cr;delete []mo55ci;

   delete []mo66sr;delete []mo66si;
   delete []mo65sr;delete []mo65si;
   delete []mo64sr;delete []mo64si;
   delete []mo63sr;delete []mo63si;
   delete []mo62sr;delete []mo62si;
   delete []mo61sr;delete []mo61si;
   delete []mo60cr;delete []mo60ci;
   delete []mo61cr;delete []mo61ci;
   delete []mo62cr;delete []mo62ci;
   delete []mo63cr;delete []mo63ci;
   delete []mo64cr;delete []mo64ci;
   delete []mo65cr;delete []mo65ci;
   delete []mo66cr;delete []mo66ci;
   
   
//ATTENTION FOR NDCU2 the AXES xyz are parallel to cab
Matrix dummy(1,dimj,1,dimj);
dummy=Jb;Jb=Jc;Jc=Ja;Ja=dummy;
ComplexMatrix dummyc(1,dimj,1,dimj);
dummyc=Jbb;Jbb=Jcc;Jcc=Jaa;Jaa=dummyc;

if (pr==1) {printf("#Axis Convention using cfield as a module:  a||y b||z  c||x\n");
printf("#xyz .... Coordinate system of the crystal field parameters used in cfield\n");
printf("#abc .... Crystal axes\n");
printf("#The interactions are described by the  PKQ Operators defined in cfield\n");
printf("#O11(s) .... Ja=Jy\n");
printf("#O10(c) .... Jb=Jz\n");
printf("#O11(c) .... Jc=Jx\n");
printf("#O22(s) .... Jd\n");
printf("#O21(s) .... Je\n");
printf("#O20(c) .... Jf\n");
printf("#O21(c) .... Jg\n");
printf("#O22(c) .... Jh\n");
printf("#O33(s) .... Ji\n");
printf("#O32(s) .... Jj\n");
printf("#O31(s) .... Jk\n");
printf("#O30(c) .... Jl\n");
printf("#O31(c) .... Jm\n");
printf("# etc ... 45 moments up to l<=6\n");
printf("#\n");
}

pr=0;
}
