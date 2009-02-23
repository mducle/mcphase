// ionpars: class to load and store matrices for internal module cfield

#include "ionpars.hpp"
#include "martin.h"
#include <cstring>
#include <cstdlib>
#include "ionpars.h"

#define PI 3.1415926535
#define NOF_OLM_MATRICES 45
#define MAXNOFCHARINLINE 1024
#define K_B  0.0862
#define SMALL 1e-6   //!!! must match SMALL in mcdisp.c and ionpars.cpp !!!
                     // because it is used to decide wether for small transition
		     // energy the matrix Mijkl contains wn-wn' or wn/kT

void myPrintComplexMat(FILE * file,ComplexMatrix & M)
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


 ionpars::ionpars (const ionpars & p) //copy constructor
 {J=p.J;
  Ja=p.Ja; Jb=p.Jb; Jc=p.Jc;Hcf=p.Hcf;
  Jaa=p.Jaa; Jbb=p.Jbb; Jcc=p.Jcc;
  gJ=p.gJ;
  alpha=p.alpha;beta=p.beta;gamma=p.gamma;
  r2=p.r2;r4=p.r4;r6=p.r6;
  Blm=p.Blm; // vector of crystal field parameters
  Llm=p.Llm; // vector of crystal field parameters
 
  Np=p.Np; Xip=p.Xip;Cp=p.Cp;
  
   int i;
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
   iontype = new char [strlen(p.iontype)+1];
   strcpy(iontype,p.iontype);
  

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

   Blm=Vector(1,45);Blm=0; // vector of crystal field parameters
   Llm=Vector(1,45);Llm=0; // vector of crystal field parameters

   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;
   alpha=0;beta=0;gamma=0;r2=0;r4=0;r6=0;

   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
   iontype = new char [MAXNOFCHARINLINE];


 int i;   
 for(i=1;i<=NOF_OLM_MATRICES;++i)
 { Olm [i]= new Matrix(1,dimj,1,dimj); 
 // define first matrix 
   OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
 }  

}

 
ionpars::ionpars (char * ion) // constructor from iontype (mind:no matrices filled with values !)
 {int dimj;
  getpar(ion, &dimj, &alpha, &beta, &gamma, &gJ,&r2, &r4,&r6 );
   iontype = new char [strlen(ion)+1];
   strcpy(iontype,ion);

  J=((double)dimj-1)/2;
  Ja=Matrix(1,dimj,1,dimj);
  Jb=Matrix(1,dimj,1,dimj);
  Jc=Matrix(1,dimj,1,dimj);
  Hcf=Matrix(1,dimj,1,dimj);
  Jaa=ComplexMatrix(1,dimj,1,dimj);
  Jbb=ComplexMatrix(1,dimj,1,dimj);
  Jcc=ComplexMatrix(1,dimj,1,dimj);

   Blm=Vector(1,45);Blm=0; // vector of crystal field parameters
   Llm=Vector(1,45);Llm=0; // vector of crystal field parameters

   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;

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
 delete []iontype;
 for (i=1;i<=NOF_OLM_MATRICES;++i)
 {delete Olm[i];delete OOlm[i];}
  delete []Olm;
  delete []OOlm;
 } //destructor


ionpars::ionpars(FILE * cf_file) 
//constructor with commands from file handle (filename of cf parameters etc)
{ 
static int pr=1;
//  FILE * tryfile;
  int dimj;complex<double> im(0,1);
  int i,j,l,dj=30; //30 ... maximum number of 2j+1
  char instr[MAXNOFCHARINLINE];
  iontype= new char[MAXNOFCHARINLINE];
  
   Blm=Vector(1,45);Blm=0; // vector of crystal field parameters
   Llm=Vector(1,45);Llm=0; // vector of crystal field parameters

   Np=Vector(1,9);Np=0; // vectors of radial wave function parameters
   Xip=Vector(1,9);Xip=0;
   Cp=Vector(1,9);Cp=0;
   alpha=0;beta=0;gamma=0;r2=0;r4=0;r6=0;

  // read in lines and get IONTYPE=  and CF parameters Blm
   while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
        extract(instr,"IONTYPE",iontype,(size_t)MAXNOFCHARINLINE);
        
        extract(instr,"N1",Np(1));extract(instr,"XI1",Xip(1));extract(instr,"C1",Cp(1));
        extract(instr,"N2",Np(2));extract(instr,"XI2",Xip(2));extract(instr,"C2",Cp(2));
        extract(instr,"N3",Np(3));extract(instr,"XI3",Xip(3));extract(instr,"C3",Cp(3));
        extract(instr,"N4",Np(4));extract(instr,"XI4",Xip(4));extract(instr,"C4",Cp(4));
        extract(instr,"N5",Np(5));extract(instr,"XI5",Xip(5));extract(instr,"C5",Cp(5));
        extract(instr,"N6",Np(6));extract(instr,"XI6",Xip(6));extract(instr,"C6",Cp(6));
        extract(instr,"N7",Np(7));extract(instr,"XI7",Xip(7));extract(instr,"C7",Cp(7));
        extract(instr,"N8",Np(8));extract(instr,"XI8",Xip(8));extract(instr,"C8",Cp(8));
        extract(instr,"N9",Np(9));extract(instr,"XI9",Xip(9));extract(instr,"C9",Cp(9));

        extract(instr,"ALPHA",alpha);
        extract(instr,"BETA",beta);
        extract(instr,"GAMMA",gamma);

        extract(instr,"R2",r2);
        extract(instr,"R4",r4);
        extract(instr,"R6",r6);

        extract(instr,"B22S",Blm(1));
        extract(instr,"B21S",Blm(2));
	extract(instr,"B20",Blm(3));
        extract(instr,"B21",Blm(4));
	extract(instr,"B22",Blm(5));
   
	extract(instr,"B33S",Blm(6));
	extract(instr,"B32S",Blm(7));
	extract(instr,"B31S",Blm(8));
	extract(instr,"B30",Blm(9));
   extract(instr,"B31",Blm(10));
   extract(instr,"B32",Blm(11));
   extract(instr,"B32",Blm(12));

   extract(instr,"B44S",Blm(13));
   extract(instr,"B43S",Blm(14));
   extract(instr,"B42S",Blm(15));
   extract(instr,"B41S",Blm(16));
   extract(instr,"B40",Blm(17));
   extract(instr,"B41",Blm(18));
   extract(instr,"B42",Blm(19));
   extract(instr,"B43",Blm(20));
   extract(instr,"B44",Blm(21));
  
   extract(instr,"B55S",Blm(22));
   extract(instr,"B54S",Blm(23));
   extract(instr,"B53S",Blm(24));
   extract(instr,"B52S",Blm(25));
   extract(instr,"B51S",Blm(26));
   extract(instr,"B50",Blm(27));
   extract(instr,"B51",Blm(28));
   extract(instr,"B52",Blm(29));
   extract(instr,"B53",Blm(30));
   extract(instr,"B54",Blm(31));
   extract(instr,"B55",Blm(32));
 
   extract(instr,"B66S",Blm(33));
   extract(instr,"B65S",Blm(34));
   extract(instr,"B64S",Blm(35));
   extract(instr,"B63S",Blm(36));
   extract(instr,"B62S",Blm(37));
   extract(instr,"B61S",Blm(38));
   extract(instr,"B60",Blm(39));
   extract(instr,"B61",Blm(40));
   extract(instr,"B62",Blm(41));
   extract(instr,"B63",Blm(42));
   extract(instr,"B64",Blm(43));
   extract(instr,"B65",Blm(44));
   extract(instr,"B66",Blm(45));

        extract(instr,"L22S",Llm(1));
        extract(instr,"L21S",Llm(2));
	extract(instr,"L20",Llm(3));
        extract(instr,"L21",Llm(4));
	extract(instr,"L22",Llm(5));
   
	extract(instr,"L33S",Llm(6));
	extract(instr,"L32S",Llm(7));
	extract(instr,"L31S",Llm(8));
	extract(instr,"L30",Llm(9));
   extract(instr,"L31",Llm(10));
   extract(instr,"L32",Llm(11));
   extract(instr,"L32",Llm(12));

   extract(instr,"L44S",Llm(13));
   extract(instr,"L43S",Llm(14));
   extract(instr,"L42S",Llm(15));
   extract(instr,"L41S",Llm(16));
   extract(instr,"L40",Llm(17));
   extract(instr,"L41",Llm(18));
   extract(instr,"L42",Llm(19));
   extract(instr,"L43",Llm(20));
   extract(instr,"L44",Llm(21));
  
   extract(instr,"L55S",Llm(22));
   extract(instr,"L54S",Llm(23));
   extract(instr,"L53S",Llm(24));
   extract(instr,"L52S",Llm(25));
   extract(instr,"L51S",Llm(26));
   extract(instr,"L50",Llm(27));
   extract(instr,"L51",Llm(28));
   extract(instr,"L52",Llm(29));
   extract(instr,"L53",Llm(30));
   extract(instr,"L54",Llm(31));
   extract(instr,"L55",Llm(32));
 
   extract(instr,"L66S",Llm(33));
   extract(instr,"L65S",Llm(34));
   extract(instr,"L64S",Llm(35));
   extract(instr,"L63S",Llm(36));
   extract(instr,"L62S",Llm(37));
   extract(instr,"L61S",Llm(38));
   extract(instr,"L60",Llm(39));
   extract(instr,"L61",Llm(40));
   extract(instr,"L62",Llm(41));
   extract(instr,"L63",Llm(42));
   extract(instr,"L64",Llm(43));
   extract(instr,"L65",Llm(44));
   extract(instr,"L66",Llm(45));


	}}
  
  
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

  &dimj,&alpha,&beta,&gamma,&gJ,&r2,&r4,&r6);



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
   
  Hcf= Matrix(1,dimj,1,dimj); 
  Hcf=0;

   if(Hcf==(double)0.0){
   // calculation of the cf matrix according 
   fprintf(stderr,"crystal field parameters\n");  
   for(l=1;l<=45;++l){Hcf+=Blm(l)*(*Olm[l]);
                   if(Blm(l)!=0){if(l<24){fprintf(stderr,"B%c=%g   ",l+99,Blm(l));}
		                     else{fprintf(stderr,"B(z+%i)=%g   ",l-23,Blm(l));}
		                }
                  }
   }
   
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




void ionpars::savBlm(FILE * outfile)
{fprintf(outfile,"units=meV\n");
   if(Blm(1)!=0){fprintf(outfile,"B22S=%g\n",Blm(1));}
   if(Blm(2)!=0){fprintf(outfile,"B21S=%g\n",Blm(2));}
   if(Blm(3)!=0){fprintf(outfile,"B20=%g\n",Blm(3));}
   if(Blm(4)!=0){fprintf(outfile,"B21=%g\n",Blm(4));}
   if(Blm(5)!=0){fprintf(outfile,"B22=%g\n",Blm(5));}
   
   if(Blm(6)!=0){fprintf(outfile,"B33S=%g\n",Blm(6));}
   if(Blm(7)!=0){fprintf(outfile,"B32S=%g\n",Blm(7));}
   if(Blm(8)!=0){fprintf(outfile,"B31S=%g\n",Blm(8));}
   if(Blm(9)!=0){fprintf(outfile,"B30=%g\n",Blm(9));}
   if(Blm(10)!=0){fprintf(outfile,"B31=%g\n",Blm(10));}
   if(Blm(11)!=0){fprintf(outfile,"B32=%g\n",Blm(11));}
   if(Blm(12)!=0){fprintf(outfile,"B32=%g\n",Blm(12));}

   if(Blm(13)!=0){fprintf(outfile,"B44S=%g\n",Blm(13));}
   if(Blm(14)!=0){fprintf(outfile,"B43S=%g\n",Blm(14));}
   if(Blm(15)!=0){fprintf(outfile,"B42S=%g\n",Blm(15));}
   if(Blm(16)!=0){fprintf(outfile,"B41S=%g\n",Blm(16));}
   if(Blm(17)!=0){fprintf(outfile,"B40=%g\n",Blm(17));}
   if(Blm(18)!=0){fprintf(outfile,"B41=%g\n",Blm(18));}
   if(Blm(19)!=0){fprintf(outfile,"B42=%g\n",Blm(19));}
   if(Blm(20)!=0){fprintf(outfile,"B43=%g\n",Blm(20));}
   if(Blm(21)!=0){fprintf(outfile,"B44=%g\n",Blm(21));}
  
   if(Blm(22)!=0){fprintf(outfile,"B55S=%g\n",Blm(22));}
   if(Blm(23)!=0){fprintf(outfile,"B54S=%g\n",Blm(23));}
   if(Blm(24)!=0){fprintf(outfile,"B53S=%g\n",Blm(24));}
   if(Blm(25)!=0){fprintf(outfile,"B52S=%g\n",Blm(25));}
   if(Blm(26)!=0){fprintf(outfile,"B51S=%g\n",Blm(26));}
   if(Blm(27)!=0){fprintf(outfile,"B50=%g\n",Blm(27));}
   if(Blm(28)!=0){fprintf(outfile,"B51=%g\n",Blm(28));}
   if(Blm(29)!=0){fprintf(outfile,"B52=%g\n",Blm(29));}
   if(Blm(30)!=0){fprintf(outfile,"B53=%g\n",Blm(30));}
   if(Blm(31)!=0){fprintf(outfile,"B54=%g\n",Blm(31));}
   if(Blm(32)!=0){fprintf(outfile,"B55=%g\n",Blm(32));}

   if(Blm(33)!=0){fprintf(outfile,"B66S=%g\n",Blm(33));}
   if(Blm(34)!=0){fprintf(outfile,"B65S=%g\n",Blm(34));}
   if(Blm(35)!=0){fprintf(outfile,"B64S=%g\n",Blm(35));}
   if(Blm(36)!=0){fprintf(outfile,"B63S=%g\n",Blm(36));}
   if(Blm(37)!=0){fprintf(outfile,"B62S=%g\n",Blm(37));}
   if(Blm(38)!=0){fprintf(outfile,"B61S=%g\n",Blm(38));}
   if(Blm(39)!=0){fprintf(outfile,"B60=%g\n",Blm(39));}
   if(Blm(40)!=0){fprintf(outfile,"B61=%g\n",Blm(40));}
   if(Blm(41)!=0){fprintf(outfile,"B62=%g\n",Blm(41));}
   if(Blm(42)!=0){fprintf(outfile,"B63=%g\n",Blm(42));}
   if(Blm(43)!=0){fprintf(outfile,"B64=%g\n",Blm(43));}
   if(Blm(44)!=0){fprintf(outfile,"B65=%g\n",Blm(44));}
   if(Blm(45)!=0){fprintf(outfile,"B66=%g\n",Blm(45));}

}

void ionpars::savLlm(FILE * outfile)
{fprintf(outfile,"units=meV\n");
   if(Llm(1)!=0){fprintf(outfile,"L22S=%g\n",Llm(1));}
   if(Llm(2)!=0){fprintf(outfile,"L21S=%g\n",Llm(2));}
   if(Llm(3)!=0){fprintf(outfile,"L20=%g\n",Llm(3));}
   if(Llm(4)!=0){fprintf(outfile,"L21=%g\n",Llm(4));}
   if(Llm(5)!=0){fprintf(outfile,"L22=%g\n",Llm(5));}
  
   if(Llm(6)!=0){fprintf(outfile,"L33S=%g\n",Llm(6));}
   if(Llm(7)!=0){fprintf(outfile,"L32S=%g\n",Llm(7));}
   if(Llm(8)!=0){fprintf(outfile,"L31S=%g\n",Llm(8));}
   if(Llm(9)!=0){fprintf(outfile,"L30=%g\n",Llm(9));}
   if(Llm(10)!=0){fprintf(outfile,"L31=%g\n",Llm(10));}
   if(Llm(11)!=0){fprintf(outfile,"L32=%g\n",Llm(11));}
   if(Llm(12)!=0){fprintf(outfile,"L32=%g\n",Llm(12));}

   if(Llm(13)!=0){fprintf(outfile,"L44S=%g\n",Llm(13));}
   if(Llm(14)!=0){fprintf(outfile,"L43S=%g\n",Llm(14));}
   if(Llm(15)!=0){fprintf(outfile,"L42S=%g\n",Llm(15));}
   if(Llm(16)!=0){fprintf(outfile,"L41S=%g\n",Llm(16));}
   if(Llm(17)!=0){fprintf(outfile,"L40=%g\n",Llm(17));}
   if(Llm(18)!=0){fprintf(outfile,"L41=%g\n",Llm(18));}
   if(Llm(19)!=0){fprintf(outfile,"L42=%g\n",Llm(19));}
   if(Llm(20)!=0){fprintf(outfile,"L43=%g\n",Llm(20));}
   if(Llm(21)!=0){fprintf(outfile,"L44=%g\n",Llm(21));}
 
   if(Llm(22)!=0){fprintf(outfile,"L55S=%g\n",Llm(22));}
   if(Llm(23)!=0){fprintf(outfile,"L54S=%g\n",Llm(23));}
   if(Llm(24)!=0){fprintf(outfile,"L53S=%g\n",Llm(24));}
   if(Llm(25)!=0){fprintf(outfile,"L52S=%g\n",Llm(25));}
   if(Llm(26)!=0){fprintf(outfile,"L51S=%g\n",Llm(26));}
   if(Llm(27)!=0){fprintf(outfile,"L50=%g\n",Llm(27));}
   if(Llm(28)!=0){fprintf(outfile,"L51=%g\n",Llm(28));}
   if(Llm(29)!=0){fprintf(outfile,"L52=%g\n",Llm(29));}
   if(Llm(30)!=0){fprintf(outfile,"L53=%g\n",Llm(30));}
   if(Llm(31)!=0){fprintf(outfile,"L54=%g\n",Llm(31));}
   if(Llm(32)!=0){fprintf(outfile,"L55=%g\n",Llm(32));}
 
   if(Llm(33)!=0){fprintf(outfile,"L66S=%g\n",Llm(33));}
   if(Llm(34)!=0){fprintf(outfile,"L65S=%g\n",Llm(34));}
   if(Llm(35)!=0){fprintf(outfile,"L64S=%g\n",Llm(35));}
   if(Llm(36)!=0){fprintf(outfile,"L63S=%g\n",Llm(36));}
   if(Llm(37)!=0){fprintf(outfile,"L62S=%g\n",Llm(37));}
   if(Llm(38)!=0){fprintf(outfile,"L61S=%g\n",Llm(38));}
   if(Llm(39)!=0){fprintf(outfile,"L60=%g\n",Llm(39));}
   if(Llm(40)!=0){fprintf(outfile,"L61=%g\n",Llm(40));}
   if(Llm(41)!=0){fprintf(outfile,"L62=%g\n",Llm(41));}
   if(Llm(42)!=0){fprintf(outfile,"L63=%g\n",Llm(42));}
   if(Llm(43)!=0){fprintf(outfile,"L64=%g\n",Llm(43));}
   if(Llm(44)!=0){fprintf(outfile,"L65=%g\n",Llm(44));}
   if(Llm(45)!=0){fprintf(outfile,"L66=%g\n",Llm(45));}

}
   // evaluate radial wave function
   double ionpars::radial_wavefunction(double rr) // rr given in Angstroems, returns R(r) in units of 1/A^1.5
   {double R=0;int p;double a0=0.5292;
    double r=rr/a0;// r is the distance in units of a0
    for(p=1;p<=9;++p){if(Np(p)!=0){
                                   R+=exp(-Xip(p)*r)*pow(r,Np(p)-1)*Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of radial wave function R(r=%g) failed due to Xi%i<=0 - continuing with R(r=%g)=0\n",r,p,r);return 0;}
                     }            }    
    // now we have R in units of 1/a0^1.5
    R/=sqrt(a0*a0*a0);
    // now we have R in units of 1/A^1.5
    return R;
   }

   //functions to calculate radial matrix elements <r^n> from radial wave function in units of a0=0.5292 A
   double ionpars::rk_from_radial_wavefunction(int k)
   {int p,q, pmax=0;
    Vector coeff(1,9);
    
    for(p=1;p<=9;++p){if(Np(p)!=0){pmax=p;
                                   coeff(p)=Cp(p)*pow(2.0*Xip(p),Np(p)+0.5)/sqrt((double)factorial(2*(int)Np(p)));
                                   if(Xip(p)<=0){fprintf (stderr,"Warning: calculation of <r^%i> failed due to Xi%i<=0 - continuing with <r^%i>=0\n",k,p,k);return 0;}
                     }            }
    if(pmax==0){fprintf (stderr,"Warning: calculation of <r^%i> failed - continuing with <r^%i>=0\n",k);return 0;}
    double rk=0;
    for(p=1;p<=pmax;++p){
    for(q=1;q<=pmax;++q){
                         rk+=coeff(p)*coeff(q)*factorial((int)Np(p)+(int)Np(q)+k)/pow(Xip(p)+Xip(q),Np(p)+Np(q)+k+1);
    }}
   return rk;
   }

   int ionpars::r2_from_radial_wavefunction() {r2=rk_from_radial_wavefunction(2);}
   int ionpars::r4_from_radial_wavefunction() {r4=rk_from_radial_wavefunction(4);}
   int ionpars::r6_from_radial_wavefunction() {r6=rk_from_radial_wavefunction(6);}

void ionpars::save_radial_wavefunction(const char * filename)
   {double r=0.1;
    FILE * fout;
    if (radial_wavefunction(r)==0){fprintf(stderr,"Warning: save_radial_wavefunction not possible\n");return;}
    fout=fopen_errchk(filename,"w");
    fprintf(fout,"# radial wave function for %s\n",iontype);
    fprintf(fout,"# the radial wave function is expanded as \n");
    fprintf(fout,"# R(r)=sum_p C_p R_Np,XIp(r)\n");
    fprintf(fout,"# R_Np,XIp(r)=r^(Np-1).exp(-xi r).(2 XIp)^(Np+0.5)/sqrt(2Np!)\n");
    fprintf(fout,"# radial wave function parameters Np XIp Cp values are\n");
    fprintf(fout,"# tabulated in clementi & roetti Atomic data and \n");
    fprintf(fout,"# nuclear data tables 14 (1974) 177-478\n");
    fprintf(fout,"# the parameters used are: \n");
    int p;    
    for(p=1;p<=9;++p){if(Np(p)!=0){fprintf(fout,"# N%i=%g XI%i=%g C%i=%g\n",p,Np(p),p,Xip(p),p,Cp(p));}}
    fprintf(fout,"# r[A]  vs R(r)[1/A^1.5]\n");
    for(r=0.01;r<=10;r*=1.05){fprintf(fout,"%8.8g  %8.8g\n",r,radial_wavefunction(r));}
    fclose(fout);
   }
//------------------------------------------------------------------------------------------------
// ROUTINE CFIELD mcalc for full crystal field + higher order interactions
//------------------------------------------------------------------------------------------------
Vector & ionpars::cfield(double & T, Vector & gjmbH, double & lnZs, double & U, ComplexMatrix & ests)
{//ABC not used !!!
    /*on input
    T		temperature[K]
    gJmbH	vector of effective field [meV]
    gJ          Lande factor
    ABC         single ion parameter values (A, B, C corresponding to <+|Ja|->,<-|Jb|->,<+|Jc|->/i
  on output    
    J		single ion momentum vector <J> (if T>0 thermal exp value <J>T 
                                                if T<0 the program asks for w_n and calculates
						       exp value <J>=sum_n w_n <n|J|n>
						       
    Z		single ion partition function
    U		single ion magnetic energy
*/
// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error internal module cfield: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

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
static Vector JJ(1,gjmbH.Hi());
   // setup hamiltonian
   int dj,i,j,k,l;
   double hkl,mukl;
   dj=Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
//   Matrix Tam(1,dj,1,dj); // transformed Hamiltonian
   ComplexMatrix z(1,dj,1,dj);
   ComplexMatrix za(1,dj,1,dj);
   ComplexMatrix zb(1,dj,1,dj);
   ComplexMatrix zc(1,dj,1,dj);
   ComplexMatrix zolm(1,dj,1,dj);    

   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;

   for(j=4;j<=JJ.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);}

/*   int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",(*Olm[j])(i1,j1));
    printf ("\n");
    }*/
   // use old eigenstates ests to transform matrix to nearly diagonal form ... however we deleted this because it needs more time to transform than to solve the eigenvalue problem
/*   Tam=0;
   for(i=1;i<=dj;++i){for(j=1;j<=dj;++j){
   if(i<j){for(k=1;k<=dj;++k){for(l=1;l<=dj;++l){
           if(k<l){hkl=Ham(l,k);mukl=-Ham(k,l);}else{hkl=Ham(k,l);if(k==l){mukl=0;}else{mukl=Ham(l,k);}}           
           Tam(i,j)-=-imag(ests(k,i))*hkl*real(ests(l,j))+imag(ests(k,i))*mukl*imag(ests(l,j))+real(ests(k,i))*mukl*real(ests(l,j))+real(ests(k,i))*hkl*imag(ests(l,j));    
           }}
          }
   else   {for(k=1;k<=dj;++k){for(l=1;l<=dj;++l){
           if(k<l){hkl=Ham(l,k);mukl=-Ham(k,l);}else{hkl=Ham(k,l);if(k==l){mukl=0;}else{mukl=Ham(l,k);}}            
           Tam(i,j)+=real(ests(k,i))*hkl*real(ests(l,j))-real(ests(k,i))*mukl*imag(ests(l,j))+imag(ests(k,i))*mukl*real(ests(l,j))+imag(ests(k,i))*hkl*imag(ests(l,j));
           }}
          }
   }}         
  */  
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=0;int maxiter=1000000;
   if (T<0) sort=1;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);

   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);
     double x,y;
     x=Min(En);
     double Zs;

     if (T>0)
     { for (i=1;i<=dj;++i)
       {if ((y=(En(i)-x)/K_B/T)<600) wn[i]=exp(-y); 
        else wn[i]=0.0;
       }
       Zs=Sum(wn);wn/=Zs;
 
       lnZs=log(Zs)-x/K_B/T;
     } 
     else
     { printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Energy     Excitation Energy\n");
     for (i=1;i<=dj;++i) printf ("%i    %4.4g meV   %4.4g meV\n",i,En(i),En(i)-x);
     char instr[MAXNOFCHARINLINE];
     for (i=1;i<=dj;++i)
      {printf("eigenstate %i: %4.4g meV %4.4g meV  - please enter probability w(%i):",i,En(i),En(i)-x,i);
       fgets(instr, MAXNOFCHARINLINE, stdin);
 
       wn(i)=strtod(instr,NULL);
      }
       Zs=Sum(wn);wn/=Zs;
 
       lnZs=log(Zs);
                         printf ("\n\nNumber   Energy     Excitation Energy   Probability\n");
     for (i=1;i<=dj;++i) printf ("%i    %4.4g meV   %4.4g meV %4.4g  \n",i,En(i),En(i)-x,wn(i));
     }

   // calculate U
     U=En*wn;
   // calculate Ja,Jb,Jc
     z=ComplexMatrix(zr,zi);
//     z=ests(1,dj,1,dj)*z; // transform to original eigenstates ... however we deleted this because it needs more time to transform than to solve the eigenvalue problem
//     ests(1,dj,1,dj)=z;
//     for (i=1;i<=dj;++i) {ests(0,i)=complex <double> (En(i),wn(i));}
//     myPrintComplexMat(stdout,ests);     
//     myPrintComplexMat(stdout,z);     
     
     za=Jaa*z;
     zb=Jbb*z;
     zc=Jcc*z;

    
     JJ=0;
//    ComplexVector ddd;
    for (i=1;i<=dj;++i)
    {
     JJ[1]+=wn(i)*real(z.Column(i)*za.Column(i));
     JJ[2]+=wn(i)*real(z.Column(i)*zb.Column(i));
     JJ[3]+=wn(i)*real(z.Column(i)*zc.Column(i));
    }
     
   for(j=4;j<=JJ.Hi();++j)
   {
    zolm=(*OOlm[j-3])*z;
    for (i=1;i<=dj;++i) JJ[j]+=wn(i)*real(z.Column(i)*zolm.Column(i));
   };
  

return JJ;
}
/**************************************************************************/
ComplexMatrix & ionpars::cfeigenstates(Vector & gjmbH, double & T)
{   /*on input
    gJmbH	vector of effective field [meV]
      on output
    Matrix containing the eigenvalues and eigenfunctions of the crystalfield problem
    eigenvalues ares stored as real part of row zero
    boltzmann population numbers are stored as imaginary part of row zero
*/

// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error internal module cfield: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

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
   // setup hamiltonian
   int dj,i,j;
//   complex <double> imag(0,1);
   dj=Hcf.Rhi();
static ComplexMatrix eigenstates(0,dj,1,dj);
   Matrix Ham(1,dj,1,dj);
   ComplexMatrix z(1,dj,1,dj);
   ComplexMatrix za(1,dj,1,dj);
   ComplexMatrix zb(1,dj,1,dj);
   ComplexMatrix zc(1,dj,1,dj);
   ComplexMatrix zolm(1,dj,1,dj);    

 
   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;

/* myPrintComplexMat(stdout,Jaa);
 myPrintComplexMat(stdout,Jbb);
 myPrintComplexMat(stdout,Jcc);*/

   for(j=4;j<=gjmbH.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);}

/*   int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",(*Olm[j])(i1,j1));
    printf ("\n");
    }*/
      
    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=1;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);

   for(i=1;i<=dj;++i){//eigenstates(0,i)=complex <double> (En(i),0);
     for(j=1;j<=dj;++j){eigenstates(i,j)=complex <double> (zr(i,j),zi(i,j));
   }}
    //calculate partition sum
     double zz=0;double KBT,E0;KBT=T*K_B;E0=En(1);
      for(j=1;j<=dj;++j){zz+=exp(-((En(j)-E0)/KBT));}
        // put boltzmann population into row 0 of eigenstates...
        for(j=1;j<=dj;++j)
         {eigenstates(0,j)=complex<double>(En(j),exp(-(En(j)-E0)/KBT)/zz);}
   
 return eigenstates;
}

/**************************************************************************/
// for mcdisp this routine is needed
int ionpars::cfielddm(int & tn,double & T,Vector & gjmbH,ComplexMatrix & mat,float & delta,ComplexMatrix & ests)
{  /*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta-+	energy of transition [meV]
    mat(i,j)	<-|Ji|+><+|Jj|-> (n+-n-),  n+,n-
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/


// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error loadable module cfield.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

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

static Vector JJ(1,gjmbH.Hi());
double lnz,u;
JJ=cfield(T,gjmbH,lnz,u,ests);  //expectation values <J>
  int pr;
  pr=1;
  if (tn<0) {pr=0;tn*=-1;}

   // setup hamiltonian
   int dj,j;
   dj=Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
    
   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;
 for(j=4;j<=gjmbH.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);}

/*   int i1,j1; //printout matrix
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    real(((*OOlm[5])-Jcc*Jcc+Jaa*Jaa)(i1,j1)));}
//    real((Jcc*Jaa+Jaa*Jcc)(i1,j1)));}
//    real((*OOlm[1])(i1,j1)));}
    printf ("\n");
    }
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    imag(((*OOlm[5])-Jcc*Jcc+Jaa*Jaa)(i1,j1)));}
//   imag((Jcc*Jaa+Jaa*Jcc)(i1,j1)));}
//   imag((*OOlm[1])(i1,j1)));}
    printf ("\n");
    }
exit(0);      
*/    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=1;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
   
   
   
   
   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);double Zs;
     double x,y;int i,k,l,m;
     x=Min(En);
     for (i=1;i<=dj;++i)
     {if ((y=(En(i)-x)/K_B/T)<700) wn[i]=exp(-y); 
      else wn[i]=0.0;
//      printf("%4.4g\n",En(i));
      }
     Zs=Sum(wn);wn/=Zs;  
     Zs*=exp(-x/K_B/T);


   // calculate Ja,Jb,Jc
     ComplexMatrix z(1,dj,1,dj);
     ComplexMatrix * zp[gjmbH.Hi()+1];
     for(l=1;l<=gjmbH.Hi();++l)
      {zp[l]= new ComplexMatrix(1,dj,1,dj);}
     z=ComplexMatrix(zr,zi);
     
     (*zp[1])=Jaa*z;
     (*zp[2])=Jbb*z;
     (*zp[3])=Jcc*z;

     
 for(j=4;j<=gjmbH.Hi();++j)
    {(*zp[j])=(*OOlm[j-3])*z;}
     
// calculate mat and delta for transition number tn
// 1. get i and j from tn
k=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=En(j)-En(i);

if (delta<-0.000001){fprintf(stderr,"ERROR module cfield.so - dmcalc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

// 3. set mat
for(l=1;l<=gjmbH.Hi();++l)for(m=1;m<=gjmbH.Hi();++m)
{if(i==j){//take into account thermal expectation values <Jl>
          mat(l,m)=((z.Column(i)*(*zp[l]).Column(j))-JJ(l))*((z.Column(j)*(*zp[m]).Column(i))-JJ(m));}
 else    {mat(l,m)=(z.Column(i)*(*zp[l]).Column(j))*(z.Column(j)*(*zp[m]).Column(i));}}



if (delta>SMALL)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Ja|%i>|^2=%4.4g |<%i|Jb|%i>|^2=%4.4g |<%i|Jc|%i>|^2=%4.4g",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
      printf(" n%i-n%i=%4.4g\n",i,j,wn(i)-wn(j));}
    mat*=(wn(i)-wn(j)); // occupation factor    
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Ja-<Ja>|%i>|^2=%4.4g |<%i|Jb-<Jb>|%i>|^2=%4.4g |<%i|Jc-<Jc>|%i>|^2=%4.4g",i,j,real(mat(1,1)),i,j,real(mat(2,2)),i,j,real(mat(3,3)));
      printf(" n%i=%4.4g\n",i,wn(i));}
    mat*=(wn(i)/K_B/T);
   }

//clean up memory
     for(l=1;l<=gjmbH.Hi();++l)
      {delete zp[l];}
     
// return number of all transitions     
 return (int)((J+1)*(2*J+1)); 
}

int ionpars::cfielddn(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexMatrix & nat)
{/*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    est		matrix with eigenstates, eigenvalues [meV], population numbers
    gjmbH	vector of effective field [meV]
  on output    
    int   	total number of transitions
    N(i,j)	<-|Q|+><+|Q|-> (n+-n-),  n+,n-
     // note that  <M(Q)>=-2x<Q>_TH in units of mb
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/
  int pr;pr=1;if (tn<0) {pr=0;tn*=-1;}
  int i,j,k,l,m;
  int dj=(int)(2*J+1);
  double delta;
// calculate nat for transition number tn
// 1. get i and j from tn (as in dmcalc
k=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=real(est(0,j))-real(est(0,i));

if (delta<-0.000001){fprintf(stderr,"ERROR module cfield.so - dncalc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

	 ComplexMatrix * MQMi[4];
         MQMi[1]=new ComplexMatrix(1,dj,1,dj);
         MQMi[2]=new ComplexMatrix(1,dj,1,dj);
         MQMi[3]=new ComplexMatrix(1,dj,1,dj);
         MQM((*MQMi[3]),(*MQMi[1]),(*MQMi[2]),th,ph,J0,J2,J4,J6,Zc);
        //      x           y         z
        //      c           a         b

// 3. set nat
         int K,M,Md;
         ComplexVector Malpha(1,3);Malpha=0;
          for(K=1;K<=3;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){
             Malpha(K)+=conj(est(M,i))*(*MQMi[K])(M,Md)*est(Md,j); 
            }}} 
if(i==j){//take into account thermal expectation values <Jl>
         ComplexVector mm(1,3); mm=0;
         for(K=1;K<=dj;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){          
           mm(1)+=imag(est(0,K))*conj(est(M,K))*(*MQMi[1])(M,Md)*est(Md,K); 
           mm(2)+=imag(est(0,K))*conj(est(M,K))*(*MQMi[2])(M,Md)*est(Md,K); 
           mm(3)+=imag(est(0,K))*conj(est(M,K))*(*MQMi[3])(M,Md)*est(Md,K); 
         }}} // --> mm(1,..3)  thermal expextation values of M
          Malpha-=mm;// subtract thermal expectation values
         
         }
         delete MQMi[1],MQMi[2],MQMi[3];


       // set matrix <i|Ml|j><j|Mm|i>
       nat=0;
          for(l=1;l<=3;++l)for(m=1;m<=3;++m)
          {nat(l,m)=Malpha(l)*conj(Malpha(m));}

// multiply by occupation number difference ...

if (delta>SMALL)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Qa|%i>|^2=%4.4g |<%i|Qb|%i>|^2=%4.4g |<%i|Qc|%i>|^2=%4.4g",i,j,real(nat(1,1)),i,j,real(nat(2,2)),i,j,real(nat(3,3)));
      printf(" n%i-n%i=%4.4g\n",i,j,imag(est(0,i))-imag(est(0,j)));}
    nat*=(imag(est(0,i))-imag(est(0,j))); // occupation factor    
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Qa-<Qa>|%i>|^2=%4.4g |<%i|Qb-<Qb>|%i>|^2=%4.4g |<%i|Qc-<Qc>|%i>|^2=%4.4g",i,j,real(nat(1,1)),i,j,real(nat(2,2)),i,j,real(nat(3,3)));
      printf(" n%i=%4.4g\n",i,imag(est(0,i)));}
    nat*=(imag(est(0,i))/K_B/T);
   }


// return number of all transitions     
 return (int)((J+1)*(2*J+1)); 

}

//**********************************************************************/
// routine to calculate the scattering operator to go beyond dip approx
// *********************************************************************

// just another routine to calculakte Z(K)
double Z(int K, float J0, float J2, float J4, float J6, Vector Zc)
{// calculate Z(K)
 if (K==1) return Zc(1)*J0+Zc(2)*J2;
 if (K==3) return Zc(3)*J2+Zc(4)*J4;
 if (K==5) return Zc(5)*J4+Zc(6)*J6;
 if (K==7) return Zc(7)*J6;
 
 return 0;
}

void ionpars::MQM(ComplexMatrix & MQXM,ComplexMatrix & MQYM,ComplexMatrix & MQZM, double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc)
{complex <double> im(0,1);
         // .... calculate scattering operator ...(formula 11.141-143 in lovesey)
	    int K,Qd,M,Md;double MJ,MJd,PKdQd,thj,factor;
	    complex <double>bracketx,brackety,bracketz;
	    MQXM=0;MQYM=0;MQZM=0;
	    for(K=1;K<=7;K+=2){
			     factor=sqrt(4.0*PI)*Z(K,J0,J2,J4,J6,Zc)/K;
                               if (factor!=0){
					     thj=threej((float)K,J,J,0,J,-J);
					     for(Qd=-K;Qd<=K;Qd+=1){
                                                 bracketx=0;brackety=0;
					       if(K-1>=Qd+1&&K-1>=-Qd-1)
					       {bracketx+=SphericalHarmonicY (K-1,Qd+1,th,ph)*sqrt((double)(K-Qd)*(K-Qd-1));
					        brackety+=SphericalHarmonicY (K-1,Qd+1,th,ph)*sqrt((double)(K-Qd)*(K-Qd-1));
					       }
					       if(K-1>=Qd-1&&K-1>=-Qd+1)
                                                 {bracketx-=SphericalHarmonicY (K-1,Qd-1,th,ph)*sqrt((double)(K+Qd)*(K+Qd-1));
						        brackety+=SphericalHarmonicY (K-1,Qd-1,th,ph)*sqrt((double)(K+Qd)*(K+Qd-1));
					       }
					       if(K-1>=Qd&&K-1>=-Qd)
                                                 {bracketz=SphericalHarmonicY (K-1,Qd,th,ph)*sqrt((double)(K-Qd)*(K+Qd));
                                                 }else {bracketz=0;}

//.(1)..USE !		     ThreeJSymbolM	(J1,J2,J3,M1,&M2min,&M2max,*thrcof,ndim,errflag);
                                                 double thrj[30];int ndim=30; double MJdmin,MJdmax; int errflag;
                                                                      
                                                 ThreeJSymbolM ((float)K,J,J,-(float)Qd,MJdmin,MJdmax,thrj,ndim,errflag);
                                                 if (errflag!=0){fprintf(stderr,"ERROR mcdiff: threejsymbol error %i\n",errflag);exit(EXIT_FAILURE);}           
                                                 for (Md=int(MJdmin+1+J);Md<=int(MJdmax+1+J);++Md){
						                 MJd=(float)Md-1-J;
								 MJ=-Qd+MJd;M=int(MJ+1+J);
							         PKdQd=thrj[Md-int(MJdmin+1+J)]/thj; 
							         PKdQd*=odd(int(J-MJd)) ? -1 : 1;
							         MQXM(M,Md)+=0.5*factor*PKdQd*bracketx;
							         MQYM(M,Md)+=-im*brackety*0.5*factor*PKdQd;
							         MQZM(M,Md)+=factor*PKdQd*bracketz;
                                                                 }

/*
//.(2)..                     3jsymb=threej(J1,J2,J3,M1,M2,M3) 
//							       for(M=1;M<=dj;++M){MJ=(float)M-1-J;  
//							        for(Md=1;Md<=dj;++Md){MJd=(float)Md-1-J; 
//							         // according to 11.140 lovesey book        
//							         PKdQd=threej((float)K,J,J,-(float)Qd,MJd,-MJ)/thj; 
//							         PKdQd*=odd(int(J-MJd)) ? -1 : 1;
//							         MQXM(M,Md)+=im*0.5*factor*PKdQd*bracketx;
//							         MQYM(M,Md)+=-0.5*factor*PKdQd*brackety;
//							         MQZM(M,Md)+=factor*PKdQd*bracketz;
//							        }
//							       }
*/
                                             }
							      }
                                                             }

}
// calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
// according to stored eigenstate matrix est
// calculates the scattering operator given the polar angles th, ph (with respect to the CEF coordinate 
// system xyz and the <jl(qr)> and the eigenstate matrix with eigenstates and thermal population numbers
ComplexVector & ionpars::MQ(double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc,ComplexMatrix & est)
    {  
       int dj=(int)(2*J+1);
       
	 ComplexMatrix MQXM(1,dj,1,dj),MQYM(1,dj,1,dj),MQZM(1,dj,1,dj);
         MQM(MQXM,MQYM,MQZM,th,ph,J0,J2,J4,J6,Zc);
							     // ... calculate thermal expectation values
							     // using the eigenstates and T 
							     // mom(1) = ....
          						     //
       static ComplexVector mm(1,3); mm=0;
       int K,M,Md;
       for(K=1;K<=dj;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){
         mm(1)+=imag(est(0,K))*conj(est(M,K))*MQXM(M,Md)*est(Md,K); 
         mm(2)+=imag(est(0,K))*conj(est(M,K))*MQYM(M,Md)*est(Md,K); 
         mm(3)+=imag(est(0,K))*conj(est(M,K))*MQZM(M,Md)*est(Md,K); 
       }}}
// myPrintComplexMatrix(stdout,MQXM);
// myPrintComplexMatrix(stdout,MQYM);
// myPrintComplexMatrix(stdout,MQZM);
// 					       myPrintComplexMatrix(stdout,est);}
					       
                     mm*=2; // this is now <M(Q)>=-2x<Q>_TH in units of mb
// myPrintComplexVector(stdout,mm);//equivalent to moment ...
    return mm;
    }
