 // *************************************************************************
 // ************************ spincf *************************************
 // *************************************************************************
// methods for class spincf 
#include <cerrno>
#include <cstdio>
#include <cmath>
#include <martin.h>
#include <vector.h>
#include <complex>
#include "spincf.hpp"
#include "chargedensity.hpp"
#include "jjjpar.hpp"

#define MAXNOFSPINS  200
#define SMALL 0.03
#define PI  3.1415
int spincf::spequal(Vector a,Vector b)
{// in this function we look if two spins are equal or not (used in order to compare or reduce
 // spinconfigurations
 //if(a*b/Norm(a)/Norm(b)<COSINE){return false;} // this is by the inner product

 //this is by looking if any of the components of the 2 spins have differen signs
 int i;
 for (i=a.Lo();i<=a.Hi();++i)
 {if (fabs(a(i))>SMALL||fabs(b(i))>SMALL)
  {if(copysign(1.0,a(i))!=copysign(1.0,b(i))){return false;}
  }
 }
return true;
}

// returns momentum of ion [i=na,j=nb,k=nc] 
Vector & spincf::m(int na, int nb, int nc)
{ return mom[in(na,nb,nc)];  // return <J(ijk)>
}

// the same but for spin number "i"
Vector & spincf::mi(int i)
{ return mom[i];
}

Vector spincf::moment(int i,int j,int k,int l) // returns moment of atom l (1,nofcomponents)
{Vector xyz(1,nofcomponents);
 int m;
 for(m=1;m<=nofcomponents;++m){xyz(m)=mom[in(i,j,k)](nofcomponents*(l-1)+m);}
 return xyz;
}
Vector spincf::magmom(int i,int j,int k,int l,Vector & gJ) // returns magnetic moment 
{       Vector xyz(1,3);int m,maxm;
       xyz=0;
       if (gJ(l)==0)  //load magnetic moment into vector xyz
       {             // intermediate coupling  <M>=2*<S>+<L>
        if(nofcomponents>6){maxm=6;}else{maxm=nofcomponents;}
        for(m=1;m<=maxm;++m){if(m==2||m==4||m==6){xyz((m+1)/2)+=mom[in(i,j,k)](nofcomponents*(l-1)+m);}
                             else                {xyz((m+1)/2)+=2.0*mom[in(i,j,k)](nofcomponents*(l-1)+m);}
                            }
       }
       else // LS coupling
       {if(nofcomponents>3){maxm=3;}else{maxm=nofcomponents;}
        for(m=1;m<=maxm;++m){xyz(m)=mom[in(i,j,k)](nofcomponents*(l-1)+m)*gJ(l);}
       }
       return xyz;
}

// (Slow) Fourier Transform of momentumconfiguration <J> ... i.e. for magnetic moments multiply by lande factor !!!
void  spincf::FT(ComplexVector * mq)
{int i,j,k,l;
 int qh,qk,ql;
    complex<double> piq(0,2*3.1415926535),g;
  for(qh=0;qh<na();++qh){for(qk=0;qk<nb();++qk){for(ql=0;ql<nc();++ql)
     {mq[in(qh,qk,ql)]=0;
     for(i=0;i<na();++i){for(j=0;j<nb();++j){for(k=0;k<nc();++k)
      { g=exp(piq*((double)qh*i/na()+(double)qk*j/nb()+(double)ql*k/nc()));
	for (l=1;l<=nofcomponents*nofatoms;++l) {mq[in(qh,qk,ql)](l)+=g*m(i+1,j+1,k+1)(l);}
      }}}
//      fprintf(stdout,"%g %g %g %g %g \n",qh,qk,ql,imag(mq[in(qh,qk,ql)](1)),real(mq[in(qh,qk,ql)](1)));
     }}}
 for (qk=1;qk<=nb();++qk)
 {for (ql=1;ql<=nc();++ql) mq[in(na(),qk,ql)]=mq[in(0,qk,ql)];}
 for (qh=1;qh<=na();++qh)
 {for (ql=1;ql<=nc();++ql) mq[in(qh,nb(),ql)]=mq[in(qh,0,ql)];}
 for (qh=1;qh<=na();++qh)
 {for (qk=1;qk<=nb();++qk) mq[in(qh,qk,nc())]=mq[in(qh,qk,0)];}
 mq[in(na(),nb(),nc())]=mq[in(0,0,0)];
 return;
}

// get index ijk=iv(1-3)  of spinconfiguration number in 
int * spincf::ijk(int in)
{div_t result; result=div(in,mxb*mxc); 
 iv[1]= result.quot;
 result=div(result.rem,mxc);
 iv[2]= result.quot;
 iv[3]= result.rem;
 return iv;}

// the inverse: get number of spin from indizes i,j,k
int spincf::in(int i, int j, int k)
{return ((i*mxb+j)*mxc+k);}

// this subtracts n2 if n1>n2
int spincf::mod(int n1,int n2)
{if (n1>n2) return n1-n2;
 else return n1;
}

// return number of spins
int spincf::n()
{return (nofa*nofb*nofc);
}
int spincf::na()
{return nofa;
}
int spincf::nb()
{return nofb;
}
int spincf::nc()
{return nofc;
}

Vector spincf::nettomagmom(Vector & gJ) // returns nettomagnetic moment [mu_b]
{
 int i,j,k,l,m,n; 
 Vector ret(1,3);
 ret=0;
 for (i=1;i<=nofa;++i)
 { for (j=1;j<=nofb;++j)
   {for (k=1;k<=nofc;++k)
    {for (l=1;l<=nofatoms;++l)
     {ret+=magmom(i,j,k,l,gJ);
     }    
    }
   }
  }     
 return ret;
}

Vector spincf::totalJ() // returns total moment <J>
{Vector ret(1,nofcomponents);
 ret=0;
 int i,j,k,l,m; 
 for (i=1;i<=nofa;++i)
 { for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {for (m=1;m<=nofcomponents;++m)
      {for (l=1;l<=nofatoms;++l)
       { ret(m)+=mom[in(i,j,k)](nofcomponents*(l-1)+m);
       }      
      }    
     }
    }
 }     
 ret/=n()*nofatoms;
 return ret;
}


// invert all spins (AND higher order moments)
void spincf::invert()
{
 int i,j,k; 
 for (i=1;i<=nofa;++i)
 { for (j=1;j<=nofb;++j)
   {for (k=1;k<=nofc;++k)
    {mom[in(i,j,k)]= -mom[in(i,j,k)];}
   }
 }  
}

// reduce spinconfiguration if it is periodic
void spincf::reduce()
{int perprob,nn,pp,j,k;
 div_t result;

//reduce nofa
pp=(int)integer(nofa/2);
for (perprob=1;perprob<=pp;++perprob)
  {result=div(nofa,perprob);
   if (result.rem==0)
    {
     for (j=1;j<=nofb;++j)
      {for (k=1;k<=nofc;++k)
       {for (nn=1;nn<=nofa-perprob;++nn)
        {if (!spequal(mom[in(nn,j,k)],mom[in(nn+perprob,j,k)])){goto nexta;}
	 if(nn==nofa-perprob&&j==nofb&&k==nofc)
        	{this->nofa=perprob;
	// here we could reduce the memory of this confuration
	// however for the moment we leave it for it does not disturb
	         goto reduceda;
	        }  
	}
       }
      } 	
    }
nexta:;
  }
reduceda:


//reduce nofb
pp=(int)integer(nofb/2);
for (perprob=1;perprob<=pp;++perprob)
  {result=div(nofb,perprob);
   if (result.rem==0)
    {
     for (j=1;j<=nofa;++j)
      {for (k=1;k<=nofc;++k)
       {for (nn=1;nn<=nofb-perprob;++nn)
        {if(!spequal(mom[in(j,nn,k)],mom[in(j,nn+perprob,k)])){goto nextb;}
	 if(nn==nofb-perprob&&j==nofa&&k==nofc)
        	{this->nofb=perprob;
	// here we could reduce the memory of this confuration
	// however for the moment we leave it for it does not disturb
	         goto reducedb;
	        }  
	}
       } 	
      }
    }
nextb:;
  }    
reducedb:
//reduce nofc
pp=(int)integer(nofc/2);
for (perprob=1;perprob<=pp;++perprob)
  {result=div(nofc,perprob);
   if (result.rem==0)
    {
     for (j=1;j<=nofa;++j)
      {for (k=1;k<=nofb;++k)
       {for (nn=1;nn<=nofc-perprob;++nn)
        {if(!spequal(mom[in(j,k,nn)],mom[in(j,k,nn+perprob)])){goto nextc;}
	 if(nn==nofc-perprob&&j==nofa&&k==nofb)
        	{this->nofc=perprob;
	// here we could reduce the memory of this confuration
	// however for the moment we leave it for it does not disturb
	         return;
	        }  
	}
       }
      } 	
    }
nextc:;
  }

return;
}


  
// create spinconfig from qvector, momentq0
// n1 n2 n3 .. periodicity of supercell
// nettom .... saturation moment (positive)
// qvector ... wave vector in units of reciprocal lattice
// momentq0 .. ferromagnetic component (between 0 and 1)
// phi ....... phase (for each component)
void spincf::spinfromq (int n1,int n2, int n3,Vector & qvector,Vector & nettom, 
       Vector & momentq0, Vector & phi)
{ int rra,rrb,rrc,qv;
  qv=1; // evtl koennte hier noch ein vorzeichen uebergeben werden (sodass nettom positiv)
  Vector rr(1,3);
  delete []mom;
  nofa=n1; nofb=n2; nofc=n3;
  mxa=n1+1; mxb=n2+1; mxc=n3+1;
  int l;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(l=0;l<=mxa*mxb*mxc;++l){mom[l]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 

// fill spinconfiguration with values
   for (rra=1;rra<=nofa;++rra)
     { for (rrb=1;rrb<=nofb;++rrb)
        {for (rrc=1;rrc<=nofc;++rrc)
             	{rr(1)=rra-1;rr(2)=rrb-1;rr(3)=rrc-1;
		 for(l=1;l<=nofcomponents*nofatoms;++l)
	          {
//		  mom[in(rra,rrb,rrc)](l)=copysign(nettom(l),qv*(0.01+
//		  copysign(1.0,momentq0(l)+cos(phi(l)+qvector*2.0*3.141529*rr)))
//				                    )                    ;	
		  mom[in(rra,rrb,rrc)](l)=nettom(l)*qv*(0.01+
		  momentq0(l)+cos(phi(l)+qvector*2.0*3.141529*rr))
				                                        ;	
                  }
	         }
	}
     }
        
}




// load spinconfiguration from file
int spincf::load(FILE * fin_coq)	
{ int i,j,k,l,nn1,nn2;  
  float na[MAXNOFSPINS];
  float nb[MAXNOFSPINS];
  float nc[MAXNOFSPINS];
  na[0]=MAXNOFSPINS;
  nb[0]=MAXNOFSPINS;
  nc[0]=MAXNOFSPINS;
  long int pos;
  
  pos=ftell(fin_coq);if (pos==-1) return 0;

nn1=0;
j=0;
while ((i=inputline (fin_coq, na))!=0)  
   {for (l=1;l<=nofcomponents*nofatoms-2;++l){inputline (fin_coq, nb);}
    nn1=i;i=inputline (fin_coq, nc);++j;}
// j (nofspins in b direction)determined
if (nn1==0) return 0; // no block to read any more 
nn2=j;

//how much is nofc ??
for (k=0;j==nn2;++k)
  {for (j=0;(i=inputline(fin_coq,na))==nn1;++j)
   {for (l=1;l<=nofcomponents*nofatoms-2;++l){inputline (fin_coq, nb);}
    i=inputline (fin_coq, nc);}
  }

nofa=nn1;
nofb=nn2;    
nofc=k;

//  printf("n1=%i n2=%i n3=%i\n",nofa,nofb,nofc);
  
  j=fseek(fin_coq,pos,SEEK_SET); if (j!=0) return 0;

  delete []mom;
  mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(l=0;l<=mxa*mxb*mxc;++l){mom[l]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL)
    {fprintf (stderr, "Out of memory\n");
     exit (EXIT_FAILURE);} 

 for (k=1;k<=nofc;++k)
 {for (j=1;j<=nofb;++j)
  {for (l=1;l<=nofcomponents*nofatoms;++l)
   {inputline (fin_coq, na);
    for (i=1;i<=nofa;++i)
       {mom[in(i,j,k)](l)=na[i];}
    }
   }
  if(inputline(fin_coq, na))
  {fprintf(stderr,"load spincf error - no  empty line between [001] planes");exit(EXIT_FAILURE);}  
  }       
  if(inputline(fin_coq, na))
  {fprintf(stderr,"load spincf error - no  2 empty lines after spinconfiguration");exit(EXIT_FAILURE);}  
         
 return 1;
}

//-----------------------------------------------------------------------
// create eps output of spinconfiguration
void spincf::eps(FILE * fout) //print spinconfiguration to stream
{eps(fout,"no title");}

void spincf::eps(FILE * fout,const char * text ) //print spinconfiguration to stream
{//viewport [-1,1,-1,1] ... distribute spins on that
  int i,j,k,l,m;
  float x0,y0,compoffset,atomoffset;
  Vector a(1,2);
  Vector b(1,2);
  float scale,d;
 
  scale=0;
  for (i=1;i<=nofa;++i)
    for (j=1;j<=nofb;++j)
     for (k=1;k<=nofc;++k)
      for (l=1;l<=nofcomponents*nofatoms;++l)
      {if ((d=fabs(mom[in(i,j,k)](l)))>scale)scale=d;
      } 
  scale=0.2/(scale+0.01)/(double)nofc;

  fprintf(fout,"%s!PS-Adobe-2.0 EPSF-2.0\n","%");
  fprintf(fout,"%sBoundingBox:0 0 549 %i\n","%%",150*nofcomponents*nofatoms);
  fprintf(fout,"%sTitle: %s\n","%%",text);
  fprintf(fout,"%sEndComments\n","%%");
  fprintf(fout,"/mm {72 mul 25.4 div} bind def\n");
  fprintf(fout,"/mx {1.2 add 70 mul mm} bind def\n");
  fprintf(fout,"/my {1.2 add 70 mul mm} bind def\n");
   fprintf(fout,"/Helvetica findfont\n15 scalefont setfont\n");
   fprintf(fout,"-0.89 mx 0.38 my moveto \n (r1) show \n");
   fprintf(fout,"-0.93 mx 0.46 my moveto \n (r2) show \n");
   fprintf(fout,"-0.99 mx 0.52 my moveto \n (r3) show \n");
   fprintf(fout,"2 mm 2 mm moveto \n (%s) show \n",text);
   
   a(1)=-0.98;a(2)=0.4;
   b(1)=a(1);b(2)=0.5;epsarrow(fout,a,b);
   b(1)=-0.94;b(2)=0.45;epsarrow(fout,a,b);
   b(1)=-0.92;b(2)=0.4;epsarrow(fout,a,b);

//  char ss='a';
  for(l=1;l<=nofatoms;++l)
  {atomoffset=(l-0)*2*nofcomponents/3;
   for(m=1;m<=nofcomponents;++m)
    {compoffset=atomoffset-1.4-0.6*(m-1); 
     fprintf(fout,"/Helvetica findfont\n15 scalefont setfont\n");
     fprintf(fout,"-0.98 mx %g my moveto \n (J%c%i) show \n",compoffset,'a'-1+m,l);

     for (i=1;i<=nofa;++i) 
      for (j=1;j<=nofb;++j)
       for (k=1;k<=nofc;++k)
       { x0=(2.0*(i-1)/nofa-1+1.7*j/nofb/nofa)*1.1+0.3;
         y0=(2.0*(k-1)/nofc-1+1.2*j/nofb/nofc)*0.2+compoffset;
         a(1)=x0;a(2)=y0;
         b(1)=x0;b(2)=mom[in(i,j,k)](nofcomponents*(l-1)+m)*scale+y0;
         epsarrow(fout,a,b);
       }
    }
  }
fprintf(fout,"showpage\n");
}
void spincf::epsarrow(FILE * fout,Vector x,Vector y)
 {  double l=0.15*Norm(y-x);
    Vector y1(1,2),y2(1,2),unn(1,2),upn(1,2);
    if (y==x){y(2)=x(2)+0.0001;}
    upn=(y-x)/Norm(y-x);unn(1)=upn(2);unn(2)=-upn(1);
    y1=y-l*upn-0.3*l*unn;
    y2=y-l*upn+0.3*l*unn;

  fprintf(fout,"0.7 setlinewidth\n");
  fprintf(fout,"%g mx  %g my moveto\n",x(1),x(2));
  fprintf(fout,"%g mx  %g my lineto\n",y(1),y(2));
  fprintf(fout,"%g mx  %g my lineto\n",y1(1),y1(2));
  fprintf(fout,"%g mx  %g my lineto\n",y(1),y(2));
  fprintf(fout,"%g mx  %g my lineto\n",y2(1),y2(2));
  fprintf(fout,"stroke\n");

  }


void spincf::calc_prim_mag_unitcell(Matrix & p,Vector & abc, Matrix & r)
{ int i,j;
  Vector dd(1,3),nofabc(1,3);
  nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  for (i=1;i<=3;++i)
  {for(j=1;j<=3;++j) {dd(j)=nofabc(j)*r(i,j)*abc(i);p(i,j)=dd(j);}
  }
 // pa=p.Column(1);  //primitive magnetic unit cell
 // pb=p.Column(2);
 // pc=p.Column(3);
}

void spincf::calc_minmax(Vector & min,Vector & max,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc)
{calc_minmax_scale(min,max,ijkmin,ijkmax,p,abc,1.0,1.0,1.0);
}

void spincf::calc_minmax_scale(Vector & min,Vector & max,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc,double scale_view_1,double scale_view_2,double scale_view_3)
{// determine max(1,2,3) min(1,2,3) (vector in Angstroem describing a quader) for viewing magnetic unit cell
  Vector ddd(1,8),dd0(1,3),dd(1,3);
  int i;
 double t;
  for (i=1;i<=3;++i)
  {ddd(1)=p.Column(1)(i);
   ddd(2)=p.Column(2)(i);
   ddd(3)=p.Column(3)(i);
   ddd(4)=p.Column(1)(i)+p.Column(2)(i);
   ddd(5)=p.Column(1)(i)+p.Column(3)(i);
   ddd(6)=p.Column(2)(i)+p.Column(3)(i);
   ddd(7)=0;
   ddd(8)=p.Column(1)(i)+p.Column(2)(i)+p.Column(3)(i);
   min(i)=Min(ddd);max(i)=Max(ddd);
   t=min(i)/abc(i);if(abs(t-int(t))>0.0001){min(i)=(int(t)-1.0)*abc(i);}
   t=max(i)/abc(i);if(abs(t-int(t))>0.0001){max(i)=(int(t)+1.0)*abc(i);}
  }
  max(1)=min(1)+(max(1)-min(1))*scale_view_1;
  max(2)=min(2)+(max(2)-min(2))*scale_view_2;
  max(3)=min(3)+(max(3)-min(3))*scale_view_3;
  // determine ijkmin ijkmax by calculating the 8 corners of the  quader
  // in terms of primitive lattice 
  // i*p.Column(1)+j*p.Column(2)+k*p.Column(3)=cornerpointvector ... i,j,k =?
  // ijk=p^-1*corerpointvector
  for (i=1;i<=3;++i)
  {dd0=min;              dd=p.Inverse()*dd0;ddd(1)=dd(i);
   dd0=min;dd0(1)=max(1);dd=p.Inverse()*dd0;ddd(2)=dd(i);
   dd0=min;dd0(2)=max(2);dd=p.Inverse()*dd0;ddd(3)=dd(i);
   dd0=min;dd0(3)=max(3);dd=p.Inverse()*dd0;ddd(4)=dd(i);
   dd0=max;              dd=p.Inverse()*dd0;ddd(5)=dd(i);
   dd0=max;dd0(1)=min(1);dd=p.Inverse()*dd0;ddd(6)=dd(i);
   dd0=max;dd0(2)=min(2);dd=p.Inverse()*dd0;ddd(7)=dd(i);
   dd0=max;dd0(3)=min(3);dd=p.Inverse()*dd0;ddd(8)=dd(i);
   ijkmin(i)=Min(ddd);ijkmax(i)=Max(ddd);
  }  
}

Vector spincf::xy(Vector xyz,int orientation,Vector min,Vector max,float bbwidth,float bbheight)
 {Vector p(1,2);
  switch(orientation)
  {case 1: p(1)=(xyz(1)-min(1))/(max(1)-min(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(2)-min(2))/(max(2)-min(2))*bbheight*0.8+bbheight*0.15;
   break;
   case 2: p(1)=(xyz(1)-min(1))/(max(1)-min(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)-min(3))/(max(3)-min(3))*bbheight*0.8+bbheight*0.15;
    break;
   case 3: p(1)=(xyz(2)-min(2))/(max(2)-min(2))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)-min(3))/(max(3)-min(3))*bbheight*0.8+bbheight*0.15;
    break;
   case 4: p(1)=(xyz(1)+(xyz(3)-min(3))*0.1-min(1))/(max(1)+(max(3)-min(3))*0.1-min(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(2)+(xyz(3)-min(3))*0.15-min(2))/(max(2)+(max(3)-min(3))*0.15-min(2))*bbheight*0.8+bbheight*0.15;
    break;
   case 5: p(1)=(xyz(1)+(xyz(2)-min(2))*0.1-min(1))/(max(1)+(max(2)-min(2))*0.1-min(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)+(xyz(2)-min(2))*0.15-min(3))/(max(3)+(max(2)-min(2))*0.15-min(3))*bbheight*0.8+bbheight*0.15;
    break;
   case 6: p(1)=(xyz(2)+(xyz(1)-min(1))*0.1-min(2))/(max(2)+(max(1)-min(1))*0.1-min(2))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)+(xyz(1)-min(1))*0.15-min(3))/(max(3)+(max(1)-min(1))*0.15-min(3))*bbheight*0.8+bbheight*0.15;
    break;
   default: p=0;
  }
  return p;
 }

Vector spincf::pos(int i, int j, int k, int l,Vector & abc,Matrix & r,float * x,float *y,float*z)
{//returns position of atom l at lattice site (i j k) (Angstrom)
Vector dd(1,3),dd0(1,3);
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  
         dd(1)=x[l]*abc(1);
         dd(2)=y[l]*abc(2);
         dd(3)=z[l]*abc(3);
         dd+=p.Column(1)*((double)(i-1)/nofa)+p.Column(2)*((double)(j-1)/nofb)+p.Column(3)*((double)(k-1)/nofc);
return dd;
}

void spincf::eps3d(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z,int orientation, Vector & gJ)
 {// function to plot spins in a 3d manner
  // orientation:1 ab 2 ac 3 bc projection
  //             4 ab 5 ac 6 bc side view
  int i,j,k,l,m,maxm;char r1,r2,r3;
  
  Vector a(1,2);
  Vector b(1,2),c(1,3);
  double scale,d,bbheight,bbwidth,t;
 
 // determine scale factor of moments
  scale=0;
  for (i=1;i<=nofa;++i)
    for (j=1;j<=nofb;++j)
     for (k=1;k<=nofc;++k)
      for(l=1;l<=nofatoms;++l)
      {c=magmom(i,j,k,l,gJ);
       if ((d=Norm(c))>scale)scale=d;
      } 
  scale=0.5/(scale+0.01);

 

  // determine max(1,2,3) min(1,2,3) (vector in Angstroem describing a quader) for viewing magnetic unit cell
  Vector max(1,3),min(1,3),dd(1,3),max_min(1,3);
  Vector ddd(1,8),xyz(1,3),dd0(1,3),ijkmax(1,3),ijkmin(1,3);

  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  calc_minmax(min,max,ijkmin,ijkmax,p,abc);
  max_min=max-min;    

 //determine bounding box for  specific view
  bbwidth=700;
  switch(orientation)
       { case 1 : 
    		 {bbheight=max_min(2)/max_min(1)*bbwidth;r1='a';r3='b';}
	         break;
	 case 2 : 
	         {bbheight=max_min(3)/max_min(1)*bbwidth;r1='a';r3='c';}
                 break;   
	 case 3 : 
	         {bbheight=max_min(3)/max_min(2)*bbwidth;r1='b';r3='c';}
                 break;
	 case 4 : 
	         {bbheight=(max_min(2)+max_min(3)*0.1)/(max_min(3)*0.15+max_min(1))*bbwidth;r1='a';r3='b';r2='c';}
                  break;
	 case 5 : 
	         {bbheight=(max_min(3)+max_min(2)*0.1)/(max_min(2)*0.15+max_min(1))*bbwidth;r1='a';r3='c';r2='b';}
                  break;
	 case 6 : 
	         {bbheight=(max_min(3)+max_min(1)*0.1)/(max_min(1)*0.15+max_min(2))*bbwidth;r1='b';r3='c';r2='a';}
	          break;
	 default:  return;
	}
  
  fprintf(fout,"%s!PS-Adobe-2.0 EPSF-2.0\n","%");
  fprintf(fout,"%sBoundingBox:0 0 %i %i\n","%%",(int)bbwidth,(int)bbheight);
  fprintf(fout,"%sTitle: %s\n","%%",text);
  fprintf(fout,"%sEndComments\n","%%");
  fprintf(fout,"/mm {72 mul 25.4 div} bind def\n");
  fprintf(fout,"/mx {1.2 add 50 mul mm} bind def\n");
  fprintf(fout,"/my {-0.3 add 50 mul mm} bind def\n");
 

  // draw abc coordinate label 
   fprintf(fout,"/Helvetica findfont\n15 scalefont setfont\n");
   fprintf(fout,"-0.89 mx 0.38 my moveto \n (%c) show \n",r1);
   fprintf(fout,"-0.99 mx 0.52 my moveto \n (%c) show \n",r3);
   a(1)=-0.98;a(2)=0.4;
   b(1)=a(1);b(2)=0.5;epsarrow(fout,a,b);
   b(1)=-0.92;b(2)=0.4;epsarrow(fout,a,b);
   if (orientation>3){ fprintf(fout,"-0.93 mx 0.46 my moveto \n (%c) show \n",r2);
                      b(1)=-0.94;b(2)=0.45;epsarrow(fout,a,b);
		     }
  fprintf(fout,"-0.7 mx 0.38 my moveto \n (%s) show \n",text);

  fprintf(fout,"/mx {} bind def\n");
  fprintf(fout,"/my {} bind def\n");

  // draw frame around min vs max   (quader)
  fprintf(fout,"0.3 setlinewidth\n");
   a=xy(min,orientation, min, max,bbwidth,bbheight);
   dd=min;dd(1)=max(1);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=min;dd(2)=max(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=min;dd(3)=max(3);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(max,orientation, min, max,bbwidth,bbheight);
   dd=max;dd(1)=min(1);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=max;dd(2)=min(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=max;dd(3)=min(3);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=min;dd(2)=max(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=min;dd(1)=max(1);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=max;dd(2)=min(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=min;dd(3)=max(3);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=max;dd(1)=min(1);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=min;dd(2)=max(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
 
   
  // draw frame around primitive unit cell
  fprintf(fout,"1 setlinewidth\n");
   dd=0;
   a=xy(dd,orientation, min, max,bbwidth,bbheight);
   b=xy(p.Column(1),orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   b=xy(p.Column(2),orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   b=xy(p.Column(3),orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(2)+p.Column(3);a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=p.Column(1)+p.Column(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(3);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(2)+p.Column(3);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(2);a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=p.Column(1);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(3);a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=p.Column(1);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(3);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(2)+p.Column(3);a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=p.Column(2);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(3);b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
  
  
  // plot atoms and moments in region xmin to xmax (quader)
int i1,j1,k1;
//i2,k2,j2,i1true,j1true,k1true;
   fprintf(fout,"/Helvetica findfont\n %i scalefont setfont\n",(int)(1000/nofa/nofb/nofc+1));

//these lines do not work if primitive lattice angles are > 90 deg ...
//i1true=1;for (i1=0;i1true==1;++i1){i1true=0;for(i2=-1;i2<=1;i2+=2){if (i1==0){i2=2;}
//j1true=1;for (j1=0;j1true==1;++j1){j1true=0;for(j2=-1;j2<=1;j2+=2){if (j1==0){j2=2;}
//k1true=1;for (k1=0;k1true==1;++k1){k1true=0;for(k2=-1;k2<=1;k2+=2){if (k1==0){k2=2;}
//   dd0=p.Column(1)*(double)(i2*i1)+p.Column(2)*(double)(j2*j1)+p.Column(3)*(double)(k2*k1);
for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
//printf("%i %i %i %i %i %i %i %i %i\n",i1,j1,k1,(int)ijkmin(1),(int)ijkmin(2),(int)ijkmin(3),(int)ijkmax(1),(int)ijkmax(2),(int)ijkmax(3));
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd+=dd0;
	    if(dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)
            {c=magmom(i,j,k,l,gJ);
              xyz(1)=dd(1)+scale*c(1);
              xyz(2)=dd(2)+scale*c(2);
              xyz(3)=dd(3)+scale*c(3);
	      a=xy(xyz,orientation, min, max,bbwidth,bbheight);
              xyz(1)=dd(1)-scale*c(1);
              xyz(2)=dd(2)-scale*c(2);
              xyz(3)=dd(3)-scale*c(3);
              b=xy(xyz,orientation, min, max,bbwidth,bbheight);
              epsarrow(fout,a,b);   
	     }
	  }
       }}}           
 }}}
  
  
fprintf(fout,"showpage\n");
  

 }
double spincf::nndist(float * x, float * y, float * z,Vector & abc,Matrix & p,Vector &dd)
{int i,j,k,l;double d,mindist=1e10;Vector ddl(1,3);
         
 for(i=-1;i<=1;++i)for(j=-1;j<=1;++j)for(k=-1;k<=1;++k)for(l=1;l<=nofatoms;++l)
 {ddl(1)=x[l]*abc(1);
  ddl(2)=y[l]*abc(2);
  ddl(3)=z[l]*abc(3);
  ddl+=p.Column(1)*((double)(i-1)/nofa)+p.Column(2)*((double)(j-1)/nofb)+p.Column(3)*((double)(k-1)/nofc);
  d=Norm(dd-ddl);if(d>0.0001&&d<mindist)mindist=d;
 }
return mindist;
}


// output for javaview
void spincf::jvx(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ,
                 double show_abc_unitcell,double show_primitive_crystal_unitcell,double show_magnetic_unitcell,double show_atoms,double scale_view_1,double scale_view_2,double scale_view_3,
                 int showprim,double phase,spincf & savev_real,spincf & savev_imag,double amplitude,Vector & hkl,
                 double spins_show_ellipses,double spins_show_direction_of_static_moment) 
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
{char *cffilenames[0];
 jvx_cd(fout,text,abc,r,x,y,z,gJ,
                 show_abc_unitcell,show_primitive_crystal_unitcell,show_magnetic_unitcell,
                 show_atoms,scale_view_1,scale_view_2,scale_view_3,
                 showprim,phase,savev_real,savev_imag,amplitude,hkl,
                 spins_show_ellipses,spins_show_direction_of_static_moment,cffilenames,0.0,1.0);
}


void spincf::jvx_cd(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ,
                 double show_abc_unitcell,double show_primitive_crystal_unitcell,double show_magnetic_unitcell,double show_atoms,double scale_view_1,double scale_view_2,double scale_view_3,
                 int showprim,double phase,spincf & savev_real,spincf & savev_imag,double amplitude,Vector & hkl,
                 double spins_show_ellipses,double spins_show_direction_of_static_moment,char ** cffilenames,double show_chargedensity,double show_spindensity)
{ int i,j,k,l,ctr=0,maxm,m;int i1,j1,k1;
 // some checks
 if(nofatoms!=savev_real.nofatoms||nofa!=savev_real.na()||nofb!=savev_real.nb()||nofc!=savev_real.nc()||
    nofatoms!=savev_imag.nofatoms||nofa!=savev_imag.na()||nofb!=savev_imag.nb()||nofc!=savev_imag.nc()||
    nofcomponents<savev_real.nofcomponents||savev_real.nofcomponents!=savev_imag.nofcomponents)
    {fprintf(stderr,"Error creating jvx movie files: eigenvector read from .qev file does not match dimension of spins structure read from sps file\n");exit(1);}

  Vector max(1,3),min(1,3),ijkmax(1,3),ijkmin(1,3),max_min(1,3),dd(1,3),dd0(1,3),c(1,3),xyz(1,3);
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  calc_minmax_scale(min,max,ijkmin,ijkmax,p,abc,scale_view_1,scale_view_2,scale_view_3);
   if(showprim==1){ijkmin(1)=1;ijkmin(2)=1;ijkmin(3)=1;ijkmax(1)=-2+(int)(scale_view_1);ijkmax(2)=-2+(int)(scale_view_2);ijkmax(3)=-2+(int)(scale_view_3);} // show only primitive magnetic unit cell
  max_min=max-min;    
  double scale=0,d,mindist=1e10;

fprintf(fout,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n"); 
fprintf(fout,"<jvx-model>\n"); 
fprintf(fout,"  <title>%s</title>\n",text); 
fprintf(fout,"  <geometries>\n"); 

if(show_abc_unitcell>0)
 {  // plot frame around crystallographic unit cell
fprintf(fout,"    <geometry name=\"crystallographic unit cell\">\n"); 
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n"); 
fprintf(fout,"        <points>\n"); 
fprintf(fout,"          <p>  %g       %g       %g </p>\n",0.0,0.0,0.0); 
fprintf(fout,"          <p name=\"a\">  %g       %g       %g </p>\n",abc(1),0.0,0.0); 
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",abc(1),abc(2),0.0); 
fprintf(fout,"          <p name=\"b\"> %g       %g       %g </p>\n",0.0,abc(2),0.0); 
fprintf(fout,"          <p name=\"c\"> %g       %g       %g </p>\n",0.0,0.0,abc(3)); 
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",abc(1),0.0,abc(3)); 
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",abc(1),abc(2),abc(3)); 
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",0.0,abc(2),abc(3)); 
fprintf(fout,"          <thickness>0.0</thickness>\n"); 
fprintf(fout,"          <colorTag type=\"rgb\">255 0 0</colorTag>\n");
fprintf(fout,"			<labelAtt horAlign=\"head\" visible=\"show\" font=\"fixed\" verAlign=\"top\">\n"); 
fprintf(fout,"				<xOffset>0</xOffset>\n"); 
fprintf(fout,"				<yOffset>0</yOffset>\n"); 
fprintf(fout,"			</labelAtt>\n"); 
fprintf(fout,"        </points>\n"); 
fprintf(fout,"      </pointSet>\n"); 
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\">\n"); 
fprintf(fout,"        <lines>\n"); 
fprintf(fout,"          <l>0 1</l>\n"); 
fprintf(fout,"          <l>1 2</l>\n"); 
fprintf(fout,"          <l>2 3</l>\n"); 
fprintf(fout,"          <l>3 0</l>\n"); 
fprintf(fout,"          <l>0 4</l>\n"); 
fprintf(fout,"          <l>1 5</l>\n"); 
fprintf(fout,"          <l>2 6</l>\n"); 
fprintf(fout,"          <l>3 7</l>\n"); 
fprintf(fout,"          <l>4 5</l>\n"); 
fprintf(fout,"          <l>5 6</l>\n"); 
fprintf(fout,"          <l>6 7</l>\n"); 
fprintf(fout,"          <l>7 4</l>\n"); 
fprintf(fout,"          <thickness>1.0</thickness>\n"); 
fprintf(fout,"        <color type=\"rgb\">%i 0 0</color>\n",(int)(255*show_abc_unitcell));
fprintf(fout,"        </lines>\n"); 
fprintf(fout,"      </lineSet>\n"); 
fprintf(fout,"    </geometry>\n"); 
 }
if(show_primitive_crystal_unitcell>0)
 {
 // plot frame around primitive crystallographic unit cell
fprintf(fout,"    <geometry name=\"primitive crystallographic unit cell\">\n"); 
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n"); 
fprintf(fout,"        <points>\n"); 
dd=0;           fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\"r1\">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(2)/(double)nofb;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd-=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\"r2\">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd =p.Column(3)/(double)nofc;fprintf(fout,"          <p name=\"r3\">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(2)/(double)nofb;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd-=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
fprintf(fout,"			<labelAtt horAlign=\"head\" visible=\"show\" font=\"fixed\" verAlign=\"bottom\">\n"); 
fprintf(fout,"				<xOffset>0</xOffset>\n"); 
fprintf(fout,"				<yOffset>0</yOffset>\n"); 
fprintf(fout,"                          <colorTag type=\"rgb\">0 255 0</colorTag>\n");
fprintf(fout,"			</labelAtt>\n"); 
fprintf(fout,"          <thickness>0.0</thickness>\n"); 
fprintf(fout,"        </points>\n"); 
fprintf(fout,"      </pointSet>\n"); 
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n"); 
fprintf(fout,"        <lines>\n"); 
fprintf(fout,"          <l>0 1</l>\n"); 
fprintf(fout,"          <l>1 2</l>\n"); 
fprintf(fout,"          <l>2 3</l>\n"); 
fprintf(fout,"          <l>3 0</l>\n"); 
fprintf(fout,"          <l>0 4</l>\n"); 
fprintf(fout,"          <l>1 5</l>\n"); 
fprintf(fout,"          <l>2 6</l>\n"); 
fprintf(fout,"          <l>3 7</l>\n"); 
fprintf(fout,"          <l>4 5</l>\n"); 
fprintf(fout,"          <l>5 6</l>\n"); 
fprintf(fout,"          <l>6 7</l>\n"); 
fprintf(fout,"          <l>7 4</l>\n"); 
fprintf(fout,"          <thickness>1.0</thickness>\n"); 
fprintf(fout,"        <color type=\"rgb\">0 %i 0</color>\n",(int)(255*show_primitive_crystal_unitcell));
fprintf(fout,"        </lines>\n"); 
fprintf(fout,"      </lineSet>\n"); 
fprintf(fout,"    </geometry>\n"); 
}
if(show_magnetic_unitcell>0)
 {
  // plot frame around primitive magnetic unit cell
fprintf(fout,"    <geometry name=\"magnetic unit cell\">\n"); 
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n"); 
fprintf(fout,"        <points>\n"); 
dd=0;           fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(2);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd-=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd =p.Column(3);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd+=p.Column(2);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
dd-=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
fprintf(fout,"        </points>\n"); 
fprintf(fout,"      </pointSet>\n"); 
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n"); 
fprintf(fout,"        <lines>\n"); 
fprintf(fout,"          <l>0 1</l>\n"); 
fprintf(fout,"          <l>1 2</l>\n"); 
fprintf(fout,"          <l>2 3</l>\n"); 
fprintf(fout,"          <l>3 0</l>\n"); 
fprintf(fout,"          <l>0 4</l>\n"); 
fprintf(fout,"          <l>1 5</l>\n"); 
fprintf(fout,"          <l>2 6</l>\n"); 
fprintf(fout,"          <l>3 7</l>\n"); 
fprintf(fout,"          <l>4 5</l>\n"); 
fprintf(fout,"          <l>5 6</l>\n"); 
fprintf(fout,"          <l>6 7</l>\n"); 
fprintf(fout,"          <l>7 4</l>\n"); 
fprintf(fout,"          <thickness>%g</thickness>\n",show_magnetic_unitcell); 
fprintf(fout,"        </lines>\n"); 
fprintf(fout,"      </lineSet>\n"); 
fprintf(fout,"    </geometry>\n"); 
}
  // plot atoms in region xmin to xmax (quader)
fprintf(fout,"    <geometry name=\"ions\">\n"); 
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n"); 
fprintf(fout,"        <points>\n"); 
  for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd+=dd0;
	    if((dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)||showprim==1&&scale_view_1>(double)(i1*nofa+i)/nofa&&scale_view_2>(double)(j1*nofb+j)/nofb&&scale_view_3>(double)(k1*nofc+k)/nofc)
            {// determine scale factor of moments
               c=magmom(i,j,k,l,gJ);if ((d=Norm(c))>scale)scale=d;      
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));                                                                                                                              
     ++ctr;
	     }
	  }
       }}}           
  }}}
fprintf(fout,"          <thickness>3.0</thickness>\n"); 
fprintf(fout,"        </points>\n"); 
fprintf(fout,"        <colors type=\"rgb\">\n"); 
  for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if((dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)||showprim==1&&scale_view_1>(double)(i1*nofa+i)/nofa&&scale_view_2>(double)(j1*nofb+j)/nofb&&scale_view_3>(double)(k1*nofc+k)/nofc)
            {// determine scale factor of moments
               c=magmom(i,j,k,l,gJ);if ((d=Norm(c))>scale)scale=d;      
fprintf(fout,"          <c>  %i       %i       %i </c>\n",(int)(255*show_atoms),(int)(show_atoms*((l*97)%256)),0);
	     }
	  }
       }}}           
  }}}
fprintf(fout,"        </colors>\n"); 
fprintf(fout,"      </pointSet>\n"); 
fprintf(fout,"    </geometry>\n"); 

  // plot magnetic moments
         for(l=1;l<=nofatoms;++l) // determine mindistance to neighbors
	 {dd(1)=x[l]*abc(1);dd(2)=y[l]*abc(2);dd(3)=z[l]*abc(3);
          if(mindist>(d=nndist(x,y,z,abc,p,dd)))mindist=d;
         }
         scale=0.4*mindist/(scale+0.001); // get a good scale factor 
if(show_spindensity>0){
fprintf(fout,"    <geometry name=\"magnetic moments\">\n"); 
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"show\">\n"); 
fprintf(fout,"        <points>\n"); 
 ctr=0;
 for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
 for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
 for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if((dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)||showprim==1&&scale_view_1>(double)(i1*nofa+i)/nofa&&scale_view_2>(double)(j1*nofb+j)/nofb&&scale_view_3>(double)(k1*nofc+k)/nofc)
            {double QR; QR=hkl(1)*dd(1)/abc(1)+hkl(2)*dd(2)/abc(2)+hkl(3)*dd(3)/abc(3);QR*=2*PI;
             xyz=magmom(i,j,k,l,gJ)+amplitude*(cos(-phase+QR)*savev_real.magmom(i,j,k,l,gJ)+sin(phase-QR)*savev_imag.magmom(i,j,k,l,gJ));
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
//fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)-xyz(1)*scale,dd(2)-xyz(2)*scale,dd(3)-xyz(3)*scale); 
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)+xyz(1)*scale,dd(2)+xyz(2)*scale,dd(3)+xyz(3)*scale);                                                                                                                              
	     ++ctr;

	     }
	  }
       }}}           
  }}}
fprintf(fout,"          <thickness>6.0</thickness>\n"); 
fprintf(fout,"        </points>\n"); 
fprintf(fout,"      </pointSet>\n"); 
fprintf(fout,"      <lineSet  arrow=\"show\" line=\"show\">\n"); 
fprintf(fout,"        <lines>\n"); 
  for(i=0;i<ctr;++i)fprintf(fout,"          <l>%i %i</l>\n",2*i,2*i+1); 
fprintf(fout,"        <color type=\"rgb\">0 0 255</color>\n");
fprintf(fout,"          <thickness>2.0</thickness>\n"); 
fprintf(fout,"        </lines>\n"); 
fprintf(fout,"      </lineSet>\n"); 
fprintf(fout,"    </geometry>\n"); 
                     }

if(spins_show_direction_of_static_moment>0)
 {
// plot a line along static magnetic moments for comparison
fprintf(fout,"    <geometry name=\"static magnetic moments\">\n"); 
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n"); 
fprintf(fout,"        <points>\n"); 
 ctr=0;
 for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
 for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
 for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if((dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)||showprim==1&&scale_view_1>(double)(i1*nofa+i)/nofa&&scale_view_2>(double)(j1*nofb+j)/nofb&&scale_view_3>(double)(k1*nofc+k)/nofc)
            {double QR; QR=hkl(1)*dd(1)/abc(1)+hkl(2)*dd(2)/abc(2)+hkl(3)*dd(3)/abc(3);QR*=2*PI;
             xyz=magmom(i,j,k,l,gJ);
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3)); 
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)+xyz(1)*scale,dd(2)+xyz(2)*scale,dd(3)+xyz(3)*scale);                                                                                                                              
	     ++ctr;

	     }
	  }
       }}}           
  }}}
fprintf(fout,"          <thickness>6.0</thickness>\n"); 
fprintf(fout,"        </points>\n"); 
fprintf(fout,"      </pointSet>\n"); 
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n"); 
fprintf(fout,"        <lines>\n"); 
  for(i=0;i<ctr;++i)fprintf(fout,"          <l>%i %i</l>\n",2*i,2*i+1); 
fprintf(fout,"          <thickness>%g</thickness>\n",spins_show_direction_of_static_moment); 
fprintf(fout,"        </lines>\n"); 
fprintf(fout,"      </lineSet>\n"); 
fprintf(fout,"    </geometry>\n"); 
}
if(spins_show_ellipses>0)
 {
// plot an ellipse along path of moment
fprintf(fout,"    <geometry name=\"ellipses\">\n"); 
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n"); 
fprintf(fout,"        <points>\n"); 
 ctr=0;
 for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
 for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
 for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if((dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)||showprim==1&&scale_view_1>(double)(i1*nofa+i)/nofa&&scale_view_2>(double)(j1*nofb+j)/nofb&&scale_view_3>(double)(k1*nofc+k)/nofc)
            {double QR; QR=hkl(1)*dd(1)/abc(1)+hkl(2)*dd(2)/abc(2)+hkl(3)*dd(3)/abc(3);QR*=2*PI;
             int phi; 
             for(phi=0;phi<=16;phi++)
             {
             xyz=magmom(i,j,k,l,gJ)+amplitude*(cos(-(double)phi*2*3.1415/16+QR)*savev_real.magmom(i,j,k,l,gJ)+sin((double)phi*2*3.1415/16-QR)*savev_imag.magmom(i,j,k,l,gJ));
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)+xyz(1)*scale,dd(2)+xyz(2)*scale,dd(3)+xyz(3)*scale);                                                                                                                              
	     }++ctr;

	     }
	  }
       }}}           
  }}}
fprintf(fout,"          <thickness>6.0</thickness>\n"); 
fprintf(fout,"        </points>\n"); 
fprintf(fout,"      </pointSet>\n"); 
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n"); 
fprintf(fout,"        <lines>\n"); 
  for(i=0;i<ctr;++i)for(j=0;j<16;++j)fprintf(fout,"          <l>%i %i</l>\n",17*i+j,17*i+j+1); 
fprintf(fout,"          <thickness>1.0</thickness>\n"); 
fprintf(fout,"        <color type=\"rgb\">0 %i %i</color>\n",(int)(100*spins_show_ellipses),(int)(100*spins_show_ellipses));
fprintf(fout,"        </lines>\n"); 
fprintf(fout,"      </lineSet>\n"); 
fprintf(fout,"    </geometry>\n"); 
}
if(show_chargedensity>0)
{int ii,tt,ff;
  double dtheta=0.2; //stepwidth to step surface
  double dfi=0.2;
  
for(l=1;l<=nofatoms;++l)
 {fprintf(fout,"    <geometry name=\"chargedensities in primitive magnetic unit cell - atom %i\">\n",l); 
  fprintf(fout,"<pointSet color=\"hide\" point=\"show\" dim=\"1\">\n");
  fprintf(fout,"<points >\n");
  double radius=0;double dx,dy,dz,R,fi,theta;
  extract(cffilenames[l],"radius",radius);  
  if(radius!=0)
  {     double rp=0.5*cbrt(abs(radius));
        for (i=1;i<=(1+(nofa-1)*scale_view_1);++i){for(j=1;j<=(1+(nofb-1)*scale_view_2);++j){for(k=1;k<=(1+(nofc-1)*scale_view_3);++k){
        dd=pos(i,j,k,l, abc, r,x,y,z);
        for(tt=0;tt<=3.1415/dtheta;++tt){for(ff=0;ff<=2*3.1415/dfi;++ff){
             theta=(double)tt*dtheta;fi=(double)ff*dfi;
             dx=rp*sin(theta)*cos(fi)+dd(3);dy=rp*sin(theta)*sin(fi)+dd(1);dz=rp*cos(theta)+dd(2);
             fprintf(fout,"<p>%4g %4g %4g</p>\n",dy,dz,dx);
             if(tt==0){ff=(int)(2*3.1415/dfi+1);}
             }}
        }}}
  }
  else
  {jjjpar ionpar(x[l],y[l],z[l],cffilenames[l]);
   chargedensity cd(dtheta,dfi);int ndd;
   for (i=1;i<=1+(nofa-1)*scale_view_1;++i){for(j=1;j<=1+(nofb-1)*scale_view_2;++j){for(k=1;k<=1+(nofc-1)*scale_view_2;++k){
   dd=pos(i,j,k,l, abc, r,x,y,z);
   Vector moments(1,nofcomponents); 
   double QR; QR=hkl(1)*dd(1)/abc(1)+hkl(2)*dd(2)/abc(2)+hkl(3)*dd(3)/abc(3);QR*=2*PI;
   for(ndd=1;ndd<=savev_real.nofcomponents;++ndd)
   {moments(ndd)=moment(i,j,k,l)(ndd)+amplitude*(cos(-phase+QR)*savev_real.moment(i,j,k,l)(ndd)+sin(phase-QR)*savev_imag.moment(i,j,k,l)(ndd));}
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
   cd.calc_cd_surface(moments,ionpar,0.05);
   for(ii=1;ii<=cd.nofpoints();++ii)
     {// here we calculate the chargedensity of ion 
     R=cd.rtf(ii)(1);theta=cd.rtf(ii)(2);fi=cd.rtf(ii)(3);
     if(ionpar.module_type==2){// mind abc||yzx in module cfield
     dz=R*sin(theta)*cos(fi)+dd(3);dx=R*sin(theta)*sin(fi)+dd(1);dy=R*cos(theta)+dd(2);
                              }
     else 
                              {// mind abc||xyz in other cases ...
     dz=R*cos(theta)+dd(3);dx=R*sin(theta)*cos(fi)+dd(1);dy=R*sin(theta)*sin(fi)+dd(2);
                              }
     fprintf(fout,"<p>%4g %4g %4g</p>\n",dx,dy,dz);
     }
   }}}
  }
    fprintf(fout,"<thickness>0.0</thickness><color type=\"rgb\">255 0 0</color><colorTag type=\"rgb\">255 0 255</colorTag>\n");
    fprintf(fout,"</points>			</pointSet>\n");
    fprintf(fout,"<faceSet face=\"show\" edge=\"show\">\n");
    fprintf(fout,"<faces >\n");
    int offset=0;
    for(i=1;i<=(1+(nofa-1)*scale_view_1)*(1+(nofb-1)*scale_view_2)*(1+(nofc-1)*scale_view_3);++i)
    {int ntt,nff,pointnr,ffnr,p1,p2,p3,p4;
    ntt=(int)(3.1415/dtheta);
    nff=(int)(2*3.1415/dfi);
    pointnr=ntt*(nff+1);
    ffnr=nff+1;
    for(tt=1;tt<=ntt;++tt){for(ff=0;ff<=nff;++ff){
    p1 = ff + 1 + (tt - 2) * ffnr+offset;
    p2 = ff + 2 + (tt - 2) * ffnr+offset;
    p3 = ff + 2 + (tt - 1) * ffnr+offset;
    p4 = ff + 1 + (tt - 1) * ffnr+offset;
    if (ff==nff){p3 = p3 - ffnr; p2 = p2 - ffnr;}
    if (tt==1) {p1 = offset; p2 = offset;}
    fprintf(fout,"<f> %i %i %i %i </f>\n",p1,p2,p3,p4);
    }}
    offset+=pointnr+1;
 }     
 //fprintf(fout,"<color type=\"rgb\">100 230 255</color>\n");
 if(radius>0)
 {fprintf(fout,"<color type=\"rgb\"> 255 0 0</color>\n");}
 else if (radius<0)
 {fprintf(fout,"<color type=\"rgb\">0  0 255</color>\n");}
 else
 {fprintf(fout,"<color type=\"rgb\">%i %i %i </color>\n",(int)(255*show_chargedensity),(int)(show_chargedensity*((l*97)%256)),0);
 }
 fprintf(fout,"<colorTag type=\"rgb\">255 0 255</colorTag>\n");
 fprintf(fout,"</faces></faceSet>\n");
 fprintf(fout,"    </geometry>\n"); 
 }

}


fprintf(fout,"  </geometries>\n"); 
fprintf(fout,"</jvx-model>\n"); 
}

//***********************************************************************************************************************************

// output for fullprof studio
void spincf::fst(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ) //print std file to stream
{int i,j,k,l,ctr=1;

int maxm,m;  
  Vector max(1,3),min(1,3),dd(1,3),max_min(1,3);
  Vector xyz(1,3),dd0(1,3),ijkmax(1,3),ijkmin(1,3);
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  calc_minmax(min,max,ijkmin,ijkmax,p,abc);

  max_min=max-min;    


fprintf(fout,"!   FILE for FullProf Studio: generated automatically by McPhase\n"); 
fprintf(fout,"!Title: %s \n",text);                                                                                         
fprintf(fout,"SPACEG P 1           \n");
fprintf(fout,"CELL     %g    %g    %g  90.0000  90.0000 90.0000   DISPLAY MULTIPLE\n",max_min(1),max_min(2),max_min(3));
fprintf(fout,"BOX   -0.15  1.15   -0.15  1.15    -0.15  1.15 \n");

  // plot atoms in region xmin to xmax (quader)
int i1,j1,k1;
for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){

   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);

      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if(dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)
            {dd(1)/=max_min(1);dd(2)/=max_min(2);dd(3)/=max_min(3);

fprintf(fout,"ATOM DY%i    RE       %g       %g       %g        \n",ctr,dd(1),dd(2),dd(3)); 
                                                                                                                             
	     ++ctr;

	     }
	  }
       }}}           
 }}}

fprintf(fout," \n");
fprintf(fout,"{\n");
fprintf(fout,"LATTICE P\n");
fprintf(fout,"K     0.00000   0.00000   0.00000\n");
fprintf(fout,"SYMM  x,y,z\n");
fprintf(fout,"MSYM  u,v,w,0.0\n");

// plot moments in region xmin to xmax (quader)
for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){

   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);

      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if(dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)
            {dd(1)/=max_min(1);dd(2)/=max_min(2);dd(3)/=max_min(3);
//             i1true=1;j1true=1;k1true=1;       

//	    a=xy(dd,orientation, min, max,bbwidth,bbheight);
            xyz=magmom(i,j,k,l,gJ);
fprintf(fout,"MATOM DY%i    DY      %g       %g       %g   GROUP\n",ctr,dd(1),dd(2),dd(3));
fprintf(fout,"SKP           1  1  %g       %g       %g       0.00000  0.00000  0.00000    0.00000\n",xyz(1),xyz(2),xyz(3));
	     ++ctr;

	     }
	  }
       }}}           
 }}}
//}}}
fprintf(fout,"}\n");
} 


void spincf::fstprim(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ) //print std file to stream
{int i,j,k,l,ctr=1;
int maxm,m; 
double alpha,beta,gamma;
  // determine max(1,2,3) min(1,2,3) (vector in Angstroem describing a quader) for viewing magnetic unit cell
  Vector ddd(1,8),xyz(1,3),xyz0(1,3),dd0(1,3),dd(1,3);

  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);

gamma=180/3.1415926*acos(p.Column(1)*p.Column(2)/Norm(p.Column(1))/Norm(p.Column(2)));
beta=180/3.1415926*acos(p.Column(1)*p.Column(3)/Norm(p.Column(1))/Norm(p.Column(3)));
alpha=180/3.1415926*acos(p.Column(2)*p.Column(3)/Norm(p.Column(2))/Norm(p.Column(3)));


fprintf(fout,"!   FILE for FullProf Studio: generated automatically by McPhase\n"); 
fprintf(fout,"!Title: %s \n",text);                                                                                         
fprintf(fout,"SPACEG P 1           \n");
fprintf(fout,"CELL     %g    %g    %g  %g %g %g   DISPLAY MULTIPLE\n",Norm(p.Column(1)),Norm(p.Column(2)),Norm(p.Column(3)),alpha,beta,gamma);
fprintf(fout,"BOX   -0.15  1.15   -0.15  1.15    -0.15  1.15 \n");

  // plot atoms 
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd0=p.Inverse()*dd;
fprintf(fout,"ATOM DY%i    RE       %g       %g       %g        \n",ctr,dd0(1),dd0(2),dd0(3));
	     ++ctr;

	     }
	  }
       }}      

fprintf(fout," \n");
fprintf(fout,"{\n");
fprintf(fout,"LATTICE P\n");
fprintf(fout,"K     0.00000   0.00000   0.00000\n");
fprintf(fout,"SYMM  x,y,z\n");
fprintf(fout,"MSYM  u,v,w,0.0\n");

// plot moments 
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd0=p.Inverse()*dd;
          xyz=magmom(i,j,k,l,gJ);
          xyz0=p.Inverse()*xyz; xyz0(1)*=Norm(p.Column(1));xyz0(2)*=Norm(p.Column(2));xyz0(3)*=Norm(p.Column(3));

fprintf(fout,"MATOM DY%i    DY      %g       %g       %g   GROUP\n",ctr,dd0(1),dd0(2),dd0(3));
fprintf(fout,"SKP           1  1  %g       %g       %g       0.00000  0.00000  0.00000    0.00000\n",xyz0(1),xyz0(2),xyz0(3));
	     ++ctr;

	     }
	  }
       }}           
fprintf(fout,"}\n");
} 


//-----------------------------------------------------------------------
//  numeric output of spinconfiguration to file
void spincf::print(FILE * fout) //print spinconfiguration to stream
{int i,j,k,l;
 for (k=1;k<=nofc;++k)
 {for (j=1;j<=nofb;++j)
  {for (l=1;l<=nofcomponents*nofatoms;++l)
   {for (i=1;i<=nofa;++i)
      {fprintf(fout," %4.4f",mom[in(i,j,k)](l));
       }
    fprintf(fout,"\n");
    }
   }
 fprintf(fout,"\n"); //new line to separate ab planes
 }
// fprintf(fout,"\n"); //new line to end spinconfiguration - removed aug 07
}               

void spincf::printall(FILE * fout,Vector & abc,Matrix & r,float * x,float *y,float*z, char ** cffilenames,float * gJ) //print spinconfiguration to stream
{ int i,j,k,l,lc,m,maxm;

 // determine primitive magnetic unit cell
  Vector dd(1,3);
  Vector xyz(1,3),dd0(1,3),mmm(1,3);
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  

 fprintf(fout,"#nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",nofa,nofb,nofc,nofatoms*nofa*nofb*nofc);
 fprintf(fout,"#{atom file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3]  <Ma> <Mb> <Mc> [mb] <Ja> <Jb> <Jc> <Jd> <Je> ...\n");
 
   // output atoms and moments in primitive unit cell
  for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd0=p.Inverse()*dd;dd0(1)*=nofa;dd0(2)*=nofb;dd0(3)*=nofc;
              fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
	              cffilenames[l],dd(1)/abc(1),dd(2)/abc(2),dd(3)/abc(3),dd0(1),dd0(2),dd0(3));
             if(gJ[l]!=0)
              {fprintf(fout," %4.4f",gJ[l]*mom[in(i,j,k)](1+nofcomponents*(l-1)));
               if(nofcomponents>=2){fprintf(fout," %4.4f",gJ[l]*mom[in(i,j,k)](2+nofcomponents*(l-1)));}else{fprintf(fout," %4.4f",0.0);}
               if(nofcomponents>=2){fprintf(fout," %4.4f",gJ[l]*mom[in(i,j,k)](3+nofcomponents*(l-1)));}else{fprintf(fout," %4.4f",0.0);}
              }
             else   // if gJ=0 it means we have so print out total moment
              { //load magnetic moment into vector mmm
               if(nofcomponents>6){maxm=6;}else{maxm=nofcomponents;}
                mmm=0;
                for(m=1;m<=maxm;++m){if(m==2||m==4||m==6){mmm((m+1)/2)+=mom[in(i,j,k)](nofcomponents*(l-1)+m);}
                                     else                {mmm((m+1)/2)+=2*mom[in(i,j,k)](nofcomponents*(l-1)+m);}
                                     }

               fprintf(fout," %4.4f %4.4f %4.4f",mmm(1),mmm(2),mmm(3));
              }
             {for (lc=1;lc<=nofcomponents;++lc)
              {fprintf(fout," %4.4f",mom[in(i,j,k)](lc+nofcomponents*(l-1)));}
              fprintf(fout,"\n");
	     }
	    
	 }
  }}}          

}               

/**************************************************************************/

//zuweisung
spincf & spincf::operator= (const spincf & op2)
{int i,j,k;
 nofa=op2.nofa; nofb=op2.nofb; nofc=op2.nofc;
 mxa=op2.mxa; mxb=op2.mxb; mxc=op2.mxc;
 nofatoms=op2.nofatoms;
 wasstable=op2.wasstable;
 nofcomponents=op2.nofcomponents;
  delete []mom;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(k=0;k<=mxa*mxb*mxc;++k){mom[k]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL)
    {fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);}
  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mom[in(i,j,k)]=op2.mom[in(i,j,k)];} 
    }
  }           
  return *this;
}

//addition
spincf & spincf::operator + (const spincf & op2)
{int i,j,k;
 if (nofa!=op2.nofa||nofb!=op2.nofb||nofc!=op2.nofc||nofatoms!=op2.nofatoms||nofcomponents!=op2.nofcomponents)
 {fprintf (stderr,"Error in adding spincfonfigurations - not equal dimension\n");
 exit (EXIT_FAILURE);}
 static spincf op1((*this)); 
 for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {op1.mom[in(i,j,k)]=mom[in(i,j,k)]+op2.mom[in(i,j,k)];} 
    }
  }           
  return op1;
}

// multiplication of spinconfiguration with constant
spincf & spincf::operator * (const double factor)
{int i,j,k;
 static spincf op1((*this)); 
  for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {op1.mom[in(i,j,k)]=mom[in(i,j,k)]*factor;} 
    }
  }           
  return op1;
}

//vergleichsoperator
int spincf::operator==(spincf & op1)
{//untersuchen der periodizitaet
int nna,nnb,nnc,chka,chkb,chkc,ma,mb,mc;
if (nofatoms!=op1.nofatoms||nofa!=op1.nofa||nofb!=op1.nofb||nofc!=op1.nofc) return 0;
// vergleichen der spins unter zuhilfename von spequal
for (chka=1;chka<=nofa;++chka)
 {for (chkb=1;chkb<=nofb;++chkb)
    {for (chkc=1;chkc<=nofc;++chkc)
    
     {for (nna=1;nna<=nofa;++nna)
         {for (nnb=1;nnb<=nofb;++nnb)
	     {for (nnc=1;nnc<=nofc;++nnc)
              {ma=mod(nna+chka,nofa);
	       mb=mod(nnb+chkb,nofb);
	       mc=mod(nnc+chkc,nofc);
               if (!spequal(mom[in(nna,nnb,nnc)],op1.mom[op1.in(ma,mb,mc)]))
                  {goto nextchk;}
	       if (nna==nofa&&nnb==nofb&&nnc==nofc) return 1;
               }
	      }
	  }     	  
nextchk: ;
     }
    }
  }    
   
return 0; //not equal
}

//constructors
//from n1xn2xn3 unit cells with na number in the cryst. basis and nc number of components of the spin of each atom
spincf::spincf (int n1,int n2,int n3,int na,int nc)
{ wasstable=0;
  int l;
  nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  nofatoms=na;
  nofcomponents=nc;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];
  for(l=0;l<=mxa*mxb*mxc;++l){mom[l]=Vector(1,nofcomponents*nofatoms);mom[l]=0;}
  if (mom == NULL)
    { fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);} 
}

//kopier-konstruktor
spincf::spincf (const spincf & p)
{ int i,j,k;
  nofa=p.nofa;nofb=p.nofb;nofc=p.nofc;
  mxa=p.mxa; mxb=p.mxb; mxc=p.mxc;
  wasstable=p.wasstable;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1];for(k=0;k<=mxa*mxb*mxc;++k){mom[k]=Vector(1,nofcomponents*nofatoms);}
  if (mom == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    } 
 for (i=1;i<=nofa;++i)
  {for (j=1;j<=nofb;++j)
    {for (k=1;k<=nofc;++k)
     {mom[in(i,j,k)]=p.mom[in(i,j,k)];} 
    }
  }           

}


//destruktor
spincf::~spincf ()
{
  delete []mom;
}


