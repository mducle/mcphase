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

#define MAXNOFSPINS  200
#define SMALL 0.03
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
 if (nofcomponents>=3){m=3;}else{m=nofcomponents;}
 Vector ret(1,3);
 ret=0;
 for (i=1;i<=nofa;++i)
 { for (j=1;j<=nofb;++j)
   {for (k=1;k<=nofc;++k)
    {for (l=1;l<=nofatoms;++l)
     { for (n=1;n<=m;++n){ret(n)+=mom[in(i,j,k)](nofcomponents*(l-1)+n)*gJ(l);}
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
  mom = new Vector[mxa*mxb*mxc+1](1,nofcomponents*nofatoms);
  if (mom == NULL){fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);} 

// fill spinconfiguration with values
   for (rra=1;rra<=nofa;++rra)
     { for (rrb=1;rrb<=nofb;++rrb)
        {for (rrc=1;rrc<=nofc;++rrc)
             	{rr(1)=rra-1;rr(2)=rrb-1;rr(3)=rrc-1;
		 for(l=1;l<=nofcomponents*nofatoms;++l)
	          {mom[in(rra,rrb,rrc)](l)=copysign(nettom(l),qv*(0.1+
		  copysign(1.0,momentq0(l)+cos(phi(l)+qvector*2.0*3.141529*rr)))
				                    )                    ;	
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

  printf("n1=%i n2=%i n3=%i\n",nofa,nofb,nofc);
  
  j=fseek(fin_coq,pos,SEEK_SET); if (j!=0) return 0;

  delete []mom;
  mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1](1,nofcomponents*nofatoms);
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

void spincf::eps(FILE * fout,char * text ) //print spinconfiguration to stream
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

Vector spincf::xy(Vector & xyz,int orientation,Vector min,Vector max,float bbwidth,float bbheight)
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
  Vector nofabc(1,3),pa(1,3),pb(1,3),pc(1,3);
  Matrix p(1,3,1,3);
  int ii,jj;
  nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  for (ii=1;ii<=3;++ii)
  {for(jj=1;jj<=3;++jj) {dd(jj)=nofabc(jj)*r(ii,jj)*abc(ii);p(ii,jj)=dd(jj);}
  }
  pa=p.Column(1);  //primitive magnetic unit cell
  pb=p.Column(2);
  pc=p.Column(3);
         dd(1)=x[l]*abc(1);
         dd(2)=y[l]*abc(2);
         dd(3)=z[l]*abc(3);
         dd+=pa*(double)(i-1)/nofabc(1)+pb*(double)(j-1)/nofabc(2)+pc*(double)(k-1)/nofabc(3);
return dd;
}

void spincf::eps3d(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z,int orientation)
 {// function to plot spins in a 3d manner
  // orientation:1 ab 2 ac 3 bc projection
  //             4 ab 5 ac 6 bc side view
  int i,j,k,l;char r1,r2,r3;
  Vector a(1,2);
  Vector b(1,2),c(1,3);
  float scale,d,bbheight,bbwidth;
 
 // determine scale factor of moments
  scale=0;
  for (i=1;i<=nofa;++i)
    for (j=1;j<=nofb;++j)
     for (k=1;k<=nofc;++k)
      for(l=1;l<=nofatoms;++l)
      {c(1)=mom[in(i,j,k)](nofcomponents*(l-1)+1);
       c(2)=mom[in(i,j,k)](nofcomponents*(l-1)+2); 
       c(3)=mom[in(i,j,k)](nofcomponents*(l-1)+3);
       if ((d=Norm(c))>scale)scale=d;
      } 
  scale=0.5/(scale+0.01);


  // determine max(1,2,3) min(1,2,3) (vector in Angstroem describing a quader) for viewing magnetic unit cell
  Vector max(1,3),min(1,3),nofabc(1,3),dd(1,3),max_min(1,3),pa(1,3),pb(1,3),pc(1,3);
  Matrix p(1,3,1,3);Vector ddd(1,8),xyz(1,3),dd0(1,3);
  nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  for (i=1;i<=3;++i)
  {for(j=1;j<=3;++j) {dd(j)=nofabc(j)*r(i,j)*abc(i);p(i,j)=dd(j);}
  }
  pa=p.Column(1);  //primitive magnetic unit cell
  pb=p.Column(2);
  pc=p.Column(3);
  for (i=1;i<=3;++i)
  {ddd(1)=pa(i);
   ddd(2)=pb(i);
   ddd(3)=pc(i);
   ddd(4)=pa(i)+pb(i);
   ddd(5)=pa(i)+pc(i);
   ddd(6)=pb(i)+pc(i);
   ddd(7)=0;
   ddd(8)=pa(i)+pb(i)+pc(i);
   min(i)=Min(ddd);max(i)=Max(ddd);
  }
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
   b=xy(pa,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   b=xy(pb,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   b=xy(pc,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pa+pb+pc;a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=pa+pb;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pa+pc;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pb+pc;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pa+pb;a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=pa;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pb;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pa+pc;a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=pa;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pc;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pb+pc;a=xy(dd,orientation, min, max,bbwidth,bbheight);
   dd=pb;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=pc;b=xy(dd,orientation, min, max,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
  
  
  // plot atoms and moments in region xmin to xmax (quader)
int i1,j1,k1,i2,k2,j2,i1true,j1true,k1true;
   fprintf(fout,"/Helvetica findfont\n %i scalefont setfont\n",(int)(1000/nofa/nofb/nofc+1));
i1true=1;for (i1=0;i1true==1;++i1){i1true=0;for(i2=-1;i2<=1;i2+=2){
j1true=1;for (j1=0;j1true==1;++j1){j1true=0;for(j2=-1;j2<=1;j2+=2){
k1true=1;for (k1=0;k1true==1;++k1){k1true=0;for(k2=-1;k2<=1;k2+=2){
   dd0=pa*(double)(i2*i1)+pb*(double)(j2*j1)+pc*(double)(k2*k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {//         r1=l+'0';
         dd(1)=x[l]*abc(1);
         dd(2)=y[l]*abc(2);
         dd(3)=z[l]*abc(3);
         dd+=pa*(double)(i-1)/nofabc(1)+pb*(double)(j-1)/nofabc(2)+pc*(double)(k-1)/nofabc(3);
         dd+=dd0;
	    if(dd(1)<=max(1)+0.0001&&dd(1)>=min(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=max(2)+0.0001&&dd(2)>=min(2)-0.0001&&
            dd(3)<=max(3)+0.0001&&dd(3)>=min(3)-0.0001)
            {
             i1true=1;j1true=1;k1true=1;       

//	    a=xy(dd,orientation, min, max,bbwidth,bbheight);
//   fprintf(fout,"%g %g moveto \n (%c) show \n",a(1),a(2),r1);
//   fprintf(fout,"%g %g moveto \n (O) show \n",a(1),a(2));
              xyz(1)=dd(1)+scale*mom[in(i,j,k)](1+nofcomponents*(l-1));
              xyz(2)=dd(2)+scale*mom[in(i,j,k)](2+nofcomponents*(l-1));
              xyz(3)=dd(3)+scale*mom[in(i,j,k)](3+nofcomponents*(l-1));
	      a=xy(xyz,orientation, min, max,bbwidth,bbheight);
              xyz(1)=dd(1)-scale*mom[in(i,j,k)](1+nofcomponents*(l-1));
              xyz(2)=dd(2)-scale*mom[in(i,j,k)](2+nofcomponents*(l-1));
              xyz(3)=dd(3)-scale*mom[in(i,j,k)](3+nofcomponents*(l-1));
              b=xy(xyz,orientation, min, max,bbwidth,bbheight);
              epsarrow(fout,a,b);   
	     }
	  }
       }}}           
 }}}}}}
  
  
  
fprintf(fout,"showpage\n");
  

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
 fprintf(fout,"\n"); //new line to end spinconfiguration
}               

void spincf::printall(FILE * fout,Vector & abc,Matrix & r,float * x,float *y,float*z) //print spinconfiguration to stream
{ int i,j,k,l,lc;

 // determine primitive magnetic unit cell
  Vector nofabc(1,3),dd(1,3),pa(1,3),pb(1,3),pc(1,3);
  Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3);
  nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  for (i=1;i<=3;++i)
  {for(j=1;j<=3;++j) {dd(j)=nofabc(j)*r(i,j)*abc(i);p(i,j)=dd(j);}
  }
  pa=p.Column(1);  //primitive magnetic unit cell
  pb=p.Column(2);
  pc=p.Column(3);


 fprintf(fout,"#nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",nofa,nofb,nofc,nofatoms*nofa*nofb*nofc);
 fprintf(fout,"#[atom number] x[a] y[b] z[c] dr1[r1] dr2[r2] dr3[r3]  <Ja> <Jb> <Jc> ...\n");

   // output atoms and moments in primitive unit cell
  for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {//         r1=l+'0';
         dd(1)=x[l]*abc(1);
         dd(2)=y[l]*abc(2);
         dd(3)=z[l]*abc(3);
         dd+=pa*(double)(i-1)/nofabc(1)+pb*(double)(j-1)/nofabc(2)+pc*(double)(k-1)/nofabc(3);
         dd0=p.Inverse()*dd;dd0(1)*=nofa;dd0(2)*=nofb;dd0(3)*=nofc;
            {
              fprintf(fout,"[%i] %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
	              l,dd(1)/abc(1),dd(2)/abc(2),dd(3)/abc(3),dd0(1),dd0(2),dd0(3));
             {for (lc=1;lc<=nofcomponents;++lc)
              {fprintf(fout," %4.4f",mom[in(i,j,k)](lc+nofcomponents*(l-1)));}
              fprintf(fout,"\n");
	     }
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
  mom = new Vector[mxa*mxb*mxc+1](1,nofcomponents*nofatoms);
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
  nofa=n1;nofb=n2;nofc=n3;
   mxa=nofa+1; mxb=nofb+1; mxc=nofc+1;
  nofatoms=na;
  nofcomponents=nc;
//dimension arrays
  mom = new Vector[mxa*mxb*mxc+1](1,nofcomponents*nofatoms);
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
  mom = new Vector[mxa*mxb*mxc+1](1,nofcomponents*nofatoms);
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


