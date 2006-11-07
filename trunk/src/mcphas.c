/***********************************************************************
 *
 * mcphas.c - program to calculate the magnetic properties for external fields
 *                    of arbitrary direction
 *
 ***********************************************************************/

#include<mcphas.h>

// maximum number of spinconfigurations allowed in phasediagramm
#define MAXNOFSPINCF 100
#define MU_B 0.05788
#define K_B  0.0862

inipar ini("mcphas.ini");
int verbose=0;

#include "htcalc.h"
#include "htcalc.c"

// main program
int main (int argc, char **argv)
{ FILE * fin=NULL;
  int im,j,l;
  int nofstapoints=0;
  float x,y,dumm;
  double z,u;
  double T;
  float nn[20];nn[0]=19;
  double sta=0;
  Vector xv(0,3),xvsav(0,3);
  Vector yv(0,3),yvsav(0,3);
  Vector h(1,3);Vector mmax1(1,3);
  xvsav=0;yvsav=0;
  
  if (ini.exit_mcphas!=0)
  {ini.exit_mcphas=0;ini.print();} // if exit was 1 - save parameters and set exit=0
  
// check command line
  for (im=1;im<=argc-1;++im)
  {if (strcmp(argv[im],"-v")==0) verbose=1;     // set verbose mode on
   if (strcmp(argv[im],"-h")==0) ini.errexit(); // display help message
  }

// as class par load  parameters from file
 par inputpars("./mcphas.j");
 Vector mmax(1,inputpars.nofatoms*inputpars.nofcomponents);
  Vector h1(1,inputpars.nofcomponents);
 
//determine saturation momentum (used for scaling the plots, generation of qvectors)
T=1.0;for(l=1;l<=inputpars.nofatoms;++l){
      for (im=1;im<=inputpars.nofcomponents;++im){h1=0;h1(im)=10*MU_B*(*inputpars.jjj[l]).gJ;
                            mmax(inputpars.nofcomponents*(l-1)+im)=(*inputpars.jjj[l]).mcalc(T,h1,z,u)(im);
			   }
                                        }
for (im=1;im<=inputpars.nofcomponents&&im<=3;++im){mmax1(im)=mmax(im);}

T=0.0;h=0;
// load testspinconfigurations (nooftstspindconfigurations,init-file,sav-file)
   testspincf testspins (MAXNOFSPINCF,"./mcphas.tst","./results/mcphas.phs",inputpars.nofatoms,inputpars.nofcomponents);
   qvectors testqs (ini,inputpars.rez,mmax,"./results/mcphas.qvc",inputpars.nofatoms,inputpars.nofcomponents,verbose);

// declare variable physprop (typa class physproperties)
   physproperties physprop(ini.nofspincorrs,ini.maxnofhkls,inputpars.nofatoms,inputpars.nofcomponents);
   
// transform mmax to contain saturation moments [muB] 
for(l=1;l<=inputpars.nofatoms;++l){for (im=1;im<=inputpars.nofcomponents&&im<=3;++im){mmax(3*(l-1)+im)*=inputpars.gJ(l);}}

if (argc>1&&strncmp(argv[argc-1],"-",1)!=0){ini.xv=0;ini.yv=0;fin=fopen_errchk (argv[argc-1],"rb");}   //input from file
// loop different H /T points in phase diagram
for (x=ini.xmin;x<=ini.xmax;x+=ini.xstep)
 { //begin initialize display file
   FILE * fout;fout = fopen_errchk ("./results/.mcphas.fum","w");
   fprintf (fout, " %4.4g %6.6g  %4.4g %4.4g %4.4g %4.4g %4.4g %8.8g %8.8g  %4.4g %4.4g %4.4g %4.4g\n",
            0.0,ini.ymin,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-Max(mmax1),0.0,0.0);
   fprintf (fout, " %4.4g %6.6g  %4.4g %4.4g %4.4g %4.4g %4.4g %8.8g %8.8g  %4.4g %4.4g %4.4g %4.4g\n",
            0.0,ini.ymax+1e-4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,Max(mmax1),0.0,0.0);
   fclose(fout); //end initialize display file
  
  for (y=ini.ymin;y<=ini.ymax;y+=ini.ystep)
  {//correct some of input errors in parameters ini (if there are any)
   if (ini.xstep<0){ini.xstep=-ini.xstep;}
   if (ini.ystep<0){ini.ystep=-ini.ystep;}
   if (ini.xmin>ini.xmax){dumm=ini.xmin;ini.xmin=ini.xmax;ini.xmax=dumm;}
   if (ini.ymin>ini.ymax){dumm=ini.ymin;ini.ymin=ini.ymax;ini.ymax=dumm;}
   
   if(argc>1&&strncmp(argv[argc-1],"-",1)!=0)  //should T-H values be read from file ?
   {while (feof(fin)==0&&0==inputline(fin,nn));  // if yes -> input them
    if (feof(fin)!=0) goto endproper;
    x=nn[1];y=nn[2];T=nn[3];h(1)=nn[5];h(2)=nn[6];h(3)=nn[7];}
   else
   {//if parameters outside specified region then put them into it ...
    xv=ini.xv;yv=ini.yv; xv/=Norm(xv);yv/=Norm(yv);
    if (x<ini.xmin) x=ini.xmin;
    if (y<ini.ymin) y=ini.ymin;    

    //if xv or yv has changed, then set new ini.zero to current vector - dos not make sense !
//    if(xv!=xvsav||yv!=yvsav)
//    {yvsav=yv;xvsav=xv;ini.zero(0)=T;ini.zero(1)=h(1);ini.zero(2)=h(2);ini.zero(3)=h(3);
//    } does not make sense !

     T=ini.zero(0)+x*xv(0)+y*yv(0);
     h(1)=ini.zero(1)+x*xv(1)+y*yv(1);
     h(2)=ini.zero(2)+x*xv(2)+y*yv(2);
     h(3)=ini.zero(3)+x*xv(3)+y*yv(3);
   } 
     
      physprop.x=x;physprop.y=y;
      physprop.T=T;
      physprop.H=h;

//calculate physical properties at HT- point
      j=htcalc(physprop.H,T,inputpars,testqs,testspins,physprop);
       switch (j)
       {case 0:
            //save physical properties of HT-point
	    sta=(sta*nofstapoints+physprop.save (verbose,j,inputpars))/(nofstapoints+1);
	    ++nofstapoints;
	      break; 
	 case 1: goto endproper;
	      break;
         case 2: //ht calculation leads to no results- save dummy line
	         physprop.save (verbose,j,inputpars);
		 sta+=1.0; // increment sta because within manifold of spincf no good solution could be found
	      break;	 
	 default:  ;
	}
 if(fin!=NULL){y=ini.ymin-ini.ystep;} // this is to switch off xy loop if xy points are read from file

    }
 }  
endproper:
  testspins.save();testqs.save();
   if(argc>1&&strncmp(argv[argc-1],"-",1)!=0) fclose(fin);
fprintf(stdout,"sta=%g\n",sta);

#ifdef linux 
return 0;
#else
exit(0);
#endif

}
