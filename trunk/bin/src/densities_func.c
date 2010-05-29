/*****************************************************************/
// functions for programs charges, spindensity, orbmomdensity momdensity
/*****************************************************************/

void check_for_best(FILE *fin_coq,double Tin, double hain,double hbin, double hcin, spincf & savmf, double & T,double & ha, double& hb,double&hc,char*outstr)
{// load mfconfigurations and check which one is nearest -------------------------------

int n;
   double ddT,ddHa,ddHb,ddHc,dd,delta;
 float numbers[13];numbers[9]=1;numbers[10]=3;
 numbers[0]=13;

 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      spins.load(fin_coq);
      ddT=Tin-numbers[3];ddT*=ddT;
      ddHa=hain-numbers[5];ddHa*=ddHa;
      ddHb=hbin-numbers[6];ddHb*=ddHb;
      ddHc=hcin-numbers[7];ddHc*=ddHc;
      dd=sqrt(ddT+ddHa+ddHb+ddHc+0.000001);
      if (dd<delta)
       {delta=dd;
        sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g n=%g spins nofatoms=%i in primitive basis nofcomponents=%i",numbers[3],numbers[5],numbers[6],numbers[7],numbers[8],(int)numbers[9],(int)numbers[10]);
        savmf=spins;T=numbers[3];ha=numbers[5];hb=numbers[6];hc=numbers[7];
       }
    }
}

int headerinput(FILE * fin_coq,FILE* fout,graphic_parameters & gp,Vector & abc,Matrix & r,float * x, float * y, float*z,float*gJ,char**cffilenames,int & nofatoms,int & nofcomponents)
{ char instr[MAXNOFCHARINLINE];
 long int pos=0,j;int n=0;
nofatoms=0;nofcomponents=3;
char *token;abc=0;
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring'
  { pos=ftell(fin_coq);
   if (pos==-1)
       {fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq);
   // strip /r (dos line feed) from line if necessary
    while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

   if (instr[strspn(instr," \t")]=='#'){fprintf(fout,"%s",instr);}
   if(abc[1]==0){extract(instr,"a",abc[1]);extract(instr,"b",abc[2]); extract(instr,"c",abc[3]);
                 extract(instr,"alpha",abc[4]);  extract(instr,"beta",abc[5]);extract(instr,"gamma",abc[6]);
   }
   extract(instr,"show_abc_unitcell",gp.show_abc_unitcell);
   extract(instr,"show_primitive_crystal_unitcell",gp.show_primitive_crystal_unitcell);
   extract(instr,"show_magnetic_unitcell",gp.show_magnetic_unitcell);
   extract(instr,"show_atoms",gp.show_atoms);
   extract(instr,"spins_scale_moment",gp.spins_scale_moment);
   extract(instr,"show_chargedensity",gp.show_chargedensity);

   extract(instr,"scale_view_1",gp.scale_view_1);
   extract(instr,"scale_view_2",gp.scale_view_2);
   extract(instr,"scale_view_3",gp.scale_view_3);

   extract(instr,"r1x",r[1][1]);extract(instr,"r2x",r[1][2]); extract(instr,"r3x",r[1][3]);
   extract(instr,"r1y",r[2][1]); extract(instr,"r2y",r[2][2]); extract(instr,"r3y",r[2][3]);
   extract(instr,"r1z",r[3][1]); extract(instr,"r2z",r[3][2]); extract(instr,"r3z",r[3][3]);
   extract(instr,"r1a",r[1][1]);extract(instr,"r2a",r[1][2]); extract(instr,"r3a",r[1][3]);
   extract(instr,"r1b",r[2][1]); extract(instr,"r2b",r[2][2]); extract(instr,"r3b",r[2][3]);
   extract(instr,"r1c",r[3][1]); extract(instr,"r2c",r[3][2]); extract(instr,"r3c",r[3][3]);
   extract(instr,"nofatoms",nofatoms);    extract(instr,"nofcomponents",nofcomponents);
   if (nofatoms>0&&(extract(instr,"x",x[n+1])+
                   extract(instr,"y",y[n+1])+
  		       extract(instr,"z",z[n+1])==0)||
		       (extract(instr,"da",x[n+1])+
                   extract(instr,"db",y[n+1])+
		       extract(instr,"dc",z[n+1])==0))
		  {++n;if(n>nofatoms||nofatoms>MAXNOFATOMS)
                    {fprintf(stderr,"ERROR charges.c reading file:maximum number of atoms in unit cell exceeded\n");exit(EXIT_FAILURE);}
                   cffilenames[n]=new char[MAXNOFCHARINLINE];
                   extract(instr,"cffilename",cffilenames[n],(size_t)MAXNOFCHARINLINE);
                   extract(instr,"gJ",gJ[n]);
//		   printf("%s\n",cffilenames[n]);
                  }
  }
    j=fseek(fin_coq,pos,SEEK_SET);
    if (j!=0){fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
return n;
}
