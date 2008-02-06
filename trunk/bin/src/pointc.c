/***********************************************************************
 *
 * point.c   - new program substituting basic program to calculate pointcharge
 *             model
 *
 ***********************************************************************/

#define MU_B  5.788378E-02 // Bohrmagneton in meV/tesla
#define PI   3.141592654
#define MAXNOFCHARINLINE 1000

#include "ionpars.hpp"
#include "martin.h"
#include<cstdio>
#include<cerrno>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector.h>

/**********************************************************************/
// main program
int main (int argc, char **argv)
{// check command line
  if (argc < 3)
    { printf ("\nProgram to calculate Crystalfield Parameters from Point Charges \n\
                Usage: pointc Ce3+ 0.2 4 1 5.3\n\
                 ... meaning calculate Blms for one pointcharge \n\
                 of +0.2|e| in distance\n\
                 x=4 A y=1 A z=5.3 A from a Ce3+ ion.\n\
                Alternative Usage: pointc Ce3+ filename\n\
                 ... meaning read several charges+coordinates from file,\n\
                 file format: column 1=charge, column 2-4 = x y z coordinate.\n\
                results are written to stdout\n\
		\n");
      exit (1);
    }

FILE * table_file;
char instr[MAXNOFCHARINLINE];
int n=0;
float invalues[100];invalues[0]=99;
  double q,x,y,z;
  printf ("#!cfield\n#<!--mcphase.sipf-->\n");
  printf ("IONTYPE=%s\n",argv[1]);

// read create class object ionpars from iontype
 ionpars * iops;
 iops=new ionpars(argv[1]);  
 Vector tetan(1,6), rl(1,6);
 double gJ,J;
// returns stevens parameters and landefactor of ion
  tetan(2)=(*iops).alpha;tetan(4)=(*iops).beta;tetan(6)=(*iops).gamma;
  gJ=(*iops).gJ; 
  J=(*iops).J;
  rl(2)=(*iops).r2;  rl(4)=(*iops).r4;  rl(6)=(*iops).r6;

  printf ("#J=%4g\n",J);
  printf("#Lande Factor: gJ = %4g\n",gJ);
  printf("#Stevens factors: alpha beta gamma = %4g %4g %4g \n",tetan(2),tetan(4),tetan(6));
  printf ("#Expextation values of radial wave function:\n");
  printf ("#<r^2>=%4g a0^2  <r^4>=%4g a0^4  <r^6>=%4g a0^6    a0=0.5292 Angstroem\n#\n",rl(2),rl(4),rl(6));

if (argc<5)
{table_file=fopen_errchk(argv[2],"r");
 while(n==0&&feof(table_file)==false)n=inputline(table_file, invalues);
  q=invalues[1];
  x=invalues[2];
  y=invalues[3];
  z=invalues[4]; 
} else 
{ n=4;
  q=strtod(argv[2],NULL);
  x=strtod(argv[3],NULL);
  y=strtod(argv[4],NULL);
  z=strtod(argv[5],NULL);
}

  printf ("#pointcharges charge[|e|]  x[A] y[A] z[A]\n",q,x,y,z);
while(n>0)
{

  printf ("pointcharge= %4g         %4g %4g %4g\n",q,x,y,z);
  Vector B(1,45); B=0;
// calculate Blm's
double r,ct,ct2,st,st2,sfi,cfi;
r = sqrt(x * x + y * y + z * z);
ct = z/r;                 //z
ct2 = ct * ct;      
st = sqrt(x*x+y*y)/r;
st2 = st * st;
if((x*x+y*y)==0){sfi=0;cfi=1;}
else
{sfi =  y/sqrt(x*x+y*y);
cfi =  x/sqrt(x*x+y*y);}

 int l,m;   
// cnst is the Zlm constants - put them into the matrix
Matrix cnst(0,6,-6,6);

cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;
cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;
cnst(6,0)=  0.06357014;
cnst(6,1)=  1.032669;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}

// this is squaring of the coefficients of Zlm, a technical trick in
// order to save a multiplication later
for(l=2;l<=6;l+=2){for(m=-l;m<=l;++m)cnst(l,m)*=cnst(l,m);}

//ro = a(0, 0) / sqrt(4.0 * 3.1415);
B(1)= cnst(2, -2)  * 2 * st2 * sfi * cfi;
B(2)= cnst(2, -1)  * st * sfi * ct;
B(3)= cnst(2, 0)  * (3 * ct2 - 1);
B(4)= cnst(2, 1)  * st * cfi * ct;
B(5)= cnst(2, 2)  * st2 * (cfi * cfi - sfi * sfi);

B(13)= cnst(4, -4) * st2 * st2 * 4 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
B(14)= cnst(4, -3) * ct * st * st2 * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
B(15)= cnst(4, -2) * (7 * ct2 - 1) * 2 * st2 * cfi * sfi;
B(16)= cnst(4, -1) * st * sfi * ct * (7 * ct2 - 3);
B(17)= cnst(4, 0) * (35 * ct2 * ct2 - 30 * ct2 + 3);
B(18)= cnst(4, 1)  * st * cfi * ct * (7 * ct2 - 3);
B(19)= cnst(4, 2)  * (7 * ct2 - 1) * st2 * (cfi * cfi - sfi * sfi);
B(20)= cnst(4, 3)  * ct * st * st2 * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
B(21)= cnst(4, 4)  * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);

B(33)= cnst(6, -6) * st2 * st2 * st2 * (6 * cfi * cfi * cfi * cfi * cfi * sfi - 20 * cfi * cfi * cfi * sfi * sfi * sfi + 6 * cfi * sfi * sfi * sfi * sfi * sfi);
B(34)= cnst(6, -5) * ct * st * st2 * st2 * (5 * cfi * cfi * cfi * cfi * sfi - 10 * cfi * cfi * sfi * sfi * sfi + sfi * sfi * sfi * sfi * sfi);
B(35)= cnst(6, -4) * (11 * ct2 - 1) * 4 * st2 * st2 * (cfi * cfi * cfi * sfi - cfi * sfi * sfi * sfi);
B(36)= cnst(6, -3) * (11 * ct * ct2 - 3 * ct) * st2 * st * (3 * cfi * cfi * sfi - sfi * sfi * sfi);
B(37)= cnst(6, -2) * 2 * st2 * sfi * cfi * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2);
B(38)= cnst(6, -1) * ct * st * sfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
B(39)= cnst(6, 0)  * (231 * ct2 * ct2 * ct2 - 315 * ct2 * ct2 + 105 * ct2 - 5);
B(40)= cnst(6, 1)  * ct * st * cfi * (33 * ct2 * ct2 - 30 * ct2 + 5);
B(41)= cnst(6, 2)  * (16 * ct2 * ct2 - 16 * ct2 * st2 + st2 * st2) * st2 * (cfi * cfi - sfi * sfi);
B(42)= cnst(6, 3)  * (11 * ct * ct2 - 3 * ct) * st2 * st * (cfi * cfi * cfi - 3 * cfi * sfi * sfi);
B(43)= cnst(6, 4)  * (11 * ct2 - 1) * st2 * st2 * (cfi * cfi * cfi * cfi - 6 * cfi * cfi * sfi * sfi + sfi * sfi * sfi * sfi);
B(44)= cnst(6, 5)  * ct * st * st2 * st2 * (cfi * cfi * cfi * cfi * cfi - 10 * cfi * cfi * cfi * sfi * sfi + 5 * cfi * sfi * sfi * sfi * sfi);
B(45)= cnst(6, 6) * st2 * st2 * st2 * (cfi * cfi * cfi * cfi * cfi * cfi - 15 * cfi * cfi * cfi * cfi * sfi * sfi + 15 * cfi * cfi * sfi * sfi * sfi * sfi - sfi * sfi * sfi * sfi * sfi * sfi);


// now calculation of the coefficients gammaLM  in cgs
int i;
//gamma2M
for (i=1;i<=5;++i)B(i)*=q/r/r/r*4*PI/5; 
//gamma4M
for (i=13;i<=21;++i)B(i)*=q/r/r/r/r/r*4*PI/9; 
//gamma6M
for (i=33;i<=45;++i)B(i)*=q/r/r/r/r/r/r/r*4*PI/13; 

double e,a0,umr,ehv2,ehv4,ehv6;
     e = 4.80325E-10; // elementarladung
// einheit von r in <r^n> ist Bohrradius^n = a0^n in angstroem^n
     a0 = .5292;//(Angstroem)
//   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in mJ = 10^4
//   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in THz =1.509166084e22
//   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in meV =0.624146e23
//   REM umr von (esu)^2*a0^n/Angstroem^(n+1) in K =0.72429024e24
     umr = 6.24146E+22;
     ehv2 = a0*a0 * umr;
     ehv4 = a0*a0*a0*a0 * umr;
     ehv6 = a0*a0*a0*a0*a0*a0 * umr;


// now calculation of the B_LM  in meV
//gamma2M
for (i=1;i<=5;++i)(*iops).Blm(i)+=-B(i)*e*e*rl(2)*tetan(2)*ehv2; 
//gamma4M
for (i=13;i<=21;++i)(*iops).Blm(i)+=-B(i)*e*e*rl(4)*tetan(4)*ehv4; 
//gamma6M
for (i=33;i<=45;++i)(*iops).Blm(i)+=-B(i)*e*e*rl(6)*tetan(6)*ehv6;
n=0;
if (argc<5)
 { while(n==0&feof(table_file)==false)n=inputline(table_file, invalues);
  q=invalues[1];
  x=invalues[2];
  y=invalues[3];
  z=invalues[4]; 
 }
}

if (argc<5){fclose(table_file);}

(*iops).savBlm(stdout);
}

/*
REM q     ...........parameter - reduced charges [|e|]
'   xNN,yNN,zNN ..... position of charge [A]
REM******berechnung der b#() aus den reduced charges q..*******

sv20 = 0: sv22 = 0: sv40 = 0: sv42 = 0: sv43 = 0: sv44 = 0: sv60 = 0: sv62 = 0: sv63 = 0: sv64 = 0: sv66 = 0

REM berechnung der Entwicklungskoeffizienten svlm der Zlm
REM svlmadd
sv20 = sv20 + FNV20(z, r) / r ^ 3 * 4 * pi# / 5 * q     '= gamma20 in cgs
sv22 = sv22 + FNV22(x, y, z) / r ^ 3 * 4 * pi# / 5 * q
sv40 = sv40 + FNV40(z, r) / r ^ 5 * 4 * pi# / 9 * q
sv42 = sv42 + FNV42(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
sv43 = sv43 + FNV43(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
sv44 = sv44 + FNV44(x, y, z) / r ^ 5 * 4 * pi# / 9 * q
sv60 = sv60 + FNV60(z, r) / r ^ 7 * 4 * pi# / 13 * q
sv62 = sv62 + FNV62(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
sv63 = sv63 + FNV63(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
sv64 = sv64 + FNV64(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
sv66 = sv66 + FNV66(x, y, z) / r ^ 7 * 4 * pi# / 13 * q
// sub for calculation of charge density given a radiu R and polar angles teta, fi and expansion coeff. alm
// some functions 
DEF FNV20 (z, r) = .25 * SQR(5 / pi#) * (3 * z * z - r * r) / r / r
  DEF FNV22 (x, y, z) = 1 / 4 * SQR(15 / pi#) * (x * x - y * y) / (x * x + y * y + z * z)
  DEF FNV40 (z, r) = 3 / 16 * SQR(1 / pi#) * (35 * z * z * z * z - 30 * z * z * r * r + 3 * r * r * r * r) / r / r / r / r
  DEF FNV42 (x, y, z) = 3 / 8 * SQR(5 / pi#) * (7 * z * z - (x * x + y * y + z * z)) * (x * x - y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV43 (x, y, z) = 3 / 8 * SQR(70 / pi#) * z * x * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV44 (x, y, z) = 3 / 16 * SQR(35 / pi#) * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV60 (z, r) = 1 / 32 * SQR(13 / pi#) * (231 * z * z * z * z * z * z - 315 * z * z * z * z * r * r + 105 * z * z * r * r * r * r - 5 * r * r * r * r * r * r) / r / r / r / r / r / r
  DEF FNV62 (x, y, z) = 1 / 64 * SQR(2730 / pi#) * (16 * z * z * z * z - 16 * (x * x + y * y) * z * z + (x * x + y * y) * (x * x + y * y)) * (x * x - y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV63 (x, y, z) = 1 / 32 * SQR(2730 / pi#) * z * x * (11 * z * z - 3 * (x * x + y * y + z * z)) * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV64 (x, y, z) = 21 / 32 * SQR(13 / 7 / pi#) * (11 * z * z - x * x - y * y - z * z) * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV66 (x, y, z) = 231 / 64 * SQR(26 / 231 / pi#) * (x * x * x * x * x * x - 15 * x * x * x * x * y * y + 15 * x * x * y * y * y * y - y * y * y * y * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)
  DEF FNV66S (x, y, z) = 231 / 32 * SQR(26 / 231 / pi#) * y * x * (3 * x * x - y * y) * (x * x - 3 * y * y) / (x * x + y * y + z * z) / (x * x + y * y + z * z) / (x * x + y * y + z * z)


REM umrechnung in Koeffizienten von Olm
     REM <r^n>-werte in bohrradius^n

IF site% = 1 THEN r2 = 1.2: r4 = 3.455: r6 = 21.226: REM ce3+
PRINT #1, "<r^2>="; r2; " a0^2  <r^4>="; r4; " a0^4  <r^6>="; r6; " a0^6    a0=0.5292 Angstroem}"

     REM stevensfaktoren
IF site% = 1 THEN alpha = -2 / (5 * 7): beta = 2 / 9 / 5 / 7: gamma = 0: REM ce3+


  b#(2, 0) = -e ^ 2 * r2 * alpha * sv20 * SQR(5 / 16 / pi#) * ehv2
  b#(2, 2) = -e ^ 2 * r2 * alpha * sv22 * SQR(15 / 16 / pi#) * ehv2

  b#(4, 0) = -e ^ 2 * r4 * beta * sv40 * 3 / 16 * SQR(1 / pi#) * ehv4
  b#(4, 2) = -e ^ 2 * r4 * beta * sv42 * 3 / 8 * SQR(5 / pi#) * ehv4
  b#(4, 3) = -e ^ 2 * r4 * beta * sv43 * 3 / 8 * SQR(70 / pi#) * ehv4
  b#(4, 4) = -e ^ 2 * r4 * beta * sv44 * 3 / 16 * SQR(35 / pi#) * ehv4

  b#(6, 0) = -e ^ 2 * r6 * gamma * sv60 / 32 * SQR(13 / pi#) * ehv6
  b#(6, 2) = -e ^ 2 * r6 * gamma * sv62 / 64 * SQR(2730 / pi#) * ehv6
  b#(6, 3) = -e ^ 2 * r6 * gamma * sv63 / 32 * SQR(2730 / pi#) * ehv6
  b#(6, 4) = -e ^ 2 * r6 * gamma * sv64 * 21 / 32 * SQR(13 / 7 / pi#) * ehv6
  b#(6, 6) = -e ^ 2 * r6 * gamma * sv66 * 231 / 64 * SQR(26 / 231 / pi#) * ehv6
*/
