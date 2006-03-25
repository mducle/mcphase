/*-----------------------------------------------------------------------------
 
                                M  O  D  U  L
 
                                   THETA C
 
-------------------------------------------------------------------------------
 
Aufgabe               : alpha_J , beta_J , gamma_J  und
 
                          n
                        <r > definieren  ( n = 2 , 4 , 6 )
 
-------------------------------------------------------------------------------
 
Literatur zur Berechnung der alpha_J , beta_J , gamma_J :
---------------------------------------------------------
[0] Diplomarbeit Peter Hoffmann : Dort wird gezeigt werden ,dass
                                  die alpha_J,beta_J und gamma_J
                                  fuer die allgemeinen Operatoren
                                  die Gleichen sind ,wie fuer die
                                  Stevensoperatoren ,welche in [1]
                                  tabelliert sind.
 
[1] Abragam und Bleaney "Electronic Paramagnetic Resonance of Transition Ions",
                        1970 ,Table 20
 
                               2      4      6
Literatur zur Berechnung der <r > , <r > , <r > :
-------------------------------------------------
[2] A.J. Freeman and J.P. Desclaux: Journal of Magnetism and Magnetic Materials
                                    12 (1979) 11
 
 
-------------------------------------------------------------------------------
 
Bemerkungen           :
 
-------------------------------------------------------------------------------
 
Definierte Funktionen :
 
-----------------------
 
-----------------------------------------------------------------------------*/
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdio.h>          /* damit FILE definiert wird               */
#include <stdlib.h>
#include <math.h>           /* damit sqrt in define_j.c definiert wird */
#define pi (4.0*atan(1.0))  /* atan() braucht <math.h>                 */
#include "types.c"          /* benutze Datentypen laden                */
/*----------------------------------------------------------------------------
                                define
-----------------------------------------------------------------------------*/
 
#define NR_RE 15
 
/*----------------------------------------------------------------------------
Extern definierte Funktionen
-----------------------------------------------------------------------------*/
 
/*----------------------------------------------------------------------------
                                alpha_J
-----------------------------------------------------------------------------*/
 
DOUBLE alpha_J[NR_RE]={
/* 4f_0        */    1.0 * 0                ,
/* 4f_1 : Ce3+ */   -1.0 * 2/5/7            ,
/* 4f_2 : Pr3+ */   -1.0 * 2*2*13/3/3/5/5/11,
/* 4f_3 : Nd3+ */   -1.0 * 7/3/3/11/11      ,
/* 4f_4 : Pm3+ */    1.0 * 2*7/3/5/11/11    ,
/* 4f_5 : Sm3+ */    1.0 * 13/3/3/5/7       ,
/* 4f_6 : Eu3+ */    1.0 * 0                ,
/* 4f_7 : Gd3+ */    1.0 * 0                ,
/* 4f_8 : Tb3+ */   -1.0 * 1/3/3/11         ,
/* 4f_9 : Dy3+ */   -1.0 * 2/3/3/5/7        ,
/* 4f_10: Ho3+ */   -1.0 * 1/2/3/3/5/5      ,
/* 4f_11: Er3+ */    1.0 * 2*2/3/3/5/5/7    ,
/* 4f_12: Tm3+ */    1.0 * 1/3/3/11         ,
/* 4f_13: Yb3+ */    1.0 * 2/3/3/7          ,
/* 4f_14       */    1.0 * 0
};
/*----------------------------------------------------------------------------
                                 beta_J
-----------------------------------------------------------------------------*/
 
DOUBLE beta_J[NR_RE]={
/* 4f_0        */    1.0 * 0                              ,
/* 4f_1 : Ce3+ */    1.0 * 2/3/3/5/7                      ,
/* 4f_2 : Pr3+ */   -1.0 * 2*2/3/3/5/11/11                ,
/* 4f_3 : Nd3+ */   -1.0 * 2*2*2*17/3/3/3/11/11/11/13     ,
/* 4f_4 : Pm3+ */    1.0 * 2*2*2*7*17/3/3/3/5/11/11/11/13 ,
/* 4f_5 : Sm3+ */    1.0 * 2*13/3/3/3/5/7/11              ,
/* 4f_6 : Eu3+ */    1.0 * 0                              ,
/* 4f_7 : Gd3+ */    1.0 * 0                              ,
/* 4f_8 : Tb3+ */    1.0 * 2/3/3/3/5/11/11                ,
/* 4f_9 : Dy3+ */   -1.0 * 2*2*2/3/3/3/5/7/11/13          ,
/* 4f_10: Ho3+ */   -1.0 * 1/2/3/5/7/11/13                ,
/* 4f_11: Er3+ */    1.0 * 2/3/3/5/7/11/13                ,
/* 4f_12: Tm3+ */    1.0 * 2*2*2/3/3/3/3/5/11/11          ,
/* 4f_13: Yb3+ */   -1.0 * 2/3/5/7/11                     ,
/* 4f_14       */    1.0 * 0
};
/*----------------------------------------------------------------------------
                                gamma_J
-----------------------------------------------------------------------------*/
 
DOUBLE gamma_J[NR_RE]={
/* 4f_0        */    1.0 * 0                                 ,
/* 4f_1 : Ce3+ */    1.0 * 0                                 ,
/* 4f_2 : Pr3+ */    1.0 * 2*2*2*2*17/3/3/3/3/5/7/11/11/13   ,
/* 4f_3 : Nd3+ */   -1.0 * 5*17*19/3/3/3/7/11/11/11/13/13    ,
/* 4f_4 : Pm3+ */    1.0 * 2*2*2*17*19/3/3/3/7/11/11/11/13/13,
/* 4f_5 : Sm3+ */    1.0 * 0                                 ,
/* 4f_6 : Eu3+ */    1.0 * 0                                 ,
/* 4f_7 : Gd3+ */    1.0 * 0                                 ,
/* 4f_8 : Tb3+ */   -1.0 * 1/3/3/3/3/7/11/11/13              ,
/* 4f_9 : Dy3+ */    1.0 * 2*2/3/3/3/7/11/11/13/13           ,
/* 4f_10: Ho3+ */   -1.0 * 5/3/3/3/7/11/11/13/13             ,
/* 4f_11: Er3+ */    1.0 * 2*2*2/3/3/3/7/11/11/13/13         ,
/* 4f_12: Tm3+ */   -1.0 * 5/3/3/3/3/7/11/11/13              ,
/* 4f_13: Yb3+ */    1.0 * 2*2/3/3/3/7/11/13                 ,
/* 4f_14       */    1.0 * 0
};
/*----------------------------------------------------------------------------
                               info_thetakq()
-----------------------------------------------------------------------------*/
info_thetakq()     /* Liste der Stevensfaktoren alpha_J,beta_J,gamma_J zeigen*/
{
    CHAR *name = "thetakq.info";
 
    FILE *fopen(),*fp;
    CHAR *s,*text=" %2s     %10.4f     %10.4f      %10.4f\n";
    DOUBLE a,b,c;
    DOUBLE x=100.0;
    DOUBLE y=10000.0;
    DOUBLE z=1000000.0;
 
    fp = fopen(name,"w");
    clearscreen;
    printf("Information contained in the File %s .\n",name);
 
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"|Table of the Stevens factors alpha_J ,beta_J und gamma_J.|\n");
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"In the disseration of  Peter Hoffmann it is shown that     \n");
    fprintf(fp,"                                    2S+1                   \n");
    fprintf(fp,"the quantisation constant  theta   (     L   ) (k=2,4,6)   \n");
    fprintf(fp,"                                 k        J                \n");
    fprintf(fp,"of the general Stevensoperatoren with alpha_J , beta_J     \n");
    fprintf(fp,"und gamma_J can be identified with the usual               \n");
    fprintf(fp,"Stevensoperatoren.                                         \n");
    fprintf(fp,"                                                           \n");
    fprintf(fp,"The Stevens factors are rational numbers , in which        \n");
    fprintf(fp,"Abragam und Bleaney   \"Electronic Paramagnetic Resonance  \n");
    fprintf(fp,"                     of Transition Ions\",1970,Table  20   \n");
    fprintf(fp,"calculated exactly ,but here on 4 decimal places           \n");
    fprintf(fp,"are shown .                                                \n");
    fprintf(fp,"\n\n");
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"|                  2             4                6       |\n");
    fprintf(fp,"|RE3+    alpha_J*10     beta_J*10       gamma_J*10        |\n");
    fprintf(fp,"===========================================================\n");
s="Ce";a=x*alpha_J[ 1];b=y*beta_J[ 1];c=z*gamma_J[ 1];fprintf(fp,text,s,a,b,c);
s="Pr";a=x*alpha_J[ 2];b=y*beta_J[ 2];c=z*gamma_J[ 2];fprintf(fp,text,s,a,b,c);
s="Nd";a=x*alpha_J[ 3];b=y*beta_J[ 3];c=z*gamma_J[ 3];fprintf(fp,text,s,a,b,c);
s="Pm";a=x*alpha_J[ 4];b=y*beta_J[ 4];c=z*gamma_J[ 4];fprintf(fp,text,s,a,b,c);
s="Sm";a=x*alpha_J[ 5];b=y*beta_J[ 5];c=z*gamma_J[ 5];fprintf(fp,text,s,a,b,c);
s="Eu";a=x*alpha_J[ 6];b=y*beta_J[ 6];c=z*gamma_J[ 6];fprintf(fp,text,s,a,b,c);
s="Gd";a=x*alpha_J[ 7];b=y*beta_J[ 7];c=z*gamma_J[ 7];fprintf(fp,text,s,a,b,c);
s="Tb";a=x*alpha_J[ 8];b=y*beta_J[ 8];c=z*gamma_J[ 8];fprintf(fp,text,s,a,b,c);
s="Dy";a=x*alpha_J[ 9];b=y*beta_J[ 9];c=z*gamma_J[ 9];fprintf(fp,text,s,a,b,c);
s="Ho";a=x*alpha_J[10];b=y*beta_J[10];c=z*gamma_J[10];fprintf(fp,text,s,a,b,c);
s="Er";a=x*alpha_J[11];b=y*beta_J[11];c=z*gamma_J[11];fprintf(fp,text,s,a,b,c);
s="Tm";a=x*alpha_J[12];b=y*beta_J[12];c=z*gamma_J[12];fprintf(fp,text,s,a,b,c);
s="Yb";a=x*alpha_J[13];b=y*beta_J[13];c=z*gamma_J[13];fprintf(fp,text,s,a,b,c);
 
    fclose(fp);
}
/*----------------------------------------------------------------------------
                                      2
                                    <r >
-----------------------------------------------------------------------------*/
 
DOUBLE r2[NR_RE]={
/* 4f_0        */    0.0                ,
/* 4f_1 : Ce3+ */   1.309               ,
/* 4f_2 : Pr3+ */   1.1963              , /* von U.Walter  Diss.   */
/* 4f_3 : Nd3+ */   1.114               ,
/* 4f_4 : Pm3+ */   1.0353              , /*        -"-            */
/* 4f_5 : Sm3+ */   0.9743              ,
/* 4f_6 : Eu3+ */   0.9175              ,
/* 4f_7 : Gd3+ */   0.8671              ,
/* 4f_8 : Tb3+ */   0.8220              ,
/* 4f_9 : Dy3+ */   0.7814              ,
/* 4f_10: Ho3+ */   0.7446              ,
/* 4f_11: Er3+ */   0.7111              ,
/* 4f_12: Tm3+ */   0.6804              ,
/* 4f_13: Yb3+ */   0.6522              ,
/* 4f_14       */   0.0
};
/*----------------------------------------------------------------------------
                                      4
                                    <r >
-----------------------------------------------------------------------------*/
 
DOUBLE r4[NR_RE]={
/* 4f_0        */    0.0                ,
/* 4f_1 : Ce3+ */   3.964               ,
/* 4f_2 : Pr3+ */   3.3335              , /* von U. Walter Diss.     */
/* 4f_3 : Nd3+ */   2.910               ,
/* 4f_4 : Pm3+ */   2.5390              , /*        -"-              */
/* 4f_5 : Sm3+ */   2.260               ,
/* 4f_6 : Eu3+ */   2.020               ,
/* 4f_7 : Gd3+ */   1.820               ,
/* 4f_8 : Tb3+ */   1.651               ,
/* 4f_9 : Dy3+ */   1.505               ,
/* 4f_10: Ho3+ */   1.379               ,
/* 4f_11: Er3+ */   1.270               ,
/* 4f_12: Tm3+ */   1.174               ,
/* 4f_13: Yb3+ */   1.089               ,
/* 4f_14       */    0.0
};
/*----------------------------------------------------------------------------
                                      6
                                    <r >
-----------------------------------------------------------------------------*/
 
DOUBLE r6[NR_RE]={
/* 4f_0        */    0.0                ,
/* 4f_1 : Ce3+ */  23.31                ,
/* 4f_2 : Pr3+ */  18.353               , /* U. Walter Diss.         */
/* 4f_3 : Nd3+ */  15.03                ,
/* 4f_4 : Pm3+ */  12.546               , /*          -"-            */
/* 4f_5 : Sm3+ */  10.55                ,
/* 4f_6 : Eu3+ */   9.039               ,
/* 4f_7 : Gd3+ */   7.831               ,
/* 4f_8 : Tb3+ */   6.852               ,
/* 4f_9 : Dy3+ */   6.048               ,
/* 4f_10: Ho3+ */   5.379               ,
/* 4f_11: Er3+ */   4.816               ,
/* 4f_12: Tm3+ */   4.340               ,
/* 4f_13: Yb3+ */   3.932               ,
/* 4f_14       */    0.0
};
/*----------------------------------------------------------------------------
                                  info_rn()
-----------------------------------------------------------------------------*/
                   /*             n                   */
info_rn()          /* Liste der <r > n=2,4,6 ausgeben */
{
    CHAR *name = "rn.info";
 
    FILE   *fopen(),*fp;
    DOUBLE a,b,c,a2,a4,a6;
 
    CHAR *s;
    CHAR *texti =" %2s       %10.4f    %10.4f    %10.4f nach U. Walter\n";
    CHAR *textfw=" %2s      %10.3f    %10.3f    %10.3f              \n";
    CHAR *sexti =" %2s      %10.4f    %10.4f    %10.4f nach U. Walter\n";
    CHAR *sextfw=" %2s      %10.4f    %10.4f    %10.4f              \n";
 
    fp = fopen(name,"w");
    clearscreen;
    printf("Information contained in the file %s.\n",name);
 
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|             n                                           |\n");
    fprintf(fp,"|Table of <r > for n = 2 , 4 , 6 : a0=Bohr Radius         |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"That over the wavefunction  R (r) of the 4f-Elektrons    \n");
    fprintf(fp,"                          n      4f                       \n");
    fprintf(fp,"average                 <r > after Freeman and            \n");
    fprintf(fp,"Desclaux Journal of Magnetism and Magnetic Materials       \n");
    fprintf(fp,"12 (1979) 11.                                             \n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|              2    2        4    4        6    6         |\n");
    fprintf(fp,"|RE3+        <r >/a0       <r >/a0       <r >/a0          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    s="Ce";a=r2[ 1];b=r4[ 1];c=r6[ 1];fprintf(fp,textfw,s,a,b,c);
    s="Pr";a=r2[ 2];b=r4[ 2];c=r6[ 2];fprintf(fp,texti ,s,a,b,c);
    s="Nd";a=r2[ 3];b=r4[ 3];c=r6[ 3];fprintf(fp,textfw,s,a,b,c);
    s="Pm";a=r2[ 4];b=r4[ 4];c=r6[ 4];fprintf(fp,texti ,s,a,b,c);
    s="Sm";a=r2[ 5];b=r4[ 5];c=r6[ 5];fprintf(fp,textfw,s,a,b,c);
    s="Eu";a=r2[ 6];b=r4[ 6];c=r6[ 6];fprintf(fp,textfw,s,a,b,c);
    s="Gd";a=r2[ 7];b=r4[ 7];c=r6[ 7];fprintf(fp,textfw,s,a,b,c);
    s="Tb";a=r2[ 8];b=r4[ 8];c=r6[ 8];fprintf(fp,textfw,s,a,b,c);
    s="Dy";a=r2[ 9];b=r4[ 9];c=r6[ 9];fprintf(fp,textfw,s,a,b,c);
    s="Ho";a=r2[10];b=r4[10];c=r6[10];fprintf(fp,textfw,s,a,b,c);
    s="Er";a=r2[11];b=r4[11];c=r6[11];fprintf(fp,textfw,s,a,b,c);
    s="Tm";a=r2[12];b=r4[12];c=r6[12];fprintf(fp,textfw,s,a,b,c);
    s="Yb";a=r2[13];b=r4[13];c=r6[13];fprintf(fp,textfw,s,a,b,c);
    fprintf(fp,"\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|              2   o 2       4   o 4       6   o 6        |\n");
    fprintf(fp,"|RE3+        <r >/ A       <r >/ A       <r >/ A          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    a2 = A0_BOHR*A0_BOHR;
    a4 = a2*a2;
    a6 = a2*a4;
    s="Ce";a=r2[ 1]*a2;b=r4[ 1]*a4;c=r6[ 1]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Pr";a=r2[ 2]*a2;b=r4[ 2]*a4;c=r6[ 2]*a6;fprintf(fp,sexti ,s,a,b,c);
    s="Nd";a=r2[ 3]*a2;b=r4[ 3]*a4;c=r6[ 3]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Pm";a=r2[ 4]*a2;b=r4[ 4]*a4;c=r6[ 4]*a6;fprintf(fp,sexti ,s,a,b,c);
    s="Sm";a=r2[ 5]*a2;b=r4[ 5]*a4;c=r6[ 5]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Eu";a=r2[ 6]*a2;b=r4[ 6]*a4;c=r6[ 6]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Gd";a=r2[ 7]*a2;b=r4[ 7]*a4;c=r6[ 7]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Tb";a=r2[ 8]*a2;b=r4[ 8]*a4;c=r6[ 8]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Dy";a=r2[ 9]*a2;b=r4[ 9]*a4;c=r6[ 9]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Ho";a=r2[10]*a2;b=r4[10]*a4;c=r6[10]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Er";a=r2[11]*a2;b=r4[11]*a4;c=r6[11]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Tm";a=r2[12]*a2;b=r4[12]*a4;c=r6[12]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Yb";a=r2[13]*a2;b=r4[13]*a4;c=r6[13]*a6;fprintf(fp,sextfw,s,a,b,c);
    fprintf(fp,"\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|                 2                 4                  6  |\n");
    fprintf(fp,"|  a := alpha_J*10    b := beta_J*10    c := gamma_J*10   |\n");
    fprintf(fp,"|                                                         |\n");
    fprintf(fp,"|              2   o 2       4   o 4       6   o 6        |\n");
    fprintf(fp,"|RE3+      a*<r >/ A     b*<r >/ A     c*<r >/ A          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    a2 = A0_BOHR*A0_BOHR*100;
    a4 = a2*a2;
    a6 = a2*a4;
    s="Ce";a=r2[ 1]*a2*alpha_J[ 1];
           b=r4[ 1]*a4*beta_J[  1];
           c=r6[ 1]*a6*gamma_J[ 1];
           fprintf(fp,sextfw,s,a,b,c);
    s="Pr";a=r2[ 2]*a2*alpha_J[ 2];
           b=r4[ 2]*a4*beta_J[  2];
           c=r6[ 2]*a6*gamma_J[ 2];
           fprintf(fp,sexti ,s,a,b,c);
    s="Nd";a=r2[ 3]*a2*alpha_J[ 3];
           b=r4[ 3]*a4*beta_J[  3];
           c=r6[ 3]*a6*gamma_J[ 3];
           fprintf(fp,sextfw,s,a,b,c);
    s="Pm";a=r2[ 4]*a2*alpha_J[ 4];
           b=r4[ 4]*a4*beta_J[  4];
           c=r6[ 4]*a6*gamma_J[ 4];
           fprintf(fp,sexti ,s,a,b,c);
    s="Sm";a=r2[ 5]*a2*alpha_J[ 5];
           b=r4[ 5]*a4*beta_J[  5];
           c=r6[ 5]*a6*gamma_J[ 5];
           fprintf(fp,sextfw,s,a,b,c);
    s="Eu";a=r2[ 6]*a2*alpha_J[ 6];
           b=r4[ 6]*a4*beta_J[  6];
           c=r6[ 6]*a6*gamma_J[ 6];
           fprintf(fp,sextfw,s,a,b,c);
    s="Gd";a=r2[ 7]*a2*alpha_J[ 7];
           b=r4[ 7]*a4*beta_J[  7];
           c=r6[ 7]*a6*gamma_J[ 7];
           fprintf(fp,sextfw,s,a,b,c);
    s="Tb";a=r2[ 8]*a2*alpha_J[ 8];
           b=r4[ 8]*a4*beta_J[  8];
           c=r6[ 8]*a6*gamma_J[ 8];
           fprintf(fp,sextfw,s,a,b,c);
    s="Dy";a=r2[ 9]*a2*alpha_J[ 9];
           b=r4[ 9]*a4*beta_J[  9];
           c=r6[ 9]*a6*gamma_J[ 9];
           fprintf(fp,sextfw,s,a,b,c);
    s="Ho";a=r2[10]*a2*alpha_J[10];
           b=r4[10]*a4*beta_J[ 10];
           c=r6[10]*a6*gamma_J[10];
           fprintf(fp,sextfw,s,a,b,c);
    s="Er";a=r2[11]*a2*alpha_J[11];
           b=r4[11]*a4*beta_J[ 11];
           c=r6[11]*a6*gamma_J[11];
           fprintf(fp,sextfw,s,a,b,c);
    s="Tm";a=r2[12]*a2*alpha_J[12];
           b=r4[12]*a4*beta_J[ 12];
           c=r6[12]*a6*gamma_J[12];
           fprintf(fp,sextfw,s,a,b,c);
    s="Yb";a=r2[13]*a2*alpha_J[13];
           b=r4[13]*a4*beta_J[ 13];
           c=r6[13]*a6*gamma_J[13];
           fprintf(fp,sextfw,s,a,b,c);
    fclose(fp);
}
/*------------------------------------------------------------------------------
ENDEMODUL    T H E T A     C
------------------------------------------------------------------------------*/
