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
[3] W. B. Lewis et al.: The Journal of Chemical Physics 53 (1970) p 809
 
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
 
#define NR_RE 29
 
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
/* 4f_14       */    1.0 * 0                ,
/* 5f_2 : U4+  */   -1.0 * 2*2*13/3/3/5/5/11,
/* 5f_3 : U3+  */   -1.0 * 7/3/3/11/11      ,
/* 5f_3 : Np4+ */   -1.0 * 7/3/3/11/11      ,
/* 4f_5 : Nd2+ */  1.0 * 13/3/3/5/7       ,
/* 4f_6 : Sm2+ */   1.0 * 0                ,
/* 4f_7 : Eu2+ */   1.0 * 0                ,
/* 4f_8 : Gd2+ */   -1.0 * 1/3/3/11         ,
/* 4f_9 : Tb2+ */   -1.0 * 2/3/3/5/7        ,
/* 4f_10: Dy2+ */   -1.0 * 1/2/3/3/5/5      ,
/* 4f_11: Ho2+ */    1.0 * 2*2/3/3/5/5/7    ,
/* 4f_12: Er2+ */    1.0 * 1/3/3/11         ,
/* 4f_13: Tm2+ */    1.0 * 2/3/3/7          ,
/* 4f_14: Yb2+ */    1.0 * 0                ,
/* 3d_2 : V3+  */    1.0 * 0  
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
/* 4f_14       */    1.0 * 0                              ,
/* 4f_2 : U4+  */   -1.0 * 2*2/3/3/5/11/11,
/* 5f_3 : U3+  */   -1.0 * 2*2*2*17/3/3/3/11/11/11/13     ,
/* 5f_3 : Np4+ */   -1.0 * 2*2*2*17/3/3/3/11/11/11/13     ,
/* 4f_5 : Nd2+ */  1.0 * 2*13/3/3/3/5/7/11              ,
/* 4f_6 : Sm2+ */    1.0 * 0                              ,
/* 4f_7 : Eu2+ */    1.0 * 0                              ,
/* 4f_8 : Gd2+ */  1.0 * 2/3/3/3/5/11/11                ,
/* 4f_9 : Tb2+ */   -1.0 * 2*2*2/3/3/3/5/7/11/13          , 
/* 4f_10: Dy2+ */  -1.0 * 1/2/3/5/7/11/13                ,
/* 4f_11: Ho2+ */  1.0 * 2/3/3/5/7/11/13                ,
/* 4f_12: Er2+ */   1.0 * 2*2*2/3/3/3/3/5/11/11          ,
/* 4f_13: Tm2+ */  -1.0 * 2/3/5/7/11                     ,
/* 4f_14: Yb2+ */  1.0 * 0                              ,
/* 3d_2 : V3+  */    1.0 * 0  
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
/* 4f_14       */    1.0 * 0                                 ,
/* 5f_2 : U4+  */    1.0 * 2*2*2*2*17/3/3/3/3/5/7/11/11/13,
/* 5f_3 : U3+  */   -1.0 * 5*17*19/3/3/3/7/11/11/11/13/13    ,
/* 5f_3 : Np4+ */   -1.0 * 5*17*19/3/3/3/7/11/11/11/13/13    ,
/* 4f_5 : Nd2+ */    1.0 * 0                                 ,
/* 4f_6 : Sm2+ */    1.0 * 0                                 ,
/* 4f_7 : Eu2+ */    1.0 * 0                                 ,
/* 4f_8 : Gd2+ */  -1.0 * 1/3/3/3/3/7/11/11/13              ,
/* 4f_9 : Tb2+ */ 1.0 * 2*2/3/3/3/7/11/11/13/13           ,
/* 4f_10: Dy2+ */ -1.0 * 5/3/3/3/7/11/11/13/13             ,
/* 4f_11: Ho2+ */ 1.0 * 2*2*2/3/3/3/7/11/11/13/13         ,
/* 4f_12: Er2+ */ -1.0 * 5/3/3/3/3/7/11/11/13              ,
/* 4f_13: Tm2+ */ 1.0 * 2*2/3/3/3/7/11/13                 ,
/* 4f_14: Yb2+ */ 1.0 * 0                                 ,
/* 3d_2 : V3+  */    1.0 * 0  
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
    fprintf(fp,"|RE      alpha_J*10     beta_J*10       gamma_J*10        |\n");
    fprintf(fp,"===========================================================\n");
s="Ce3+";a=x*alpha_J[ 1];b=y*beta_J[ 1];c=z*gamma_J[ 1];fprintf(fp,text,s,a,b,c);
s="Pr3+";a=x*alpha_J[ 2];b=y*beta_J[ 2];c=z*gamma_J[ 2];fprintf(fp,text,s,a,b,c);
s="Nd3+";a=x*alpha_J[ 3];b=y*beta_J[ 3];c=z*gamma_J[ 3];fprintf(fp,text,s,a,b,c);
s="Pm3+";a=x*alpha_J[ 4];b=y*beta_J[ 4];c=z*gamma_J[ 4];fprintf(fp,text,s,a,b,c);
s="Sm3+";a=x*alpha_J[ 5];b=y*beta_J[ 5];c=z*gamma_J[ 5];fprintf(fp,text,s,a,b,c);
s="Eu3+";a=x*alpha_J[ 6];b=y*beta_J[ 6];c=z*gamma_J[ 6];fprintf(fp,text,s,a,b,c);
s="Gd3+";a=x*alpha_J[ 7];b=y*beta_J[ 7];c=z*gamma_J[ 7];fprintf(fp,text,s,a,b,c);
s="Tb3+";a=x*alpha_J[ 8];b=y*beta_J[ 8];c=z*gamma_J[ 8];fprintf(fp,text,s,a,b,c);
s="Dy3+";a=x*alpha_J[ 9];b=y*beta_J[ 9];c=z*gamma_J[ 9];fprintf(fp,text,s,a,b,c);
s="Ho3+";a=x*alpha_J[10];b=y*beta_J[10];c=z*gamma_J[10];fprintf(fp,text,s,a,b,c);
s="Er3+";a=x*alpha_J[11];b=y*beta_J[11];c=z*gamma_J[11];fprintf(fp,text,s,a,b,c);
s="Tm3+";a=x*alpha_J[12];b=y*beta_J[12];c=z*gamma_J[12];fprintf(fp,text,s,a,b,c);
s="Yb3+";a=x*alpha_J[13];b=y*beta_J[13];c=z*gamma_J[13];fprintf(fp,text,s,a,b,c);
s="U4+ ";a=x*alpha_J[15];b=y*beta_J[15];c=z*gamma_J[15];fprintf(fp,text,s,a,b,c);
s="U3+ ";a=x*alpha_J[16];b=y*beta_J[16];c=z*gamma_J[16];fprintf(fp,text,s,a,b,c);
s="Np4+";a=x*alpha_J[17];b=y*beta_J[17];c=z*gamma_J[17];fprintf(fp,text,s,a,b,c);
s="Nd2+";a=x*alpha_J[18];b=y*beta_J[18];c=z*gamma_J[18];fprintf(fp,text,s,a,b,c);
s="Sm2+";a=x*alpha_J[19];b=y*beta_J[19];c=z*gamma_J[19];fprintf(fp,text,s,a,b,c);
s="Eu2+";a=x*alpha_J[20];b=y*beta_J[20];c=z*gamma_J[20];fprintf(fp,text,s,a,b,c);
s="Gd2+";a=x*alpha_J[21];b=y*beta_J[21];c=z*gamma_J[21];fprintf(fp,text,s,a,b,c);
s="Tb2+";a=x*alpha_J[22];b=y*beta_J[22];c=z*gamma_J[22];fprintf(fp,text,s,a,b,c);
s="Dy2+";a=x*alpha_J[23];b=y*beta_J[23];c=z*gamma_J[23];fprintf(fp,text,s,a,b,c);
s="Ho2+";a=x*alpha_J[24];b=y*beta_J[24];c=z*gamma_J[24];fprintf(fp,text,s,a,b,c);
s="Er2+";a=x*alpha_J[25];b=y*beta_J[25];c=z*gamma_J[25];fprintf(fp,text,s,a,b,c);
s="Tm2+";a=x*alpha_J[26];b=y*beta_J[26];c=z*gamma_J[26];fprintf(fp,text,s,a,b,c);
s="Yb2+";a=x*alpha_J[27];b=y*beta_J[27];c=z*gamma_J[27];fprintf(fp,text,s,a,b,c);
s="V3+";a=x*alpha_J[28];b=y*beta_J[28];c=z*gamma_J[28];fprintf(fp,text,s,a,b,c);
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
/* 4f_14       */   0.0                 ,
/* 5f_2 : U4+  */   2.042               , /* Freeman et al. PRB 13 (1976) 1168 */
/* 5f_3 : U3+  */   2.346               , /* Freeman et al. PRB 13 (1976) 1168 */
/* 5f_3 : Np4+ */   1.884               , /* Lewis et al. J. Chem Phys. 53 (1970) 809*/
/* 4f_5 : Nd2+ */   1.392               ,
/* 4f_6 : Sm2+ */   1.197               ,
/* 4f_7 : Eu2+ */   1.098               ,
/* 4f_8 : Gd2+ */   1.028               ,
/* 4f_9 : Tb2+ */   0.968               ,
/* 4f_10: Dy2+ */   0.913               ,
/* 4f_11: Ho2+ */   0.866               ,
/* 4f_12: Er2+ */   0.824               ,
/* 4f_13: Tm2+ */   0.785               ,
/* 4f_14: Yb2+ */   0.750               ,
/* 3d_2 : V3+  */   4                     /* just estimate !*/ 
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
/* 4f_14       */    0.0                ,
/* 5f_2 : U4+  */   7.632               , /* Freeman et al. PRB 13 (1976) 1168 */
/* 5f_3 : U3+  */  10.906               , /* Freeman et al. PRB 13 (1976) 1168 */
/* 5f_3 : Np4+ */   6.504               , /* Lewis et al. J. Chem Phys. 53 (1970) 809*/
/* 4f_5 : Nd2+ */   5.344               ,
/* 4f_6 : Sm2+ */   3.861               ,
/* 4f_7 : Eu2+ */   3.368               ,
/* 4f_8 : Gd2+ */   2.975               ,
/* 4f_9 : Tb2+ */   2.655               ,
/* 4f_10: Dy2+ */   2.391               ,
/* 4f_11: Ho2+ */   2.169               ,
/* 4f_12: Er2+ */   1.979               ,
/* 4f_13: Tm2+ */   1.819               ,
/* 4f_14: Yb2+ */   1.677               ,
/* 3d_2 : V3+  */   16                     /* just estimate !*/ 
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
/* 4f_14       */    0.0                ,
/* 5f_2 : U4+  */  47.774               , /* Freeman et al. PRB 13 (1976) 1168 */
/* 5f_3 : U3+  */  90.544               , /* Freeman et al. PRB 13 (1976) 1168 */
/* 5f_3 : Np4+ */  37.80                , /* Lewis et al. J. Chem Phys. 53 (1970) 809*/
/* 4f_5 : Nd2+ */  45.450               ,
/* 4f_6 : Sm2+ */  28.560               ,
/* 4f_7 : Eu2+ */  23.580               ,
/* 4f_8 : Gd2+ */  19.850               ,
/* 4f_9 : Tb2+ */  16.980               ,
/* 4f_10: Dy2+ */  14.730               ,
/* 4f_11: Ho2+ */  12.920               ,
/* 4f_12: Er2+ */  11.450               ,
/* 4f_13: Tm2+ */  10.240               ,
/* 4f_14: Yb2+ */   9.232               ,
/* 3d_2 : V3+  */  64                     /* just estimate !*/ 
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
    fprintf(fp,"|RE          <r >/a0       <r >/a0       <r >/a0          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    s="Ce3+";a=r2[ 1];b=r4[ 1];c=r6[ 1];fprintf(fp,textfw,s,a,b,c);
    s="Pr3+";a=r2[ 2];b=r4[ 2];c=r6[ 2];fprintf(fp,texti ,s,a,b,c);
    s="Nd3+";a=r2[ 3];b=r4[ 3];c=r6[ 3];fprintf(fp,textfw,s,a,b,c);
    s="Pm3+";a=r2[ 4];b=r4[ 4];c=r6[ 4];fprintf(fp,texti ,s,a,b,c);
    s="Sm3+";a=r2[ 5];b=r4[ 5];c=r6[ 5];fprintf(fp,textfw,s,a,b,c);
    s="Eu3+";a=r2[ 6];b=r4[ 6];c=r6[ 6];fprintf(fp,textfw,s,a,b,c);
    s="Gd3+";a=r2[ 7];b=r4[ 7];c=r6[ 7];fprintf(fp,textfw,s,a,b,c);
    s="Tb3+";a=r2[ 8];b=r4[ 8];c=r6[ 8];fprintf(fp,textfw,s,a,b,c);
    s="Dy3+";a=r2[ 9];b=r4[ 9];c=r6[ 9];fprintf(fp,textfw,s,a,b,c);
    s="Ho3+";a=r2[10];b=r4[10];c=r6[10];fprintf(fp,textfw,s,a,b,c);
    s="Er3+";a=r2[11];b=r4[11];c=r6[11];fprintf(fp,textfw,s,a,b,c);
    s="Tm3+";a=r2[12];b=r4[12];c=r6[12];fprintf(fp,textfw,s,a,b,c);
    s="Yb3+";a=r2[13];b=r4[13];c=r6[13];fprintf(fp,textfw,s,a,b,c);
    s="U4+ ";a=r2[15];b=r4[15];c=r6[15];fprintf(fp,textfw,s,a,b,c);
    s="U3+ ";a=r2[16];b=r4[16];c=r6[16];fprintf(fp,textfw,s,a,b,c);
    s="Np4+";a=r2[17];b=r4[17];c=r6[17];fprintf(fp,textfw,s,a,b,c);
    s="Nd2+";a=r2[18];b=r4[18];c=r6[18];fprintf(fp,textfw,s,a,b,c);
    s="Sm2+";a=r2[19];b=r4[19];c=r6[19];fprintf(fp,textfw,s,a,b,c);
    s="Eu2+";a=r2[20];b=r4[20];c=r6[20];fprintf(fp,textfw,s,a,b,c);
    s="Gd2+";a=r2[21];b=r4[21];c=r6[21];fprintf(fp,textfw,s,a,b,c);
    s="Tb2+";a=r2[22];b=r4[22];c=r6[22];fprintf(fp,textfw,s,a,b,c);
    s="Dy2+";a=r2[23];b=r4[23];c=r6[23];fprintf(fp,textfw,s,a,b,c);
    s="Ho2+";a=r2[24];b=r4[24];c=r6[24];fprintf(fp,textfw,s,a,b,c);
    s="Er2+";a=r2[25];b=r4[25];c=r6[25];fprintf(fp,textfw,s,a,b,c);
    s="Tm2+";a=r2[26];b=r4[26];c=r6[26];fprintf(fp,textfw,s,a,b,c);
    s="Yb2+";a=r2[27];b=r4[27];c=r6[27];fprintf(fp,textfw,s,a,b,c);
    s="V3+ ";a=r2[28];b=r4[28];c=r6[28];fprintf(fp,textfw,s,a,b,c);
    fprintf(fp,"\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|              2   o 2       4   o 4       6   o 6        |\n");
    fprintf(fp,"|RE3+        <r >/ A       <r >/ A       <r >/ A          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    a2 = A0_BOHR*A0_BOHR;
    a4 = a2*a2;
    a6 = a2*a4;
    s="Ce3+";a=r2[ 1]*a2;b=r4[ 1]*a4;c=r6[ 1]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Pr3+";a=r2[ 2]*a2;b=r4[ 2]*a4;c=r6[ 2]*a6;fprintf(fp,sexti ,s,a,b,c);
    s="Nd3+";a=r2[ 3]*a2;b=r4[ 3]*a4;c=r6[ 3]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Pm3+";a=r2[ 4]*a2;b=r4[ 4]*a4;c=r6[ 4]*a6;fprintf(fp,sexti ,s,a,b,c);
    s="Sm3+";a=r2[ 5]*a2;b=r4[ 5]*a4;c=r6[ 5]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Eu3+";a=r2[ 6]*a2;b=r4[ 6]*a4;c=r6[ 6]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Gd3+";a=r2[ 7]*a2;b=r4[ 7]*a4;c=r6[ 7]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Tb3+";a=r2[ 8]*a2;b=r4[ 8]*a4;c=r6[ 8]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Dy3+";a=r2[ 9]*a2;b=r4[ 9]*a4;c=r6[ 9]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Ho3+";a=r2[10]*a2;b=r4[10]*a4;c=r6[10]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Er3+";a=r2[11]*a2;b=r4[11]*a4;c=r6[11]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Tm3+";a=r2[12]*a2;b=r4[12]*a4;c=r6[12]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Yb3+";a=r2[13]*a2;b=r4[13]*a4;c=r6[13]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="U4+ ";a=r2[15]*a2;b=r4[15]*a4;c=r6[15]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="U3+ ";a=r2[16]*a2;b=r4[16]*a4;c=r6[16]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Np4+";a=r2[17]*a2;b=r4[17]*a4;c=r6[17]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Nd2+";a=r2[18]*a2;b=r4[18]*a4;c=r6[18]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Sm2+";a=r2[19]*a2;b=r4[19]*a4;c=r6[19]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Eu2+";a=r2[20]*a2;b=r4[20]*a4;c=r6[20]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Gd2+";a=r2[21]*a2;b=r4[21]*a4;c=r6[21]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Tb2+";a=r2[22]*a2;b=r4[22]*a4;c=r6[22]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Dy2+";a=r2[23]*a2;b=r4[23]*a4;c=r6[23]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Ho2+";a=r2[24]*a2;b=r4[24]*a4;c=r6[24]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Er2+";a=r2[25]*a2;b=r4[25]*a4;c=r6[25]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Tm2+";a=r2[26]*a2;b=r4[26]*a4;c=r6[26]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="Yb2+";a=r2[27]*a2;b=r4[27]*a4;c=r6[26]*a6;fprintf(fp,sextfw,s,a,b,c);
    s="V3+ ";a=r2[28]*a2;b=r4[28]*a4;c=r6[26]*a6;fprintf(fp,sextfw,s,a,b,c);
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
    s="Ce3+";a=r2[ 1]*a2*alpha_J[ 1];
           b=r4[ 1]*a4*beta_J[  1];
           c=r6[ 1]*a6*gamma_J[ 1];
           fprintf(fp,sextfw,s,a,b,c);
    s="Pr3+";a=r2[ 2]*a2*alpha_J[ 2];
           b=r4[ 2]*a4*beta_J[  2];
           c=r6[ 2]*a6*gamma_J[ 2];
           fprintf(fp,sexti ,s,a,b,c);
    s="Nd3+";a=r2[ 3]*a2*alpha_J[ 3];
           b=r4[ 3]*a4*beta_J[  3];
           c=r6[ 3]*a6*gamma_J[ 3];
           fprintf(fp,sextfw,s,a,b,c);
    s="Pm3+";a=r2[ 4]*a2*alpha_J[ 4];
           b=r4[ 4]*a4*beta_J[  4];
           c=r6[ 4]*a6*gamma_J[ 4];
           fprintf(fp,sexti ,s,a,b,c);
    s="Sm3+";a=r2[ 5]*a2*alpha_J[ 5];
           b=r4[ 5]*a4*beta_J[  5];
           c=r6[ 5]*a6*gamma_J[ 5];
           fprintf(fp,sextfw,s,a,b,c);
    s="Eu3+";a=r2[ 6]*a2*alpha_J[ 6];
           b=r4[ 6]*a4*beta_J[  6];
           c=r6[ 6]*a6*gamma_J[ 6];
           fprintf(fp,sextfw,s,a,b,c);
    s="Gd3+";a=r2[ 7]*a2*alpha_J[ 7];
           b=r4[ 7]*a4*beta_J[  7];
           c=r6[ 7]*a6*gamma_J[ 7];
           fprintf(fp,sextfw,s,a,b,c);
    s="Tb3+";a=r2[ 8]*a2*alpha_J[ 8];
           b=r4[ 8]*a4*beta_J[  8];
           c=r6[ 8]*a6*gamma_J[ 8];
           fprintf(fp,sextfw,s,a,b,c);
    s="Dy3+";a=r2[ 9]*a2*alpha_J[ 9];
           b=r4[ 9]*a4*beta_J[  9];
           c=r6[ 9]*a6*gamma_J[ 9];
           fprintf(fp,sextfw,s,a,b,c);
    s="Ho3+";a=r2[10]*a2*alpha_J[10];
           b=r4[10]*a4*beta_J[ 10];
           c=r6[10]*a6*gamma_J[10];
           fprintf(fp,sextfw,s,a,b,c);
    s="Er3+";a=r2[11]*a2*alpha_J[11];
           b=r4[11]*a4*beta_J[ 11];
           c=r6[11]*a6*gamma_J[11];
           fprintf(fp,sextfw,s,a,b,c);
    s="Tm3+";a=r2[12]*a2*alpha_J[12];
           b=r4[12]*a4*beta_J[ 12];
           c=r6[12]*a6*gamma_J[12];
           fprintf(fp,sextfw,s,a,b,c);
    s="Yb3+";a=r2[13]*a2*alpha_J[13];
           b=r4[13]*a4*beta_J[ 13];
           c=r6[13]*a6*gamma_J[13];
           fprintf(fp,sextfw,s,a,b,c);
    s="U4+ ";a=r2[15]*a2*alpha_J[15];
           b=r4[15]*a4*beta_J[ 15];
           c=r6[15]*a6*gamma_J[15];
           fprintf(fp,sextfw ,s,a,b,c);
    s="U3+ ";a=r2[16]*a2*alpha_J[16];b=r4[16]*a4*beta_J[ 16];c=r6[16]*a6*gamma_J[16];fprintf(fp,sextfw ,s,a,b,c);
    s="Np4+";a=r2[17]*a2*alpha_J[17];b=r4[17]*a4*beta_J[ 17];c=r6[17]*a6*gamma_J[17];fprintf(fp,sextfw ,s,a,b,c);
    s="Nd2+";a=r2[18]*a2*alpha_J[18];b=r4[18]*a4*beta_J[ 18];c=r6[18]*a6*gamma_J[18];fprintf(fp,sextfw ,s,a,b,c);
    s="Sm2+";a=r2[19]*a2*alpha_J[19];b=r4[19]*a4*beta_J[ 19];c=r6[19]*a6*gamma_J[19];fprintf(fp,sextfw ,s,a,b,c);
    s="Eu2+";a=r2[20]*a2*alpha_J[20];b=r4[20]*a4*beta_J[ 20];c=r6[20]*a6*gamma_J[20];fprintf(fp,sextfw ,s,a,b,c);
    s="Gd2+";a=r2[21]*a2*alpha_J[21];b=r4[21]*a4*beta_J[ 21];c=r6[21]*a6*gamma_J[21];fprintf(fp,sextfw ,s,a,b,c);
    s="Tb2+";a=r2[22]*a2*alpha_J[22];b=r4[22]*a4*beta_J[ 22];c=r6[22]*a6*gamma_J[22];fprintf(fp,sextfw ,s,a,b,c);
    s="Dy2+";a=r2[23]*a2*alpha_J[23];b=r4[23]*a4*beta_J[ 23];c=r6[23]*a6*gamma_J[23];fprintf(fp,sextfw ,s,a,b,c);
    s="Ho2+";a=r2[24]*a2*alpha_J[24];b=r4[24]*a4*beta_J[ 24];c=r6[24]*a6*gamma_J[24];fprintf(fp,sextfw ,s,a,b,c);
    s="Er2+";a=r2[25]*a2*alpha_J[25];b=r4[25]*a4*beta_J[ 25];c=r6[25]*a6*gamma_J[25];fprintf(fp,sextfw ,s,a,b,c);
    s="Tm2+";a=r2[26]*a2*alpha_J[26];b=r4[26]*a4*beta_J[ 26];c=r6[26]*a6*gamma_J[26];fprintf(fp,sextfw ,s,a,b,c);
    s="YB2+";a=r2[27]*a2*alpha_J[27];b=r4[27]*a4*beta_J[ 27];c=r6[27]*a6*gamma_J[27];fprintf(fp,sextfw ,s,a,b,c);
    s="V3+ ";a=r2[28]*a2*alpha_J[28];b=r4[28]*a4*beta_J[ 28];c=r6[28]*a6*gamma_J[28];fprintf(fp,sextfw ,s,a,b,c);
    fclose(fp);
}
/*------------------------------------------------------------------------------
ENDEMODUL    T H E T A     C
------------------------------------------------------------------------------*/
