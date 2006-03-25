//File:stdinc.h
// $Id: stdinc.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: stdinc.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef STDINC_H
#define STDINC_H 1

#define TRUE 1
#define FALSE 0

#define ON 1
#define OFF 0

#define INCH 25.4

#define MAXPATH 256
#define MAX_LINELENGTH 1024

#define SQUARE(X) ((X)*(X))

#define SWAP(TYP,X,Y) 	{TYP D=X;X=Y;Y=D;}

#define STRSWAP(A,B,H) 	{strcpy(H,A);strcpy(A,B);strcpy(B,H);}

#define GDATA int
 typedef struct graf_point {GDATA x,y,v;} GPoint, *pGPoint;
 typedef struct graf_wind {GDATA l,t,r,b;} GWind, *pGWind;

#define MDATA double
 typedef struct m_point {MDATA x,y,v;} MPoint, *pMPoint;

#define WINDOW(W) {window(W.l,W.t,W.r,W.b);}
/* left,top,right,bottom window coorinates */

// typedef struct plot_point
//	  {int x,y0,y1;} PLOTP;

#endif //STDINC_H
