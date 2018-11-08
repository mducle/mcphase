// File grafic.h
// $Id: grafic.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: grafic.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef GRAFIC_H
#define GRAFIC_H 1

#include <qlabel.h>
#include <qlineedit.h>
#include <qevent.h>
#include <qdialog.h>
#include <qpushbutton.h>
#include <qpainter.h>
#include <qpixmap.h>
#include <qstack.h>
#include <qmultilinedit.h>

#include "cdata.h"
#include "dfile.h"

#define     CIRCLE 0
#define       QUAD 1
#define   TRIANGLE 2
#define    DIAMOND 3
#define      CROSS 4
#define LEFT_CROSS 8
#define     V_LINE 9
#define  ERROR_BAR -9
#define    CONNECT 10
#define        BAR 11

#define TOP_FILENAME 0x0001
#define BOT_FILELIST 0x0002
#define BOT_FILEHEAD 0x0004
#define TOP_IDN      0x0008
#define ZERO_X_LINE  0x0010
#define TOP_INFOLINE 0x0020

#define NORM_FORM 1
#define ARCH_FORM 2

#define MAX_PLOTS 8

 int MinDist(QPointArray *p);
void DrawSymbols(const int iSymbol, const int iSize,
                   QPointArray *p, QPainter *qp, int iStart=0, int iEnd=-1);

//void DrawSymbol(const int iSymbol, const int iSize, QPoint p,QPainter *qp);
void DrawSymbol(const int iSymbol, const int iSize, QPainter *qp,
                QPoint p, QPoint pp=QPoint(-1,-1));
		
//void PaintSymbol(const int iSymbol, const int iSize,const QColor col, QPoint pt,
//                QWidget *PaintW);
void PaintSymbol(const int iSymbol, const int iSize, const QColor col, 
                 QWidget *PaintW, QPoint pt, QPoint ppt=QPoint(-1,-1));

void DrawData(QPainter &p, DataFile *DF, CRData *CD, struct PlotPars &PP);
 int GetMinMax(MPoint &Min, MPoint &Max);

void DrawFormText(QPainter &p, QFont F, const int x0, const int y0,
                  const int flag,  const char *szT);

void PrintText(QMultiLineEdit *w, FILE *f, const int iAppend=0);

inline void PrintRect(QRect &R,const char* szN=0, FILE *f=stderr)
     {fprintf(f,"%s: x:%d y:%d w:%d h:%d\n",(szN?szN:" "),R.x(),R.y(),R.width(),R.height());}

struct PlotPars
   { QFont PFont;
     QRect PRect;
    char * xLab;
    char * yLab;
       int pFlag;
       int iStart;
       int iEnd;
//       int iSize;
//       int iErrBar;
    MPoint DMin,DMax;
    
    int nPlots;
    int pSize[MAX_PLOTS];
    int pSymbol[MAX_PLOTS];
    int pColor[MAX_PLOTS];
    
    PlotPars(void);
    inline void SetSize(const int i, const int v)
                   {if(i>=0 && i<MAX_PLOTS)pSize[i]=v;}
    inline void SetSymbol(const int i, const int v)
                   {if(i>=0 && i<MAX_PLOTS)pSymbol[i]=v;}
    inline void SetColor(const int i, const int v)
                   {if(i>=0 && i<MAX_PLOTS)pColor[i]=v;}
		   
    void Init(void);
    void Copy(const struct PlotPars & P);
    PlotPars(const struct PlotPars & P);
    struct PlotPars& operator=(const struct PlotPars & P);


   };
 
 
// ******************************************************************
class MultiPlot : public QWidget
{
//  Q_OBJECT
private:
  struct PlotPars PLP;

        void paintEvent(QPaintEvent *);
  static int PFontSz;

public:
  MultiPlot(QWidget *parent=0, const char *name=0);
    void SetPlotPars(struct PlotPars &P){PLP=P;}
  struct PlotPars GetPlotPars(void) {return PLP;}
    void DrawIt(QPainter &p);
};
// ***************************************************************
#define END_CLOSE 1
#define END_QUIT 2

class DisplayData: public QWidget
{
  Q_OBJECT
private:
  int nPlots;
  int iFormat;        // Normal Format or Archive Format

  char szXLabel[80];
  char szYLabel[80];

  QLineEdit * xLab;
  QLineEdit * yLab;

  QPushButton * eBut;  // Error Bars
  QPushButton * aBut;  // Archive
  QPushButton * pBut;  // Print
  QPushButton * cBut;  // Close
  int iQuit;           // iQuit=END_CLOSE -> hide Window
                       // iQuit=END_QUIT -> Quit application
 
  MultiPlot *PlotA;
     QFrame *LegA;

  static int w,h;
  static QRect     B0;
  static QRect WMain;
  static QRect   WxL;
  static QRect   WxE;
  static QRect   WyL;
  static QRect   WyE;
  static QRect   WeB;
  static QRect   WaB;
  static QRect   WpB;
  static QRect   WcB;

  static QRect  WPlot;
  static QRect WFrame;

  void keyPressEvent(QKeyEvent *e);

public slots:
  void DataChanged();
  void PrintIt();
  void ErrorBars();
  void ArchiveFormat();
  void Close();  

public:
  DisplayData(QWidget *parent=0, const char *name=0);
           int DataIsValid(void) const {return nPlots!=0;}
           int AddData(DataFile *df);
           int InitPlot(void);
          void SetEndClose(void){iQuit=END_CLOSE;}
          void SetEndExit(void){iQuit=END_QUIT;}
  static QRect GetPlotW(void){return WPlot;}

};
// ***************************************************************
class MTrackWid: public QFrame
{
 Q_OBJECT
 bool RDown; 


protected:
 void mousePressEvent(QMouseEvent *);
 void mouseReleaseEvent(QMouseEvent *);
 void mouseMoveEvent(QMouseEvent *);

 signals:
 void RCoords(QPoint p);
 void LCoords(QPoint p);
 
 public:
 MTrackWid(QWidget *parent=0, const char *name=0);

};
// ****************************************************
class ManipData: public QWidget
{
  Q_OBJECT

  DataFile * DF;
     CRData * TD;
  
  QLabel * InfoLT;
  QLabel * InfoLB;

  MTrackWid * PlotF;
        QRect EffP;
       struct PlotPars PLP;

  QPushButton *Bts[10];  //F1-F10
  int iQuit;           // iQuit=END_CLOSE -> hide Window
                       // iQuit=END_QUIT -> Quit application

  int iStep;
 
  QPixmap CursPM;
   QPoint PM0;
     bool First;
   QPoint CurCurs;

  QStack<QPoint> Range;
            int iRCount;
           bool bRSelect;

  static   int w,h;
  static   int iFS;
  static QRect WMain;
//  static QRect WIT;
//  static QRect WIB;
  static QSize WPl;
//  static QRect B0;

 protected:
  void paintEvent(QPaintEvent *);
  void keyPressEvent(QKeyEvent *e);
  void MoveCursor(const  QColor c);

 public slots:
  void ShowCoords(QPoint p);
  void ClosestPoint(QPoint p);
  void ErrorBars();
  void MeanVal();
  void RemovePoint();
  void NewRange();
  void RangeBack();
  void Integrate();
  void SaveRange();
  void Close();

 public:
  ManipData(DataFile *df, QWidget *parent=0, const char *name=0); 
  int IsValid(void){return TD!=NULL;}
 void SetEndClose(void){iQuit=END_CLOSE;}
 void SetEndExit(void){iQuit=END_QUIT;}

};
// ***************************************************************
class ShowWid : public QWidget
{
  Q_OBJECT
private:
        DataFile * TD;
 struct PlotPars PLP;
    int iQuit;
   char szXLabel[80];
   char szYLabel[80];
    int iColX,iColY;
       
  static int w,h;
   
  void paintEvent(QPaintEvent *);
  void keyPressEvent(QKeyEvent *e);
  void InitPlot(void);
  void Close(void);
  void PrintIt(void);
public:
  ShowWid(DataFile *df, QWidget *parent=0, const char *name=0);
  int DataIsValid(void){return TD!=NULL;}
 void SetEndClose(void){iQuit=END_CLOSE;}
 void SetEndExit(void){iQuit=END_QUIT;}
};
// ***************************************************************
class TheView : public QWidget
{
  Q_OBJECT
private:
        DataFile * TD;
        DataFile * AD;
        DataFile * Tmp;

 struct PlotPars PLP;
    int iQuit;
   char szXLabel[80];
   char szYLabel[80];
    int iColX[8],iColY[8];
       
  static int w,h;
  static QRect WFrame;
  static QRect WPlot;
  QFrame *PlotF;

  void paintEvent(QPaintEvent *);
  void keyPressEvent(QKeyEvent *e);
  void InitPlot(void);
  void Close(void);
  void PrintIt(void);
public:
  TheView(DataFile *df, DataFile *af, QWidget *parent=0, const char *name=0);
  int DataIsValid(void){return (TD!=NULL && AD!=0);}
 void SetEndClose(void){iQuit=END_CLOSE;}
 void SetEndExit(void){iQuit=END_QUIT;}
};
// ***************************************************************
class ShowTxt : public QWidget
{
  Q_OBJECT
private:
  char szFile[MAXPATH+1];
  QMultiLineEdit *Txt;
  int iQuit;           // iQuit=END_CLOSE -> hide Window
                       // iQuit=END_QUIT -> Quit application
      
  static int w,h;
   
  void keyPressEvent(QKeyEvent *e);
  void resizeEvent(QResizeEvent *);
  void Close();

signals: 
  void Removed();

public:
  ShowTxt(const char *szF, QRect *Pos=0,  QWidget *parent=0, const char *name=0);
  int FileValid(void) {return (*szFile!=0);}
  void SetEndClose(void){iQuit=END_CLOSE;}
  void SetEndExit(void){iQuit=END_QUIT;}
  const char* textLine (int line){return (Txt ? Txt->textLine(line) : 0);}
  void setText (const char *szT){if(Txt)Txt->setText(szT);}
  void setCursorPosition(int line)
                         {if(Txt)Txt->setCursorPosition(line,0,FALSE);}
  int numLines (void) {return (Txt? Txt->numLines() : 0);} 
};
// ***************************************************************
#endif




















