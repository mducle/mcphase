// File sfit.h
// $Id: sfit.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: sfit.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef SFIT_H
#define SFIT_H 1

#include <qlistbox.h>
#include <qlabel.h>
#include <qpushbutton.h>
#include <qmultilinedit.h>
#include <qevent.h>
#include <qprogressbar.h>
#include <qcheckbox.h>

#include "simplex.h"
#include "grafic.h"

#define DEF_LIMIT 1e-6
 
// ******************************************************************
class MLineEd : public QMultiLineEdit
{
  Q_OBJECT
  
 int iHeight;
 private slots:
  void CheckLength();

 public:
  MLineEd(QWidget *parent=0, const char *name=0);
  void setCellHeight(int i) {QTableView::setCellHeight(i);}
  void setFixedHeight(int H,int nL)
              {QWidget::setFixedHeight(H); iHeight=nL;}
   int getHeight(void) {return iHeight;}
  void setHeight(const int i) {iHeight=i;}			     
   int IsNumFilled(void);

//    protected:
//    void keyPressEvent(QKeyEvent *e);
 

};
// ***************************************************************
class SFitWid : public QWidget
{
  Q_OBJECT

private slots:
    void Close();
    void FunctionSel(int);
    void NoCurvChg();
    void ChkLength();
    void Show();
    void Fit();
    void PrintIt();
    void StopFit();    
    void Iterate();
private:
    Simplex *SI; 

    //CRDataInfo * TDI;

         QListBox * FitFns;
                int nFns;
                int iSelFn;
		
           QLabel * ParText;
                int nPars;
		int ncPars;

          MLineEd * ParMEd;
           double * lfP;
           double lfC[1]; //Sonderfall fuer Lorentz a12

          MLineEd * LimitMEd;
           double * lfL;

          MLineEd * StepMEd;
           double * lfSt;

    QLineEdit * NoCurEd;
            int nCurv; 
	    
    QLineEdit * MaxItEd;
            int nMaxIt;

    QPushButton * FitB;
    QPushButton * ShB;
    QPushButton * ClB;
    QPushButton * PrB;
    QPushButton * StB;
              int bStop;
     QProgressBar * Prg;
         QLabel * PrgL;
      QCheckBox * Anim;
	 
         QTimer * t;
	      int FitEnd;   
         double *lfH;

        QFrame  * FitData;
      struct PlotPars PLP,PLE;
          QFont WF,BF,PF; //Fonts for Windows Buttons and Plots
        int iQuit;
  static int PFontSz;
  static int iChx,iChy;
  static int nMaxPars;

        char szOut[2048];
  
     DataFile *DF;  //Inputdata
    AsciiFile *OF;  //Fitted Curve
        CRData *TD;  

    AsciiFile *EF;  //Error Data
        CRData *ED;  
	
  static int w,h;
  static QRect MainW;
  static QRect FitFnsW;
  static QRect ParTextW;
  static QRect ClBW;
  static QRect ShBW;
  static QRect FitBW;
  static QRect PrBW;
  static QRect StBW;
  static QRect PrgW;
  static QRect ResultW;
  static QRect FitW;
  static QRect ErrW;
  static QRect NoCurW;     
  static QRect MaxItW;     
  static QRect ParMEdW;     
  static QRect StepMEdW;     
  static QRect LimitMEdW;     

protected:
//    void DrawIt(QPainter &p);
      void paintEvent(QPaintEvent *);
      void keyPressEvent(QKeyEvent *e);
      void ShowParList(void);
      void InitPars(void);      
       int CollectPars(void);
       int CollectSteps(void);
public:
    SFitWid(DataFile *df,QWidget *parent=0, const char *name=0);
    AsciiFile * GetData(void){return OF;}   

    void SetEndClose(void){iQuit=END_CLOSE;}
    void SetEndExit(void){iQuit=END_QUIT;}
};
// ******************************************************************
class FunWid : public QWidget
{
  Q_OBJECT

private slots:
    void Close();
    void FunctionSel(int);
    void NoCurvChg();
    void ChkLength();
    void Calculate();
    void PrintIt();
private:

         QListBox * FitFns;
                int nFns;
                int iSelFn;
		
           QLabel * ParText;
                int nPars;
		
          MLineEd * ParMEd;
           double * lfP;
	      
    QLineEdit * NoCurEd;
            int nCurv; 
	    
    QLineEdit * xStartEd;
    QLineEdit * xEndEd;
    QLineEdit * nStepsEd;
 
    QPushButton * CalcB;
    QPushButton * ClB;
    QPushButton * PrB;
	 
        QFrame  * FitData;
      struct PlotPars PLP;
       QFont WF,BF,PF; // Fonts for Windows Buttons and Plots
        CRData  * OutData;

  static int PFontSz;
  static int iChx,iChy;
  static int nMaxPars;

  static int w,h;
  static QRect MainW;
  static QRect FitFnsW;
  static QRect ParTextW;
  static QRect CalcBW;
  static QRect ClBW;
  static QRect PrBW;
  static QRect ResultW;
  static QRect NoCurW;     
  static QRect ParMEdW;     
  static QRect xStEdW;     
  static QRect xEnEdW;     
  static QRect nStEdW;     

protected:
//    void DrawIt(QPainter &p);
      void paintEvent(QPaintEvent *);
      void keyPressEvent(QKeyEvent *e);
      void ShowParList(void);
      void InitPars(void);      
       int CollectPars(void);
public:
    FunWid(QWidget *parent=0, const char *name=0);

};
// ******************************************************************

#endif
