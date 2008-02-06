// File: rgb.h
// $Id: rgb.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: rgb.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef RGB_H
#define RGB_H 1

#include <qapplication.h>
#include <qpushbutton.h>
#include <qscrollbar.h>
#include <qlcdnumber.h>
#include <qfont.h>
#include <qframe.h>
#include <qlabel.h>

void SetWidgBackColor(QWidget *w, QColor c);
void SetWidgColor(QWidget *w, QColor c);

// *********************************************************
class LCDBar : public QWidget
{
  Q_OBJECT

   QScrollBar  *sBar;
   QLCDNumber  *lcd;
       QLabel  *lab;
signals:
    void ValueChanged(int);
private slots:
    void Display(int);
  
 public:
    LCDBar(const char *szText, QColor c, QWidget *parent=0, const char *name=0 );
};
// *******************************************************
class RGBWidget : public QWidget
{
  Q_OBJECT

public:
    RGBWidget( QWidget *parent=0, const char *name=0 );

private slots:
void NewR(int iR);
void NewG(int iR);
void NewB(int iR);


private:

    QPushButton *quit;
    
         LCDBar *RedLCD;
         LCDBar *GreenLCD;
         LCDBar *BlueLCD;

         QFrame *Sample;
	 
    unsigned uRed;
    unsigned uGreen;
    unsigned uBlue;
    
  void SetRGB(void);    
};
// ********************************************************
#endif
