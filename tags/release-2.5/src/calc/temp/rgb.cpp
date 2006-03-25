// File: rgb.cpp
// $Log: rgb.cpp,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//
// Revision 1.3  1999/03/15 09:08:37  herbie
// *** empty log message ***
//

#include <qwidget.h>
#include <qapplication.h>
#include <qpushbutton.h>
#include <qscrollbar.h>
#include <qlcdnumber.h>
#include <qfont.h>

#include <stdio.h>

#include "rgb.h"

const char *szRCSID={"$Id: rgb.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};

// **********************************************************
void SetWidgColor(QWidget *w, QColor c)
{QColorGroup cg=w->colorGroup();

 QColorGroup nc( c, cg.background(),
                 cg.light(), cg.dark(), cg.mid(),
                 c, cg.base() );
 QPalette p=w->palette();
 p.setNormal(nc);
 p.setActive(nc);
 w->setPalette(p);
 w->show();

}
// ************************************************************
void SetWidgBackColor(QWidget *w, QColor c)
{QColorGroup cg=w->colorGroup();

 QColorGroup nc( cg.foreground(), c,
                 cg.light(), cg.dark(), cg.mid(),
                 cg.text(), c );
 QPalette p=w->palette();
 p.setNormal(nc);
 p.setActive(nc);
 w->setPalette(p);
 w->setBackgroundMode(QWidget::PaletteBackground);
 w->show();

}
// **********************************************************

// *******************************************************
RGBWidget::RGBWidget( QWidget *parent, const char *name ) : 
                                                       QWidget( parent, name )
{
    setGeometry(20,20,300,320);

    quit = new QPushButton( "Quit", this, "quit" );
    quit->setGeometry( 130, 280, 75, 30 );
    quit->setFont( QFont( "Helvetica", 14, QFont::Bold ) );

    connect( quit, SIGNAL(clicked()), qApp, SLOT(quit()) );

    RedLCD = new LCDBar("Red",red,this,"red");
    connect(RedLCD, SIGNAL(ValueChanged(int)) , SLOT(NewR(int)) );

    GreenLCD = new LCDBar("Green",green,this,"green");
    GreenLCD->move(0,105);
    connect(GreenLCD, SIGNAL(ValueChanged(int)) , SLOT(NewG(int)) );

    BlueLCD = new LCDBar("Blue",blue,this,"blue");
    BlueLCD->move(0,210);
    connect(BlueLCD, SIGNAL(ValueChanged(int)) , SLOT(NewB(int)) );

    Sample = new QFrame(this,"sample");
    Sample->setFrameStyle( QFrame::Panel | QFrame::Sunken );
    Sample->setGeometry( 150, 105, 100, 100 );
    SetWidgBackColor(Sample, black);
                
    show();
    
    uRed=0;
    uGreen=0;
    uBlue=0;
  
}
// ********************************************************
void RGBWidget::NewR(int iR)
{ uRed=iR;
 //fprintf(stderr,"New R: %02x\n",iR);
 SetRGB();
}
// ********************************************************
void RGBWidget::NewG(int iR)
{ uGreen=iR;
 //fprintf(stderr,"New G: %02x\n",iR);
 SetRGB();
}
// ********************************************************
void RGBWidget::NewB(int iR)
{uBlue=iR;
 //fprintf(stderr,"New B: %02x\n",iR);
 SetRGB();
}
// ********************************************************
void RGBWidget::SetRGB(void)
{char szB[20];
 sprintf(szB,"#%02X%02X%02X",uRed,uGreen,uBlue);
 QColor bc;
 bc.setNamedColor(szB);
 SetWidgBackColor(Sample, bc );
}    
// ********************************************************
LCDBar::LCDBar(const char *szText, QColor c, QWidget *parent, const char *name ):
                                   QWidget(parent,name)
{
    setGeometry(0,0,120,100);
    lab = new QLabel(this);
    lab->setFrameStyle( QFrame::Panel | QFrame::Sunken );
    lab->setGeometry(3,3,117,97);
    lab->setAlignment(AlignLeft|AlignTop);
    lab->setText(szText);
   
    lcd  = new QLCDNumber( 2, this, "lcd" );
    lcd->setHexMode();
    lcd->setGeometry(30,25,50,50);
    lcd->setSegmentStyle(QLCDNumber::Filled);
    sBar = new QScrollBar( 0, 255,                       // range
                           1, 10,                       // line/page steps
                           0,                           // inital value
                           QScrollBar::Horizontal,      // orientation
                           this, "scrollbar" );
    SetWidgColor(lcd,c);
    sBar->setGeometry(8,80,110,20); 
    connect( sBar, SIGNAL(valueChanged(int)), SLOT(Display(int)) );
//    connect( this, SIGNAL(ValueChanged(int)), lcd, SLOT(display(int)) ); 
}
// *********************************************************
void LCDBar::Display(int iD)
{lcd->display(iD);
 emit ValueChanged(iD);
}
// *********************************************************
int main( int argc, char **argv )
{
    QApplication a( argc, argv );

    RGBWidget w;
    a.setMainWidget( &w );
    w.show();
    return a.exec();
}
// ***********************************************************
