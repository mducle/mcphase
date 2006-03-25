#include <qapplication.h>
#include <qstring.h>
#include <qlabel.h>
#include <qlistbox.h>
#include <qlineedit.h>
#include <qevent.h>
#include <qmultilinedit.h>
#include <qtooltip.h>
#include <qfile.h>
#include <qmessagebox.h>
#include <qdatastream.h>
#include <qlayout.h>
#include <qpaintdevice.h>
#include <qpainter.h>
#include <qpointarray.h>
#include <qfontmetrics.h>
#include <qkeycode.h>

#include <stdio.h>
#include <math.h>

#include "WidGen.h"

#include "grafic.h"
#include "stdinc.h"
#include "stdfunc.h"
#include "cdata.h"
#include "dfile.h"
#include "strings.h"

//extern char szAnswer[2048];
extern QPoint SCREEN;

// ***************************************************************
TDataInfo::TDataInfo(DataFile *df, const int x, const  int y,
                     QWidget *parent, const char *name):
           QMultiLineEdit(parent,name)
{
 const int InfoW=550;
 const int InfoH=350;
 setBackgroundColor(lightGray);
 setReadOnly(TRUE);
 setGeometry(x,y,InfoW,InfoH);
 String S;
 if(df!=0)df->SPrintInfo(S,PR_ALL);
 else S.Set("Illegal Data Block");
 setText((const char *)S);
 QToolTip::add(this,"Statistics of Data File");
}
// ****************************************************************
TDataInfo::TDataInfo(CRData *td,
                     const int x, const  int y, const int w, const int h,
                     QWidget *parent, const char *name):
           QMultiLineEdit(parent,name)
{
 setBackgroundColor(lightGray);
 setReadOnly(TRUE);
 setGeometry(x,y,w,h);
 String S;
 if(td!=0)td->SPrintStatistic(S);
 else S.Set("Illegal Data Block");
 setText((const char *)S);
 QToolTip::add(this,"Statistics of Data File");
}
// ****************************************************************
TDataInfo::TDataInfo(CRData *td, const int x, const  int y,
                     QWidget *parent, const char *name):
           QMultiLineEdit(parent,name)
{
 const int InfoW=550;
 const int InfoH=350;
 setBackgroundColor(lightGray);
 setReadOnly(TRUE);
 setGeometry(x,y,InfoW,InfoH);
 String S;
 if(td!=0)td->SPrintStatistic(S);
 else S.Set("Illegal Data Block");
 setText((const char *)S);
 QToolTip::add(this,"Statistics of Data File");
}
// *********************************************************
// **********************************************************
//
// Construct SaveFile
//
// *********************************************************

SaveFile::SaveFile(const char * szFileName,DataFile *df,
                         QWidget *parent, const char *name):
             QDialog( parent, name,TRUE )
{GWind WMain= {100,         50,                   600,        450};
 GWind WEdit= {10,          10,                   WMain.r-20, WMain.b-100};
 GWind WTextL={10,          WEdit.t+WEdit.b+10,   100,        20};
 GWind WTextE={WTextL.r+10, WTextL.t,             430,        20};
// GWind WColL= {10,          WTextL.t+WTextL.b+10, 50,         20};
// GWind WColE= {WColL.r+20,  WColL.t,              70,         20};        
 GWind WOk=   {10,          WTextE.t+WTextE.b+10,   90,         20};
 GWind WCa=   {WOk.r+10,    WOk.t,                90,         20};

 GWind EditW={10,10,WMain.r/2-10,WEdit.b};
 GWind WInfo={EditW.r+10,EditW.t,EditW.r-10,WEdit.b};

// int Offs=WColL.r+WColE.r+30;


 DF=df;
 strncpy(szFile,szFileName,MAXPATH);
 szFile[MAXPATH]=0;

 //nCols=(TD->GetCols()>2?3:2);
 
 String S;
 S.Setf("Save %s",szFileName);
 setCaption((const char *)S);

 setGeometry(WMain.l,WMain.t,WMain.r,WMain.b);

 TextW=0;

 MsgText=0;
 TextEd=0;

 int iTxtL=DF->GetNoTxtL();
 if(iTxtL>1)
   {TextW = new QMultiLineEdit(this,"textwrite");
    TextW->setGeometry(EditW.l,EditW.t,EditW.r,EditW.b);
    int i;
    for(i=0;i<iTxtL;i++)TextW->insertLine(DF->GetInfoText(i,0));
    TDI=new TDataInfo(DF->GetCRData(),
                       WInfo.l,WInfo.t,WInfo.r, WInfo.b,this,"DefInfo");
   } 
 
 else
 //-----------------------MULTI LINE EDIT INFO DEFAULT BLOCK:
 {  
  TDI=new TDataInfo(DF,WEdit.t,WEdit.l,this,"DefInfo");
  //FileW=new QMultiLineEdit(this,"DefInfo");
  //FileW->setBackgroundColor(lightGray);
  //FileW->setReadOnly(TRUE);
  //FileW->setGeometry(WEdit.l,WEdit.t,WEdit.r,WEdit.b);
  //TD->SPrintStatistic(szAnswer);
  //FileW->setText(szAnswer);
  //QToolTip::add(FileW,"Statistics of Data File");
 
   S.Setf("Enter Text line:");
   MsgText=new QLabel((const char *)S, this, "TextLabel");
   MsgText->setAlignment(AlignLeft);
   MsgText->setFrameStyle( QFrame::Panel | QFrame::Sunken );
   MsgText->setGeometry(WTextL.l,WTextL.t,WTextL.r,WTextL.b);

   TextEd = new QLineEdit( this, "TextEdit" );
   TextEd->setGeometry(WTextE.l,WTextE.t,WTextE.r,WTextE.b);
   TextEd->setMaxLength(80);
   TextEd->setText(DF->GetInfoText());
   connect( TextEd, SIGNAL(returnPressed()), SLOT(Accept()) );
   QToolTip::add( TextEd, "Enter text line" );

   //   MsgCID = new QLabel* [nCols];
   //   CIDEd  = new QLineEdit* [nCols];

   //   for(int i=0; i<nCols;i++)
   //      {sprintf(szAnswer, "Col ID%1d:",i+1);
   //    MsgCID[i]=new QLabel(szAnswer,this);

   //    MsgCID[i]->setAlignment(AlignLeft);
   //    MsgCID[i]->setFrameStyle( QFrame::Panel | QFrame::Sunken );
   //    MsgCID[i]->setGeometry(WColL.l+i*Offs,WColL.t,WColL.r,WColL.b);
       
   //    CIDEd[i]=new QLineEdit(this);
   //    CIDEd[i]->setGeometry(WColE.l+i*Offs,WColE.t,WColE.r,WColE.b);
   //    CIDEd[i]->setText(TD->GetColID(i+1));
   //   }
 }
  OkB=new QPushButton("Ok",this);
  OkB->setGeometry(WOk.l,WOk.t,WOk.r,WOk.b);
  OkB->setFixedHeight( OkB->sizeHint().height() );
  QToolTip::add( OkB, "Accept input" );
  connect( OkB, SIGNAL(clicked()), SLOT(Accept()));

  CaB=new QPushButton("Cancel",this);
  CaB->setGeometry(WCa.l,WCa.t,WCa.r,WCa.b);
  CaB->setFixedHeight( CaB->sizeHint().height() );
  QToolTip::add( CaB, "Cancel dialog" );
  connect( CaB, SIGNAL(clicked()), SLOT(reject()));

}
// *********************************************************
void SaveFile::Accept()
{
 if(TextW==0)
   { DF->SetInfoText(TextEd->text()); }
 else
   {int iL=DF->GetNoTxtL();
    int i;
    const char *t;
    for(i=0;i<TextW->numLines() && i<iL ;i++)
       {t=TextW->textLine(i);
        if(t[0])DF->SetInfoText(t,i); 
       }
   }
 //DF->SetDataFile(szFile);
 //char szB[MAX_CID_LEN+1];
 //Col IDs
 //for(int j=0;j<nCols;j++)
 //   {strncpy(szB,CIDEd[j]->text(),MAX_CID_LEN);
 //    szB[MAX_CID_LEN]=0;
 //    TD->SetColID(szB,j+1);
 //   }
 emit accept(); 
}
// ***************************************************
//
// Construct EnterNumber
//
// *********************************************************

EnterNumber::EnterNumber(QWidget *parent,  const char *text,
                         const double min, const double max,const int iT,
                         const char *szDefText):
             QDialog( parent,NULL,TRUE )
{GWind WMain= {100,         50,                   300,     120 };
 GWind WLab={10,            10,      280,        60};
 GWind WEdi={10,  WLab.t+WLab.b+10,   80,        20};

if(min>=max || !(iT==T_INT || iT==T_DBL))
  {String S;
   S.Setf("Illegal boundaries or type flag\n Min:%f Max:%f",min,max);
   QMessageBox::message("ERROR",(const char *)S);
   emit reject();
   return;
  }
 
 lfMin=min;
 lfMax=max;
 iType=iT;

 setCaption("Enter Value");

 setGeometry(WMain.l,WMain.t,WMain.r,WMain.b);

 String S;
 if(iType==T_INT)
   S.Setf("Enter a number (%d...%d)\n <ESC>:cancel <CR>:accept\n",
                     (int)min,(int)max);

 if(iType==T_DBL)
   S.Setf("Enter a number (%f...%f)\n <ESC>:cancel <CR>:accept\n",
                     min,max);

 S.Add(text);
 Msg=new QLabel((const char *)S, this);
 Msg->setAlignment(AlignLeft);
 Msg->setFrameStyle( QFrame::Panel | QFrame::Sunken );
 Msg->setGeometry(WLab.l,WLab.t,WLab.r,WLab.b);

 Ed = new QLineEdit( this );
 Ed->setGeometry(WEdi.l,WEdi.t,WEdi.r,WEdi.b);
 if(szDefText)
   {Ed->setText(szDefText);
    Ed->selectAll();
   }
 Ed->setFocus();
}

// *********************************************************
void EnterNumber::Accept()
{ 
 double a=atof(Ed->text());
 int i=atoi(Ed->text());
 
 if(iType==T_DBL)
   {if( a>=lfMin && a<=lfMax ){ emit accept(); return;} }
 
 if(iType==T_INT)
   {if( i>=(int)lfMin && i<=(int)lfMax ){ emit accept(); return;}}
 
 String S;
 S.Setf("Number out of range\n (%f...%f)",lfMin,lfMax);
 QMessageBox::message("ERROR",(const char *)S);
 
}
// *********************************************************
void EnterNumber::keyPressEvent(QKeyEvent *e)
{if(Ed->hasFocus())
  {
   //fprintf(stderr,"EnterNumber: Key event: %x\n",e->key());
   switch(e->key()) 
     {case Key_Enter:
     case Key_Return:
                     Accept();
                     break; 
     case Key_Escape:emit reject();
                     break;
	  default:
                   e->ignore();
	  }//switch
  }//if
}
// ***********************************************************

// ********************************************************
// **********************************************************
//
// Construct EnterPassword
//
// *********************************************************

EnterPassword::EnterPassword(QWidget *parent,  const char *text):
             QDialog( parent,NULL,TRUE )
{GWind WMain= {100,         50,                   300,     120 };
 GWind WLab={10,            10,      280,        60};
 GWind WEdi={10,  WLab.t+WLab.b+10,   80,        20};

 setCaption("Enter Password");

 setGeometry(WMain.l,WMain.t,WMain.r,WMain.b);

 String S;
 S.Setf("Enter Password\n <ESC>:cancel <CR>:accept\n");

 S.Add(text);
 Msg=new QLabel((const char *)S, this);
 Msg->setAlignment(AlignLeft);
 Msg->setFrameStyle( QFrame::Panel | QFrame::Sunken );
 Msg->setGeometry(WLab.l,WLab.t,WLab.r,WLab.b);

 Ed = new QLineEdit( this );
 Ed->setGeometry(WEdi.l,WEdi.t,WEdi.r,WEdi.b);
 Ed->setEchoMode(QLineEdit::Password);
 Ed->setFocus();
}

// *********************************************************
void EnterPassword::keyPressEvent(QKeyEvent *e)
{if(Ed->hasFocus())
  {switch(e->key()) 
    {case Key_Enter:
     case Key_Return:
                     emit accept();
                     break; 
     case Key_Escape:emit reject();
                     break;
	  default:
                   e->ignore();
	  }//switch
  }//if
}
// ***********************************************************




























