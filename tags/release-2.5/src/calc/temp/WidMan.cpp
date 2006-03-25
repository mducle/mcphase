// File: WidMan.cpp
//$Log: WidMan.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.3  1999/03/15 09:08:37  herbie
//*** empty log message ***
//
// Widget Manager WidMan.cpp

#include <qapplication.h>
#include <qpushbutton.h>
#include <qframe.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qtooltip.h>
#include <qfile.h>
#include <qdatastream.h>
#include <qmessagebox.h>
#include <qlabel.h>
#include <qevent.h>
#include <qkeycode.h>
#include <qmultilinedit.h>
#include <qfiledialog.h>
#include <qwidget.h>
#include <qsocketnotifier.h>

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#ifndef WidMan_H
#include "WidMan.h"
#endif

#ifndef STDFUNC_H
#include "stdfunc.h"
#endif

#ifndef BASICIO_H
#include "basicio.h"
#endif

const char *szRCSID={"$Id: WidMan.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"};
const char *szConfPath=getenv("AUSW_CONF_PATH");

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

//  QColorGroup nc( cg.foreground(), c,
//                  cg.light(), cg.dark(), cg.mid(),
//                  cg.text(), cg.base() );
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
// **********************************************************
//
// Class WidMan functions
//
// *********************************************************
// **********************************************************
// Global gives screen size
 QPoint SCREEN=QPoint(0,0);

// statics of WidgetAusw
 QRect WidgetAusw::WTot;
 QRect WidgetAusw::WCmdL;
 QRect WidgetAusw::WCmdLB;
 QRect WidgetAusw::WCmdLBL;
 QRect WidgetAusw::WCmdEd;  

 QRect WidgetAusw::WQuitB;
 QRect WidgetAusw::WFileB;  
 QRect WidgetAusw::WHelpB;
  
 QRect WidgetAusw::WStdOL;
 QRect WidgetAusw::WStdOut;
 QRect WidgetAusw::WStdEL;
 QRect WidgetAusw::WStdErr;  

 QRect WidgetAusw::TypeBW;

 QRect WidgetAusw::WHelpL;
 QRect WidgetAusw::WHelp;

 QRect WidgetAusw::WHistL;
 QRect WidgetAusw::WHistLB;

 const char *WidgetAusw::szType[]={" -t default", " -t dos "," -t unix "};
 
 QRect WidgetAusw::VerifyBW;
 
WidgetAusw::WidgetAusw(const char *szU,  QWidget *parent, const char *name )
    : QWidget( parent, name )
{
  // Position of windows
  // Total

  const int w=SCREEN.x()*95/100;
  const int h=SCREEN.y()*90/100;

  const int XOfs0=10;
  const int w3=(w-4*XOfs0)/3;

  const int XOfs1=2*XOfs0+w3;
  const int XOfs2=XOfs1+XOfs0+w3;

  const int YOfs0=5;
//  const int h3=(h-4*YOfs0)/3;    
  const int ButH=h/30; 
  const int ButW=60;
  
     WTot.setRect(5,5,w,h); 
   WCmdEd.setRect(XOfs0,  WTot.height()-ButH-10, 2*w3+XOfs0, 3*ButH/2);
    WCmdL.setRect(XOfs0,  WCmdEd.y()-ButH, w3, ButH);

   WQuitB.setRect(XOfs0,  WCmdL.y()-ButH-3*YOfs0, ButW, ButH);    
   WFileB.setRect(2*XOfs0+ButW, WQuitB.y(), ButW, ButH);    
//    WHelpB.setRect(3*XOfs0+2*ButW, WQuitB.y(), ButW, ButH);    

    TypeBW.setRect(XOfs0,  WQuitB.y()-4*ButH-2*YOfs0, w/8, 4*ButH);
 VerifyBW.setRect(TypeBW.x()+TypeBW.width()+5, TypeBW.y(), w/8, 2*ButH);

  WCmdLBL.setRect(XOfs0,  YOfs0/2, w/4, ButH);
   WCmdLB.setRect(XOfs0,  WCmdLBL.y()+WCmdLBL.height()+YOfs0/2,
                  w/4, TypeBW.y()-WCmdLBL.height()-3*YOfs0);

//    WStdOL.setRect(XOfs2,  YOfs0, w3, ButH);
//   WStdOut.setRect(XOfs2,  WStdOL.y()+WStdOL.height()+YOfs0,
//                   w3, 2*h3-WStdOL.height());

//    WStdEL.setRect(XOfs2,  WStdOut.y()+WStdOut.height()+YOfs0, w3, ButH);
//   WStdErr.setRect(XOfs2,  WStdEL.y()+WStdEL.height()+YOfs0,
//                   w3, h3-WStdEL.height());
//   WHelpL.setRect(WCmdLBL.x()+WCmdLBL.width()+XOfs0, WCmdLBL.y(),
//                  w-w3-w/4-4*XOfs0,ButH);
//    WHelp.setRect(WHelpL.x(),WCmdLB.y(),WHelpL.width(),2*WStdOut.height()/3);
// ************
   WHelpL.setRect(WCmdLBL.x()+WCmdLBL.width()+XOfs0, WCmdLBL.y(),
                  w-w/4-3*XOfs0,ButH);
    WHelp.setRect(WHelpL.x(),WCmdLB.y(),WHelpL.width(),2*WCmdLB.height()/3);
// ************

    WStdOL.setRect(XOfs2,  WHelp.y()+WHelp.height()+YOfs0/2, w3, ButH);
   WStdOut.setRect(XOfs2,  WStdOL.y()+WStdOL.height()+YOfs0/2,
                   w3, h-WHelp.y()-WHelp.height()-WStdOL.height());

//   WHistL.setRect(WHelp.x(), WHelp.y()+WHelp.height()+2*YOfs0,
//                  WHelp.width(),ButH);
//  WHistLB.setRect(WHistL.x(),WHistL.y()+WHistL.height()+2*YOfs0,
//                  WHelpL.width(),WCmdLB.y()+WCmdLB.height()-WHistL.y()-WHistL.height()-2*YOfs0);
// *************
   WHistL.setRect(WHelp.x(), WHelp.y()+WHelp.height()+YOfs0/2,
                  w-w3-w/4-4*XOfs0,ButH);
  WHistLB.setRect(WHistL.x(),WHistL.y()+WHistL.height()+YOfs0/2,
                  WHistL.width(),WCmdLB.y()+WCmdLB.height()-WHistL.y()-WHistL.height()-YOfs0);
// **************
    WStdEL.setRect(WHistL.x(),  WHistLB.y()+WHistLB.height()+YOfs0/2,
                   WHistLB.width(), ButH);
   WStdErr.setRect(WStdEL.x(),  WStdEL.y()+WStdEL.height()+YOfs0/2,
                   WStdEL.width(), h-WStdEL.y()-WStdEL.height()-2*YOfs0-WCmdEd.height());


//    TypeBW.setRect(XOfs0,  WCmdLB.y()+WCmdLB.height()+YOfs0, w/7, 4*ButH);
//  VerifyBW.setRect(TypeBW.x()+TypeBW.width()+5, TypeBW.y(), w/7, 2*ButH);

   const int iButFontSz=(SCREEN.y()==600 ? 12 : 14);
   const int iLabFontSz=12;
   const int iWriFontSz=14;

   FSel=0;

  szManPath=getenv(MAN_PATH);
 //Set the window caption/title

   char szB[MAX_LINELENGTH+1];
   sprintf(szB,"HMs XAusw 2.0 %s",__DATE__);
   setCaption(szB);
   setGeometry(WTot);

   AppendPath(szConfPath,PARAMETER_FILE,szB,MAXPATH);
   FP=new FilePar(szB);
   if(FP==0 || FP->GetStatus()!=FOUND)
     {PRINT_DEBUG("Error opening %s\a\n",PARAMETER_FILE)
      QMessageBox::critical(0,"XAusw",QString("File ausw.ini not found\n")+
                                           "XAusw exits now!"); 
     }

   FP->GetPar("Files","DefDir","%s",szDefPath);
   FP->GetPar("Files","HelpFile","%s",szHelpFile);

   szUser[0]=0;
   szHistory[0]=0;

   if(szU)strncpy(szUser,szU,MAXUSER);
   else sprintf(szUser,"XAusw");
   sprintf(szHistory,"%s.history",szUser);

// Create Buttons
//-----------------------READ BUTTON:
    QPushButton *QuitB;
    QuitB = new QPushButton(this, "quit" );
    QuitB->setFont(QFont("Helvetica", iButFontSz,QFont::DemiBold));
    QuitB->setText("Quit");
    QuitB->setFixedHeight( QuitB->sizeHint().height() );
    QuitB->setGeometry(WQuitB);
    connect( QuitB, SIGNAL(clicked()),SLOT(QuitBClicked()) );
    QToolTip::add( QuitB, "Exit XAusw" );

//-----------------------FILE BUTTON:
    QPushButton *FileB;
    FileB = new QPushButton("File", this, "help" );
    FileB->setFont(QFont("Helvetica", iButFontSz,QFont::DemiBold));
    FileB->setFixedHeight( FileB->sizeHint().height() );
    FileB->setGeometry(WFileB);
    connect( FileB, SIGNAL(clicked()), SLOT(FileBClicked()) );
    QToolTip::add( FileB, "Display File screen" );

//-----------------------HELP BUTTON:
//     QPushButton *HelpB;
//     HelpB = new QPushButton("Help", this, "help" );
//     HelpB->setFont(QFont("Helvetica", iButFontSz,QFont::DemiBold));
//     HelpB->setFixedHeight( HelpB->sizeHint().height() );
//     HelpB->setGeometry(WHelpB);
//     connect( HelpB, SIGNAL(clicked()), SLOT(HelpBClicked()) );
//     QToolTip::add( HelpB, "Display Help screen" );

//-----------------------LABEL CommandLineEdit:
    QLabel * CmdL = new QLabel( this, "CmdLabel" );
    CmdL->setText( "Command line:" );
    CmdL->setFont(QFont("Helvetica",iLabFontSz,QFont::DemiBold));
    CmdL->setGeometry(WCmdL);

//----------------------LINE EDIT COMMAND:
// Create a single line edit for commands
    CmdEd = new QLineEdit( this, "lineEdit" );
    CmdEd->setFont(QFont("Helvetica",iWriFontSz,QFont::Normal));
    CmdEd->setGeometry(WCmdEd);
    connect( CmdEd, SIGNAL(returnPressed()), SLOT(CmdReady()) );
    QToolTip::add( CmdEd, "Command line" );
    CmdEd->setMaxLength(MAX_LINELENGTH);
    
//-----------------------LABEL availableCommandListBox:
    QLabel * CmdLBL = new QLabel( this, "ACmdLabel" );
    CmdLBL->setText( "Available commands:" );
    CmdLBL->setFont(QFont("Helvetica",iLabFontSz,QFont::DemiBold));
    CmdLBL->setGeometry(WCmdLBL);

//-----------------------LISTBOX of available COMMANDS:
    CmdLB = new QListBox( this, "CmdLB" );
    CmdLB->setGeometry(WCmdLB);
    CmdLB->setFont(QFont("Helvetica",iLabFontSz,QFont::Normal));
    connect( CmdLB, SIGNAL(selected(const char *)),
                    SLOT(CmdSelected(const char *)) );
    connect( CmdLB, SIGNAL(highlighted(const char *)),
                    SLOT(CmdHelper(const char *)) );
    QToolTip::add( CmdLB, "Available commands" );
    if(!GetCommands())
      {PRINT_DEBUG("Error opening %s\a\n",szCmdFile)
      QMessageBox::critical(0,"XAusw",QString("Command list not found\n")+
                                           "XAusw exits now!"); 
     }

//-----------------------LABEL StdOut:
    QLabel * StdOL = new QLabel( this, "StdOut" );
    StdOL->setText( "Output of commands:" );
    StdOL->setFont(QFont("Helvetica",iLabFontSz,QFont::DemiBold));
    StdOL->setGeometry(WStdOL);

    QColor bc;
//---------------------MULTI LINE EDIT STDOUT:
    StdOut=new QMultiLineEdit(this,"CmdOut");
    StdOut->setBackgroundColor(white);
    StdOut->setReadOnly(TRUE);
    // CmdOut->setFrameStyle(QFrame::Panel|QFrame::Plain);
    StdOut->setGeometry(WStdOut);
    StdOut->setFont(QFont("Courier",12,QFont::Normal));
    bc.setNamedColor("#FFFFE1");
    SetWidgBackColor(StdOut, bc );
    QToolTip::add(StdOut,"Output of command");
    connect( StdOut, SIGNAL(textChanged()), SLOT(CheckOutLines()) );
//-----------------------LABEL StdErr:
    QLabel * StdEL = new QLabel( this, "StdErr" );
    StdEL->setText( "Errors:" );
    StdEL->setFont(QFont("Helvetica",iLabFontSz,QFont::DemiBold));
    StdEL->setGeometry(WStdEL);


//---------------------MULTI LINE EDIT STDerr:
    StdErr=new QMultiLineEdit(this,"CmdErr");
    StdErr->setReadOnly(TRUE);
    StdErr->setGeometry(WStdErr);
    QToolTip::add(StdErr,"Error messages");
    //bc.setNamedColor("LightGoldenrod1");
    bc.setNamedColor("#FFE6D5");
    SetWidgBackColor(StdErr, bc );
    connect( StdErr, SIGNAL(textChanged()), SLOT(CheckErrLines()) );

//-----------------------LABEL Helper:
    QLabel * HelpL = new QLabel( this, "Help" );
    HelpL->setText( "Help:" );
    HelpL->setFont(QFont("Helvetica",iLabFontSz,QFont::DemiBold));
    HelpL->setGeometry(WHelpL);

//---------------------MULTI LINE EDIT HELPER:
    Helper=new QMultiLineEdit(this,"HlpCmd");
    Helper->setBackgroundColor(white);
    Helper->setReadOnly(TRUE);
    Helper->setGeometry(WHelp);
    bc.setNamedColor("#DFE8DF");
    SetWidgBackColor(Helper, bc );
    Helper->setFont(QFont("Courier",12,QFont::Normal));
    QToolTip::add(Helper,"Help message");

//-----------------------LABEL CommandHistoryListBox:
    QLabel * HistL = new QLabel( this, "HistLab" );
    HistL->setText( "Command history:" );
    HistL->setFont(QFont("Helvetica",iLabFontSz,QFont::DemiBold));
    HistL->setGeometry(WHistL);

//-----------------------LISTBOX of Command History:
    HistLB = new QListBox( this, "CmdHist" );
    HistLB->setGeometry(WHistLB);
    HistLB->setFont(QFont("Helvetica",iLabFontSz,QFont::Normal));
    connect( HistLB, SIGNAL(selected(int)), SLOT(PrevCmdSelected(int)) );
    QToolTip::add( HistLB, "Previous Commands" );
    if(!GetHistory())
      {PRINT_DEBUG("History file %s not found\n",szHistory)
       QMessageBox::message("XAusw Warning", "History not available"); 
      }
 
// -------------------Radio Buttons for type (-t)

    QButtonGroup *TypeBox = new QButtonGroup( "Type (-t)", this, "TypeBox" );

    TypeBox->setGeometry(TypeBW);
    connect( TypeBox, SIGNAL(clicked(int)), SLOT(NewType(int)) );
    int i;
    for( i = 0 ; i < 3 ; i++ )
       {sprintf(szB,"TypeB%i",i);
        TypeRB[i] = new QRadioButton( TypeBox, szB );
        TypeRB[i]->setGeometry( XOfs0, 3*YOfs0+i*ButH ,4*ButW/3, ButH );
        TypeRB[i]->setText( WidgetAusw::szType[i]+1 );
    }
    TypeRB[0]->setChecked( TRUE );    
    iType=1;

// -------------------Radio Buttons for verify (-t)
    QButtonGroup *VerifyBox = new QButtonGroup( "Verify (-v)", this, "VerifyBox" );

    VerifyBox->setGeometry(VerifyBW);

    VerifyB = new QCheckBox( VerifyBox, "set -v" );
    VerifyB->setGeometry( XOfs0, 3*YOfs0 , 4*ButW/3, ButH );

//    QSocketNotifier *sn=new QSocketNotifier(1,QSocketNotifier::Read, parent );
//    QObject::connect( sn, SIGNAL(activated(int)), SLOT(DataReceived()) );

 show();
  
 CmdEd->setFocus();    

 FSel=new QFileDialog(".");
 connect( FSel, SIGNAL(fileSelected(const char *)), SLOT(Filename(const char *)) );

}
// *************************************************************
// void WidgetAusw::HelpBClicked()
// {
//  fprintf(stdout,"Help clicked\n");
// }
// ************************************************************
// void WidgetAusw::DataReceived()
// {char szC[1025];
//  while(fgets(szC,1024,stdout))
//  fprintf(stderr,"Data:%s",szC);
// }
// ************************************************************
void WidgetAusw::FileBClicked()
{
 FSel->show();
}
// *************************************************************
void WidgetAusw::CmdSelected(const char *szC)
{
 if(strstr(szC,"***"))
   {QApplication::beep();
    Append(StdErr,"XAusw: Illegal selection %s",szC);
   } 
 else
 { char szB[MAX_LINELENGTH+1];
   strncpy(szB,szC,MAX_LINELENGTH);
   if(szB[strlen(szB)-1]=='\n')szB[strlen(szB)-1]=0;
   CmdEd->setText(szB);
   CmdEd->setFocus(); 
 }
}
// *************************************************************
void WidgetAusw::CmdHelper(const char *szC)
{
 if(!strstr(szC,"***"))
   {char szB[MAX_LINELENGTH+1];
    strncpy(szB,szC,MAX_LINELENGTH);
    if(szB[strlen(szB)-1]=='\n')szB[strlen(szB)-1]=0;
    char *p=strchr(szB,' ');
    if(p)*p=0;
    String B;
    Usage(szManPath,szB,B);
    Helper->setText(B.GetBuf());
   }
}
// *************************************************************
void WidgetAusw::Filename(const char *szF)
{
 const char *p=CmdEd->text();
 char *s=new char [strlen(szF) + strlen(p) +1];
 CHECK_POINTER_RET(s);
 *s=0;
 strcat(s,p);
 strcat(s,szF);

 //StdErr->insertLine(s);
 CmdEd->setText(s);
 delete s;
}
// *************************************************************
void WidgetAusw::CheckOutLines()
{if(StdOut->numLines()>MAX_OUTLINES)
   {int i;
    int iL=StdOut->numLines();
    for (i=1; i<iL-MAX_OUTLINES+50; i++)StdOut->removeLine(i);
    Append(StdErr,"XAusw :%d of %d lines removed from StdOut",i-1,iL);
   }
}
// *************************************************************
void WidgetAusw::CheckErrLines()
{if(StdErr->numLines()>MAX_ERRLINES)
   {int i;
    int iL=StdErr->numLines();
    for (i=1; i<iL-MAX_ERRLINES+10; i++) StdErr->removeLine(i);
    Append(StdErr,"XAusw: %d of %d lines removed from StdErr",i-1,iL);
   }
}
// **************************************************************
void WidgetAusw::CheckHistLines()
{if(HistLB->count()>MAX_HISTLENGTH)
   {HistLB->removeItem(0);
    Append(StdErr,"XAusw: Remove item from HistLB");
   }
}
// **************************************************************
void WidgetAusw::QuitBClicked()
{
//    QMessageBox::message("WINDOWS Error 4711",
//                         "Type any key to stop\nor any other key to continue ...");
 if( !SaveHistory(szHistory) )
   {PRINT_DEBUG("Can not save command history %s\n",szHistory)
   }
 qApp->quit();
}
// ***********************************************************
void WidgetAusw::CmdReady( )
{   char szC[MAX_LINELENGTH+1];
    strncpy(szC,CmdEd->text(),MAX_LINELENGTH);
    szC[MAX_LINELENGTH]=0;

    if(iType>1)
       {char *p=strchr(szC,' ');
        if(p)*p=0;
        if(IsXauswCmd(szC) && !HasOpt(szC,"-t"))
          {String cb(szC);
	   cb.Add(szType[iType-1]);
	   if(p)cb.Add(p+1);
	   strcpy(szC,cb.GetBuf());
          }
	 if(p)*p=' ';  
         CmdEd->setText(szC);
        }

    if(VerifyB->isChecked())
      {char *p=strchr(szC,' ');
        if(p)*p=0;
        if(IsXauswCmd(szC) && !HasOpt(szC,"-v"))
          {String cb(szC);
	   cb.Add(" -v ");
	   if(p)cb.Add(p+1);
	   strcpy(szC,cb.GetBuf());
          }
	 if(p)*p=' ';  
       CmdEd->setText(szC);
      }

    if(szC[strlen(szC)-1]=='\n')szC[strlen(szC)-1]=0;
    if(szC && szC[0]!=0)
      {StdOut->append(szC);
       struct FILE2 f=popen2("/bin/sh","/bin/sh","-c",szC,NULL);
       if(!f.sto || !f.ste)
         {Append(StdOut,"XAusw: Error open pipe");
          PRINT_DEBUG("Error open pipe")
	  return;
	 }
    setCursor(waitCursor);
    String cB;
    while(fgets(szC,MAX_LINELENGTH,f.sto))cB.Add(szC);
    StdOut->append(cB.GetBuf()); 	  
    StdOut->setCursorPosition(StdOut->numLines(),0);	    

    cB.Set("");
    while(fgets(szC,MAX_LINELENGTH,f.ste))cB.Add(szC);
    StdErr->append(cB.GetBuf()); 	  
    StdErr->setCursorPosition(StdErr->numLines(),0);	    

    pclose2(f);
    HistLB->insertItem(CmdEd->text());
    HistLB->setTopItem(HistLB->count()-HistLB->numItemsVisible()+1);
    CheckHistLines();
    setCursor(arrowCursor);
   }
}
// **********************************************************
void WidgetAusw::NewType(int iT)
{iType=iT+1;
 Append(StdErr,"XAusw: Type set to %s",szType[iT] );
}
// **********************************************************
void WidgetAusw::PrevCmdSelected( int index )
{
    const char *szItem=HistLB->text(index);
    HistLB->insertItem(szItem);
    CmdEd->setText(szItem);
    CmdEd->setFocus();
}
// ***********************************************************
void WidgetAusw::keyPressEvent(QKeyEvent *e)
{if(CmdEd->hasFocus())
  {
   switch(e->key()) 
	  {
// 	   case Key_Up:if(nCmds>0 && iCmd>0)
// 	                 {iCmd--;
//                           CmdEd->setText(CmdLB->text(iCmd));
//                          }//if
//                        else{//CmdOut->append("Up out of range");
// 		            QApplication::beep();
// 			   }
//                        break;
//           case Key_Down:if(iCmd<nCmds)
// 	                  {iCmd++;
//                            CmdEd->setText(CmdLB->text(iCmd));
//                           }//if
//                        else{//CmdOut->append("Down out of range");
// 		            QApplication::beep();
// 			   }  
//                        break;
	  case Key_PageUp:{char szB[2048];
                           int i=GetFiles(CmdEd->text(),szB,2047);
                           //fprintf(stderr,"%d: %s\n",i,szB);
                           if(i==1){CmdEd->setText(szB);
			            StdErr->append(" ");
			            }
                           else StdErr->append(szB);
                           StdErr->setCursorPosition(StdErr->numLines(),0);	    
                           break; 
                           }
// 	  case Key_F1:if(bCmdInfo)
// 	                {bCmdInfo=0;
// 			 CmdH->hide();
// 			} 
// 	               else
// 		        {bCmdInfo=1;
//                          CmdH = new CmdHelp(WTot);
//                          CmdH->show();
// 			 //setFocusPolicy(QWidget::StrongFocus);
// 			 //fprintf(stderr,"Focus enabled: %d\n",isFocusEnabled());
// 			 //setActiveWindow();
// 			 CmdEd->setFocus();
// 			}
//                        break; 
            default://fprintf(stderr,"Key event: %x\n",e->key());
                   e->ignore();
	  }//switch
  }//if
}
// **********************************************************
int WidgetAusw::GetCommands(void)
{if(FP->GetPar("Files","CommandList","%s",szCmdFile)!=FOUND)return 0;
 FILE *f=fopen(szCmdFile,"r");
 if(!f)return 0;
 char szB[MAX_LINELENGTH+1];
 while(fgets(szB,MAX_LINELENGTH,f))
      {szB[MAX_LINELENGTH]=0;
       CmdLB->insertItem(szB);
      }
 fclose(f);     
 if(CmdLB->count()==0)return 0;
 return 1;      
}
// **********************************************************
int WidgetAusw::GetHistory(void)
{FILE *f=fopen(szHistory,"r");
 if(!f)return 0;
 char szB[MAX_LINELENGTH+1];
 while(fgets(szB,MAX_LINELENGTH,f))
      {szB[MAX_LINELENGTH]=0;
       HistLB->insertItem(szB);
      }
 fclose(f);     
 if(HistLB->count()==0)return 0;
 return 1;      
}
// **********************************************************************
int WidgetAusw::SaveHistory(const char *szFile)
{
 if(HistLB->count()==0)return 0;
 FILE *f;

 f=fopen(szFile,"w");

 CHECK_FILE_POINTER_RETURN(f,szFile,0)
 char szB[MAX_LINELENGTH+1];
 for(unsigned i=0; i<HistLB->count(); i++)
    {strncpy(szB,HistLB->text(i),MAX_LINELENGTH);
     RemoveCR_LF(szB);
     fprintf(f,"%s\n",szB);
    }
 fclose(f);
 return 1;
}
// **********************************************************************
void WidgetAusw::Append(QMultiLineEdit *w, const char *szF,...)
{char szB[MAX_LINELENGTH+1];
 va_list arg_ptr;
 va_start(arg_ptr, szF);
 vsnprintf(szB,MAX_LINELENGTH,szF,arg_ptr);
 szB[MAX_LINELENGTH]=0;
 va_end(arg_ptr);
 w->append(szB);
}
// ***********************************************************
int WidgetAusw::IsXauswCmd(const char * szT)
{unsigned i;
 const char *p;
 for(i=0; i<CmdLB->count();i++)
    {p=CmdLB->text(i);
     if(!strchr(p,'*'))
       {if(strstr(p,szT)==p)return 1;
       }
     }  
 return 0;
}
// **********************************************************
int WidgetAusw::HasOpt(const char *szT, const char *szFind)
{const char *p=strchr(szT,'|');
 const char *f=strstr(szT,szFind);
 if(!f)return 0;

 if(p){if(f>p)return 0;}

 return 1; 
}
// **********************************************************
//
// Create and display WidgetAusw.
//

int main( int argc, char **argv )
{
  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication a( argc, argv );
  QApplication::setFont( QFont("Helvetica") );
  
  char *b=NULL;
  if(a.argc()==2)b=a.argv()[1];
  else
    {if(a.argc()!=1)
      {fprintf(stderr,"Illegal Parameter\nUsage: %s [username]\n",a.argv()[0]);
       exit(EXIT_FAILURE);
      }
    }
  if(b)
    {if(strlen(b)>MAXUSER)b[MAXUSER]=0;}

  QWidget *d=QApplication::desktop();
  SCREEN.setX(d->width());
  SCREEN.setY(d->height());

    WidgetAusw w(b);
    a.setMainWidget( &w );
    w.show();
    return a.exec();
}




