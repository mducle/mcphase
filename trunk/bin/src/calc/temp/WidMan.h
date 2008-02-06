//Widget Manager: WidMan.h
// $Id: WidMan.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: WidMan.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//


#ifndef WidMan_H
#define WidMan_H

#include <qapplication.h>
#include <qstring.h>
#include <qlabel.h>
#include <qlistbox.h>
#include <qlineedit.h>
#include <qevent.h>
#include <qmultilinedit.h>
#include <qpoint.h>
#include <qfiledialog.h>
#include <qradiobutton.h>
#include <qbuttongroup.h>
#include <qcheckbox.h>
#include <qwidget.h>

#include "filepar.h"

#include "ausw.h"
#include "stdinc.h"

void  SetWidgColor(QWidget *w, QColor c);
void  SetWidgBackColor(QWidget *w, QColor c);

#define MAX_OUTLINES 2000
#define MAX_ERRLINES 200
#define MAX_HISTLENGTH 100
#define MAN_PATH "AUSW_MANPATH"

class WidgetAusw : public QWidget
{
    Q_OBJECT
public:
      WidgetAusw(const char *szU,  QWidget *parent=0, const char *name=0 );
      int IsXauswCmd(const char *szT);
private slots:
    void QuitBClicked();
    void FileBClicked();
//    void HelpBClicked();
    void Filename(const char *szF);
    void CheckOutLines();
    void CheckErrLines();           
    void CheckHistLines();
    void NewType(int);
      
    void PrevCmdSelected(int);
    void CmdSelected(const char * szC); 
    void CmdReady();
    void CmdHelper(const char * szC);
//    void DataReceived();
private:
     QListBox * CmdLB;
    QLineEdit * CmdEd;
         char * szCmd;

    QMultiLineEdit * StdOut;
    QMultiLineEdit * StdErr;

    QMultiLineEdit * Helper;

    QListBox * HistLB;
    
    QRadioButton *TypeRB[3];
    int iType;

    QCheckBox *VerifyB;

    FilePar *FP;
 
    QFileDialog * FSel;
    
    char szDefPath[MAXPATH+1];
    char szHelpFile[MAXPATH+1];
    char szCmdFile[MAXPATH+1];
    char szUser[MAXUSER+1];
    char szHistory[MAXPATH+1];
    char *szManPath;
   
//    static QPoint Screen;

     void  keyPressEvent(QKeyEvent*);
      int  SaveHistory(const char * szFile);
      int  GetHistory(void);
      int  GetCommands(void);
     void  Append(QMultiLineEdit *w, const char *szF,...);
      int  HasOpt(const char *szT, const char *szFind);
      
     static QRect WTot;
     static QRect WCmdL;
     static QRect WCmdLB;
     static QRect WCmdLBL;
     static QRect WCmdEd;

     static QRect WQuitB;
     static QRect WFileB;
     static QRect WHelpB;
     
     static QRect WStdOL;
     static QRect WStdOut;
     static QRect WStdEL;
     static QRect WStdErr;

     static QRect WHelpL;
     static QRect WHelp;

     static QRect WHistL;
     static QRect WHistLB;

     static QRect TypeBW;          
     static const char *szType[];

     static QRect VerifyBW;
};


#endif //WidMan_H





