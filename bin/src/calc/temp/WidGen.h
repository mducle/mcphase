#ifndef WidGen_H
#define WidGen_H

#include <qlabel.h>
#include <qlineedit.h>
#include <qevent.h>
#include <qmultilinedit.h>
#include <qdialog.h>
#include <qpushbutton.h>
#include <qpainter.h>
#include <qfiledialog.h>

#ifndef CDATA_H
#include "cdata.h"
#endif

#ifndef DFILE_H
#include "dfile.h"
#endif

// $Id: WidGen.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: WidGen.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

// ***************************************************************
class TDataInfo : public QMultiLineEdit
{    
  public:
    TDataInfo(DataFile *df, const int x,const  int y,
              QWidget *parent=0, const char *name=0); 
    TDataInfo(CRData *td, const int x, const  int y,
              QWidget *parent, const char *name=0);
    TDataInfo(CRData *td,const int x, const  int y, const int w, const int h,
              QWidget *parent, const char *name);
};
// *****************************************************************
class SaveFile : public QDialog
{
    Q_OBJECT

private slots:
    void Accept();

private:
    TDataInfo * TDI;

       QLabel * MsgText;
    QLineEdit * TextEd;

    QMultiLineEdit *TextW;

//    QLabel **MsgCID;
//    QLineEdit **CIDEd;

    QPushButton *OkB;
    QPushButton *CaB;

    DataFile *DF;
    char szFile[MAXPATH+1];
  
public:
    SaveFile(const char *szFileName, DataFile *df,
             QWidget *parent=0, const char *name=0);

};
// ******************************************************************
#define T_INT 1
#define T_DBL 2

class EnterNumber : public QDialog
{
    Q_OBJECT

private slots:
    void Accept();

private:

    
       QLabel * Msg;
    QLineEdit * Ed;
    double lfMin,lfMax;
       int iType;

    void keyPressEvent(QKeyEvent *e);


public:
    EnterNumber(QWidget *parent, const char *text, 
                const double min, const double max, const int iF,
                const char *szDefText=NULL);
       int GetIntNo(void){return atoi(Ed->text());}
    double GetDoubleNo(void){return atof(Ed->text());}
};
// ******************************************************************
class EnterPassword : public QDialog
{
private:

    
       QLabel * Msg;
    QLineEdit * Ed;

    void keyPressEvent(QKeyEvent *e);


public:
    EnterPassword(QWidget *parent, const char *text);
    const char * GetText(void){return Ed->text();}
};
// ******************************************************************


#endif







