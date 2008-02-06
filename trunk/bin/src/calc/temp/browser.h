// File: browser.h
// $Id: browser.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: browser.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef BROWSER_H
#define BROWSER_H 1

#include <qlist.h>
#include <qdialog.h>
#include <qlineedit.h>
#include <qdialog.h>
#include <qpushbutton.h>
#include <qlistbox.h>

#include "stdinc.h"
#include "stdfunc.h"

#define MAX_NETBIOSN 60
#define    MAX_HOSTS 256

#define SERVER_LIST 0x0001u
#define  SHARE_LIST 0x0002u


class BrowseList
{
 char szHost[MAX_NETBIOSN];
 char szP[MAX_NETBIOSN];
 char szU[MAX_NETBIOSN];

 char NameList[MAX_HOSTS][MAX_NETBIOSN];
 int iNHosts;

public:
 
 BrowseList(const char *szH, const char *szu, const char *szp); 
 
 int GetNHosts(void) const {return iNHosts;}
 int GetError(void) const {return iNHosts==-1;}
 int ReqList(const unsigned What);
 int GetHostList( char *szL) const ;
 
 const char * GetHost(const int iH) const ;

};
// ******************************************************************
struct ShareItem
   {char szMount[MAXPATH+1];
    char szServer[MAXPATH+1];
    char szUser[MAXPATH+1];
    char szPass[MAXPATH+1];
   };

class MountShare : public QWidget
{
    Q_OBJECT

private slots:
    void Accept();
    void ServerSelected(int);
    void ShareSelected(int);
    void UnMount();

private:

    QLineEdit * Drive;
    QLineEdit * MountDir;
    QLineEdit * UserName;
    QLineEdit * PassWord;

    QListBox * Server;
    QListBox * Share;
    QListBox * CurMount;

    QPushButton * OkB;
    QPushButton * CaB;
    QPushButton * UnmountB;

    struct ShareItem Current;

    static QList<struct ShareItem> OpenShares; 
       int GetShareList(void);
      void UpdateOpenShareLB(void);
public:
    MountShare(QWidget *parent=0, const char *name=0);
           const struct ShareItem *GetCurrentItem(void){return &Current;}
    static const struct ShareItem* GetShare(const int iS);
    static          int GetNofShares(void);
    static         void Append(struct ShareItem *SI){OpenShares.append(SI);}
                    int CreateMount(struct ShareItem* SI);
};

// *******************************************************
int ReadTimeOut(int fd, String *CB,const int iTimeOut); 

#endif
