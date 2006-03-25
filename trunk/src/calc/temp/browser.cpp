//File: browser.cpp
//$Id: browser.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
//$Log: browser.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.2  1999/03/15 09:08:37  herbie
//*** empty log message ***
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <sys/time.h>

#include <qlist.h>
#include <qlabel.h>
#include <qtooltip.h>
#include <qmessagebox.h>
#include <qlistbox.h>
#include <qlineedit.h>
#include <qapplication.h>

#include "browser.h"
#include "basicio.h"
#include "stdfunc.h"
#include "strings.h"

//#define BROWSER_TEST

#ifndef BROWSER_TEST

 QList<struct ShareItem> MountShare::OpenShares; 

// extern char szAnswer[2048];
#endif

char NoHost[]={"Browser: Index out of range!"};


int ReadTimeOut(int fd,String *CB,const int iTimeOut)
{const int BSize=2048;
 char szCmd[BSize+1];

 int iCnt,iTim=0;
  while(1)
       {iCnt=read(fd,szCmd,BSize);
        if(iCnt==0)
          {sleep(1);
           iTim++;
           if(iTim<iTimeOut)continue;
           else break;
          }
        if(iCnt<0)
          {PRINT_DEBUG("read: %s\n",strerror(errno))
	   break;
          }
        CB->Add(szCmd,iCnt); 
        if(iCnt<BSize)break;
       }
 return iCnt;
}
// ***************************************************
// class Browser functions
// ***************************************************
BrowseList::BrowseList(const char *szH, const char *szu, const char *szp)
{iNHosts=-1;
 szHost[0]=szU[0]=szP[0]=0;
 if(szH)strncpy(szHost,szH,MAX_NETBIOSN);
 if(szu)strncpy(szU,szu,MAX_NETBIOSN); 
 if(szp)strncpy(szP,szp,MAX_NETBIOSN); 
 if(! (szH[0] && szU[0] && szP[0]) )return;

 iNHosts=0;
}
// **************************************************
int BrowseList::ReqList(const unsigned What)
{const int BSize=2048;
 char szCmd[BSize+1];

 //String *Buf=new CharBuf();

//char *Argv[]={"smbclient","-L",szHost,"-U",szU,"-d 1",NULL};
char *Argv[]={"smbclient","-L",szHost,"-U",szU,"-Wxphys","-d1",NULL};

 int fd=fd_popen("smbclient",Argv);

 if(fd<0){PRINT_DEBUG("Error opening pipe %s\n",Argv[0])
          return 0;
         }
 
 fd_set rfds;
 //fd_set wfds;
 struct timeval tv;

 /* Watch fd to see when it has input. */
 FD_ZERO(&rfds);
 FD_SET(fd,&rfds);

 /* Wait up to 20 seconds. */
 tv.tv_sec = 20;
 tv.tv_usec = 0;


 int retval = select(FD_SETSIZE,&rfds, NULL, NULL, &tv);
  /* Don't rely on the value of tv now! */

 if (!retval)
    {PRINT_DEBUG("No data within timeout.\n")
     return 0;
    }
// sleep(5);

 FILE *f=fdopen(fd,"r");
 FILE *w=fdopen(fd,"w");
 int iCnt;
 for(iCnt=0; iCnt<BSize;iCnt++)
    {if( (szCmd[iCnt]=fgetc(f)) == EOF)
       {PRINT_DEBUG("Unexpected end of input stream\n")
        return 0;
       }
      szCmd[iCnt+1]=0; 
      if(strstr(szCmd,"Password:")){fputs("xfi\n",w);break;}
      if(strstr(szCmd,"Domain="))break;
      //putc(szCmd[iCnt],stdout);
     }

 iCnt=0;
  
 int SelServer=0, SelShare=0; 

 while(fgets(szCmd,BSize,f)!=NULL)
   {if(What==SERVER_LIST)
      {if(strstr(szCmd,"This machine has a browse list:"))
	  {for(int i=0;i<3;i++)if(fgets(szCmd,BSize,f)==NULL)break;
           SelServer=1;
           continue;
          }
       if(SelServer==1)
         {char *s=szCmd,*b;
          if(szCmd[0]=='\n' || szCmd[0]=='\r')break;
          while(strchr(" \t",*s))s++;
          b=strtok(s," \t");
          if(b && iCnt < MAX_HOSTS && 
            (!(*b=='\n' || *b=='\r' || *b==0)) )
                             strncpy(NameList[iCnt],b,MAX_NETBIOSN);
          else {iNHosts=iCnt; break;}
          iCnt++;
         }
       }
    if(What==SHARE_LIST)
      {if(strstr(szCmd,"Sharename"))
	  {if(fgets(szCmd,BSize,f)==NULL)break;
           SelShare=1;
           continue;
          }
       if(SelShare==1)
         {char *s=szCmd,*b;
          if(szCmd[0]=='\n' || szCmd[0]=='\r')break;
          if(strstr(szCmd,"Disk")==NULL)continue;
          while(strchr(" \t",*s))s++;
          b=strtok(s," \t");
          if(strchr(b,'$'))continue;
          if(b && iCnt < MAX_HOSTS && 
            (!(*b=='\n' || *b=='\r' || *b==0)) )
                             strncpy(NameList[iCnt],b,MAX_NETBIOSN);
          else {iNHosts=iCnt; break;}
          iCnt++;
         }
       }
   }
 iNHosts=iCnt;

 fclose (f);
 fclose (w);
 return iCnt;
}
// ****************************************************
const char *BrowseList::GetHost(const int iH) const
{if(iH>=0 && iH<iNHosts)return NameList[iH];
 else return NoHost;
}
// ****************************************************
int BrowseList::GetHostList( char *szL) const
{int i;

 if(iNHosts<=0) return 0;
 char *sp=szL;

 for(i=0;i<iNHosts; i++)sp+=sprintf(sp,"%s\n",NameList[i]);
 
 return iNHosts;
}

#ifndef BROWSER_TEST
// *****************************************************
// MountShare Functions
// *****************************************************
MountShare::MountShare(QWidget *parent=0, const char *name=0) :
            QWidget(parent,name)

{
  setCursor(waitCursor);

  Current.szServer[0]=0;
  Current.szMount[0]=0;
  Current.szUser[0]=0;
  Current.szPass[0]=0;

 
  const QRect WTot(20,20,650,520); 
  const int ButH=20;
  const int ButW=70;

  int WCent=WTot.width()/2;

//  const QRect WOkB(WCent,               WTot.height()-ButH-20,
//                   ButW,             ButH);  
  const QRect WCaB(WTot.width()-3*ButW/2-10, WTot.height()-ButH-20,
                   3*ButW/2,             ButH); 

  const QRect  WSeL(20,       20,          100, ButH);
  const QRect WSeLB(WSeL.x(), WSeL.y()+WSeL.height()+5, 160, 200);

  const QRect WShL(30+WSeL.x()+WSeLB.width(), WSeL.y(),
                   WSeL.width(),              WSeL.height());
  const QRect WShLB(WShL.x(),                 WSeLB.y(),
                    WSeLB.width(),            WSeLB.height()); 

  const QRect WDrL(WSeLB.x(),  WSeLB.y()+WSeLB.height()+25,
                    WCent-40,   ButH);
  const QRect WDrE(WDrL.x(),   WDrL.y()+WDrL.height()+5,
               WCent-40,        ButH);

  const QRect WCMoL(WShLB.x()+WShLB.width()+40,  WShL.y(),
                    WShL.width(),   ButH);
  const QRect WCMoLB(WCMoL.x(),   WShLB.y(),
               WShLB.width()+40,     WShLB.height());
  const QRect WUnMB(WCMoLB.x(), WCMoLB.y()+WCMoLB.height()+10,
                   3*ButW/2,             ButH); 


  const QRect WMoL(WSeL.x(),     WDrE.y()+WDrE.height()+10,
                   WDrL.width(), WDrL.height());
  const QRect WMoE(WMoL.x(),     WMoL.y()+WMoL.height()+5,
                   WDrE.width(), WDrE.height());

  const QRect WUnL(WSeLB.x(),     WMoE.y()+WMoE.height()+10,
                    WCent-40,     ButH);
  const QRect WUnE(WSeLB.x(),     WUnL.y()+WUnL.height()+5,
                   WSeLB.width(), ButH );

  const QRect WPwL(WUnL.x(),      WUnE.y()+WUnE.height()+10,
                   WSeL.width(),  ButH);
  const QRect WPwE(WUnL.x(),      WPwL.y()+WPwL.height()+5,
                   WUnE.width(),  ButH);
  const QRect WOkB(WPwE.x()+WPwE.width()+20, WUnE.y(),
                   3*ButW/2,             ButH);  

 const char szDefServer[]="e1316b";
 const char szDefUser[]="browser";
 const char szDefPw[]="xfi";
 
 setCaption("Mount shared drive");
 setGeometry(WTot);

 QLabel *CMountL;
 CurMount=0;
 OpenShares.setAutoDelete(TRUE);
 GetShareList();
 if(OpenShares.count()>0)
   {CMountL = new QLabel( this, "CMountLabel" );
    CMountL->setText( "Open mounts:" );
    CMountL->setFont(QFont("Helvetica",12,QFont::DemiBold));
    CMountL->setGeometry(WCMoL);

    CurMount = new QListBox( this, "CMount" );
    CurMount->setGeometry(WCMoLB);
    CurMount->setFont(QFont("Helvetica",12,QFont::Normal));
    UpdateOpenShareLB();
    QToolTip::add( CurMount, "List of open mounts" );

    UnmountB=new QPushButton("Unmount",this);
    UnmountB->setFont(QFont("Helvetica",14,QFont::DemiBold));
    UnmountB->setGeometry(WUnMB);
    UnmountB->setFixedHeight( UnmountB->sizeHint().height() );
    QToolTip::add( UnmountB, "Unmount selected item" );
    connect( UnmountB, SIGNAL(clicked()), SLOT(UnMount()));
   }

 QLabel *ServerL = new QLabel( this, "ServLabel" );
 ServerL->setText( "Server list:" );
 ServerL->setFont(QFont("Helvetica",12,QFont::DemiBold));
 ServerL->setGeometry(WSeL);

 Server = new QListBox( this, "server" );
 Server->setGeometry(WSeLB);
 Server->setFont(QFont("Helvetica",14,QFont::Normal));
 connect( Server, SIGNAL(selected(int)), SLOT(ServerSelected(int)) );
 QToolTip::add( Server, "List of servers available" );

 QLabel *ShareL = new QLabel( this, "ShareLabel" );
 ShareL->setText( "Shared disks:" );
 ShareL->setFont(QFont("Helvetica",12,QFont::DemiBold));
 ShareL->setGeometry(WShL);

 Share = new QListBox( this, "share" );
 Share->setGeometry(WShLB);
 Share->setFont(QFont("Helvetica",14,QFont::Normal));
 connect( Share, SIGNAL(highlighted(int)), SLOT(ShareSelected(int)) );
 QToolTip::add( Share, "List of shared disks on selected server" );

 QLabel *DriveL = new QLabel( this, "DriveLabel" );
 DriveL->setText( "Path to connect:" );
 DriveL->setFont(QFont("Helvetica",12,QFont::DemiBold));
 DriveL->setGeometry(WDrL);

 Drive = new QLineEdit( this, "driveEdit" );
 Drive->setFont(QFont("Helvetica",14,QFont::Normal));
 Drive->setGeometry(WDrE);
 QToolTip::add( Drive, "Path to connect" );

 QLabel *MountL = new QLabel( this, "MountLabel" );
 MountL->setText( "Mount name:" );
 MountL->setFont(QFont("Helvetica",12,QFont::DemiBold));
 MountL->setGeometry(WMoL);

 MountDir = new QLineEdit( this, "mountEdit" );
 MountDir->setFont(QFont("Helvetica",14,QFont::Normal));
 MountDir->setGeometry(WMoE);
 QToolTip::add( MountDir, "Directory for mounting" );

 QLabel *UserL = new QLabel( this, "UserLabel" );
 UserL->setText( "Username:" );
 UserL->setFont(QFont("Helvetica",12,QFont::DemiBold));
 UserL->setGeometry(WUnL);

 UserName = new QLineEdit( this, "UserEdit" );
 UserName->setFont(QFont("Helvetica",14,QFont::Normal));
 UserName->setGeometry(WUnE);
 QToolTip::add( UserName, "Username to connect to share" );

 QLabel *PwL = new QLabel( this, "PwLabel" );
 PwL->setText( "Password:" );
 PwL->setFont(QFont("Helvetica",12,QFont::DemiBold));
 PwL->setGeometry(WPwL);

 PassWord = new QLineEdit( this, "PwEdit" );
 PassWord->setFont(QFont("Helvetica",14,QFont::Normal));
 PassWord->setGeometry(WPwE);
 PassWord->setEchoMode(QLineEdit::Password);
 QToolTip::add( UserName, "Password" );

  OkB=new QPushButton("Mount",this);
  OkB->setFont(QFont("Helvetica",14,QFont::DemiBold));
  OkB->setGeometry(WOkB);
  OkB->setFixedHeight( OkB->sizeHint().height() );
  QToolTip::add( OkB, "Mount selected share" );
  connect( OkB, SIGNAL(clicked()), SLOT(Accept()));

  CaB=new QPushButton("Close",this);
  CaB->setFont(QFont("Helvetica",14,QFont::DemiBold));
  CaB->setGeometry(WCaB);
  CaB->setFixedHeight( CaB->sizeHint().height() );
  QToolTip::add( CaB, "Close the window" );
  connect( CaB, SIGNAL(clicked()),qApp,SLOT(quit()));

  BrowseList BL(szDefServer,szDefUser,szDefPw);
  BL.ReqList(SERVER_LIST);
  setCursor(arrowCursor);
  if(BL.GetNHosts()<=0)
    {//PRINT_DEBUG("Error getting server list\n")
     QMessageBox::message("ERROR","Error getting server list");
     qApp->quit();
     return;
    }
  int i;

  for(i=0;i<BL.GetNHosts();i++)
     Server->insertItem(BL.GetHost(i),i);
  UserName->setText(szDefUser);
  PassWord->setText(szDefPw);
  Drive->setFocus();  

  show();
}
// *****************************************************
int MountShare::GetShareList(void)
{ struct FILE2 f=popen2("mount","mount",NULL);
     int i=0;

 if(!f.sto || !f.ste)
   {PRINT_DEBUG("Error open pipe")
    return i;
   }
 if(OpenShares.count())OpenShares.clear();

 char szC[MAX_LINELENGTH+1];
 while(fgets(szC,MAX_LINELENGTH,f.sto))
      {if(strstr(szC,"//")==szC)
         {char *p=strstr(szC," on");
          if(!p)break;
          char s=*p;
          *p=0;
          struct ShareItem *S;
          S=new struct ShareItem[1];
          strcpy(S->szServer,szC);
          *p=s;
          p+=4;
          char *e=strstr(p," type");
          *e=0;
          strcpy(S->szMount,p);
          OpenShares.append(S);
          //fprintf(stderr,"Share: %s\n",szC);
          i++;
         }
       }
 pclose2(f);
 return i;
}
// *******************************************************
void MountShare::UpdateOpenShareLB(void)
{char *p;
 if(!CurMount)return;
 if(CurMount->count())CurMount->clear(); 
 for( unsigned i=0; i<OpenShares.count();i++)
    {p=0;
     int il=strlen( OpenShares.at(i)->szMount );
     String S;
     if(il>18)p=OpenShares.at(i)->szMount+il-15;
     if(p)S.Setf("%s -> ...%s",OpenShares.at(i)->szServer,p);
     else S.Setf("%s -> %s",OpenShares.at(i)->szServer,OpenShares.at(i)->szMount);
     CurMount->insertItem((const char *)S,i);
    }   
}
// *******************************************************
void MountShare::Accept()
{strcpy(Current.szServer,Drive->text());
 strcpy(Current.szMount,MountDir->text());
 strcpy(Current.szUser,UserName->text());
 strcpy(Current.szPass,PassWord->text());
//  if(!CreateMount(&Current))emit reject();
//  else emit accept();
 CreateMount(&Current);
}
// *****************************************************
int MountShare::CreateMount(struct ShareItem* SI)
 {if(SI==0)return 0;
  String S;
  if(mkdir(SI->szMount,0750)<0)
   {if(errno==EEXIST)
      {struct stat Buf;
       if(stat(SI->szMount,&Buf)<0)
         {S.Setf("Unexpected error opening %s",SI->szMount);
          QMessageBox::message("ERROR",(const char *)S);
          return 0;
         }
       //PRINT_DEBUG("Nlink of %s = %d\n",SI->szMount,Buf.st_nlink)
       if(Buf.st_nlink!=2) // not an empty directory
         {S.Setf("directory %s is not empty",SI->szMount);
          QMessageBox::message("ERROR",(const char *)S);
          return 0;
         }

       S.Setf("Directory %s already exist",SI->szMount);
       if(QMessageBox::warning(0,"WARNING",(const char *)S,
                                              "Mount", "Cancel",0,1)==1)
        {return 0;} 

      }
    else
     {S.Setf("MountShare: Can not create directory %s: %s",
                        SI->szMount,strerror(errno));
      PRINT_DEBUG("%s\n",(const char *)S)
      QMessageBox::message("ERROR",(const char *)S);
      return 0;
     }
   }
 S.Setf("smbmount %s %s -U %s -P %s -D xphys",
         SI->szServer,SI->szMount,SI->szUser,SI->szPass);

 int i=system((const char *)S);

 //fprintf(stderr,"Executing: smbmount %s %s -U %s -P Uups -D xphys = %d\n",
 //        SI->szServer,SI->szMount,SI->szUser,i);

 if(i!=0)
   {S.Setf("ERROR> executing smbmount  %s %s -U %s -P Uups -D xphys",
           SI->szServer,SI->szMount,SI->szUser);
    PRINT_DEBUG("%s\n",(const char *)S)
    QMessageBox::message("ERROR",(const char *)S);
    return 0;
   }

 GetShareList();  
 UpdateOpenShareLB();
 return 1;
}
// *****************************************************
void MountShare::ServerSelected(int iS)
{BrowseList BL(Server->text(iS),"Browser","xfi");
 setCursor(waitCursor);
 BL.ReqList(SHARE_LIST);
 setCursor(arrowCursor);
 if(BL.GetNHosts()<=0)
    {PRINT_DEBUG("Share list empty\n")
     String S;
     S.Setf("Share list of %s empty", Server->text(iS));
     QMessageBox::message("ERROR",(const char *)S);
     Share->clear();
     return;
    }
 int i;
 if(Share->count()>0)Share->clear();

 for(i=0;i<BL.GetNHosts();i++)
     Share->insertItem(BL.GetHost(i),i);
 sprintf(Current.szServer,"//%s/",Server->text(iS));
 Drive->setText(Current.szServer);
}
// *****************************************************
void MountShare::ShareSelected(int iS)
{char *s=strrchr(Current.szServer,'/');
 s++;
 *s=0;
 sprintf(s,"%s",Share->text(iS));
 Drive->setText(Current.szServer);
}
// *****************************************************
void MountShare::UnMount()
{int i=CurMount->currentItem();

 if(!CurMount || !CurMount->count() || ! OpenShares.count())
   {PRINT_DEBUG("Nothing to unmount\n")
    return;
   }
 String S;
 if(i>=0)
   {struct ShareItem *SI;
    SI=OpenShares.at(i);
    S.Setf("smbumount %s",SI->szMount);
    int j=system((const char *)S);

    //fprintf(stderr,"Executing: smbumount %s = %d\n",SI->szMount,j);

    if(j!=0)
      {S.Setf("ERROR> executing smbumount %s",SI->szMount);
       PRINT_DEBUG("%s\n",(const char *)S)
       QMessageBox::message("ERROR",(const char *)S);
       return ;
      }

    //sprintf(szAnswer,"Unmounting %s",SI->szMount);
    CurMount->removeItem(i);
    CurMount->show();
    OpenShares.remove(i);
    //delete SI;
    
   }
}
// *****************************************************
const struct ShareItem* MountShare::GetShare(const int iS)
{if(iS>=0 && iS<(int)OpenShares.count())return OpenShares.at(iS);
 else return NULL;
}
// *****************************************************
int MountShare::GetNofShares(void)
{return OpenShares.count();}
// *****************************************************

#endif //!BROWSER_TEST

#if defined (BROWSER_TEST)
main(int ArgC, char ** ArgV)
{
 if(ArgC != 2)
   {printf("Usage: %s server\n",ArgV[0]);
    exit(EXIT_FAILURE);
   }

 printf("Testing Browser:\n\n");
 BrowseList B(ArgV[1],"BROWSER","xfi");

 char szB[2048]; 

 szB[0]=0;
 printf("Return: %d \n",B.ReqList(SERVER_LIST));
 B.GetHostList(szB);
 puts(szB);
 // printf("\n*******************************\n");
 //printf("Return: %d\n",B.ReqList(SHARE_LIST));
 //B.GetHostList(szB);
 //puts(szB);
 //printf("\n*******************************\n");

}


#endif

