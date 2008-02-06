//File: grafic.cpp
//$Id: grafic.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
//$Log: grafic.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.4  1999/06/14 09:07:12  herbie
//Laout manager added
//
//Revision 1.3  1999/03/15 09:08:37  herbie
//*** empty log message ***
//

#include <math.h>

#include <qcolor.h>
#include <qlist.h>
#include <qtooltip.h>
#include <qprinter.h>
#include <qkeycode.h>
#include <qapplication.h>
#include <qlayout.h>
#include <qmultilinedit.h>

#include "grafic.h"
#include "stdinc.h"
#include "stdfunc.h"

//extern char szAnswer [2048];

QList<DataFile> PlotList;

const QColor PlColors[]={black, red, darkGreen, blue, cyan,
                         magenta, darkYellow, gray};

QRect B(10,20, 90, 20);

extern QPoint SCREEN;

// ***************************************************
// General Functions
// **************************************************
int MinDist(QPointArray *p)
{unsigned i=p->size();
 if(i==0)return -1;
 QPoint p1,p2,d;
 int s,smin;
 for(i=0; i<p->size()-1;i++)
   {p1=p->point(i);
    p2=p->point(i+1);
    d=p2-p1;
    s=d.x()*d.x()+d.y()*d.y();
    if(i==0)smin=s;
    else {if(s<smin)smin=s;}
   }
 return smin;
}
// ***************************************************
void DrawSymbols(const int iSymbol, const int iSize, 
                 QPointArray *p, QPainter *qp, int iStart, int iEnd)
{int i;

 if(iEnd==-1)iEnd=p->size();
 if(iStart>=iEnd)return;

 for(i=iStart;i<iEnd;i++)
     DrawSymbol(iSymbol,iSize,qp,p->point(i),
                (i==iStart ? p->point(i):p->point(i-1)) );	

}
// *************************************************
void DrawSymbol(const int iSymbol, const int iSize, QPainter *qp,
                QPoint p, QPoint pp)
{//pp: previos point if iSymbol==CONNECT line to previous point is drawn
 switch(iSymbol)
   {case 5:
    case CIRCLE://0: CIRCLE
	    qp->drawEllipse(p.x()-iSize,p.y()-iSize, 2*iSize,2*iSize);
           break;          
   case 6:
   case QUAD://1: QUAD (SQUARE)
	   qp->drawRect(p.x()-iSize,p.y()-iSize,2*iSize,2*iSize);
           break;          
   case 7:
   case TRIANGLE:// 2: TRIANGLE
	  qp->drawLine(p.x()-iSize,p.y()+iSize,p.x(),p.y()-iSize);
          qp->drawLine(p.x(),p.y()-iSize, p.x()+iSize,p.y()+iSize);
          qp->drawLine(p.x()+iSize,p.y()+iSize,p.x()-iSize,p.y()+iSize);
          break;
   case DIAMOND:// 3: DIAMOND
	  qp->drawLine(p.x(),p.y()+iSize,p.x()-iSize,p.y());
          qp->drawLine(p.x()-iSize,p.y(),p.x(),p.y()-iSize);
          qp->drawLine(p.x(),p.y()-iSize,p.x()+iSize,p.y());
          qp->drawLine(p.x()+iSize,p.y(),p.x(),p.y()+iSize);
           break;          
   case CROSS://4: CROSS
          qp->drawLine(p.x(),p.y()-iSize,p.x(),p.y()+iSize);
          qp->drawLine(p.x()-iSize,p.y(),p.x()+iSize,p.y());
           break;          
   case LEFT_CROSS://8: LEFT_CROSS
          qp->drawLine(p.x(),p.y()-iSize,p.x(),p.y()+iSize);
          qp->drawLine(p.x()-iSize,p.y(),p.x(),p.y());
          break;          
   case V_LINE://9: V_LINE
          qp->drawLine(p.x(),p.y()-iSize,p.x(),p.y()+iSize);
          break;        
   case CONNECT://10: CONNECT
          if(pp.x()==-1 && pp.y()==-1)break;
          qp->drawLine(p.x(),p.y(),pp.x(),pp.y());
          break;        
   case    BAR://11: BAR
          qp->drawLine(p.x(),p.y(),pp.x(),pp.y());
          break;        

    default:PRINT_DEBUG("No symbol #%d avalilable\n",iSymbol)
           break;
   }

}
// *************************************************
void PaintSymbol(const int iSymbol, const int iSize, const QColor col, 
                 QWidget *PaintW, QPoint pt, QPoint ppt)
{ QPainter p;
  p.begin(PaintW);
  p.setPen(col);
  DrawSymbol(iSymbol, iSize, &p, pt, ppt);
  p.end();
}
// *************************************************
void DrawData(QPainter &p, DataFile *DF, CRData *CD, struct PlotPars &PP)
{
  // DF == NULL plots complete DataList
  // DF != NULL plots this DF
  // CD != NULL plots this CD 
  // One of CD or DF must be zero	 
 QPen pen(black);
 
 QPointArray **pts=NULL; 

 QFontMetrics FM(PP.PFont);
 FM=p.fontMetrics();

 p.setPen(pen);
 p.setFont(PP.PFont);

 unsigned iPCount;

 const int xDig=7,yDig=6;
 int iChx=FM.width('0');
 int iChy=FM.lineSpacing();
 
 if( DF && CD)
   {PRINT_DEBUG("Unambigous data parameter selection\n");
    p.drawText(10,10,"ERROR >DISP: No data!");
    return;
   }
	 
 if( (DF==NULL && PlotList.count()==0) && (DF && DF->GetCRData()==NULL))
   {PRINT_DEBUG("No data to plot\n");
    p.drawText(10,10,"ERROR >DISP: No data!");
    return;
   }
   
 else
   {
    if((PP.pFlag & TOP_FILENAME) && DF)p.drawText(iChx,iChy,DF->GetFileName());
    if((PP.pFlag & TOP_INFOLINE) && DF)p.drawText(iChx,iChy,DF->GetInfoText());

  CRData *TD = (DF==NULL ?(CRData*) CD :(CRData *) DF->GetCRData());
  
    if(TD)
      {TD->NewMinMax(PP.iStart,PP.iEnd);
       TD->GetMinMax(PP.DMin,PP.DMax);
      }
    else 
      {if(!GetMinMax(PP.DMin,PP.DMax))return;}
 
int iTop = ((PP.pFlag & TOP_FILENAME) ||
	    (PP.pFlag & TOP_IDN)      || (PP.pFlag & TOP_INFOLINE)? 3 : 1);

//    QRect EffA(iChx*(yDig+5)+PP.PRect.x(),iTop*iChy+PP.PRect.y(),
//               PP.PRect.width()-(yDig+8)*iChx,PP.PRect.height()-(2+iTop)*iChy);
    QRect EffA(iChx*(yDig+3)+PP.PRect.x(),iTop*iChy+PP.PRect.y(),
               PP.PRect.width()-(yDig+5)*iChx,PP.PRect.height()-(2+iTop)*iChy);
    //raw y-Axis
    p.drawLine(EffA.x(),EffA.y(),
               EffA.x(),EffA.y()+EffA.height()+iChy/4);

    p.drawLine(EffA.x()+EffA.width(),EffA.y(),
               EffA.x()+EffA.width(),EffA.y()+EffA.height());

    //raw x-Axis
    p.drawLine(EffA.x()-iChx/2,EffA.height()+EffA.y(),
               EffA.x()+EffA.width(),EffA.height()+EffA.y());

    p.drawLine(EffA.x(),EffA.y(),
               EffA.x()+EffA.width(), EffA.y());
    
    
    MPoint NMin,NMax,NDelta;
    //Scale y-axis
    if(!Scale(PP.DMin.y,PP.DMax.y,EffA.height()/2/iChy,
                                  NMin.y,NMax.y,NDelta.y))
      {p.setPen(red);
       p.drawText(EffA.x(),EffA.height()/2,"Error scaling y-Axis");
       p.end();
       return;
      }

    //Scale x-axis
    if(!Scale(PP.DMin.x,PP.DMax.x,(int)(EffA.width()/(xDig+2.)/iChx+.5),
                                  NMin.x,NMax.x,NDelta.x)) //old res: (int)(EffA.width()/1.5/xDig/iChx)
      {p.setPen(red);
       p.drawText(EffA.x(),EffA.height()/2,"Error scaling y-Axis");
       p.end();
       return;
      }

    int iR=(int) ( (NMax.y-NMin.y)/NDelta.y+.5) ; 
    int i;
    int tick;
    double label;
    char szN[MAXPATH+1];
    char szT[MAXPATH+1];
    GetCurrentTime(szT);
    if(PP.pFlag & TOP_IDN)
      {sprintf(szN,"Plot from  XAusw (c) HM %s",szT);
       p.drawText(5*iChx,iChy,szN);
      }

    //y-Axis ticks
    if(PP.yLab && PP.yLab[0]) p.drawText(PP.PRect.x()+5*iChx,
                                         PP.PRect.y()+(int)((iTop+0.5)*iChy/2),PP.yLab);
    for(i=0;i<=iR;i++)
      {label=NMin.y+i*NDelta.y;
       tick=EffA.y()+EffA.height()-(i*EffA.height()/iR);

       p.drawLine(EffA.x()-iChx/2,tick,EffA.x(),tick);
       gcvt(label,yDig,szN);

       p.drawText(EffA.x()-FM.width(szN)-iChx,tick+iChy/2,szN);
        
      }

    //x-Axis ticks
    if(PP.xLab && PP.xLab[0])
       p.drawText(PP.PRect.x()+PP.PRect.width()-iChx*(strlen(PP.xLab)+2),
                  PP.PRect.y()+PP.PRect.height()-5,PP.xLab);

    iR=(int) ( (NMax.x-NMin.x)/NDelta.x+.5);
    for(i=0;i<=iR;i++)
      {label=NMin.x+i*NDelta.x;
       tick=EffA.x()+(i*EffA.width()/iR);

       p.drawLine(tick,EffA.height()+EffA.y(),
                  tick,EffA.height()+EffA.y()+iChy/4);
       gcvt(label,xDig,szN);

       p.drawText(tick-FM.width(szN)/2,EffA.y()+EffA.height()+iChy,szN);
        
      }
 
    PP.DMin=NMin;
    PP.DMax=NMax;
    p.setPen(red);

    MPoint PMin={EffA.x(),EffA.height()+EffA.y(),0};
    MPoint PMax={EffA.width()+EffA.x(),EffA.y(),0};

    if(PP.pFlag & ZERO_X_LINE)
      {int y0=(int) fLinTrafo(PMax.y,PMin.y,NMax.y,NMin.y,0);
        p.drawLine(EffA.x(),y0, EffA.x()+EffA.width(),y0);
      }

    if(TD)
      {pts=new QPointArray * [1];
       iPCount=1;
      }
    else
      {iPCount=PlotList.count();
       pts=new QPointArray * [iPCount];
      }
 // Initialize point array
    CRData *T;
    DataFile *D;
    unsigned u;
    for(u=0 ; u<iPCount; u++)
      {D=(DF ? DF : PlotList.at(u));
       T=(D ? D->GetCRData() : CD);
       pts[u]= new QPointArray(T->GetSteps());
       }
 // Perform plot(s)
    unsigned j;
    for(j=0; j < iPCount && j < (unsigned)MAX_PLOTS; j++)
      {D=(DF ? DF :  PlotList.at(j));
       T=(D ? D->GetCRData() : CD);
       for(u=0; u<(unsigned)T->GetSteps(); u++)
        {
         pts[j]->setPoint(u,
               (int) fLinTrafo(PMax.x,PMin.x,NMax.x,NMin.x,T->GetPointX(u)),
               (int) fLinTrafo(PMax.y,PMin.y,NMax.y,NMin.y,T->GetPointY(u)) );
        }//for u
      }//for j
     int im;//iSy=0;
     int nPts=(PP.iEnd==-1?PP.iEnd:PP.iEnd-PP.iStart+1);    
     QPoint ep;

     for(u=0; u<iPCount && u<MAX_PLOTS ;u++)
        {//Old struct p.setPen(PlColors[u]);
	 p.setPen(PlColors[PP.pColor[u]]);
         //Old struct im=(PP.iSize >= 0  ?  MinDist(pts[u]) : -4*PP.iSize);
         im=( (PP.pSize[u] == 0 || !PP.pSize)  ?  MinDist(pts[u]) : 4*PP.pSize[u]);
         //Old struct if(im<3 && !PP.iErrBar)
	 if(im<3 && !(PP.pSymbol[u]==ERROR_BAR || PP.pSymbol[u]==BAR))
           {p.drawPoints(*pts[u],PP.iStart,nPts);
            //Old struct if(u==iPCount-1)PP.iSize=0;
            if(u==iPCount-1 && PP.pSize)PP.pSize[u]=1;
           }
         else 
           {//if(!PP.iErrBar)
	    //if(PP.pSymbol[u]!=ERROR_BAR)
            switch(PP.pSymbol[u])
              {default:if(im>16)im=16;
                   DrawSymbols(PP.pSymbol[u],im/4,pts[u],&p,PP.iStart,PP.iEnd);
                   if(u==iPCount-1)PP.pSize[u]=im/4;
                   break;
            case ERROR_BAR: //draw error bars
                  {int ip;
	           for(ip=PP.iStart;ip<(PP.iEnd==-1?(int)pts[u]->size():PP.iEnd); ip++)
		      {double y=T->GetPointY(ip);
		       double e=(y>0?sqrt(y):0);
                       int iy=(int) fLinTrafo(PMax.y,PMin.y,NMax.y,NMin.y,y);
                       int ie=(int) fLinTrafo(PMax.y,PMin.y,NMax.y,NMin.y,y+e);
                       DrawSymbol(V_LINE,(iy==ie?2:ie-iy),&p,pts[u]->point(ip));
                      }
                   }
                   break;
             case BAR:      
                  {int ip;
	           for(ip=PP.iStart;ip<(PP.iEnd==-1?(int)pts[u]->size():PP.iEnd); ip++)
		       DrawSymbol(BAR,1,&p,pts[u]->point(ip),
                                  QPoint(pts[u]->point(ip).x(),EffA.height()+EffA.y()));
                      
                   }
                   break;

	       }//switch
           }
         if(PP.pFlag & BOT_FILELIST && ( D=(DF ? DF : PlotList.at(u)) ) )
           {//D=(DF ? DF : PlotList.at(u));
            T=D->GetCRData();
            im=T->GetIMaxX();
            sprintf(szN,"%2d",u+1);
            ep=pts[u]->point(im);
            p.drawText(EffA.width()+EffA.x()-3*iChx,ep.y(),szN);

            strcpy(szN,D->GetFileName());
            im=strlen(szN);

            char szB[MAX_LINELENGTH+1];
	    sprintf(szB,"%2d: %s (%d points)>%s",u+1,&szN[im>30?im-30:0],
                                      T->GetSteps(),
                                      D->GetInfoText());
            p.drawText(5*iChx,PP.PRect.height()+PP.PRect.y()+iChy+iChy*u,szB);
           }
         if(PP.pFlag & BOT_FILEHEAD && iPCount==1 && (D=(DF ? DF : PlotList.at(u))))
           {//D=(DF ? DF : PlotList.at(u));
            String Answer("");
            D->SPrintInfo(Answer,PR_ALL & ~PR_COLS);
            
	    QFont hf(PP.PFont);
	    hf.setPointSize(8);
            DrawFormText(p,hf,10*iChx,PP.PRect.height(),
	                 AlignTop|AlignLeft,Answer.GetBuf());
  	   }		 
	    
        }// for (u)
    PP.PRect=EffA;

   }

 for(unsigned u=0;u<iPCount;u++)if(pts[u])delete pts[u];
 if(pts)delete [] pts; 
}
// **************************************************
int GetMinMax(MPoint &Min, MPoint &Max)
{if(PlotList.count()==0)return 0;
 CRData *T;
 for(unsigned i=0; i<PlotList.count(); i++)
    {T=PlotList.at(i)->GetCRData();
     T->NewMinMax();
     if(i==0)T->GetMinMax(Min,Max);

     if(T->GetMinX()<Min.x)Min.x=T->GetMinX();
     if(T->GetMaxX()>Max.x)Max.x=T->GetMaxX();
     if(T->GetMinY()<Min.y)Min.y=T->GetMinY();
     if(T->GetMaxY()>Max.y)Max.y=T->GetMaxY();

    }
 return 1;

}
// **************************************************
void DrawFormText(QPainter &p, QFont F, const int x0, const int y0,
                  const int flag,  const char *szT)
{QFont f(p.font());
 p.setFont(F);
 //QFontMetrics FM(F);
 QFontMetrics FM=p.fontMetrics();
 int iSp=FM.lineSpacing();
 int ix=FM.width("0");
 int icr=0;
 const char *pc=szT;
 int iW=0,iWMax=0;
 while(*pc)
      {iW++;
       if (*pc=='\n')
	  {icr++;
	   if(iW>iWMax)iWMax=iW;
	   iW=0;
	  }
	pc++;
       } 
 //fprintf(stderr,"%d lines; spacing: %d -> %d lines free\n"
 //       ,icr,iSp, (600-PP.PRect.y())/iSp);		 
 p.drawText(x0,y0,ix*iWMax,icr*iSp,flag, szT);
 p.setFont(f);
}
// *****************************************************
void PrintText(QMultiLineEdit *w, FILE *f, const int iAppend)
{char szB[MAX_LINELENGTH+1];
 int iC=0;
  while(1)
      {if(!fgets(szB,MAX_LINELENGTH,f))break;
       szB[MAX_LINELENGTH]=0;
       if(szB[strlen(szB)-1]=='\n')szB[strlen(szB)-1]=0;
       if(!iAppend && iC==0) w->setText(szB);
       else w->append(szB);
       iC++;
     }    
}
// *****************************************************
// struct PlotPars
PlotPars::PlotPars(void)
{PlotPars::Init();}
// *****************************************************
PlotPars::PlotPars(const struct PlotPars & P)
{PlotPars::Copy(P);}
// *****************************************************
struct PlotPars& PlotPars::operator=(const struct PlotPars & P)
{PlotPars::Copy(P);
 return *this;
}
// ******************************************************
void PlotPars::Init(void)
{nPlots=MAX_PLOTS;
 xLab=yLab=NULL;
 pFlag=iStart=0;
 iEnd=-1;
 int i;
 for(i=0;i<nPlots;i++)
    {pSize[i]=0;
     pSymbol[i]=i;
     pColor[i]=i;
    } 
}
// *****************************************************
void PlotPars::Copy(const struct PlotPars & P)
{ PFont=P.PFont;
  PRect=P.PRect;
  xLab=P.xLab;
  yLab=P.yLab;
  pFlag=P.pFlag;
  iStart=P.iStart;
  iEnd=P.iEnd;
  //iSize=P.iSize;
  //iErrBar=P.iErrBar;
  DMin=P.DMin;
  DMax=P.DMax;
    
  nPlots=P.nPlots;
 
 int i;
 for(i=0;i<nPlots;i++)
    {pSize[i]=P.pSize[i];
     pSymbol[i]=P.pSymbol[i];
     pColor[i]=P.pColor[i];
    } 
}
// *****************************************************
// *****************************************************
//DisPlayData
// *****************************************************
 int DisplayData::w=0;
 int DisplayData::h=0;
 QRect DisplayData::B0;
 QRect DisplayData::WMain;
 QRect DisplayData::WxL;
 QRect DisplayData::WxE;
 QRect DisplayData::WyL;
 QRect DisplayData::WyE;
 QRect DisplayData::WeB;
 QRect DisplayData::WaB;
 QRect DisplayData::WpB;
 QRect DisplayData::WcB;
 QRect DisplayData::WPlot;
 QRect DisplayData::WFrame;
// **************************************************
// DisplayData Functions
// **************************************************
DisplayData:: DisplayData(QWidget *parent, const char *name):
              QWidget(parent,name)
{
 //DisplayData::w=(WidgetAusw::ScreenHeight()==600 ? 530: 700);
 //DisplayData::h=(WidgetAusw::ScreenHeight()==600 ? 530: 700);

    DisplayData::w=(SCREEN.y()==600 ? 530: 700);
    DisplayData::h=(SCREEN.y()==600 ? 530: 700);
 
    DisplayData::B0.setRect(::B.x(),::B.y(),::B.width(),::B.height());
    DisplayData::WMain.setRect(10, 10,w,h);
    DisplayData::WxL.setRect(B.x(), B.y(), B.width(), B.height());
    DisplayData::WxE.setRect(B.x(), WxL.y()+WxL.height()+10,
                             B.width(),B.height());
    DisplayData::WyL.setRect(B.x(), WxE.y()+WxE.height(),B.width(),B.height());
    DisplayData::WyE.setRect(B.x(), WyL.y()+WyL.height()+10,
                             B.width(),B.height());
    DisplayData::WeB.setRect(B.x(), WMain.height()-4*B.y()-4*B.height(),
                             B.width(),   B.height());
    DisplayData::WaB.setRect(B.x(), WMain.height()-6*B.y()-6*B.height(),
                             B.width(),   B.height());
    DisplayData::WpB.setRect(B.x(), WMain.height()-2*B.y()-2*B.height(),
                             B.width(),   B.height());
    DisplayData::WcB.setRect(B.x(), WMain.height()-B.y()-B.height(),
                             B.width(),B.height());
    DisplayData::WFrame.setRect(2*B0.x()+B0.width(),  B.y(),
                                WMain.width()-3*B.x()-B.width(),
				WMain.height()-2*B.y());

    DisplayData::WPlot.setRect(0,0,
                               WFrame.width()-5,h*3/4);


 iQuit = END_CLOSE;
 nPlots=0;
 iFormat=NORM_FORM;
 PlotList.setAutoDelete(FALSE);
 if(PlotList.count())PlotList.clear();
}
// ********************************************************
int DisplayData::AddData(DataFile *df)
{
 if(!df || df->GetError())
   {PRINT_DEBUG("Error adding DataFile-Object to plot list\n")
    return 0;
   }
 PlotList.append(df);
 nPlots++;
 return 1;

}
// *********************************************************
int DisplayData::InitPlot(void)
{
  setGeometry(WMain);
  setCaption("Plot selected data");

  //Label Edit
  QLabel *xMsg=new QLabel("x-Label:", this);
  xMsg->setAlignment(AlignLeft);
  xMsg->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  xMsg->setGeometry(WxL);
  xLab = new QLineEdit( this );
  xLab->setGeometry(WxE);
  DataFile *DF=PlotList.at(0);
  int iColX=DF->GetCRData()->GetColX();
  int iColY=DF->GetCRData()->GetColY();
  //fprintf(stderr,"DisplayData::InitPlot xCOL: %d, yCOL:%d\n",iColX,iColY);

  xLab->setText(DF->GetColID(iColX));
  QToolTip::add( xLab, "Enter new text for x-Label" );
   connect(xLab , SIGNAL(returnPressed()), SLOT(DataChanged()) );

  QLabel *yMsg=new QLabel("y-Label:", this);
  yMsg->setAlignment(AlignLeft);
  yMsg->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  yMsg->setGeometry(WyL);
  yLab = new QLineEdit( this );
  yLab->setGeometry(WyE);
  yLab->setText(DF->GetColID(iColY));
  QToolTip::add( yLab, "Enter new text for y-Label" );
  connect(yLab , SIGNAL(returnPressed()), SLOT(DataChanged()) );

  //Push buttons
  if(nPlots==1)
    {aBut = new QPushButton("Form",this);
     aBut->setFont(QFont("Helvetica",13));
     aBut->setText("Archive Form");
     aBut->setGeometry(WaB);
     aBut->setFixedHeight( aBut->sizeHint().height() );
     connect(aBut, SIGNAL(clicked()), SLOT(ArchiveFormat()) );
     QToolTip::add( aBut, "Toggle archive/normal format" );
    }
  else aBut=NULL;
      
    eBut = new QPushButton("EBars",this);
    eBut->setFont(QFont("Helvetica",13));
    eBut->setText("Error bars");
    eBut->setGeometry(WeB);
    eBut->setFixedHeight( eBut->sizeHint().height() );
    connect(eBut, SIGNAL(clicked()), SLOT(ErrorBars()) );
    QToolTip::add( eBut, "Toggle error bars/symbols" );
 
    pBut = new QPushButton("Print",this);
    pBut->setFont(QFont("Helvetica",13));
    pBut->setText("Print");
    pBut->setGeometry(WpB);
    pBut->setFixedHeight( pBut->sizeHint().height() );
    connect(pBut, SIGNAL(clicked()), SLOT(PrintIt()) );
    QToolTip::add( pBut, "Print graphic" );

    cBut = new QPushButton(this, "close" );
    cBut->setFont(QFont("Helvetica",13,QFont::DemiBold));
    cBut->setText("Close");
    cBut->setGeometry(WcB);
    cBut->setFixedHeight( cBut->sizeHint().height() );
    connect(cBut, SIGNAL(clicked()), SLOT(Close()) );

    QToolTip::add( cBut, "Close this window" );

    QFrame *PlotF = new QFrame(this,"PlotFrame");
    PlotF->setGeometry(WFrame);
    PlotF->setFrameStyle( QFrame::Panel | QFrame::Sunken);
    PlotF->setBackgroundColor(white);
    PlotF->setLineWidth(3);

    strncpy(szXLabel,DF->GetColID(iColX),79);
    strncpy(szYLabel,DF->GetColID(iColY),79);

    PlotA=new MultiPlot(this,"PlotArea");
    struct PlotPars P=PlotA->GetPlotPars();;
    //PlotA->GetPlotPars().Init();
    //PlotA->SetPlotPars(P);
    P.xLab=szXLabel;
    P.yLab=szYLabel;
    P.PRect.setRect(WPlot.x(),WPlot.y(),WPlot.width(),WPlot.height());
    //PlotA->GetPlotPars().PRect.setRect(WPlot.x(),WPlot.y(),WPlot.width(),WPlot.height());
    PlotA->SetPlotPars(P);
    PlotA->setGeometry(WFrame.x()+3,WFrame.y()+3,WFrame.width()-3,WFrame.height()-3);
    PlotA->setBackgroundColor(white);
    PlotA->show();

    //LegA = new QFrame(this,"LegFrame");
    //LegA->setGeometry(WLeg.l,WLeg.t,WLeg.r,WLeg.b);
    //LegA->setFrameStyle( QFrame::Panel | QFrame::Sunken);
    //LegA->setBackgroundColor(white);
    //LegA->setLineWidth(3);
    show();

  return 1;
 }
// ******************************************************
void DisplayData::DataChanged()
{
 strncpy(szXLabel,xLab->text(),79);
 strncpy(szYLabel,yLab->text(),79);
 PlotA->repaint();
}
// ******************************************************
void DisplayData::Close()
{if(iQuit==END_QUIT) qApp->quit();
 else hide();
}
// ******************************************************
void DisplayData::ErrorBars()
{setUpdatesEnabled(TRUE);
 struct PlotPars P=PlotA->GetPlotPars();
 int i;
 for(i=0;i<nPlots;i++)
    {if(P.pSymbol[i]!=ERROR_BAR)
       {P.pSymbol[i]=ERROR_BAR; 
        if(i==0) eBut->setText("Symbols");
       }
     else
       {P.pSymbol[i]=i;
        if(i==0)eBut->setText("Error bars");
       }
    }     
 PlotA->SetPlotPars(P);
 PlotA->update();
 repaint();
 
}
// **************************************************
void DisplayData::ArchiveFormat()
{if(nPlots!=1)
   {PRINT_DEBUG("Archive format only for a single plot possible\n")
    return;
   }
 setUpdatesEnabled(TRUE);
 struct PlotPars P=PlotA->GetPlotPars();
 if(iFormat==NORM_FORM)
   {
    P.pFlag &= ~BOT_FILELIST;
    P.pFlag |= BOT_FILEHEAD;
    aBut->setText("Normal Form");
    iFormat=ARCH_FORM;   
    //DisplayData::WPlot.setRect(2*B0.x()+B0.width(),
    //                           B.y(),WMain.width()-3*B.x()-B.width(),h*1/2);
    DisplayData::WPlot.setHeight(h*1/2);

   }
 else
   {iFormat=NORM_FORM;   
    aBut->setText("Archive Form");
    P.pFlag |= BOT_FILELIST;
    P.pFlag &= ~BOT_FILEHEAD;
    //DisplayData::WPlot.setRect(2*B0.x()+B0.width(),
    //                           B.y(),WMain.width()-3*B.x()-B.width(),h*3/4);
    DisplayData::WPlot.setHeight(h*3/4);
   }

 PlotA->SetPlotPars(P);
 PlotA->update();
 repaint();
 
}
// **************************************************
void DisplayData::keyPressEvent(QKeyEvent *e)
{struct PlotPars P=PlotA->GetPlotPars();
 static int iP=0; 
 switch(e->key()) 
       {case Key_1:iP=0;break;
        case Key_2:iP=1;break;
        case Key_3:iP=2;break;
        case Key_4:iP=3;break;
        case Key_5:iP=4;break;
        case Key_6:iP=5;break;
        case Key_7:iP=6;break;
        case Key_8:iP=7;break;
 	case Key_C:if(P.pSymbol[iP]!=CONNECT)P.pSymbol[iP]=CONNECT;
	           else P.pSymbol[iP]=iP;
    	           PlotA->SetPlotPars(P);
		   PlotA->update();
		   repaint();
		   break;
  	case Key_E:if(P.pSymbol[iP]!=ERROR_BAR)P.pSymbol[iP]=ERROR_BAR;
	           else P.pSymbol[iP]=iP;
    	           PlotA->SetPlotPars(P);
		   PlotA->update();
		   repaint();
		   break;
 	case Key_Plus:if(P.pSymbol[iP]>=0)P.pSize[iP]++;
	              //fprintf(stderr,"Sybolsize[%d] + : %d\n",iP,P.pSize[iP]);
		      PlotA->SetPlotPars(P);
		      PlotA->update();
		      repaint();
		      break;   
 	case Key_Minus:if(P.pSymbol[iP]>=0)P.pSize[iP]--;
	              //fprintf(stderr,"Sybolsize[%d] - : %d\n",iP,P.pSize[iP]);
		      PlotA->SetPlotPars(P);
		      PlotA->update();
		      repaint();
		      break;   
        case Key_Escape:Close();
	                break;
           default:e->ignore();
                   return;
	  }//switch
}
// **************************************************
void DisplayData::PrintIt()
{
  QPrinter Prt;
  if(Prt.setup(PlotA))
   {QPainter paint;
    paint.begin(&Prt);
    struct PlotPars P=PlotA->GetPlotPars();
    P.PFont=QFont("Helvetica",13);
    P.PRect=QRect(WPlot);
    P.PRect.setX(0);
    P.PRect.setY(0);
    P.xLab=szXLabel;
    P.yLab=szYLabel;
    //P.pFlag= TOP_IDN | BOT_FILELIST;
    P.pFlag |= TOP_IDN;
    P.iStart=0;
    P.iEnd=-1;

    DrawData(paint,NULL,NULL,P);
    paint.end();
   }
  //else fprintf(stderr,"Printer setup failed\n");

}
// ******************************************************
// MultiPlot
// ******************************************************
MultiPlot::MultiPlot(QWidget *parent, const char *name):QWidget(parent,name)
{ PLP.PRect=QRect(DisplayData::GetPlotW());
  PLP.pFlag= TOP_IDN | BOT_FILELIST,
  PLP.PFont=QFont("Helvetica",13);
  //if(WidgetAusw::ScreenHeight()==600)PLP.PFont.setPointSize(11);
  if(SCREEN.y()==600)PLP.PFont.setPointSize(11);
  PLP.iStart=0;
  PLP.iEnd=-1;
  //PLP.iErrBar=0;          
}
// ******************************************************
void MultiPlot::paintEvent(QPaintEvent *)
{
 QPainter p;
 p.begin(this);
 PLP.PRect=QRect(DisplayData::GetPlotW());
 DrawData(p,NULL,NULL,PLP);
 p.end();
}
// ******************************************************
// **************************************************
// ManipData
// **************************************************
  int ManipData::w=0;
  int ManipData::h=0;
 
  QRect ManipData::WMain;
//  QRect ManipData::WIT;
  QSize ManipData::WPl;
//  QRect ManipData::WIB;
//  QRect ManipData::B0;
  int ManipData::iFS;
   const int NoOfBut=10;
  const char *szButLab[NoOfBut]=
  {"Help F1", "Error bars F2", "Meanvalue F3", "Remove F4", "New Range F5",
   "Prev Range F6", "Integrate F7", "None F8", "Sav Range F9","Exit F10"};
  const char *szBTip[NoOfBut]=
  {"Display Help","Toggle error bars","Meanvalue of neigbour points","Remove actual point",
   "Select a new range", "Goto previous range","Numerical integral",
   "","Save range & exit","Exit now"}; 
  const int ButHere[NoOfBut]={1,1,1,1,1,1,1,0,1,1};
 
// ***************************************************
ManipData::ManipData(DataFile *df, QWidget *parent, const char *name):
           QWidget(parent,name)
{
 ManipData::w=(SCREEN.y()==600 ? 750: 800);
 ManipData::h=(SCREEN.y()==600 ? 510: 600);
 ManipData::iFS=(SCREEN.y()==600 ? 12 :14); 
  ManipData::WMain.setRect(10,    10,w,h);
//   ManipData::WIT.setRect(B.x()+B.width()+10,  10, 
//                          WMain.width()-B.y()-B.width()-20, B.height());

  ManipData::WPl.setWidth(WMain.width()-B.y()-B.width()-20);
  ManipData::WPl.setHeight(WMain.height()-2*B.height()-40);
    
//  ManipData::WPl.setRect(WIT.x(), WIT.y()+WIT.height()+10, WIT.width(),
//                             WMain.height()-2*WIT.height()-40) ;
//   ManipData::WIB.setRect(WIT.x()+20,  WPl.y()+WPl.height()+10, WIT.width(),
//                              WIT.height());
//   ManipData::B0.setRect(10, 10, B.width(), B.height());

  if(!df ||df->GetError())
    {TD=NULL;return;}

  DF=df;
  TD=df->GetCRData();

  if(!TD ||TD->GetSteps()<=0)
    {TD=NULL; DF=NULL; return;}
  
  iStep=0;
  First=TRUE;
  iRCount=0;
  bRSelect=FALSE;

  iQuit=END_CLOSE;
  
  PLP.iStart=0;
  PLP.iEnd=TD->GetSteps();
  //PLP.iErrBar=0;
  //PLP.iSize=1;

  setGeometry(ManipData::WMain); 
  setCaption("Change selected data");

  QBoxLayout *topLayout = new QVBoxLayout( this, 5 );

  
  char szA[MAX_LINELENGTH+1];
  QBoxLayout *TopLabel = new QHBoxLayout();
  topLayout->addLayout( TopLabel );  

  InfoLT=new QLabel(szA, this, "TopLabel");
  InfoLT->setAlignment(AlignLeft|AlignVCenter);
  InfoLT->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  InfoLT->setMinimumSize( InfoLT->sizeHint() );
  InfoLT->setFixedHeight(InfoLT->height());
  //InfoLT->setMaximumSize(SCREEN.x(),B.height());
  //InfoLT->setGeometry(WIT);
  TopLabel->addWidget(InfoLT,100);
  

  QBoxLayout *Middle = new QHBoxLayout();
  topLayout->addLayout( Middle );  

  QBoxLayout *buttons = new QVBoxLayout();
  Middle->addLayout(buttons,1);  

 int i;
  for(i=0;i<NoOfBut;i++)
   {if(ButHere[i])
     {sprintf(szA,"B%2d",i);
      Bts[i] = new QPushButton(this,szA);
      Bts[i]->setFont(QFont("Helvetica",iFS,QFont::Light));
      Bts[i]->setText(szButLab[i]);
      Bts[i]->setMinimumSize( Bts[0]->sizeHint() );
      Bts[i]->setFixedWidth(Bts[0]->sizeHint().width()*2 );
      //Bts[i]->setGeometry(B0.x(), B0.y()+i*(B0.height()+20),
      //                    B0.width(),B0.height());
      //Bts[i]->setFixedHeight( Bts[i]->sizeHint().height() );
      buttons->addWidget( Bts[i]);
      QToolTip::add( Bts[i], szBTip[i] );
     }
        
    else 
     {Bts[i]=NULL;
      QLabel *L=new QLabel(this);
      L->setFrameStyle( QFrame::NoFrame );
      L->setMinimumSize(Bts[0]->sizeHint());
      L->setFixedWidth(Bts[0]->width());
      buttons->addWidget(L);
     }
   }//for
  connect( Bts[1], SIGNAL(clicked()), SLOT(ErrorBars()) );
  connect( Bts[2], SIGNAL(clicked()), SLOT(MeanVal()) );
  connect( Bts[3], SIGNAL(clicked()), SLOT(RemovePoint()) );
  connect( Bts[4], SIGNAL(clicked()), SLOT(NewRange()) );
  connect( Bts[5], SIGNAL(clicked()), SLOT(RangeBack()) );
  connect( Bts[6], SIGNAL(clicked()), SLOT(Integrate()) );
  connect( Bts[8], SIGNAL(clicked()), SLOT(SaveRange()) );
  connect( Bts[9], SIGNAL(clicked()), SLOT(Close()) );

  
  PlotF = new MTrackWid(this,"PlotFrame");
  //PlotF->setGeometry(WPl);
  PlotF->setMinimumSize( QSize(WPl.width(),WPl.height()));
  PlotF->setFrameStyle( QFrame::Panel | QFrame::Sunken);
  PlotF->setBackgroundColor(white);
  PlotF->setLineWidth(3);
  connect(PlotF,SIGNAL(RCoords(QPoint)),SLOT(ShowCoords(QPoint)));
  connect(PlotF,SIGNAL(LCoords(QPoint)),SLOT(ClosestPoint(QPoint)));
  Middle->addWidget(PlotF,20);
 

  QBoxLayout *BotLabel = new QHBoxLayout();
  topLayout->addLayout( BotLabel );  

  InfoLB=new QLabel(szA, this, "BotLabel");
  InfoLB->setAlignment(AlignLeft|AlignVCenter);
  InfoLB->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  //InfoLB->setGeometry(WIB);
  InfoLB->setMinimumSize(InfoLB->sizeHint());
  InfoLB->setFixedHeight(InfoLB->height());
  BotLabel->addWidget(InfoLB);

  //PlotF->setEnabled(TRUE);
  //PlotF->repaint();
  //setMouseTracking(TRUE); 
  topLayout->activate();
  show();

}
// ***************************************************
void ManipData::ShowCoords(QPoint p)
{
 char szA[MAX_LINELENGTH+1];
 if(EffP.contains(p))
  {MPoint D;
   D.x=fBackTrafo(EffP.width()+EffP.x(),EffP.x(),PLP.DMax.x,PLP.DMin.x,p.x());
   D.y=fBackTrafo(EffP.height()+EffP.y(),EffP.y(),PLP.DMin.y,PLP.DMax.y,p.y());
   sprintf(szA,"X: %12.5g Y: %12.5g",D.x,D.y);
  }
   else sprintf(szA,"Out of plot area");
   InfoLB->setText(szA);
}
// **************************************************
void ManipData::ClosestPoint(QPoint p)
{
 if(EffP.contains(p))
  {MPoint D;
   D.x=fBackTrafo(EffP.width()+EffP.x(),EffP.x(),PLP.DMax.x,PLP.DMin.x,p.x());
   D.y=fBackTrafo(EffP.height()+EffP.y(),EffP.y(),PLP.DMin.y,PLP.DMax.y,p.y());
   int s;
   double d[3];
   int sig=(int) sign(TD->GetPointX(PLP.iStart)-D.x);
   for(s=PLP.iStart;sign(TD->GetPointX(s)-D.x)==sig && s<PLP.iEnd-1;s++);
   int s1=s;
   if(s>PLP.iEnd-3)s1=PLP.iEnd-3;
   int s2=s1+3;
   if(s1<0 || s2<0)return;
   int i;
   for(i=s1; i<s2; i++)
      {double dx=D.x-TD->GetPointX(i);
       double dy=D.y-TD->GetPointY(i);
       d[i-s1]=dx*dx+dy*dy;
      }
  if(d[2] < (d[0] < d[1] ? d[0] : d[1]) )iStep=s+1;
  if(d[0] < (d[1] < d[2] ? d[1] : d[2]) )iStep=s-1;
  if(D.x<=TD->GetPointX(PLP.iStart))iStep=PLP.iStart;
  if(D.x>=TD->GetPointX(PLP.iEnd-1))iStep=PLP.iEnd-1;

  MoveCursor(red);
  char szA[MAX_LINELENGTH+1];
  sprintf(szA,"Point: %d of %d X:%12.5g Y:%12.5g prev. ranges:%d",
              iStep+1,TD->GetSteps(),TD->GetPointX(iStep),
              TD->GetPointY(iStep),iRCount);
  InfoLT->setText(szA);  

  //sprintf(szA,"Mouse to %d",iStep);
  }
 // else sprintf(szA,"Out of plot area");
 //  InfoLB->setText(szA);
}
// **************************************************
void ManipData::MeanVal()
{
 int iL=iStep-1;
 if(iStep==0)iL=1;
 if(iStep==TD->GetSteps())iL=TD->GetSteps()-2;
 
 int iR=iStep+1;
 if(iStep==0)iR=2;
 if(iStep==TD->GetSteps())iR=TD->GetSteps()-1;

 double DeltaX=TD->GetPointX(iR) - TD->GetPointX(iL);
 double DeltaY=TD->GetPointY(iR) - TD->GetPointY(iL);

 double dx=TD->GetPointX(iR)-TD->GetPointX(iStep);
 double dy=DeltaY/DeltaX*dx;

 double oldY=TD->GetPointY(iStep);
 TD->SetPointY(TD->GetPointY(iStep+1)-dy,iStep);
 TD->NewMinMax(); 
  QPoint D;
  D.setX((int)fLinTrafo(EffP.width()+EffP.x(),EffP.x(),PLP.DMax.x,PLP.DMin.x,
             TD->GetPointX(iStep)) );
  D.setY((int)fLinTrafo(EffP.height()+EffP.y(),EffP.y(),PLP.DMin.y,PLP.DMax.y,
             TD->GetPointY(iStep)) );

  //fprintf(stderr,"New point: %d / %d\n",D.x(),D.y());
  PaintSymbol(DIAMOND, 3, red, PlotF, D);
  char szA[MAX_LINELENGTH+1];  
  sprintf(szA,"#%d Y(%12.5g): %12.5g set to %12.5g",iStep,
                     TD->GetPointX(iStep),oldY,TD->GetPointY(iStep));
  InfoLB->setText(szA);

}
// **************************************************
void ManipData::RemovePoint()
{

  TD->Remove(iStep);

  PLP.iEnd--;
  //QPainter p;
  //p.begin(PlotF);
  //p.eraseRect(5,5,WPl.width()-5,WPl.height()-5);
  //p.end();
  setUpdatesEnabled(TRUE);
  update();
  char szA[MAX_LINELENGTH+1];
  sprintf(szA,"Point %d removed (%d remain)",iStep+1, TD->GetSteps());
  InfoLB->setText(szA);
  
}
// **************************************************
void ManipData::NewRange()
{ static int i1=-1,i2=-1;
  char szA[MAX_LINELENGTH+1]; 
  if(i1==-1)
    {i1=iStep;
     sprintf(szA,"First point %d for new range",i1);
     InfoLB->setText(szA);
     Bts[4]->setFont(QFont("Helvetica",iFS,QFont::DemiBold));
     Bts[4]->show();
     QPoint pp=CurCurs;
     pp.setY(PLP.PRect.y()-7);
     PaintSymbol(LEFT_CROSS,10,darkGreen,PlotF,pp);
     Bts[5]->setFont(QFont("Helvetica",iFS,QFont::Light));
     Bts[5]->show();
     return;
    }    
  if(i2==-1)i2=iStep;
  if(i1>i2)SWAP(int,i1,i2)
  if(i2-i1<=2)
    {sprintf(szA,"Range to small %d to %d",i1,i2);
     i1=i2=-1;
     InfoLB->setText(szA);
     return;
    }    
   
   sprintf(szA,"New range from %d to %d",i1,i2);
   InfoLB->setText(szA);

   iRCount++;
   QPoint *p= new QPoint(PLP.iStart,PLP.iEnd);
   Range.push(p);
   PLP.iStart=i1;
   PLP.iEnd=i2;
   i1=i2=-1;
   Bts[4]->setFont(QFont("Helvetica",iFS,QFont::Normal));
   Bts[4]->show();
   Bts[5]->setFont(QFont("Helvetica",iFS,QFont::DemiBold));
   Bts[5]->show();
   iStep=i1+1;
   bRSelect=TRUE;
   MoveCursor(red);
   PlotF->setUpdatesEnabled(TRUE);
   PlotF->update();

   setUpdatesEnabled(TRUE);
   repaint();
  
}
// **************************************************
void ManipData::RangeBack()
{char szA[MAX_LINELENGTH+1]; 
  if(!iRCount || !bRSelect) 
    { sprintf(szA,"Cannot go back to previous range");
     InfoLB->setText(szA);
     return; 
    }

   iRCount--;
   QPoint *p;
   p=Range.pop();
   if(!p)
     {PRINT_DEBUG("Uups, unexpected error popping range stack\n")
     } 
   PLP.iStart=p->x();
   PLP.iEnd=p->y();
   iStep=PLP.iStart+1;
   //bRSelect=FALSE;
 
   sprintf(szA,"Back one range from %d to %d",p->x(),p->y());
   InfoLB->setText(szA);
   if(!iRCount)
     {Bts[5]->setFont(QFont("Helvetica",iFS,QFont::Light));
      Bts[5]->show();
      }
   PlotF->setUpdatesEnabled(TRUE);
   PlotF->update();

   setUpdatesEnabled(TRUE);
   repaint();
 
}
// **************************************************
void ManipData::Integrate()
{MDATA SumI=0,SumU=0,Dx,y; 
 int i;

 for(i=PLP.iStart+1;i<PLP.iEnd;i++)
   {Dx=TD->GetPointX(i)-TD->GetPointX(i-1);
    y=0.5*(TD->GetPointY(i-1)+TD->GetPointY(i))*Dx;
    SumI+=y;
    SumU+=y*TD->GetPointX(i);
   }
 y=SumU/SumI;
 char szA[MAX_LINELENGTH+1]; 
 sprintf(szA,"Integral: %-14.7g Mass center: %-14.7g",SumI,y);
 InfoLB->setText(szA);
  
}
// **************************************************
void ManipData::ErrorBars()
{setUpdatesEnabled(TRUE);
 if(PLP.pSymbol[0]!=ERROR_BAR)
   {PLP.pSymbol[0]=ERROR_BAR;
    Bts[1]->setText("Symbols F2");
    PlotF->update();
    repaint();
    return;
   }
 else
   {PLP.pSymbol[0]=0;
    Bts[1]->setText(szButLab[1]);
   }
 PlotF->update();
 repaint();
 
}
// **************************************************
void ManipData::SaveRange()
{TD->NewRowRange(PLP.iStart,PLP.iEnd);
 Close();
}
// **************************************************
void ManipData::Close()
{if(iQuit==END_QUIT)
   {if(!DF->SaveData(stdout))
        fprintf(stderr,"ERROR: Can not write to stdout\n");
    qApp->quit();
   }			  
 else hide();
}
// ******************************************************
void ManipData::paintEvent(QPaintEvent *)
{
 QPainter p;
 p.begin(PlotF);
 QFont PF("Helvetica",14);
 //if(WidgetAusw::ScreenHeight()==600)PF.setPointSize(11);
 if(SCREEN.y()==600)PF.setPointSize(11);

 char szX[80],szY[80];
 strncpy(szX,DF->GetColID(TD->GetColX()),79);
 strncpy(szY,DF->GetColID(TD->GetColY()),79);
    PLP.PFont=PF;              //QFont("Helvetica",13);
    //PLP.PRect=QRect(WPl);
    WPl.setWidth(PlotF->width());
    WPl.setHeight(PlotF->height());
    PLP.PRect=QRect(0,0,WPl.width(),WPl.height());
    PLP.xLab=szX;
    PLP.yLab=szY;
    PLP.pFlag=TOP_FILENAME,
 DrawData(p,DF,NULL,PLP);
 EffP=PLP.PRect;
 p.end();
 
 MoveCursor(red);

 
}
// *************************************************
void ManipData::keyPressEvent(QKeyEvent *e)
{ 
 switch(e->key()) 
       {case Key_Left:if(iStep>PLP.iStart)iStep--;
                      else {if(iStep==PLP.iStart)iStep=PLP.iEnd-PLP.iStart-1;}
                      MoveCursor(red); 
                      break;
       case Key_Right:if(iStep<PLP.iEnd)iStep++;
                      if(iStep==PLP.iEnd)iStep=PLP.iStart;
                      MoveCursor(red);
                      break;
	case Key_Home:iStep=PLP.iStart;
                      MoveCursor(red);
                      break;
         case Key_End:iStep=PLP.iEnd-1;
                      MoveCursor(red);
                      break; 
          case Key_F2:ErrorBars();
                      break;
          case Key_F3:MeanVal();
                      return;
          case Key_F4:RemovePoint();
                      return;
          case Key_F5:NewRange();
                      return;
          case Key_F6:RangeBack();
                      return;
          case Key_F7:Integrate();
                      return;
         case Key_F10:Close();
                      return;
 	case Key_Plus:if(PLP.pSymbol[0]>=0)PLP.pSize[0]++;
	              //fprintf(stderr,"Sybolsize[%d] + : %d\n",iP,PLP.pSize[0]);
                      PlotF->update();
		      repaint();
		      break;   
 	case Key_Minus:if(PLP.pSymbol[0]>=0)PLP.pSize[0]--;
	              //fprintf(stderr,"Sybolsize[%d] - : %d\n",iP,PLP.pSize[0]);
                      PlotF->update();
		      repaint();
		      break;   
          default:e->ignore();
                   return;
	  }//switch
  char szB[80];
  sprintf(szB,"Point: %d of %d X:%12.5g Y:%12.5g Prev. Ranges:%d",
               iStep+1,TD->GetSteps(),
               TD->GetPointX(iStep), TD->GetPointY(iStep),iRCount);
  InfoLT->setText(szB);  
}
// *************************************************
void ManipData::MoveCursor(const QColor c)
{ QPainter p;
  const int size=10;
  static QPoint POld;
  
  CurCurs.setX((int)fLinTrafo(PLP.PRect.width()+PLP.PRect.x(),PLP.PRect.x(),
                     PLP.DMax.x,PLP.DMin.x,TD->GetPointX(iStep)));

  CurCurs.setY((int)fLinTrafo(PLP.PRect.y(),PLP.PRect.height()+PLP.PRect.y(),
                     PLP.DMax.y,PLP.DMin.y,TD->GetPointY(iStep)));

  PM0.setX(CurCurs.x()-size-1);  PM0.setY(CurCurs.y()-size-1);
  
  p.begin(PlotF);
  if(!First)
    {p.setPen(PlotF->backgroundColor());
     DrawSymbol(CROSS, size, &p, POld);
     //fprintf(stderr,"%d  Pixmap to: %d / %d\n\n",iStep,POld.x(),POld.y());
     p.drawPixmap(POld,CursPM);
    }        
   else
    {First=FALSE;}

     CursPM=QPixmap::grabWindow(PlotF->winId(),PM0.x(),PM0.y(),
                                2*size+2,2*size+2);
     POld=PM0;
     //p.drawRect(PM0.x(),PM0.y(),2*size,2*size);
     //fprintf(stderr,"%d Pixmap from: %d / %d\n",iStep, PM0.x(),PM0.y());
     //p.drawPixmap(100,100,CursPM);

    p.setPen(c);
    DrawSymbol(CROSS, size, &p, CurCurs);
    
  p.end();

}

// ***************************************************
// ***************************************************
// MTrackWid functions
// ***************************************************
MTrackWid::MTrackWid(QWidget *parent, const char *name):
		     QFrame(parent,name)
{ RDown=FALSE;
  setCursor(crossCursor);
}
// ***************************************************
void MTrackWid::mousePressEvent(QMouseEvent *e)
{if(e->button()==RightButton)RDown=TRUE;

}
 // ***************************************************
void MTrackWid::mouseReleaseEvent(QMouseEvent *e)
{if(e->button()==RightButton)RDown=FALSE;
 if(e->button()==LeftButton)emit LCoords(QPoint(e->pos()));
 
}
// ***************************************************
void MTrackWid::mouseMoveEvent(QMouseEvent *e)
{
 if(RDown)
   {emit RCoords(QPoint(e->pos()));
   }
}
// ***************************************************
// *****************************************************
//ShowWid
// *****************************************************
  int ShowWid::w=0;
  int ShowWid::h=0;
// **************************************************
// DisplayData Functions
// **************************************************
ShowWid::ShowWid(DataFile *df, QWidget *parent, const char *name):
              QWidget(parent,name)
{
 ShowWid::w=(SCREEN.x()/4);
 ShowWid::h=(SCREEN.y()/4);

 iQuit = END_CLOSE;
 iColX = -1;
 iColY = -1;

 TD = df;
 if(!df)return;
 if(!TD->GetCRData())TD=0;
 InitPlot();
}
// ********************************************************
void ShowWid::InitPlot(void)
{
  setCaption(TD->GetFileName());
  setGeometry(SCREEN.x()-w-10,0,w,h);

  iColX=TD->GetCRData()->GetColX();
  iColY=TD->GetCRData()->GetColY();
  strncpy(szXLabel,TD->GetColID(iColX),79);
  strncpy(szYLabel,TD->GetColID(iColY),79);

  enum FileType FileT=TD->GetFType();
  if(FileT==FtSXS_NEW || FileT==FtSXS_OLD || FileT==FtXDIF_RAW)
                                                   PLP.pSymbol[0]=ERROR_BAR;  
  else PLP.pSymbol[0]=CIRCLE;

  PLP.xLab=szXLabel;
  PLP.yLab=szYLabel;
  PLP.PRect.setRect(3,3,width()-3,height()-3);
  PLP.pFlag = 0;
  PLP.iStart=0;
  PLP.iEnd=-1;
  show();
 }
// ******************************************************
void ShowWid::Close()
{if(iQuit==END_QUIT) qApp->quit();
 else hide();
}
// ******************************************************
void ShowWid::keyPressEvent(QKeyEvent *e)
{ switch(e->key()) 
       {case Key_C:if(PLP.pSymbol[0]!=CONNECT)PLP.pSymbol[0]=CONNECT;
	           else PLP.pSymbol[0]=CIRCLE;
		   repaint();
		   break;
  	case Key_E:if(PLP.pSymbol[0]!=ERROR_BAR)PLP.pSymbol[0]=ERROR_BAR;
	           else PLP.pSymbol[0]=CIRCLE;
		   repaint();
		   break;
 	case Key_Plus:if(PLP.pSymbol[0]>=0)PLP.pSize[0]++;
		      repaint();
		      break;   
 	case Key_Minus:if(PLP.pSymbol[0]>=0)PLP.pSize[0]--;
		      repaint();
		      break;   
        case Key_P:PrintIt();
		   break;
        case Key_Q:
        case Key_Escape:Close();
	                break;
        default:e->ignore();
               return;
       }//switch
}
// **************************************************
void ShowWid::paintEvent(QPaintEvent *)
{
 QPainter paint;
 if(!paint.begin(this))
   {  PRINT_DEBUG("Painter: begin failed\n")}
 int iF=height()/17;
 if(iF==9)iF=10;
 if(iF<8)iF=8;
 PLP.PFont=QFont("Helvetica",iF);
 //fprintf(stderr,"Fontsize: %d\n",height()/17);
 PLP.PRect.setRect(3,3,width()-3,height()-3);
 DrawData(paint,TD,NULL,PLP);
 paint.end();
}
// ******************************************************
void ShowWid::PrintIt(void)
{
  QPrinter Prt;
  if(Prt.setup(this))
   {QPainter paint;
    paint.begin(&Prt);
int iF=height()/17;
 if(iF==9)iF=10;
 if(iF<8)iF=8;
 PLP.PFont=QFont("Helvetica",iF);
 //fprintf(stderr,"Fontsize: %d\n",height()/17);
 PLP.PRect.setRect(3,3,width()-3,height()-3);
 DrawData(paint,TD,NULL,PLP);
    paint.end();
   }
  //else fprintf(stderr,"Printer setup failed\n");

}
// ******************************************************
// ***************************************************
// *****************************************************
//TheView
// *****************************************************
  int TheView::w=0;
  int TheView::h=0;
  QRect TheView::WPlot;
  QRect TheView::WFrame;

  
// **************************************************
// DisplayData Functions
// **************************************************
TheView::TheView(DataFile *df, DataFile *af, QWidget *parent, const char *name):
              QWidget(parent,name)
{
  TheView::w=SCREEN.x()/2;
  TheView::h=SCREEN.y()-50;

 TheView::WFrame.setRect(10, 0,w-5-2*5,h-4*5);
 TheView::WPlot.setRect(1,1,w/2-5,h/4-5);
 iQuit = END_CLOSE;

 TD = df;
 if(!df)return;
 if(!TD->GetCRData())TD=0;
 AD = af;
 if(!af)return;
 if(!AD->GetCRData())AD=0;

 Tmp=new TheFile(TD->GetFileName());
 if( Tmp->ReadData(1)>0 || Tmp->GetCRData()->GetSteps()==0 )
   {fprintf(stderr,"ERROR Theview: Error reading data\n");
    QApplication::exit(EXIT_FAILURE);
   }
 Tmp->GetCRData()->Deriv(3,3,0);
 Tmp->SetColID("dTs/dt",0); 
 InitPlot();
}
// ********************************************************
void TheView::InitPlot(void)
{
  setCaption(TD->GetFileName());
  setGeometry(SCREEN.x()-w-10,0,w,h);

  iColX[0]=1;
  iColY[0]=2;

  iColX[1]=4;
  iColY[1]=1;

  iColX[2]=2;
  iColY[2]=7;

  iColX[3]=4;
  iColY[3]=1;

  iColX[4]=2;
  iColY[4]=5;

  iColX[5]=2;
  iColY[5]=6;

  iColX[6]=2;
  iColY[6]=4;

  iColX[7]=2;
  iColY[7]=3;


  PLP.pSymbol[0]=CIRCLE;
  PLP.pSize[0]=0;
  PLP.pFlag = 0;
  PLP.iStart=0;
  PLP.iEnd=-1;
  //int iF=height()/17;
  //if(iF==9)iF=10;
  //if(iF<8)iF=8;
  //fprintf(stderr,"Font: %d\n",iF); 
 PLP.PFont=QFont("Helvetica",8);
  
  //PrintRect(WFrame,"Frame");
  //PrintRect(WPlot,"Plot");
  PlotF = new QFrame(this,"PlotFrame");
  PlotF->setGeometry(WFrame);
  PlotF->setFrameStyle( QFrame::Panel | QFrame::Sunken);
  PlotF->setBackgroundColor(white);
  PlotF->setLineWidth(3);
  PlotF->setUpdatesEnabled(TRUE);

  show();
 }
// ******************************************************
void TheView::Close()
{if(iQuit==END_QUIT) qApp->quit();
 else hide();
}
// ******************************************************
void TheView::keyPressEvent(QKeyEvent *e)
{ switch(e->key()) 
       {case Key_C:if(PLP.pSymbol[0]!=CONNECT)PLP.pSymbol[0]=CONNECT;
	           else PLP.pSymbol[0]=CIRCLE;
		   repaint();
		   break;
 	case Key_Plus:if(PLP.pSymbol[0]>=0)PLP.pSize[0]++;
		      repaint();
		      break;   
 	case Key_Minus:if(PLP.pSymbol[0]>=0)PLP.pSize[0]--;
		      repaint();
		      break;   
        case Key_P:PrintIt();
		   break;
        case Key_Q:
        case Key_Escape:Close();
	                break;
        default:e->ignore();
               return;
       }//switch
}
// **************************************************
void TheView::paintEvent(QPaintEvent *)
{
int i,j,s;
QPainter paint;
if(!paint.begin(PlotF))
          {PRINT_DEBUG("Painter: begin failed\n")}
for(i=0,s=0;i<4;i++)
   {for(j=0;j<2;j++)
       {DataFile *df; 
        df=AD;
        if(s==0 || s==1)df=TD;
        if(s==3)df=Tmp;
        strncpy(szXLabel,df->GetColID(iColX[s]-1),79);
        strncpy(szYLabel,df->GetColID(iColY[s]-1),79);
        PLP.xLab=szXLabel;
        PLP.yLab=szYLabel;
        PLP.PRect.setRect(3+WPlot.width()*j,3+WPlot.height()*i,WPlot.width()-3,WPlot.height()-3);
        //PrintRect(PLP.PRect);
        df->GetCRData()->SetColX(iColX[s]-1);
        df->GetCRData()->SetColY(iColY[s]-1);

        DrawData(paint,df,NULL,PLP);
        s++;
       }
   }
paint.end();
}
// ******************************************************
void TheView::PrintIt(void)
{
  QPrinter Prt;
  if(Prt.setup(this))
   {QPainter paint;
    //paint.begin(&Prt);
    if(!paint.begin(&Prt))
          {PRINT_DEBUG("Painter (print): begin failed\n")}
    int i,j,s;
   for(i=0,s=0;i<4;i++)
    {for(j=0;j<2;j++)
       {DataFile *df; 
        df=AD;
        if(s==0 || s==1)df=TD;
        if(s==3)df = Tmp;
        strncpy(szXLabel,df->GetColID(iColX[s]-1),79);
        strncpy(szYLabel,df->GetColID(iColY[s]-1),79);
        PLP.xLab=szXLabel;
        PLP.yLab=szYLabel;
        PLP.PRect.setRect(20+WPlot.width()*j,10+WPlot.height()*i,WPlot.width()-20,WPlot.height()-10);
        //PrintRect(PLP.PRect);
        df->GetCRData()->SetColX(iColX[s]-1);
        df->GetCRData()->SetColY(iColY[s]-1);

        DrawData(paint,df,NULL,PLP);
        s++;
       }
     }
   String Answer("");
   TD->SPrintInfo(Answer,PR_BASIC);
   QFont hf(PLP.PFont);
   hf.setPointSize(8);
   DrawFormText(paint,hf,50,5+WPlot.height()*i,
                AlignTop|AlignLeft,Answer.GetBuf());
   Answer.Set("");
    //paint.drawText(20+WPlot.width(),20+WPlot.height()*i,TD->GetFileName());
   TD->SPrintInfo(Answer,PR_TIME);
   DrawFormText(paint,hf,50,30+WPlot.height()*i,
                AlignTop|AlignLeft,Answer.GetBuf());

   Answer.Set("");
   TD->SPrintInfo(Answer,PR_HEAD);
   DrawFormText(paint,hf,20+WPlot.width(),5+WPlot.height()*i,
                AlignTop|AlignLeft,Answer.GetBuf());

    paint.end();
   }
  //else fprintf(stderr,"Printer setup failed\n");

}
// ******************************************************
// *****************************************************
//ShowTxt
// *****************************************************
  int ShowTxt::w=0;
  int ShowTxt::h=0;
// **************************************************
// DisplayData Functions
// **************************************************
ShowTxt::ShowTxt(const char *szF, QRect *Pos, QWidget *parent, const char *name):
                 QWidget(parent,name)
{
 ShowTxt::w=(SCREEN.x()/4);
 ShowTxt::h=(SCREEN.y()/2);

 *szFile=0;
 
 FILE *f=fopen(szF,"r");
 CHECK_FILE_POINTER_RET(f,szF)
  
 strncpy(szFile,szF,MAXPATH);		  
 setCaption(szFile);
 if(!Pos)setGeometry(SCREEN.x()-w-10,0,w,h);
 else setGeometry(*Pos);
 
 Txt=new QMultiLineEdit(this,"TxT");
 Txt->setReadOnly(TRUE);
 Txt->setGeometry(3,3,width()-6,height()-6);
 Txt->setAutoUpdate(TRUE);
 char szB[MAX_LINELENGTH+1];
 while(1)
      {if(!fgets(szB,MAX_LINELENGTH,f))break;
       RemoveCR_LF(szB);
       Txt->append(szB);
      }
 fclose(f); 
 show();
 }
// ******************************************************
void ShowTxt::keyPressEvent(QKeyEvent *e)
{ switch(e->key()) 
       { case Key_Q:
  	 case Key_Escape:Close();
	                break;
        default:e->ignore();
               return;
       }//switch
}
// ********************************************************
void ShowTxt::resizeEvent(QResizeEvent *)
{ Txt->setGeometry(3,3,width()-6,height()-6);
  Txt->show();
}
// ********************************************************
void ShowTxt::Close()
{emit Removed();
 if(iQuit==END_QUIT)qApp->quit();
 else hide();
}
// ******************************************************




