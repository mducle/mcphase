//File: sfit.cpp
//$Id: sfit.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
//$Log: sfit.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.2  1999/03/15 09:08:37  herbie
//*** empty log message ***
//

#include <qtooltip.h>
#include <qkeycode.h>
#include <qprinter.h>
#include <qapplication.h>

#include <stdio.h>
#include <stdlib.h>

#include "sfit.h"
#include "simplex.h"
#include "grafic.h"
#include "stdfunc.h"

#define PARAMETER_FILE "ausw.conf"

extern char szAnswer[2048];

extern QList<DataFile> PlotList;

extern QPoint SCREEN;

// #define MAX_FITFN 5
// const struct FitFn FIT_FUNCS[MAX_FITFN]=
//        {{"Lorentzian",1,3,Lorentz,szLorentzP}, 
//         {"Lorentz Alfa 12",2,2, Lorentz12,szLorentz12P},
// 	{"Curie-Weiss",3,0, CurieWeiss,szCurieWeissP},
// 	{"NonFermiRes",3,0, NonFermiRes,szNonFermiResP},
// 	{"Straight line",2,0, StraightLine,szStraightLineP}
//        };

const FitFunction  * FF[MAX_FITFN]=
      { new LorentzFn(),     new LorentzA12Fn(),  new CurieWeissFn(),
        new NonFermiResFn(), new StraightLineFn(), new PolynomFn()
      };


extern QRect B;  //defined in grafic.cpp
//extern char *szLorentz12P[MAX_NAME_LEN];
//extern char *szLorentzP[MAX_NAME_LEN];
//extern char *szCurieWeissP[MAX_NAME_LEN];

// *********************************************************
 int SFitWid::w=0;
 int SFitWid::h=0;
        
 QRect SFitWid::MainW;
 QRect SFitWid::FitFnsW;   //ListBox of Fit Fns
 QRect SFitWid::ParTextW;  //Label fo Parameter names
 QRect SFitWid::ClBW;      //Close button
 QRect SFitWid::ShBW;      //Show button
 QRect SFitWid::FitBW;     //FitButton
 QRect SFitWid::PrBW;      //PrintButton
 QRect SFitWid::StBW;      //StopButton
 QRect SFitWid::PrgW;      //Fit Progress bar
 QRect SFitWid::ResultW;   //Frame displaying results
 QRect SFitWid::FitW;      //Plot of fit data
 QRect SFitWid::ErrW;      //Plot of error data
 QRect SFitWid::NoCurW;    //LineEdit of # of Cures
 QRect SFitWid::MaxItW;    //LineEdit of # of Iterations
 QRect SFitWid::ParMEdW;   //MultiLineEdit of Pars 
 QRect SFitWid::StepMEdW;  //MultiLineEdit of Steps    
 QRect SFitWid::LimitMEdW; //MultiLineEdit of Limits    

   int SFitWid::iChx=0;
   int SFitWid::iChy=0;
   int SFitWid::nMaxPars=0;

// ********************************************************* 
SFitWid::SFitWid(DataFile *df, QWidget *parent, const char *name): 
                  QWidget(parent,name)
{//int sh=WidgetAusw::ScreenHeight();
 w = (SCREEN.y()==600 ? 750 : 900);
 h = (SCREEN.y()==600 ? 550 : 700);

    SI = NULL;
    OF = NULL;
    TD = NULL;
 iQuit = END_CLOSE;
 
 if(df==NULL)return;
    DF = df;
 
 if(DF->SaveData(DF->GetInfoText(),"ausw.fit")==0)
   {PRINT_DEBUG("Can not create ausw.fit !\n")
    return;
   }
   
 OF = new AsciiFile("ausw.fit");
 CHECK_POINTER_RET(OF)
 if(OF->ReadData()>0)
   {PRINT_DEBUG("Can not read ausw.fit !\n")
    return;
   }
 
 EF = new AsciiFile("ausw.fit");
 CHECK_POINTER_RET(EF)
 if(EF->ReadData()>0)
   {PRINT_DEBUG("Can not read ausw.fit !\n")
    return;
   }
 
 TD = OF->GetCRData();
 CHECK_POINTER_RET(TD)

 ED = EF->GetCRData();
 CHECK_POINTER_RET(ED) // Only to See if evrthing is OK;
 ED = NULL; // No error data available; fit has not performed
   
 const int LeftO0=10;
 BF = QFont("Helvetica",(SCREEN.y()==600 ? 12 : 14),QFont::DemiBold); // for buttons
 WF = QFont("Helvetica",(SCREEN.y()==600 ? 12 : 14),QFont::Normal); // for sub window text
 PF = QFont("Helvetica",10);

 QFontMetrics FM(WF);
 iChx = FM.width('0');
 iChy = FM.lineSpacing();

   nFns = MAX_FITFN;
 iSelFn = -1;

 MainW.setRect(10, 10,w,h);
 const int nF = (SCREEN.y()==600 ? 3 : 5);
 FitFnsW.setRect(LeftO0, 30,20*iChx,nF*3*iChy/2);
 NoCurW.setRect(LeftO0, FitFnsW.y()+FitFnsW.height()+2*B.height()+10,
                FitFnsW.width(),B.height());
 MaxItW.setRect(LeftO0, NoCurW.y()+NoCurW.height()+2*B.height()+10,
                FitFnsW.width(),B.height());

 const int LeftO1 = LeftO0+FitFnsW.width()+20;

  int iMH = MainW.height()-(MaxItW.y()+MaxItW.height()+20+10);
 nMaxPars = iMH/iChy-1;
 
 ParTextW.setRect(LeftO0,MaxItW.y()+MaxItW.height()+B.height()+15,
                  iChx*12,nMaxPars*iChy);
 ParMEdW.setRect(ParTextW.x()+ParTextW.width()+5,ParTextW.y(),
                 ParTextW.width(),ParTextW.height());
 StepMEdW.setRect(ParMEdW.x()+ParMEdW.width()+5,ParTextW.y(),
                 ParTextW.width(),ParTextW.height());
 LimitMEdW.setRect(StepMEdW.x()+StepMEdW.width()+5,ParTextW.y(),
                 ParTextW.width(),ParTextW.height());
		  
// ParTextW.setRect(LeftO1,10,iChx*15,20*iChy);


 ClBW.setRect(LeftO1,10,B.width(),B.height());
 ShBW.setRect(LeftO1,ClBW.y()+ClBW.height()+15,B.width(),B.height());
 FitBW.setRect(LeftO1,ShBW.y()+ShBW.height()+15,B.width(),B.height());
 PrBW.setRect(LeftO1,FitBW.y()+FitBW.height()+25,B.width(),B.height());
 PrgW.setRect(LeftO1, ParTextW.y()-3*B.height()-5,
              2*B.width(),B.height());

 StBW.setRect(LeftO1+B.width()+10, PrgW.y()-4*B.height(),
              B.width(),B.height());

 ResultW.setRect(LimitMEdW.x()+LimitMEdW.width()+10,ClBW.y(),
                 MainW.width()-LimitMEdW.x()-LimitMEdW.width()-20,
		 MainW.height()-2*ClBW.y());
 FitW.setRect(1,1,ResultW.width(), LimitMEdW.y()-ClBW.height()-10);

 ErrW.setRect(0,FitW.height()-5, ResultW.width(), FitW.height()/2);

  setGeometry(MainW);
  setFixedSize(w,h);
  setCaption("Fit function to data");
  // ------------------------------------  
  //Label Function List Box
  QLabel *FnLbL = new QLabel("Fit Funs:", this);
  FnLbL->setAlignment(AlignLeft);
  FnLbL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  FnLbL->setGeometry(LeftO0,10,FitFnsW.width(),B.height());

  //Listbox for fit functions
  FitFns = new QListBox( this , "FnLB");
  FitFns->setGeometry(FitFnsW);
  FitFns->setFont(WF);
  int i;
  //for(i=0; i<nFns; i++)FitFns->insertItem(FIT_FUNCS[i].FnName);
  for(i=0; i<nFns; i++)FitFns->insertItem(FF[i]->GetFnName());
  connect( FitFns, SIGNAL(selected(int)),SLOT(FunctionSel(int)) );
  QToolTip::add( FitFns, "List of available functions" );
  // ------------------------------------  
  //Label: Number of curves single line edit window
  QLabel *CurvL=new QLabel("Number of curves:", this);
  CurvL->setAlignment(AlignLeft);
  CurvL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  CurvL->setGeometry(LeftO0,FitFnsW.y()+FitFnsW.height()+20,
                     NoCurW.width(),B.height());

  //Single Line edit: No of curves
  NoCurEd = new QLineEdit( this );
  NoCurEd->setGeometry(NoCurW);
  NoCurEd->setText("1");
  connect( NoCurEd, SIGNAL(returnPressed()),SLOT(NoCurvChg()) );
  QToolTip::add( NoCurEd, "Enter # of curves" );
  nCurv=1;
  // ------------------------------------  
  nMaxIt=1000;
  //Label: Number of Max Iterations single line edit window
  QLabel *MItL = new QLabel("Max. number of iterations:", this);
  MItL->setAlignment(AlignLeft);
  MItL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  MItL->setGeometry(LeftO0,NoCurW.y()+NoCurW.height()+20,
                     NoCurW.width(),B.height());

  //Single Line edit: Max. no of iterations
  MaxItEd = new QLineEdit( this );
  MaxItEd->setGeometry(MaxItW);
  MaxItEd->setText("1000");
  QToolTip::add( MaxItEd, "Enter Max. # of iterations" );
 
  //FitProgress
  Prg = new QProgressBar( this );
  Prg->setGeometry(PrgW);
  QToolTip::add( Prg, "Fit progress" );
  Prg->hide();

  QLabel *PrgL = new QLabel("Fit progress", this);
  PrgL->setAlignment(AlignLeft);
  PrgL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  PrgL->setGeometry(PrgW.x(),PrgW.y()-B.height(),B.width(),B.height());
  //PrgL->hide();

  nPars = 0;
 ncPars = 0;
    lfP = NULL;
    lfL = NULL;
   lfSt = NULL;
  // ------------------------------------  
  // Names of paramters
  //Label: For parameter names 
  QLabel *PNaL = new QLabel("Parameter name:", this);
  PNaL->setAlignment(AlignLeft);
  PNaL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  PNaL->setGeometry(ParTextW.x(),MaxItW.y()+MaxItW.height()+10,
                     ParTextW.width(),B.height());

  ParText = new QLabel(this);
  ParText->setAlignment(AlignRight);
  ParText->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  ParText->setGeometry(ParTextW);
  ParText->setFont(WF);
  ParText->setText("Names of\nparamters\nof selected\nfunction");
  QToolTip::add( ParText, "List of pars for selected function" );

  // ------------------------------------  
  // Values of paramters
  //Label: For parameter values 
  QLabel *PVaL = new QLabel("Initial value:", this);
  PVaL->setAlignment(AlignLeft);
  PVaL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  PVaL->setGeometry(ParMEdW.x(),MaxItW.y()+MaxItW.height()+10,
                     ParTextW.width(),B.height());

  ParMEd = new MLineEd(this);
  ParMEd->setGeometry(ParMEdW);
  ParMEd->setFont(WF);
  ParMEd->setCellHeight(iChy);
  connect( ParMEd, SIGNAL(returnPressed()),SLOT(ChkLength()) );
  QToolTip::add( ParMEd, "Enter initial parameter values for fit" );

  // ------------------------------------  
  // Values of start steps
  //Label: For start steps values 
  QLabel *PStL = new QLabel("Initial steps:", this);
  PStL->setAlignment(AlignLeft);
  PStL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  PStL->setGeometry(StepMEdW.x(),MaxItW.y()+MaxItW.height()+10,
                     ParTextW.width(),B.height());

  StepMEd = new MLineEd(this);
  StepMEd->setGeometry(StepMEdW);
  StepMEd->setFont(WF);
  StepMEd->setCellHeight(iChy);
  //connect( StepMEd, SIGNAL(returnPressed()),SLOT(ChkLength()) );
  QToolTip::add( StepMEd, "Correct initial step values for fit" );

  // ------------------------------------  
  // Values of limits
  //Label: For limits
  QLabel *PLiL = new QLabel("Limits:", this);
  PLiL->setAlignment(AlignLeft);
  PLiL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  PLiL->setGeometry(LimitMEdW.x(),MaxItW.y()+MaxItW.height()+10,
                     ParTextW.width(),B.height());

  LimitMEd = new MLineEd(this);
  LimitMEd->setGeometry(LimitMEdW);
  LimitMEd->setFont(WF);
  LimitMEd->setCellHeight(iChy);
  //connect( LimitMEd, SIGNAL(returnPressed()),SLOT(ChkLength()) );
  QToolTip::add( LimitMEd, "Correct initial step values for fit" );

  // ------------------------------------  
  // Close Button
  ClB = new QPushButton("Close", this, "close" );
  ClB->setFont(BF);
  ClB->setFixedHeight(ClB->sizeHint().height() );
  ClB->setGeometry(ClBW);
  connect( ClB, SIGNAL(clicked()), SLOT(Close()) );
  QToolTip::add( ClB, "Close (exit) this window" );
  // ------------------------------------  

  // Show Button
  ShB = new QPushButton("Show", this, "show" );
  ShB->setFont(BF);
  ShB->setFixedHeight(ShB->sizeHint().height() );
  ShB->setGeometry(ShBW);
  connect( ShB, SIGNAL(clicked()), SLOT(Show()) );
  QToolTip::add( ShB, "Show function with initial pars" );
  ShB->hide();
  // ------------------------------------  

  // Fit Button
  FitB = new QPushButton("Fit", this, "fit" );
  FitB->setFont(BF);
  FitB->setFixedHeight(FitB->sizeHint().height() );
  FitB->setGeometry(FitBW);
  connect( FitB, SIGNAL(clicked()), SLOT(Fit()) );
  QToolTip::add( FitB, "Start fit with initial pars" );
  FitB->hide();
  // ------------------------------------  

  // CheckBox to animate fit
  Anim = new QCheckBox("Animate fit",this,"anim");
  Anim->setGeometry(FitBW.x()+FitBW.width()+10,FitBW.y(),B.width(),B.height());
  QToolTip::add( Anim, "Check it for animated fit" );
  Anim->hide();
  // ------------------------------------
   
  // Print Button
  PrB = new QPushButton("Print", this, "print" );
  PrB->setFont(BF);
  PrB->setFixedHeight(PrB->sizeHint().height() );
  PrB->setGeometry(PrBW);
  connect( PrB, SIGNAL(clicked()), SLOT(PrintIt()) );
  QToolTip::add( PrB, "Print grafic" );
  // ------------------------------------  

  // Stop Button
  StB = new QPushButton("Stop", this, "stop" );
  StB->setFont(BF);
  StB->setFixedHeight(StB->sizeHint().height() );
  StB->setGeometry(StBW);
  connect( StB, SIGNAL(clicked()), SLOT(StopFit()) );
  QToolTip::add( StB, "Stop fit" );
  //StB->hide();
  bStop=0;
  
  t=new QTimer(this);
  connect( t, SIGNAL(timeout()), SLOT(Iterate()) );
  FitEnd = 0;
     lfH = NULL;
  // ------------------------------------  

  // Grafic display of fitted data
  PlotList.setAutoDelete(FALSE);
  if(PlotList.count())PlotList.clear();
  PlotList.append(DF);
  PlotList.append(OF);

  FitData = new QFrame(this,"FitFrame");
  FitData->setGeometry(ResultW);
  FitData->setFrameStyle( QFrame::Panel | QFrame::Sunken);
  FitData->setBackgroundColor(white);
  FitData->setFocusPolicy(QWidget::ClickFocus);
  //repaint();
  show();

}
// ************************************************************
void SFitWid::Close()
{if(iQuit==END_QUIT)
   {fprintf(stdout,"Fit Result:\n%s",szOut);
    qApp->quit();
   } 
 else hide();
}
// ************************************************************
void SFitWid::StopFit()
{//fprintf(stderr,"Stop pressed\n");
 bStop = 1;
}
// ************************************************************
void SFitWid::Show()
{ if(!CollectPars())return;
  int nS = TD->GetSteps();
  int  i;
  for(i=0;i<nS; i++)
     {//TD->SetPointY(FIT_FUNCS[iSelFn].FitFunc(TD->GetPointX(i),nCurv,lfP),i);
      TD->SetPointY(FF[iSelFn]->GetFFStruct().FitFunc(TD->GetPointX(i),nCurv,lfP,ncPars,lfC),i);
     }
  repaint();   
}
// ************************************************************
void SFitWid::Fit()
{ if(!CollectPars())return;
  if(!CollectSteps())return;

  delete SI;
//  SI = new Simplex(DF->GetCRData(), nMaxIt,lfSt,lfL,nPars, lfP, nCurv,
//                   FIT_FUNCS[iSelFn]);
  SI = new Simplex(DF->GetCRData(), nMaxIt,lfSt,lfL,nPars, lfP, nCurv,
                   FF[iSelFn],ncPars,lfC);

  CHECK_POINTER_RET(SI)
  if(SI->GetStatus()!=STATUS_OK)
    {PRINT_DEBUG("Error creating Simplex\n")
     return;}
      
  //int i,j;

  delete [] lfH;
  lfH = new double [nPars];
  CHECK_POINTER_RET(lfH)
  ED = EF->GetCRData();  
  bStop = 0;
  Prg->show();
  FitB->hide();
  FitEnd = 0;
  t->start(0,FALSE);
}
// ************************************************************
void SFitWid::Iterate()
{
  int i;
  FitEnd = SI->Iterate(); //==0 it runs
  Prg->setProgress(100*SI->GetNIter()/nMaxIt);
       
//         fprintf(stderr,"Fit iteration %d of %d\n",SI->GetNIter(), nMaxIt);
//         for(i=0; i<nPars; i++)
// 	   {lfH[i]=SI->GetFitPar(i);
//              fprintf(stderr,"%s: %f error: %12.5g limit: %12.5g\n",
// 	             SI->GetPName(i),lfH[i],
// 		     SI->GetErr(i),SI->GetLimit(i));
//            }//for
// 	 fprintf(stderr,"Squaresum: %15.5g (%12.5g) error: %12.5g limit: %12.5g\n",
//          SI->GetFitPar(i),sqrt(SI->GetFitPar(i))/DF->GetCRData()->GetSteps(),
// 		     SI->GetErr(i),SI->GetLimit(i));
//  
  if(Anim->isChecked())
    {for(i=0; i<nPars; i++)lfH[i]=SI->GetFitPar(i);
     for(i=0;i<TD->GetSteps(); i++)
     //TD->SetPointY(FIT_FUNCS[iSelFn].FitFunc(TD->GetPointX(i),nCurv,lfH),i);
     TD->SetPointY(FF[iSelFn]->GetFFStruct().FitFunc(TD->GetPointX(i),nCurv,lfH),i);
     repaint();
    } 
  if(bStop)fprintf(stderr,"Fit stopped by user\n");

 //If fit is finished or stopped by user
  if(FitEnd || bStop)
    {t->stop();
     fprintf(stderr,"Iteration end (%d)\n",FitEnd);
     Prg->reset();
     Prg->hide();     	

     char *sp=szOut;
     szOut[0]=0;
     sp+=sprintf(sp,"Fit function %s: End at iteration %d of %d\n",
                     SI->GetFnName(),SI->GetNIter(), nMaxIt);
     for(i=0; i<nPars; i++)
        {lfH[i]=SI->GetFitPar(i);
         sp+=sprintf(sp,"%15.15s: %-14.7g Error: %-12.5g Limit: %-12.5g\n",
	             SI->GetPName(i),lfH[i],
		     SI->GetErr(i),SI->GetLimit(i));
        }//for
     double lfSoS=SI->GetFitPar(i);   
     sp+=sprintf(sp,"\nSquaresum of deviations: %-15.5g   ",lfSoS);

   
     ED=EF->GetCRData();  
     double lfM=0;
     for(i=0;i<TD->GetSteps(); i++)
        {//TD->SetPointY(FIT_FUNCS[iSelFn].FitFunc(TD->GetPointX(i),nCurv,lfH),i);
         TD->SetPointY(FF[iSelFn]->GetFFStruct().FitFunc(TD->GetPointX(i),nCurv,lfH,ncPars,lfC),i); 
        lfM+=TD->GetPointY(i);
        }
     lfM/=TD->GetSteps();
     
     double lfE,lfV=0;
     for(i=0;i<ED->GetSteps(); i++)
        {lfE=(FF[iSelFn]->GetFFStruct().FitFunc(TD->GetPointX(i),nCurv,lfH,ncPars,lfC) -
         DF->GetCRData()->GetPointY(i))/TD->GetPointY(i);
         ED->SetPointY(lfE,i);
         lfV+=SQUARE(TD->GetPointY(i)-lfM);
        }
 
     sp+=sprintf(sp,"r-value: %-15.5g\n", sqrt(fabs(1-lfSoS/lfV)) );
     sp+=sprintf(sp,"Mean deviation: %-12.5g Error: %-12.5g Limit: %-12.5g\n",
                     sqrt(SI->GetFitPar(nPars))/DF->GetCRData()->GetSteps(),
		     SI->GetErr(nPars),SI->GetLimit(nPars));
     
     //ED->AssignCol();    
     EF->SaveData("xausw.err");

     //TD->AssignCol();    
     OF->SaveData("lastfit");
     FitB->show();
     repaint();
     StB->hide();
    }
}
// ************************************************************
void SFitWid::ChkLength()
{
 if(ParMEd->IsNumFilled())
   {InitPars();
    int i;
    char *sp=szAnswer;
    *sp=0;
    for(i=0;i<nPars; i++)
       {lfP[i]=atof(ParMEd->textLine(i));
        lfL[i]=0.1*lfP[i];
	lfSt[i]=DEF_LIMIT;
        sp+=sprintf(sp,"%-14.7g\n",lfL[i]);
       }
     szAnswer[strlen(szAnswer)-1]=0; //remove last <lf>
     StepMEd->setText(szAnswer);

    sp=szAnswer;
    *sp=0;
    for(i=0;i<nPars; i++)sp+=sprintf(sp,"%-14.7g\n",lfSt[i]);
    szAnswer[strlen(szAnswer)-1]=0; //remove last <lf>
    LimitMEd->setText(szAnswer);    
    ShB->show();
    FitB->show();
    Anim->setStyle(MotifStyle);
    Anim->update();
    Anim->show();
   }
   
}
// ************************************************************
int SFitWid::CollectPars(void)
{int i;
 if(ParMEd->IsNumFilled())
   {for(i=0;i<nPars; i++)lfP[i]=atof(ParMEd->textLine(i));
     const RAW_DIF_HEAD *rdh=DF->GetDifHead();
    ncPars=0;
    if(rdh!=0 && 
       strcmp(FF[iSelFn]->GetFnName(),"Lorentz Alfa 12")==0)
        {ncPars=1;
         lfC[0]=rdh->fWLKa2/rdh->fWLKa1;
         //RWLa21=rdh->fWLKa2/rdh->fWLKa1;
        }
    if(strcmp(FF[iSelFn]->GetFnName(),"Lorentz Alfa 12")==0)
     {ncPars=1;lfC[0]=1.79278/1.78892;} //Cobalt as default	
    return 1;
   }
 else
   {QApplication::beep();
    fprintf(stderr,"Parameter window incomplete\n");
    return 0;
   } 
}
// ************************************************************
int SFitWid::CollectSteps(void)
{int i;
 if(StepMEd->IsNumFilled())
   {for(i=0;i<nPars; i++)lfSt[i]=atof(StepMEd->textLine(i));
   }
 else
   {QApplication::beep();
    fprintf(stderr,"Step window incomplete\n");
    return 0;
   } 

 if(LimitMEd->IsNumFilled())
   {for(i=0;i<nPars; i++)lfL[i]=atof(LimitMEd->textLine(i));
    lfL[nPars]=lfL[nPars-1];
   }
 else
   {QApplication::beep();
    fprintf(stderr,"Limit window incomplete\n");
    return 0;
   } 
  nMaxIt=atoi(MaxItEd->text());
  if(nMaxIt<=0)
    {fprintf(stderr,"Max # of iterations must be > than 0!\n");
     QApplication::beep();
     return 0;
    }

 return 1;
}
// *************************************************************
void SFitWid::FunctionSel(int in)
{//fprintf(stderr,"Function %d of %d selected\n",in,nFns);
 if(in>=0 && in<nFns)
   {iSelFn=in;
    if(FF[iSelFn]->GetNSumP()==0)NoCurEd->hide();
    else NoCurEd->show();
    ShowParList();
    return;
   }
 else
  ParText->setText("No valid\nfunction\nselected !");
  QApplication::beep();
  iSelFn=-1;   
}
// ************************************************************
void SFitWid::NoCurvChg()
{ if(iSelFn<0 || iSelFn>=nFns)return;

  nCurv=atoi(NoCurEd->text());
  if(nCurv<=0)
    {ParText->setText("Invalid\n# of curves\nentered !");
     QApplication::beep();
     return;
    }
   ShowParList(); 
}    
// ************************************************************
void SFitWid::ShowParList(void)
{const int iSP=FF[iSelFn]->GetNSumP();
 const int iFP=FF[iSelFn]->GetNFixP();

 nCurv=atoi(NoCurEd->text());
 
 nPars=nCurv*iSP+iFP;
 int i,j;
 char *sp=szAnswer;
 szAnswer[0]=0;

 for(i=0; i<nCurv; i++)
    {for(j=0; j<iSP; j++)
         sp+=sprintf(sp,"%s[%d]:\n",FF[iSelFn]->GetParName(j),i+1);  
    }	 

 if(iSP==0)nCurv=nPars;
 for(i=0; i<iFP; i++)
         sp+=sprintf(sp,"%s:\n",FF[iSelFn]->GetParName(iSP+i));  

 ParText->hide();
 int iH;
 if(nPars>nMaxPars)
   {sprintf(szAnswer,"Number of\ncurves: %d\ntoo big!!\n Pars: %d of %d",
                      nCurv,nPars,nMaxPars);
    iH=nMaxPars*iChy+iChy/2;
   }
 else iH=nPars*iChy+iChy/2;
 ParTextW.setHeight(iH);

 //ParText->setGeometry(ParTextW);
 //ParMEd->setGeometry(ParMEdW);

 ParText->setFixedHeight(iH);
 if(ParMEd->numLines()>nPars)
   {
    for(i=nPars;i<ParMEd->numLines();i++)ParMEd->removeLine(i);
   }
   
   ParMEd->setFixedHeight(iH,nPars);
  StepMEd->setFixedHeight(iH,nPars);
 LimitMEd->setFixedHeight(iH,nPars);

 ParText->setText(szAnswer);
 ParText->show();
 
}
// ************************************************************
void SFitWid::InitPars(void)
{
  delete lfP;
  delete lfL;
  delete lfSt;
  
  nPars=nCurv*FF[iSelFn]->GetNSumP()+FF[iSelFn]->GetNFixP();
   lfP = new double [nPars];
   CHECK_POINTER_RET(lfP)
   lfL = new double [nPars+1];
   CHECK_POINTER_RET(lfL)
  lfSt = new double [nPars];
   CHECK_POINTER_RET(lfSt)
  
}
// ************************************************************
void SFitWid::paintEvent(QPaintEvent *)
{
 QPainter p;
 p.begin(FitData);
 //if(WidgetAusw::ScreenHeight()==600)PF.setPointSize(11);
 if(SCREEN.y()==600)PF.setPointSize(11);

  static char szX[80],szY[80],szYE[80];
 strncpy(szX,DF->GetColID(DF->GetCRData()->GetColX()),79);
 strncpy(szY,DF->GetColID(DF->GetCRData()->GetColY()),79);
    PLP.PFont=PF;  
    PLP.PRect=QRect(FitW);
    PLP.xLab=szX;
    PLP.yLab=szY;
    PLP.pFlag=0;
    //PLP.iErrBar=0;
    PLP.iStart=0;
    PLP.iEnd=-1;
    //PLP.iSize=2;
    PLP.pFlag &= ~ZERO_X_LINE;
    PLP.pSymbol[1]=CONNECT;
 p.eraseRect(PLP.PRect.x()+5,PLP.PRect.y()+5,
             PLP.PRect.width()-5,PLP.PRect.height()-5);    
 DrawData(p,NULL,NULL,PLP);

 if(ED)
   {sprintf(szYE,"r. Error");
    PLE.yLab=szYE;
    PLE.PRect=QRect(ErrW);
    PLE.xLab=szX;
    PLE.pFlag=0;
    //PLE.iErrBar=0;
    PLE.iStart=0;
    PLE.iEnd=-1;
    //PLE.iSize=2;
    PLE.pFlag |= ZERO_X_LINE;
    p.eraseRect(PLE.PRect.x()+5,PLE.PRect.y()+5,
                PLE.PRect.width()-5,PLE.PRect.height()-5);    
    DrawData(p,EF,NULL,PLE);

    p.eraseRect(5,ErrW.y()+ErrW.height()+5,
                ResultW.width()-5,
                ResultW.height()-ErrW.y()-ErrW.height()-5);    
    DrawFormText(p,PF,5*iChx,ErrW.y()+ErrW.height(),
	                 AlignTop|AlignLeft,szOut);

   }
 p.end();

}
// *************************************************
void SFitWid::PrintIt()
{
  QPrinter Prt;
  if(Prt.setup(FitData))
   {QPainter paint;
    paint.begin(&Prt);
    PLP.PRect=QRect(FitW);
    DrawData(paint,NULL,NULL,PLP);
    if(ED)
      {PLE.PRect=QRect(ErrW);
       DrawData(paint,EF,NULL,PLE);
       DrawFormText(paint,PF,5*iChx,ErrW.y()+ErrW.height(),
       	                 AlignTop|AlignLeft,szOut);

      }

    paint.end();
   }
  //else fprintf(stderr,"Printer setup failed\n");

}
// *****************************************************
void SFitWid::keyPressEvent(QKeyEvent *e)
{static int iP=0; 
 switch(e->key()) 
       {case Key_1:iP=0;break;
        case Key_2:iP=1;break;
        case Key_3:iP=2;break;
 	case Key_C:if(iP==0 || iP==1)
	             {if(PLP.pSymbol[iP]!=CONNECT)PLP.pSymbol[iP]=CONNECT;
	              else PLP.pSymbol[iP]=iP;
		     } 
		   repaint();
		   break;
  	case Key_E:if(iP==0 || iP==1)
	             {if(PLP.pSymbol[iP]!=ERROR_BAR)PLP.pSymbol[iP]=ERROR_BAR;
	              else PLP.pSymbol[iP]=iP;
		     } 
		   repaint();
		   break;
     case Key_Plus:if(iP==0 || iP==1)
	             {if(PLP.pSymbol[iP]>=0)PLP.pSize[iP]++;}
		    else  
	             {if(PLE.pSymbol[0]>=0)PLE.pSize[0]++;}
		      //fprintf(stderr,"Sybolsize[%d] + : %d\n",iP,P.pSize[iP]);
		      repaint();
		      break;   
    case Key_Minus:if(iP==0 || iP==1)
	             {if(PLP.pSymbol[iP]>=0)PLP.pSize[iP]--;}
                   else  
		     {if(PLE.pSymbol[0]>=0)PLE.pSize[0]--;}
                      //fprintf(stderr,"Sybolsize[%d] - : %d\n",iP,P.pSize[iP]);
		      repaint();
		      break;   
    default:e->ignore();
            repaint();
            //fprintf(stderr,"A key pressed in SFitWid\n");
      }
}
// *******************************************************
MLineEd::MLineEd(QWidget *parent=0, const char *name=0):
        QMultiLineEdit(parent,name)
{connect( this, SIGNAL(returnPressed()), SLOT(CheckLength()) );

 iHeight=0;}
// *************************************************
void MLineEd::CheckLength()
{//fprintf(stderr,"Line %d of %d\n",numLines(),iHeight);
 if(numLines()>iHeight)
    {QApplication::beep();
     //int l,c;
     //cursorPosition(&l,&c);
     //if(l>=iHeight)removeLine(l);
     //else removeLine(l-1);
     backspace();
    }
}
// *************************************************
int MLineEd::IsNumFilled(void)
{//fprintf(stderr,"Analyzing line %d of %d\n",numLines(),iHeight);
 if(numLines()<iHeight)return 0;
 const char *cp;
 int i;
 for(i=0;i<numLines();i++)
   {cp=textLine(i);
    //fprintf(stderr,"Line: %d:%s\n",i,cp);
    while(*cp==' ')cp++;
    if(! *cp)
      {PRINT_DEBUG("Line: %d:*%s* -> empty\n",i,cp)
       return 0;
      } 
    while(*cp)
         {if(!strchr("0123456789-+eE. ",*cp))
	    {PRINT_DEBUG("Line: %d:*%c* (%x) -> non num\n",i,*cp,*cp)
	     return 0;
	    } 
          cp++;
	 }
	  
   }
 return 1;
}   
// *************************************************	
// FunWid Functions
// *********************************************************
 int FunWid::w=0;
 int FunWid::h=0;
        
 QRect FunWid::MainW;
 QRect FunWid::FitFnsW;   //ListBox of Fit Fns
 QRect FunWid::ParTextW;  //Label fo Parameter names
 QRect FunWid::ClBW;      //Close button
 QRect FunWid::CalcBW;    //FitButton
 QRect FunWid::PrBW;      //PrintButton
 QRect FunWid::ResultW;   //Frame displaying results
 QRect FunWid::NoCurW;    //LineEdit of # of Curves
 QRect FunWid::ParMEdW;   //MultiLineEdit of Pars 
 QRect FunWid::xStEdW;    //LineEdit of xStart 
 QRect FunWid::xEnEdW;    //LineEdit of xEnd 
 QRect FunWid::nStEdW;    //LineEdit of no. of steps

   int FunWid::iChx=0;
   int FunWid::iChy=0;
   int FunWid::nMaxPars=0;

// ********************************************************* 
FunWid::FunWid(QWidget *parent, const char *name): 
                  QWidget(parent,name)
{//int sh=WidgetAusw::ScreenHeight();
//  w = (SCREEN.y()==600 ? 750 : 900);
//  h = (SCREEN.y()==600 ? 550 : 700);

 w = 80*SCREEN.x()/100;
 h = 85*SCREEN.y()/100;
  
 const int LeftO0=10;
 BF = QFont("Helvetica",(SCREEN.y()==600 ? 12 : 14),QFont::DemiBold); // for buttons
 WF = QFont("Helvetica",(SCREEN.y()==600 ? 12 : 14),QFont::Normal); // for sub window text
 PF = QFont("Helvetica",10);

 QFontMetrics FM(WF);
 iChx = FM.width('0');
 iChy = FM.lineSpacing();

 nFns = MAX_FITFN;
 iSelFn = -1;
 OutData=0;

 MainW.setRect(10, 10,w,h);

 const int nF = (SCREEN.y()==600 ? 3 : 5);
 const QRect FitFLW(LeftO0,10,FitFnsW.width(),B.height());
 FitFnsW.setRect(LeftO0, 30,20*iChx,nF*3*iChy/2);
 NoCurW.setRect(LeftO0, FitFnsW.y()+FitFnsW.height()+2*B.height()+10,
                FitFnsW.width(),B.height());

 //const int LeftO1 = LeftO0+FitFnsW.width()+20;

  int iMH = MainW.height()-(NoCurW.y()+NoCurW.height()+20+10);
 nMaxPars = iMH/iChy-1;
 
 ParTextW.setRect(LeftO0,NoCurW.y()+NoCurW.height()+B.height()+15,
                  iChx*15,nMaxPars*iChy);
 ParMEdW.setRect(ParTextW.x()+ParTextW.width()+5,ParTextW.y(),
                 ParTextW.width()+3*iChx,ParTextW.height());
		  

 ClBW.setRect(ParMEdW.x()+ParMEdW.width()+20,ParMEdW.y()+ParMEdW.height()-B.height(),
              B.width(),B.height());

 ResultW.setRect(ParMEdW.x()+ParMEdW.width()+10, FitFLW.y(),
                 MainW.width()-ParMEdW.x()-ParMEdW.width()-20,
		 ParMEdW.y()+2*ClBW.height()-10);

 const QRect xStLW(ResultW.x(),ResultW.y()+ResultW.height()+B.height(),B.width(),B.height());
 xStEdW.setRect(xStLW.x()+xStLW.width()+10,xStLW.y(),B.width(),B.height());

 const QRect xEnLW(xStLW.x(),xStLW.y()+xStLW.height()+10,B.width(),B.height());
 xEnEdW.setRect(xStEdW.x(),xEnLW.y(),B.width(),B.height());
 
 const QRect nStLW(xEnLW.x(),xEnLW.y()+xEnLW.height()+20,B.width(),B.height());
 nStEdW.setRect(xStEdW.x(),nStLW.y(),B.width(),B.height());

 CalcBW.setRect(xEnEdW.x()+xEnEdW.width()+20,xEnEdW.y(),B.width(),B.height());

  setGeometry(MainW);
  setFixedSize(w,h);
  setCaption("Calculate a function");
  // ------------------------------------  
  //Label Function List Box
  QLabel *FnLbL = new QLabel("Fit Funs:", this);
  FnLbL->setAlignment(AlignLeft);
  FnLbL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  FnLbL->setGeometry(FitFLW);

  //Listbox for fit functions
  FitFns = new QListBox( this , "FnLB");
  FitFns->setGeometry(FitFnsW);
  FitFns->setFont(WF);
  int i;
  //for(i=0; i<nFns; i++)FitFns->insertItem(FIT_FUNCS[i].FnName);
  for(i=0; i<nFns; i++)FitFns->insertItem(FF[i]->GetFnName());
  connect( FitFns, SIGNAL(selected(int)),SLOT(FunctionSel(int)) );
  QToolTip::add( FitFns, "List of available functions" );
  // ------------------------------------  
  //Label: Number of curves single line edit window
  QLabel *CurvL=new QLabel("Number of curves:", this);
  CurvL->setAlignment(AlignLeft);
  CurvL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  CurvL->setGeometry(LeftO0,FitFnsW.y()+FitFnsW.height()+20,
                     NoCurW.width(),B.height());

  //Single Line edit: No of curves
  NoCurEd = new QLineEdit( this );
  NoCurEd->setGeometry(NoCurW);
  NoCurEd->setText("1");
  connect( NoCurEd, SIGNAL(returnPressed()),SLOT(NoCurvChg()) );
  QToolTip::add( NoCurEd, "Enter # of curves" );
  nCurv=1;
  // ------------------------------------  

  //Label: xStart single line edit window
  QLabel *xStL = new QLabel("x Start:", this);
  xStL->setAlignment(AlignRight);
  xStL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  xStL->setGeometry(xStLW);

  //Single Line edit: xStart
   xStartEd = new QLineEdit( this );
   xStartEd->setGeometry(xStEdW);
   QToolTip::add( xStartEd, "Enter x-start" );
  // ---------------------------------------

  //Label: xEnd single line edit window
  QLabel *xEnL = new QLabel("x End:", this);
  xEnL->setAlignment(AlignRight);
  xEnL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  xEnL->setGeometry(xEnLW);

  //Single Line edit: xEnd
   xEndEd = new QLineEdit( this );
   xEndEd->setGeometry(xEnEdW);
   QToolTip::add( xEndEd, "Enter x-end" );
  // ---------------------------------------

  //Label: nSteps single line edit window
  QLabel *nStL = new QLabel("# of Steps:", this);
  nStL->setAlignment(AlignRight);
  nStL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  nStL->setGeometry(nStLW);

  //Single Line edit: nSteps
   nStepsEd = new QLineEdit( this );
   nStepsEd->setGeometry(nStEdW);
   QToolTip::add( nStepsEd, "Enter # of steps" );


   nPars = 0;
    lfP = NULL;
 // ------------------------------------  
  // Names of paramters
  //Label: For parameter names 
  QLabel *PNaL = new QLabel("Parameter name:", this);
  PNaL->setAlignment(AlignLeft);
  PNaL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  PNaL->setGeometry(ParTextW.x(),NoCurW.y()+NoCurW.height()+10,
                     ParTextW.width(),B.height());

  ParText = new QLabel(this);
  ParText->setAlignment(AlignRight);
  ParText->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  ParText->setGeometry(ParTextW);
  ParText->setFont(WF);
  ParText->setText("Names of\nparamters\nof selected\nfunction");
  QToolTip::add( ParText, "List of pars for selected function" );

  // ------------------------------------  
  // Values of paramters
  //Label: For parameter values 
  QLabel *PVaL = new QLabel("Initial value:", this);
  PVaL->setAlignment(AlignLeft);
  PVaL->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  PVaL->setGeometry(ParMEdW.x(),NoCurW.y()+NoCurW.height()+10,
                     ParTextW.width(),B.height());

  ParMEd = new MLineEd(this);
  ParMEd->setGeometry(ParMEdW);
  ParMEd->setFont(WF);
  ParMEd->setCellHeight(iChy);
  connect( ParMEd, SIGNAL(returnPressed()),SLOT(ChkLength()) );
  QToolTip::add( ParMEd, "Enter initial parameter values for fit" );


  // ------------------------------------  
  // Close Button
  ClB = new QPushButton("Close", this, "close" );
  ClB->setFont(BF);
  ClB->setFixedHeight(ClB->sizeHint().height() );
  ClB->setGeometry(ClBW);
  connect( ClB, SIGNAL(clicked()),SLOT(Close()) );
  QToolTip::add( ClB, "Close (exit) this window" );
  // ------------------------------------  


  // Calc Button
  CalcB = new QPushButton("Calculate", this, "fit" );
  CalcB->setFont(BF);
  CalcB->setFixedHeight(CalcB->sizeHint().height() );
  CalcB->setGeometry(CalcBW);
  connect( CalcB, SIGNAL(clicked()), SLOT(Calculate()) );
  QToolTip::add( CalcB, "Start fit with initial pars" );
  //CalcB->hide();
  // ------------------------------------  

  // Print Button
  PrB = new QPushButton("Print", this, "print" );
  PrB->setFont(BF);
  PrB->setFixedHeight(PrB->sizeHint().height() );
  PrB->setGeometry(PrBW);
  connect( PrB, SIGNAL(clicked()), SLOT(PrintIt()) );
  QToolTip::add( PrB, "Print grafic" );
  // ------------------------------------  

  FitData = new QFrame(this,"FitFrame");
  FitData->setGeometry(ResultW);
  FitData->setFrameStyle( QFrame::Panel | QFrame::Sunken);
  FitData->setBackgroundColor(white);
  FitData->setFocusPolicy(QWidget::ClickFocus);
  show();

}
// ************************************************************
void FunWid::Close()
{fprintf(stderr,"Close pressed\n");
 qApp->quit();
}
// ************************************************************
void FunWid::Calculate()
{if(iSelFn<0)
    {PRINT_DEBUG("Illegal function %d\n",iSelFn)
     QApplication::beep();
     return;
    } 
 InitPars();
 if(!CollectPars())return;
  double xStart=atof(xStartEd->text());
  double xEnd=atof(xEndEd->text()); 
  int nSteps=strtol(nStepsEd->text(),(char **)NULL, 10);

  if(xStart>=xEnd || nSteps <=1)
    {PRINT_DEBUG("Incompatible steps [%f, %f] %d\n",xStart,xEnd,nSteps);
     QApplication::beep();
     return;
    } 
  delete OutData;
  OutData=new CRData(2,nSteps);
  //fprintf(stderr,"Allocating new data: xStart:%f xEnd:%f Steps:%i\n",xStart,xEnd,nSteps);
  CHECK_POINTER_RET(OutData)
  OutData->SetColX(0);
  OutData->SetColY(1);      
   int i;
   double lfM=(xEnd-xStart)/(nSteps-1);
   for(i=0;i<OutData->GetSteps(); i++)
      {OutData->SetPointX(xStart+i*lfM,i);
       OutData->SetPointY(FF[iSelFn]->GetFFStruct().FitFunc(OutData->GetPointX(i),nCurv,lfP),i); 
      }
 
     repaint();
}
// ************************************************************
void FunWid::ChkLength()
{
 if(ParMEd->IsNumFilled())
   {InitPars();
    int i;
    for(i=0;i<nPars; i++)
       {lfP[i]=atof(ParMEd->textLine(i));}
   }
   
}
// ************************************************************
int FunWid::CollectPars(void)
{int i;
 if(ParMEd->IsNumFilled())
   {for(i=0;i<nPars; i++)lfP[i]=atof(ParMEd->textLine(i));
    if(strcmp(FF[iSelFn]->GetFnName(),"Lorentz Alfa 12")==0)
        {
         RWLa21=0.5;
        }
    return 1;
   }
 else
   {QApplication::beep();
    fprintf(stderr,"Parameter window incomplete\n");
    return 0;
   } 
}
// ************************************************************
void FunWid::FunctionSel(int in)
{//fprintf(stderr,"Function %d of %d selected\n",in,nFns);
 if(in>=0 && in<nFns)
   {iSelFn=in;
    if(FF[iSelFn]->GetNSumP()==0)NoCurEd->hide();
    else NoCurEd->show();
    ShowParList();
    return;
   }
 else
  ParText->setText("No valid\nfunction\nselected !");
  QApplication::beep();
  iSelFn=-1;   
}
// ************************************************************
void FunWid::NoCurvChg()
{ if(iSelFn<0 || iSelFn>=nFns)return;

  nCurv=atoi(NoCurEd->text());
  if(nCurv<=0)
    {ParText->setText("Invalid\n# of curves\nentered !");
     QApplication::beep();
     return;
    }
   ShowParList(); 
}    
// ************************************************************
void FunWid::ShowParList(void)
{const int iSP=FF[iSelFn]->GetNSumP();
 const int iFP=FF[iSelFn]->GetNFixP();

 nCurv=atoi(NoCurEd->text());
 
 nPars=nCurv*iSP+iFP;
 int i,j;
 char *sp=szAnswer;
 szAnswer[0]=0;

 for(i=0; i<nCurv; i++)
    {for(j=0; j<iSP; j++)
         sp+=sprintf(sp,"%s[%d]:\n",FF[iSelFn]->GetParName(j),i+1);  
    }	 

 if(iSP==0)nCurv=nPars;
 for(i=0; i<iFP; i++)
         sp+=sprintf(sp,"%s:\n",FF[iSelFn]->GetParName(iSP+i));  

 ParText->hide();
 int iH;
 if(nPars>nMaxPars)
   {sprintf(szAnswer,"Number of\ncurves: %d\ntoo big!!\n Pars: %d of %d",
                      nCurv,nPars,nMaxPars);
    iH=nMaxPars*iChy+iChy/2;
   }
 else iH=nPars*iChy+iChy/2;
 ParTextW.setHeight(iH);

 //ParText->setGeometry(ParTextW);
 //ParMEd->setGeometry(ParMEdW);

 ParText->setFixedHeight(iH);
 if(ParMEd->numLines()>nPars)
   {
    for(i=nPars;i<ParMEd->numLines();i++)ParMEd->removeLine(i);
   }
   
   ParMEd->setFixedHeight(iH,nPars);
 ParText->setText(szAnswer);
 ParText->show();
 
}
// ************************************************************
void FunWid::InitPars(void)
{
  delete lfP;
  
  nPars=nCurv*FF[iSelFn]->GetNSumP()+FF[iSelFn]->GetNFixP();
   lfP = new double [nPars];
   CHECK_POINTER_RET(lfP)
}
// ************************************************************
void FunWid::paintEvent(QPaintEvent *)
{
 QPainter p;
 p.begin(FitData);
 //if(WidgetAusw::ScreenHeight()==600)PF.setPointSize(11);
 if(SCREEN.y()==600)PF.setPointSize(11);

 char *szX="x", *szY="y";
    PLP.PFont=PF;  
    PLP.PRect=QRect(0,0,ResultW.width(),ResultW.height());
    PLP.xLab=szX;
    PLP.yLab=szY;
    PLP.pFlag=0;
    //PLP.iErrBar=0;
    PLP.iStart=0;
    PLP.iEnd=-1;
    //PLP.iSize=2;
    //PLP.pFlag &= ~ZERO_X_LINE;
    PLP.pSymbol[0]=CONNECT;
 p.eraseRect(PLP.PRect.x()+5,PLP.PRect.y()+5,
             PLP.PRect.width()-5,PLP.PRect.height()-5);    

 DrawData(p,NULL,OutData,PLP);

  p.end();

}
// *************************************************
void FunWid::PrintIt()
{
  QPrinter Prt;
  if(Prt.setup(FitData))
   {QPainter paint;
    paint.begin(&Prt);
    PLP.PRect=QRect(ResultW);
    DrawData(paint,NULL,OutData,PLP);
    paint.end();
   }
  //else fprintf(stderr,"Printer setup failed\n");

}
// *****************************************************
void FunWid::keyPressEvent(QKeyEvent *e)
{const int iP=0; 
 switch(e->key()) 
       {case Key_C:if(PLP.pSymbol[iP]!=CONNECT)PLP.pSymbol[iP]=CONNECT;
                   else PLP.pSymbol[iP]=iP;
		   repaint();
		   break;
  	case Key_E:if(PLP.pSymbol[iP]!=ERROR_BAR)PLP.pSymbol[iP]=ERROR_BAR;
	           else PLP.pSymbol[iP]=iP;
		   repaint();
		   break;
     case Key_Plus:if(PLP.pSymbol[iP]>=0)PLP.pSize[iP]++;
		      //fprintf(stderr,"Sybolsize[%d] + : %d\n",iP,P.pSize[iP]);
		      repaint();
		      break;   
    case Key_Minus:if(PLP.pSymbol[iP]>=0)PLP.pSize[iP]--;
                   //fprintf(stderr,"Sybolsize[%d] - : %d\n",iP,P.pSize[iP]);
	           repaint();
		   break;   
    default:e->ignore();
            repaint();
            //fprintf(stderr,"A key pressed in FunWid\n");
      }
}
// *******************************************************
