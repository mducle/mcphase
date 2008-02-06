//file: thefunc.cpp
//$Id: thefunc.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
//$Log: thefunc.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.3  1999/03/15 09:08:37  herbie
//*** empty log message ***
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <bits/nan.h>

#include "thefunc.h"
#include "dfile.h"

// ***************************************************
//XAUSW
int GetCellPars(CAP_CELL *CC, LineString &PP, LineString &TH)
{

 CC->uFlag=0;
 CC->szZero[0]=0;
 CC->szCal[0]=0;
 CC->szCal2[0]=0;
 int icnt=0;


 char szCell[12];
 char cCell;
 if(!TH.GetPar("[FileHeader]","Cell","%c",&cCell))return 0;
 sprintf(szCell,"[Cell_%c]",cCell);
 
 if(PP.GetPar(szCell,"Diam","%lf",&(CC->fCDiam)))
              {CC->uFlag|=CC_CELL_DIAM; icnt++;}
 if(PP.GetPar(szCell,"Gap","%lf",&(CC->fCGap)))
              {CC->uFlag|=CC_CELL_GAP; icnt++;}
 if(PP.GetPar(szCell,"ZeroTemp","%lf",&(CC->fZeroTemp)))
              {CC->uFlag|=CC_CELL_ZTEM; icnt++;}
 if(PP.GetPar(szCell,"CalFile","%s",CC->szCal))
              {CC->uFlag|=CC_CAL1_FILE; icnt++;}
 if(PP.GetPar(szCell,"ZeroFile","%s",CC->szZero))
              {CC->uFlag|=CC_ZERO_FILE; icnt++;}

 char szB[MAX_LINELENGTH+1];
 if(!TH.GetPar("[FileHeader]","CalFile","%S",szB))return 0;
 if(strcmp(szB,"NULL"))
   {strcpy(CC->szZero,szB);
    CC->uFlag|=CC_CAL2_FILE;
    icnt++;
  }

 double k0;
 if(!TH.GetPar("[FileHeader]","K0","%lf",&k0))return 0;
 CC->K0=k0;
 CC->uFlag|=CC_CELL_C0;
 icnt++;

 if(!TH.GetPar("[FileHeader]","L0","%lf",&k0))return 0;
 CC->l0=k0*1.e-3;     // L0 in meters
 CC->uFlag|=CC_CELL_L0;
 icnt++;
 
 if(PP.GetPar(szCell,"CalFile2","%s",CC->szCal2))
   {//strcpy(CC->szCal2,"NULL");
    CC->uFlag|=CC_CAL2_FILE;
    icnt++;
   }

 if(PP.GetPar(szCell,"Pivot","%lf",&(CC->fCPivot)))
 {CC->uFlag|=CC_CELL_PIVO; icnt++;}
 return icnt;

}
// *******************************************************
int SPrintCC(CAP_CELL &CC, char *szB, const unsigned W)
{ char *sp=szB;

  if(szB)sp=szB;
  else return 0;

  int icnt=0;

  if( (CC.uFlag & CC_CELL_DIAM) && (W &  CC_CELL_DIAM) )
    {sp+=sprintf(sp, "%4.2f\n",CC.fCDiam); icnt++;}
  if( (CC.uFlag & CC_CELL_GAP) && ( W & CC_CELL_GAP)  )
    {sp+=sprintf(sp, "%4.2f\n",CC.fCGap); icnt++;}
  if( (CC.uFlag & CC_CELL_PIVO)  && ( W & CC_CELL_PIVO)  )
    {sp+=sprintf(sp, "%4.2f\n",CC.fCPivot); icnt++;}

  if((CC.uFlag & CC_ZERO_FILE)  && ( W & CC_ZERO_FILE) )
    {sp+=sprintf(sp, "%s\n",CC.szZero); icnt++;}
  if( (CC.uFlag & CC_CAL1_FILE) && ( W &  CC_CAL1_FILE) )
    {sp+=sprintf(sp, "%s\n",CC.szCal); icnt++;}
  if( (CC.uFlag & CC_CAL2_FILE)  && ( W &  CC_CAL2_FILE) )
    {sp+=sprintf(sp, "%s\n",CC.szCal2); icnt++;}
   
  if( (CC.uFlag & CC_CELL_ZTEM)  && ( W & CC_CELL_ZTEM)  )
    {sp+=sprintf(sp, "%4.2f\n",CC.fZeroTemp); icnt++;};

  if( (CC.uFlag & CC_CELL_C0)  && (W & CC_CELL_C0)  )
    {sp+=sprintf(sp, "%f\n",CC.K0); icnt++;}

  if( (CC.uFlag & CC_CELL_L0)  && (W  &  CC_CELL_L0) )
    {sp+=sprintf(sp, "%f4.2\n",CC.l0); icnt++;}
   
  return icnt;
}
// ******************************************************
//XAUSW
int CalcGap(DataFile *DF, LineString &FP, char * szA)
{char *sp=szA;

 if(DF->GetFType()/10*10 != (int)TheGROUP)
    {if(szA)sp+=sprintf(sp,"CalcGap: No thermal expansion data file");
     else fprintf(stderr,"CalcGap: No thermal expansion data file(%d)!=%d\n",DF->GetFType()*10/10,(int)TheGROUP);
     return 0;
    }

 CRData *TD=DF->GetCRData();

 if(TD==NULL || TD->GetSteps()==0)
   {if(szA)sp+=sprintf(sp,"Block empty");
    else fprintf(stderr,"Block empty\n");
    return 0;
   }

 CAP_CELL CC;

 
 if(!GetCellPars(&CC,FP,DF->GetHead()))
   {if(szA)sp+=sprintf(sp,"Parameter not found");
    else fprintf(stderr,"Parameter not found\n");
    return 0;
   }
 
 if(! (CC.uFlag & CC_CELL_DIAM) )return 0;

 if(   (CC.uFlag & CC_CAL2_FILE) && (CC.uFlag & CC_CAL1_FILE)
    && (CC.uFlag & CC_CELL_GAP)  && (CC.uFlag & CC_CELL_PIVO))
  {if(szA)sp+=sprintf(sp,"Ag file:%s\nSaphire file:%s\n",CC.szCal,CC.szCal2);
   else fprintf(stderr,"Ag file:%s\nSaphire file:%s\n",CC.szCal,CC.szCal2);

 TheFile *AgF=new TheFile(CC.szCal);

 if( AgF->ReadData(1)>0)
   {if(szA)sp+=sprintf(sp,"Error reading data of %s\n",CC.szCal);
    fprintf(stderr,"Error reading data of %s\n",CC.szCal);
    return 0;
   }

  CRData *Ag=AgF->GetCRData();
  if(AgF->GetError() || !Ag || Ag->GetSteps()==0)
    {if(szA)sp+=sprintf(sp,"%s not found\n",CC.szCal);
     else fprintf(stderr,"%s not found\n",CC.szCal);
     return 0;
    }

 Ag->SetColX(0);
 Ag->SetColY(1);
   
 TheFile *SapF=new TheFile(CC.szCal2);

 if( SapF->ReadData(1)>0)
   {if(szA)sp+=sprintf(sp,"Error reading data of %s\n",CC.szCal2);
    fprintf(stderr,"Error reading data of %s\n",CC.szCal2);
    return 0;
   }

  CRData *Sap=SapF->GetCRData();
  if(SapF->GetError() || !Sap || Sap->GetSteps()==0)
    {if(szA)sp+=sprintf(sp,"%s not found\n",CC.szCal2);
     else fprintf(stderr,"%s not found\n",CC.szCal2);
     return 0;
    }

   Sap->SetColX(0);
   Sap->SetColY(1);

    double ri=0.5*CC.fCGap;
    double ro=0.5*CC.fCDiam;

    TiltPlate *TP=new TiltPlate(ro,ri,CC.fCPivot,0.8, CC.K0,Ag,Sap);
    if(!TP->IsValid())
      {if(szA)sp+=sprintf(sp,"Invalid parameter\n");
       else fprintf(stderr,"Invalid parameter\n");
       return 0;
      }

    TP->CalcGap(TD);
    if(TP)delete TP;
    if(AgF)delete AgF;
    if(SapF)delete SapF;
    return 1;
  } 
 //else
 if(szA)
   sp+=sprintf(sp,"Calculate Gap");
 if( ! (CC.uFlag & CC_CELL_GAP) )CC.fCGap=-10.;
 return (*TD)[TD->GetColY()].CalculateGap(CC.fCDiam*0.1,CC.fCGap*0.1);
}
// ********************************************************************
//XAUSW
int Calibrate(DataFile *DF, LineString &PP, char * szA)
{
 char *sp=szA;

 CRData *TD=DF->GetCRData();

 if(!TD || TD->GetSteps()<2)
   {if(szA)sp+=sprintf(sp,"Block empty\n");
    else fprintf(stderr,"Block empty\n");
    return 0;
   }
 TD->SetColX(0);
 TD->SetColY(1);
   
 CAP_CELL CC;
 if(!GetCellPars(&CC,PP,DF->GetHead()))
   {if(szA)sp+=sprintf(sp,"Calibrate: Parameter not found\n");
    else fprintf(stderr,"Calibrate: Parameter not found\n");
    return 0;}

 if(! (CC.uFlag & CC_ZERO_FILE) && (CC.uFlag & CC_CAL1_FILE) )
   {if(szA)sp+=sprintf(sp,"Calibrate: Insufficient cell Parameters\n");
    else fprintf(stderr,"Calibrate: Insufficient cell Parameters\n");
    return 0;
   }
 if(szA)
   {sp+=sprintf(sp,"Calibration of THE-Data\n");
    sp+=sprintf(sp,"%s : Zero signal correction\n",CC.szZero);
    sp+=sprintf(sp,"%s : Reference correction\n",CC.szCal);
   }
else 
   {fprintf(stderr,"Calibration of THE-Data\n");
    fprintf(stderr,"%s : Zero signal correction\n",CC.szZero);   
    fprintf(stderr,"%s : Reference correction\n",CC.szCal);
   }
  if( (CC.uFlag&CC_CAL2_FILE) && (CC.uFlag&CC_CELL_C0) )
    {if(szA)
       {sp+=sprintf(sp,"%s : Reference correction\n",CC.szCal2);
        sp+=sprintf(sp,"C0 : %f\nL0 = %f\n",CC.K0,CC.l0);
       }
     else 
       {fprintf(stderr,"%s : Reference correction\n",CC.szCal2);
        fprintf(stderr,"C0 : %f\nL0 = %f\n",CC.K0,CC.l0);
       }
    }   
 //Calculate Gap
  if(szA)sp+=sprintf(sp,"1: Calculating Gap\n");
  else fprintf(stderr,"1: Calculating Gap\n");

  if(CalcGap(DF,PP,szA)==0)
    {if(szA)sp+=sprintf(sp,"Error calculating gap\n");
     else fprintf(stderr,"Error calculating gap\n");
     return -3;
    }

  if(szA)sp+=sprintf(sp,"2: Calculating Dl/l\n");
  else fprintf(stderr,"2: Calculating Dl/l\n");
  
  if( (*TD)[TD->GetColY()].CalculateDl_l(CC.l0)==0)
    {if(szA)sp+=sprintf(sp,"Error calculating Dl/l");
     else fprintf(stderr,"Error calculating Dl/l");
     return -3;
    }

 TD->ZeroShift(CC.fZeroTemp);

 //Zerosignal
 if(szA)sp+=sprintf(sp,"3:Subtract Zero Signal: %s\n",CC.szZero);
 else fprintf(stderr,"3:Subtract Zero Signal: %s\n",CC.szZero);
 
 TheFile *HF=new TheFile(CC.szZero);
 if( HF->ReadData(1)>0)
   {if(szA)sp+=sprintf(sp,"Error reading data of %s\n",CC.szZero);
    fprintf(stderr,"Error reading data of %s\n",CC.szZero);
    return -3;
   }
 CRData *HD=HF->GetCRData();

 if(HF->GetError())
   {if(szA)sp+=sprintf(sp,"Bad zero file %s\n",CC.szZero);
    else fprintf(stderr,"Bad zero file %s\n",CC.szZero);
    return -3;
   } //Illegal DataFile ?

 if(HD==NULL || HD->GetSteps()==0)
   {if(szA)sp+=sprintf(sp,"Out of memory allocating zero signal\n");
    else fprintf(stderr,"Out of memory allocating zero signal\n");
    delete HF;
    return -2;
    } //No memory

 HD->SetColX(0);
 HD->SetColY(1);

 //Calculate Gap

 if(CalcGap(HF,PP,sp)==0)
   {if(szA)sp+=sprintf(sp,"Error calculating gap from zero file\n");
    else fprintf(stderr,"Error calculating gap from zero file\n");
    delete HF;
    return -3;
   }

 //Calculate Dl/l
 if( (*HD)[HD->GetColY()].CalculateDl_l(CC.l0)==0)
   {if(szA)sp+=sprintf(sp,"Error calculating Dl/l from zero file\n");
    else fprintf(stderr,"Error calculating Dl/l from zero file\n");
    delete HF;
    return -3;
   }

 HD->ZeroShift(CC.fZeroTemp);

 TD->MathOper(HD,'-');
 delete HF;

 // Cu Calibration
 if(szA)sp+=sprintf(sp,"4:Add %s\n",CC.szCal);
 else fprintf(stderr,"4:Add %s\n",CC.szCal);

 HF=new TheFile(CC.szCal);
 if( HF->ReadData(1)>0)
   {if(szA)sp+=sprintf(sp,"Error reading data of %s\n",CC.szCal);
    fprintf(stderr,"Error reading data of %s\n",CC.szCal);
    return -3;
   }
 HD=HF->GetCRData();

 if(HF->GetError())
   {if(szA)sp+=sprintf(sp,"Bad reference file %s\n",CC.szCal);
    else fprintf(stderr,"Bad reference file %s\n",CC.szCal);
    return -3;
   } //Illegal DataFile ?

 if(HD==NULL || HD->GetSteps()==0)
   {if(szA)sp+=sprintf(sp,"Out of memory allocating reference data\n");
    else fprintf(stderr,"Out of memory allocating reference data\n");
    delete HF;
    return -2;
   } //No memory
 HD->SetColX(0);
 HD->SetColY(1);

 HD->ZeroShift(CC.fZeroTemp);

 TD->MathOper(HD,'+');

// 2. Spalte = Index 1 
 DF->SetColID("Dl/l0 [m/m]",1);

 if(szA)sp+=sprintf(sp,"Calibration complete.\n");
 else fprintf(stderr,"Calibration complete.\n");

 delete HF;

 return 1;
}
// ******************************************************
// int u2th(CRData *CR, const char * szF, const double lfH)
// {
//   if(CR->GetSteps<=0)
//    {PRINT_DEBUG("Empt data block\n");
//     return 0;
//    }
// 
//   FieldSensor F;
//   if(!F.ReadSensor(szF))
//    {PRINT_DEBUG("Bad sensor file %s\n",szF);
//     return 0;
//    }
//     
//  int i;
// // double lfV;
//  for(i=0;i<CR->GetSteps();i++)
//     {(*CR)[0][i]=F.GetT((*CR)[2][i],lfH,FROM_U);
//     }
//  return i-1;
// }
// ******************************************************
// ******************************************************
// TiltPlate Functions
// ******************************************************
TiltPlate::TiltPlate(const double Nro, const double Nri, const double Nb,
                     const double Nd, const double nk,CRData *Ag, CRData *Saph)
{ro=Nro*1.E-3; ri=Nri*1.E-3; b=Nb*1.E-3; ds=Nd*1.E-3;
 //Formula requires [m] input is in mm

 k0=EPS0/nk/1.E-12*(ro*ro*M_PI - ri*ri*M_PI);
 //Formula requires k0 in [m] input is in pF (from measurement)
 T=C=0;

 DLAg=DLSaph=0;
 if(Ag->GetSteps()>2)DLAg=Ag;
 if(Saph) {if(Saph->GetSteps()>2)DLSaph=Saph;}
 iErr=(DLAg==0 || (Saph ? DLSaph==0 : 0) );
 Ag->ZeroShift(300);
 if(Saph)Saph->ZeroShift(300);
 }
// *******************************************************
double TiltPlate::DistFunc(const double x)
{MDATA dllAg;

 if( (iErr= (DLAg->LinIntpol(T,dllAg)<=0) ))return 0;

 MDATA dllSp;
 if( (iErr= (DLSaph->LinIntpol(T,dllSp)<=0) ))return 0;

 double kT=k0+2.*ds*(dllAg-dllSp);

 double go=ro/b*(kT/x-1);
 double gi=ri/b*(kT/x-1);

 double AExp=(1+dllAg)*(1+dllAg);
 double Ai=ri*ri*M_PI*AExp;
 double Ao=ro*ro*M_PI*AExp;

 double to=Ao*(1-sqrt(1-go*go))/go/go;
 double ti=Ai*(1-sqrt(1-gi*gi))/gi/gi;

 return C - (2*EPS0/x *( to-ti));
}
// ******************************************************************
double TiltPlate::CompDistFunc(const double x)
{MDATA dllAg;

 if( (iErr= (DLAg->LinIntpol(T,dllAg)<=0) ))return 0;

//Sapphire not relevant for compensated cell
 double kT=k0*(1+dllAg);

 double go=ro/b*(kT/x-1);
 double gi=ri/b*(kT/x-1);

// double AExp=(1+dllAg)*(1+dllAg);
double AExp=(1+dllAg); // (1+dllag)nicht quadrieren, da capacitaetssensor
                       // einerseits durch vergroessern der platten mit
                       // steigender T um (1+dllag)^2 zu grosse cap werte
                       // liefert, andererseits aber wegen der thermischen
                       // verkuerzung des plattenabstands um (1-dllag) zu
                       // kleine cap werte liefert - beide effekte zusammen
                       // sollten der korrekturfaktor (1+dllag) am
                       // besten beschreiben
 double Ai=ri*ri*M_PI*AExp;
 double Ao=ro*ro*M_PI*AExp;

 double to=Ao*(1-sqrt(1-go*go))/go/go;
 double ti=Ai*(1-sqrt(1-gi*gi))/gi/gi;

 return C - (2*EPS0/x *( to-ti));
}
// ******************************************************************
void TiltPlate::CalcGap(CRData *Meas,int iFuncType)
{ int i,iS=Meas->GetSteps(),iI;
  double Eps=1.e-30, Del= 1.E-30,x0,x1,x;
  //double Ai=ri*ri*M_PI;
  //double Ao=ro*ro*M_PI;

  for(i=0;i<iS;i++)
	  {x0=0.5*k0; //2*EPS0/Meas->GetPointY(i)/1.E-12*(Ao - Ai);
		x1=2*k0;
		SetT_C(Meas->GetPointX(i),Meas->GetPointY(i));
		iI=RegulaFalsi(x0,x1,Eps,Del,x,iFuncType);
		if(iI<0)
		  PRINT_DEBUG("Cannot calculate gap point %d; error in regula falsi\n",i);
		//printf("%d: x0:%lf x1:%lf x:%lf #:%d\n",i,x0,x1,x,iI);
		//x=EPS0/Meas->GetPointY(i)/1.E-12*(Ao-Ai);
		Meas->SetPointY(x,i);
	  }
	 Meas->NewMinMax();
}
// *****************************************************************
int TiltPlate::RegulaFalsi(const double x0, const double x1, const double Eps,
					 const double Delta, double &x, int iFuncType)
{double F0=SelectFunction(x0,iFuncType);
 double F1=SelectFunction(x1,iFuncType);
 if(F0*F1>0)return -1;
 if(F0==0){x=x0;return 1;}
 if(F1==0){x=x1;return 1;}

 double xm=x1,xmm1=x0,xmp1,falt=0,fneu;
 double fm,fmm1;
 while(1)
 {fm=SelectFunction(xm,iFuncType);
  fmm1=SelectFunction(xmm1,iFuncType);
  xmp1=xm - fm /( (fm-fmm1)/(xm-xmm1) );

  //printf("xm-1=%lf  xm=%lf  xm+1=%lf\n",xmm1,xm,xmp1);
  //printf("F(xm-1)=%lf  F(xm)=%lf  F(xm+1)=%lf\n",fmm1,fm,DistFunc(xmp1));
  fneu=SelectFunction(xmp1,iFuncType);
  //printf("Eps: %12.5g  Delta:%12.5g\n",fneu-falt,DistFunc(xmp1));
  if(fabs(fneu-falt)<=Eps)break;
  if(fabs(fneu)<Delta)break;
  falt=fneu;
  double Fmp1=SelectFunction(xmp1,iFuncType);
  if(sign(Fmp1)!=sign(fm)){xmm1=xmp1;continue;}
  if(sign(Fmp1)!=sign(fmm1)){xm=xmp1;continue;}
 }
  x=xmp1;
  return 1;
 }
// ****************************************************************
double TiltPlate::SelectFunction(double x, int iFuncType)
{ switch (iFuncType)
         {  case FUNC_MICRO: return DistFunc(x);
          case FUNC_COMPENS: return CompDistFunc(x);
                    default: return  NAN;
         }

}// ****************************************************************
