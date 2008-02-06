import sys,types
from qt import *

try: locale()['MOD_ftypes']
except:
 from asciifile import *
 try: from thefile import *
 except: pass
 try: from sxsfile import *
 except: pass
#print "grafic", MOD_ftypes

from xydata import *
from axis import *
from time import *
from stdfunc import *
from math import *
from string import *

#$Log: grafic.py,v $
#Revision 1.2  2006/01/04 14:41:36  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 08:57:36  herbie
#Initial revision
#
CVS_ID="$Id: grafic.py,v 1.2 2006/01/04 14:41:36 herbie Exp herbie $"

SCREEN=QPoint(0,0)
END_CLOSE = 1
END_QUIT = 2

MAX_PLOTS = 8
TOP_FILENAME = 0x0001
BOT_FILELIST = 0x0002
BOT_FILEHEAD = 0x0004
TOP_IDN = 0x0008
ZERO_X_LINE = 0x0010
TOP_INFOLINE = 0x0020

NORM_FORM =  1
ARCH_FORM =  2

SYMBOL=["CIRCLE", "QUAD", "TRIANGLE", "DIAMOND", "CROSS", "LEFT_CROSS", "V_LINE",
        "ERROR_BAR", "CONNECT", "BAR"]

COLORS = [Qt.black,   Qt.red,        Qt.darkGreen, Qt.blue, Qt.cyan,
          Qt.magenta, Qt.darkYellow, Qt.gray]

iP=0

B=QRect(10,20, 90, 20) # button prototype
# ---------------------------------------
def PrintRect(R,N=0,f=sys.stderr):
    f.write("%s: x:%d y:%d w:%d h:%d\n" % (N,R.x(),R.y(),R.width(),R.height()))
# ---------------------------------------
def IsDFileObj(obj):
    try: fn=obj.FileName
    except: return 0
    else: return 1
# ---------------------------------------
def SqDist(p,q):
    """Square of distance of two QPoints"""
    return (p.x()-q.x())**2 + (p.y()-q.y())**2
# ---------------------------------------
def MinDist(p):
    """Calculates the smallest distance of two points in
       a QPointArray
       Parameters:
         p: QPointArray
       Return:
         minimal distance of two neighbour points
       Errors:
         ValueError: if size of QPointArray is to small
    """
    l=p.size();
    if l < 2:
       raise ValueError("Size must be >= 2 (%d)" % l)
    for i in range(0,l-1):
       s=sqrt((p.point(i+1)[0]-p.point(i)[0])**2+(p.point(i+1)[1]-p.point(i)[1])**2)
       if i==0: smin=s
       else:
        if s<smin: smin=s
   
    return smin
# ---------------------------------------
def DrawFormText(p, font, x0,  y0, flag, text):
    f=QFont(p.font())
    p.setFont(font)
    FM=p.fontMetrics()
    sp=FM.lineSpacing()
    ix=FM.width("0")
    if type(text)==types.ListType: l=list(text); t=string.join(text,'\n')
    elif type(text)==types.StringType:l=text.split('\n'); t=text
    else: raise TypeError("Bad type of parameter text")
    
    icr=len(l)
    max_w=len(l[0])
    for i in l:
       ll=len(i)
       if ll>max_w: max_w=ll
#     print "%d lines; spacing: %d -> %d lines free" % (icr,iSp,0) #(600-PP.PRect.y())/iSp);		 
    p.drawText(x0,y0,ix*max_w,icr*sp,flag,t)
    p.setFont(f)
    return (ix*max_w,icr*sp)
# ---------------------------------------
def DrawSymbols(symbol, size, pa, qp, istart=0, iend=-1):
    """Draw a list of points with the same symbol
       Parameters:
          symbol: see DrawSymbol
            size: see DrawSymbol
              pa: QPointArray of points to be drawn
              qp: QPainter for painting
          istart: Index of first point to be drawn
            iend: Index of last point to be drawn
       Errors:
         ValueError: Wrong istart/iend
    """
    if iend == -1: iend=pa.size()
    if istart>=iend or iend > pa.size() or istart < 0:
       raise ValueError("Wrong istart (%d) or iend (%d)" % (istart,iend))
    for i in range(istart,iend):
        if i==istart: pi=i
        else: pi=i-1
        DrawSymbol(symbol,size,qp,pa.point(i),pa.point(pi))
# ---------------------------------------
def DrawSymbol(symbol, size, qp, p,  pp=None):
    """Draw one of the symbols with the specifid size to a QPainter
       Available sysmbols:
       "CIRCLE", "QUAD", "TRIANGLE", "DIAMOND", "CROSS", "LEFT_CROSS", "V_LINE",
       "ERROR_BAR", "CONNECT", "BAR"
       Parameters:
        symbol: text; one of SYMBOL
          size: symbol size in pixels
            qp: QPainter for painting
             p: QPoint to paint
            pp: previos QPoint only for symbols 'CONNECT' and 'BAR' to draw a line
       Errors:
        ValuError: If symbol is not available
    """
#    print pp
    if symbol == 0 or symbol == 5 or symbol == 'CIRCLE':
       qp.drawEllipse(p[0]-size,p[1]-size, 2*size,2*size)
    elif symbol == 1 or symbol == 6 or symbol == 'QUAD':
       qp.drawRect(p[0]-size,p[1]-size,2*size,2*size);
    elif symbol == 2 or symbol == 7 or symbol == 'TRIANGLE':
       qp.drawLine(p[0]-size,p[1]+size,p[0],p[1]-size)
       qp.drawLine(p[0],p[1]-size, p[0]+size,p[1]+size)
       qp.drawLine(p[0]+size,p[1]+size,p[0]-size,p[1]+size)
    elif symbol == 3 or symbol == 'DIAMOND':
       qp.drawLine(p[0],p[1]+size,p[0]-size,p[1])
       qp.drawLine(p[0]-size,p[1],p[0],p[1]-size)
       qp.drawLine(p[0],p[1]-size,p[0]+size,p[1])
       qp.drawLine(p[0]+size,p[1],p[0],p[1]+size)
    elif symbol == 4 or symbol == 'CROSS':
       qp.drawLine(p[0],p[1]-size,p[0],p[1]+size)
       qp.drawLine(p[0]-size,p[1],p[0]+size,p[1])
    elif symbol == 8 or symbol == 'LEFT_CROSS':
       qp.drawLine(p[0],p[1]-size,p[0],p[1]+size)
       qp.drawLine(p[0]-size,p[1],p[0],p[1])
    elif symbol == 9 or symbol == 'V_LINE':
       qp.drawLine(p[0],p[1]-size,p[0],p[1]+size)
    elif symbol == 10 or symbol == 'CONNECT':
       if  not pp: return
#       if pp[0] == -1 and pp[1]==-1: return
       qp.drawLine(p[0],p[1],pp[0],pp[1])
    elif symbol == 11 or symbol == 'BAR':
       qp.drawLine(p[0],p[1],pp[0],pp[1])
    else: raise ValueError("No symbol '%s' avalilable",str(symbol))
# ---------------------------------------
def DrawData(p, plot_obj, ppars):
    """Plot a object with automatic scaling
       Parameters:
               p: QPainter where to plot
        plot_obj: Object to plot; either
                  a list of DataFile or a XYData
           ppars: A dict describing the plot
                 {   "font": QFont("Helvetica",10), ## Font for labels
                     "rect": QRect(),               ## Complete plotting area in pixels
                  "x_label":'Temperature [K],       ## Text for x-axis
                  "y_label":'R [Ohm]',              ## Text for y-axis
	             "flag":0,                      ## Defines where info will be placed
	           "istart":0,                      ## First index to be plotted
	             "iend":-1,                     ## Last index to be plotted
                     "dmin":{"x":0.0,"y":0.0},      ## min-Values of x- and y-axis
                     "dmax":{"x":0.0,"y":0.0},      ## max-Values of x- and y-axis
                     "pmin":{"x":0.0,"y":0.0},      ## min-Values of plotting area
                     "pmax":{"x":0.0,"y":0.0},      ## max-Values of plotting area
	          "n_plots":0,                      ## number of plots (may be obsolete)
	           "symbol":{"size":[],"color":[],"type":[]} ## symbols for each plot
	        }
    """
    FM=QFontMetrics(ppars["font"])
#    FM=p.fontMetrics()
    p.setPen(QPen(Qt.black))
    p.setFont(ppars["font"])
    xDig=7
    yDig=6
    iChx=FM.width('0')
    iChy=FM.lineSpacing()
#    print "cx:",iChx,"cy:",iChy
    if type(plot_obj)!=types.ListType:
       raise TypeError("plot_obj must be a list")

    l=len(plot_obj)
    if l == 0:
       raise ValueError("No data to plot")
    else:
#       pprint(ppars["symbol"])
#       print l,len(ppars["symbol"]["size"])
       if len(ppars["symbol"]["size"]) != l or len(ppars["symbol"]["color"]) != l \
          or len(ppars["symbol"]["type"]) != l:
          raise ValueError("Plot pars: Size of 'symbol' != length of plot list")


    if IsDFileObj(plot_obj[0]):
     if ppars['flag'] & TOP_FILENAME:
#         p.drawText(iChx,iChy,plot_obj[0].FileName)
         l=len(str(plot_obj[0].FileName))
         p.drawText(ppars["rect"].x()+ppars["rect"].width()-iChx*l,iChy,plot_obj[0].FileName)
     if ppars['flag'] & TOP_INFOLINE:
         p.drawText(iChx,iChy,plot_obj[0].Info)
     if ppars['flag'] & TOP_IDN :
         p.drawText(FM.width(ppars['y_label'])+10*iChx,iChy,"Plot from pat (c) HM %s" % ctime())


    xy_obj=[]
    for i in plot_obj:
        if IsDFileObj(i):
#           xy_obj.append(XYData(i.Nums.Export(i.iCx)[ppars["istart"]:ppars["iend"]-1],i.Nums.Export(i.iCy)[ppars["istart"]:ppars["iend"]-1]))
           xy_obj.append(XYData(i.Nums.Export(i.iCx),i.Nums.Export(i.iCy)))
        else: xy_obj.append(i)
    first=1
    for i in xy_obj:
        if first:
           ppars['dmin']['x']=i.Xmin(ppars["istart"],ppars["iend"])[0]
           ppars['dmin']['y']=i.Ymin(ppars["istart"],ppars["iend"])
           ppars['dmax']['x']=i.Xmax(ppars["istart"],ppars["iend"])[0]
           ppars['dmax']['y']=i.Ymax(ppars["istart"],ppars["iend"])
           first=0
        else:
           t=i.Xmin()[0]
           if t < ppars['dmin']['x']: ppars['dmin']['x']=t
           t=i.Ymin()
           if t < ppars['dmin']['y']: ppars['dmin']['y']=t
           t=i.Xmax()[0]
           if t > ppars['dmax']['x']: ppars['dmax']['x']=t
           t=i.Ymax()
           if t > ppars['dmax']['y']: ppars['dmax']['y']=t
    if (ppars['flag'] & TOP_FILENAME) or (ppars['flag'] & TOP_IDN) or (ppars['flag'] & TOP_INFOLINE):
       itop=3
    else: itop=1.5              
   
    EffA=QRect(iChx*(yDig+3)+ppars['rect'].x(),itop*iChy+ppars['rect'].y(),\
    	       ppars['rect'].width()-(yDig+5)*iChx,ppars['rect'].height()-(2+itop)*iChy)
    #EffA=QRect(iChx*(yDig+7)+ppars['rect'].x(),itop*iChy+ppars['rect'].y(),\
    #           ppars['rect'].width()-(yDig+9)*iChx,ppars['rect'].height()-(2+itop)*iChy)
#    print "y:",EffA.height()/iChy
#    print "x:",EffA.width()/iChx
#    p.drawRect(ppars['rect'])
    #raw y-Axis
    p.drawLine(EffA.x(),EffA.y(),EffA.x(),\
               EffA.y()+EffA.height()+iChy/4)

    p.drawLine(EffA.x()+EffA.width(),EffA.y(),\
               EffA.x()+EffA.width(),EffA.y()+EffA.height())

    #raw x-Axis
    p.drawLine(EffA.x()-iChx/2,EffA.height()+EffA.y(),\
               EffA.x()+EffA.width(),EffA.height()+EffA.y());

    p.drawLine(EffA.x(),EffA.y(),\
               EffA.x()+EffA.width(), EffA.y());
    # scale y-axis
    ya=axis( (ppars['dmin']['y'],ppars['dmax']['y']),7,EffA.height()/iChy,'y')

    # scale x-axis
    xa=axis( (ppars['dmin']['x'],ppars['dmax']['x']),7,EffA.width()/iChx,'x')
    # y-Axis label
#    p.drawText(ppars['rect'].x()+5*iChx,ppars['rect'].y()+(int)((itop+0.5)*iChy/4),ppars['y_label'])
    p.drawText(ppars['rect'].x()+4*iChx,ppars['rect'].y()+(int)((itop+3.)*iChy/4),ppars['y_label'])
    # y-Axis ticks
    ir=len(ya.scale)
    #print ya.scale
    for i in range(ir):
       tick=EffA.y()+EffA.height()-(i*EffA.height()/(ir-1))
       p.drawLine(EffA.x()-iChx/2,tick,EffA.x(),tick)
       p.drawText(EffA.x()-FM.width(ya.scale[i][0])-iChx,tick+iChy/4,ya.scale[i][0])

    # x-Axis label
    p.drawText(ppars['rect'].x()+ppars['rect'].width()-iChx*(len(ppars['x_label'])+2),\
               ppars['rect'].y()+ppars['rect'].height()-5,ppars['x_label'])
    # x-Axis ticks
    ir=len(xa.scale)
    #print xa.scale
    for i in range(ir):
       tick=EffA.x()+(i*EffA.width()/(ir-1))
       p.drawLine(tick,EffA.height()+EffA.y(),tick,EffA.height()+EffA.y()+iChy/4)
       p.drawText(tick-FM.width(xa.scale[i][0])/2,EffA.y()+EffA.height()+iChy,xa.scale[i][0])

    ppars['dmax']['y']=ya.max 
    ppars['dmin']['y']=ya.min 
    ppars['dmax']['x']=xa.max 
    ppars['dmin']['x']=xa.min 

    pmin=QPoint(EffA.x(),EffA.height()+EffA.y())
    pmax=QPoint(EffA.width()+EffA.x(),EffA.y())

    ppars['pmax']['y']=pmax.y()
    ppars['pmin']['y']=pmin.y() 
    ppars['pmax']['x']=pmax.x() 
    ppars['pmin']['x']=pmin.x()
    
#    print (pmin.x(),pmax.x()),(pmin.y(),pmax.y())
    if ppars['flag'] &  ZERO_X_LINE:
       y0=LinTrafo( (pmin.y(),pmax.y()), (ppars['dmin']['y'],ppars['dmax']['y']), 0)
       p.setPen(Qt.red)
       p.drawLine(EffA.x(),y0, EffA.x()+EffA.width(),y0)
       p.setPen(Qt.black)

    # Perform plot(s)
    plot_list=[]
    pcount=0
    #Set QArray
    for i in xy_obj:
        if pcount > MAX_PLOTS:
               print "Plot %d ignored; max is %d" % (pcount+1,MAX_PLOTS)
               continue
        pa=QPointArray(len(i))
        for j in range(len(i)):
            pa.setPoint(j, LinTrafo( (pmin.x(),pmax.x()), (ppars['dmin']['x'],ppars['dmax']['x']), i.x[j][0]),\
                           LinTrafo( (pmin.y(),pmax.y()), (ppars['dmin']['y'],ppars['dmax']['y']), i.y[j]))
        plot_list.append(pa)
        pcount+=1

    ep=QPoint()
    if ppars['iend']==-1: npoints=-1
    else: npoints=ppars['iend']-ppars['istart']+1
    for i in range(len(plot_list)):
        p.setPen(ppars['symbol']['color'][i])
  
        if ppars['symbol']['size'][i] == 0: im=MinDist(plot_list[i])
        else: im=4*ppars['symbol']['size'][i]
#        print im,ppars['symbol']['type'][i]
        if im<3 and not (ppars['symbol']['type'][i]=="ERROR_BAR" or ppars['symbol']['type']=="BAR"):
           p.drawPoints(plot_list[i],ppars['istart'],npoints)
#            if (u==iPCount-1 && PP.pSize)PP.pSize[u]=1;
        else: 
           if ppars['symbol']['type'][i]=="ERROR_BAR":
              end=ppars['iend']
              if ppars['iend']==-1:end=plot_list[i].size()
              for ip in range(ppars['istart'],end):
                  y=xy_obj[i].y[ip]
                  e=0
                  if y>0: e=sqrt(y)
                  iy=LinTrafo( (pmin.y(),pmax.y()), (ppars['dmin']['y'],ppars['dmax']['y']),y)
                  ie=LinTrafo((pmin.y(),pmax.y()), (ppars['dmin']['y'],ppars['dmax']['y']),y+e)
                  t=ie-iy
                  if t==0: t=2
                  DrawSymbol("V_LINE",t,p,plot_list[i].point(ip))

           elif ppars['symbol']['type'][i]=="BAR":
              end=ppars['iend']
              if ppars['iend']==-1:end=plot_list[i].size()
              for ip in range(ppars['istart'],end):
	          DrawSymbol("BAR",1,p,plot_list[i].point(ip),
                              QPoint(plot_list[i].point(ip).x(),EffA.height()+EffA.y()))
#
           else:
              if im>16: im=16
              DrawSymbols(ppars['symbol']['type'][i],im/4,plot_list[i],p,ppars['istart'],ppars['iend'])
#                   if(u==iPCount-1)PP.pSize[u]=im/4;
        if ppars['flag'] & BOT_FILELIST:
	   w=list(xy_obj[i].x)
	   w.sort()
	   ep=plot_list[i].point(w[-1][1]-1)
	   pp_x=EffA.width()+EffA.x()
	   pp_y=ep[1]
           p.drawText(pp_x-2*iChx,pp_y,"%d" % (i+1))
	   if ppars['symbol']['type'][i]=='ERROR_BAR': sy='V_LINE'
	   else: sy=ppars['symbol']['type'][i]
           DrawSymbol(sy,3,p,(pp_x-0.5*iChx,pp_y-3))
	   if IsDFileObj(plot_obj[i]):
	      fn=plot_obj[i].FileName
	      text="%2d: %s (%d points)>%s" % (i+1, fn[-30:],len(w),plot_obj[i].Info)
	   else:
	      fn='No filename available'
	      text="%2d: %s (%d points)>No info" % (i+1, fn[-30:],len(w))
#	   print fn
           l=len(fn)
           p.drawText(5*iChx,ppars['rect'].height()+ppars['rect'].y()+iChy+iChy*i,text)

        if ppars['flag'] & BOT_FILEHEAD and len(plot_obj)==1:
	   hf=QFont(ppars['font'])
#	   hf.setPointSize(4)
	   hf.setPointSize(6)
	   
           if IsDFileObj(plot_obj[i]):
	      text1=plot_obj[i].InfoText()
	      text2=plot_obj[i].FText
	   else: text1='No info available'; text2='No info available'
	   r=DrawFormText(p,hf,10*iChx,ppars['rect'].height(),Qt.AlignTop|Qt.AlignLeft,text1)
#	   DrawFormText(p,hf,5*iChx+r[0],ppars['rect'].height(),Qt.AlignTop|Qt.AlignLeft,text2)
	   DrawFormText(p,hf,10*iChx,ppars['rect'].height()+r[1],Qt.AlignTop|Qt.AlignLeft,text2)
    return list(plot_list)
#  	   }
#	    
#        }// for (u)
#    PP.PRect=EffA;
#
#   }
# -------------------------------
# ShowWid
# -------------------------------
class ShowWid(QWidget):
  w=0
  h=0
  def __init__(self,df,parent=None,name=None):
     QWidget.__init__(self,parent,name)
#   Change here if initial window size of Show.py does not fit
#   self.w: Width of initial window in parts of the full screen width
#   self.h: Height of initial window in parts of the full screen height
     self.w=SCREEN.x()/4
     self.h=SCREEN.y()/4

     self.Quit = END_QUIT
     self.iColX = -1
     self.iColY = -1
     self.FontSize = self.height()/24
     
     self.TD = [df]
     if self.TD[0].nNLines == 0:
       raise ValueError("Data File empty")
     self.InitPlot()

  def InitPlot(self):
     self.setCaption(str(self.TD[0].FileName))
     self.setGeometry(SCREEN.x()-self.w-30,30,self.w,self.h);

     self.iColX=self.TD[0].iCx
     self.iColY=self.TD[0].iCy
 
#   Change here if font size of Show.py does not fit
#   Fontsize in parts of the window height
     self.FontSize = self.height()/24
     
     self.PPar={   "font": QFont("Helvetica",10),
                   "rect": QRect(),
                "x_label":self.TD[0].ColId[self.iColX],
                "y_label":self.TD[0].ColId[self.iColY],
	           "flag":0,
	         "istart":0,
	           "iend":-1,
                   "dmin":{"x":0.0,"y":0.0},
                   "dmax":{"x":0.0,"y":0.0},
                   "pmin":{"x":0.0,"y":0.0},
                   "pmax":{"x":0.0,"y":0.0},
	        "n_plots":0,
	         "symbol":{"size":[0],"color":[Qt.black],"type":['CIRCLE']}
	      }
	 		
     if self.TD[0].FileType == FT_SXSMulti and MOD_ftypes['sxs']: self.PPar["symbol"]["type"][0]="ERROR_BAR"  

     self.PPar["rect"].setRect(3,3,self.width()-3,self.height()-3)
#     pprint(self.PPar)
     self.show()

  def Close(self):
     if self.Quit == END_QUIT: qApp.quit()
     else: self.hide()

  def keyPressEvent(self,e):
      if e.key() == Qt.Key_Q: self.Close()
      elif e.key() == Qt.Key_Escape: self.Close()
      elif e.key() == Qt.Key_C:
          if self.PPar["symbol"]["type"][0] != "CONNECT":
             self.PPar["symbol"]["type"][0] = "CONNECT"
          else: self.PPar["symbol"]["type"][0] = "CIRCLE"
          self.repaint()
      elif e.key() == Qt.Key_E:
          if self.PPar["symbol"]["type"][0] != "ERROR_BAR":
             self.PPar["symbol"]["type"][0] = "ERROR_BAR"
          else: self.PPar["symbol"]["type"][0] = "CIRCLE"
          self.repaint()
      elif e.key() == Qt.Key_Plus:
          if self.PPar["symbol"]["size"][0] >= 0:
             self.PPar["symbol"]["size"][0]+=1
          self.repaint()
      elif e.key() == Qt.Key_Minus:
          if self.PPar["symbol"]["size"][0] > 0:
             self.PPar["symbol"]["size"][0]-=1
          self.repaint()
      elif e.key() == Qt.Key_P:
          self.PrintIt()
      elif e.key() == Qt.Key_0:
          if self.PPar['flag'] & ZERO_X_LINE: self.PPar['flag'] &= ~ZERO_X_LINE
          else: self.PPar['flag'] |= ZERO_X_LINE
          self.repaint()
      elif e.key() == Qt.Key_Up:
          self.FontSize+=1
	  self.repaint()
      elif e.key() == Qt.Key_Down: 
          if self.FontSize>4: self.FontSize-=1
	  self.repaint()
      else: e.ignore()

  def paintEvent(self,p):
      paint = QPainter()
      if not paint.begin(self):
         print "Painter: begin failed"
#      iF=self.height()/17
#      iF=self.height()/24
##     if iF==9: iF=10
##      iF=iF/2*2-2
#      if iF<8: iF=8
#      if iF>30: iF=30
#      self.PPar["font"]=QFont("sansserif",iF);
      self.PPar["font"]=QFont("sansserif",self.FontSize);
#      print "Fontsize:",iF
      self.PPar["rect"].setRect(3,3,self.width()-3,self.height()-3)
      DrawData(paint,self.TD,self.PPar)

      paint.setFont(QFont("sanserif",8))
      paint.drawText(3,self.height()-3,str(self.TD[0].FileName))
      paint.end()

  def PrintIt(self):
      Prt = QPrinter()
      if Prt.setup(self):
         paint = QPainter()
         paint.begin(Prt)

#         iF=self.height()/17
         iF=self.height()/24
#         if iF==9: iF=10
#         if iF<8: iF=8
         self.PPar["font"]=QFont("Helvetica",iF);
#         print "Fontsize:",iF
         self.PPar["rect"].setRect(3,3,self.width()-3,self.height()-3);
         DrawData(paint,self.TD,self.PPar);

         paint.setFont(QFont("Helvetica",8))
         paint.drawText(3,self.height()-3,str(self.TD[0].FileName))
         paint.end()
# -------------------------------
# PlotData
# -------------------------------
B=QRect(10,20, 90, 20)

class PlotData(QWidget):
 def __init__(self,df,parent=None,name=None):
    QWidget.__init__(self,parent,name)
#    pdm=QPaintDeviceMetrics( QApplication.desktop() )
#    print pdm.widthMM(),pdm.heightMM()
#   self.w: Width of window in pixels
#   self.h: Height of window in pixels
#   Change here if window size of Plot.py does not fit
    if SCREEN.y() == 600:
       self.w=530
       self.h=530
    elif SCREEN.y() == 768:
       self.w=650
       self.h=700
    elif SCREEN.y() == 1024:
       self.w=650*1.3
       self.h=700*1.3
    WMain=QRect(10,10,self.w,self.h)
    self.setGeometry(WMain)
    self.setCaption("Plot data")

#   Change here if font size of Plot.py does not fit
    self.PPar={   "font": QFont("sanserif",10),
    		  "rect": QRect(),
    	       "x_label":df[0].ColId[df[0].iCx],
    	       "y_label":df[0].ColId[df[0].iCy],
        	  "flag": TOP_IDN | BOT_FILELIST,
        	"istart":0,
        	  "iend":-1,
    		  "dmin":{"x":0.0,"y":0.0},
    		  "dmax":{"x":0.0,"y":0.0},
    		  "pmin":{"x":0.0,"y":0.0},
    		  "pmax":{"x":0.0,"y":0.0},
               "n_plots":len(df),
        	"symbol":{"size":[],"color":[],"type":[]}
             }
    self.df=df          
    for i in range(0,len(df)):
       if df[i].FileType == FT_SXSMulti and MOD_ftypes['sxs']:
         self.PPar["symbol"]["type"].append("ERROR_BAR")  
       else: self.PPar['symbol']['type'].append(i)
       self.PPar['symbol']['size'].append(0)
       self.PPar['symbol']['color'].append(COLORS[i])
#    pprint(self.PPar)
#   Label Edit
#   x-Label
    WxL=QRect(B)
    WxE=QRect(B.x(), WxL.y()+WxL.height()+10,B.width(),B.height())
 
    xMsg=QLabel("x-Label:",self)
    xMsg.setAlignment(Qt.AlignLeft)
    xMsg.setFrameStyle( QFrame.Panel | QFrame.Sunken )
    xMsg.setGeometry(WxL)
    self.xLab = QLineEdit(self)
    self.xLab.setGeometry(WxE)
    self.xLab.setText(self.PPar["x_label"])
    QToolTip.add( self.xLab, "Enter new text for x-Label" )
    self.connect(self.xLab,SIGNAL("returnPressed()"),self.DataChanged )

#   y-Label
    WyL=QRect(B.x(), WxE.y()+WxE.height(),B.width(),B.height())
    WyE=QRect(B.x(), WyL.y()+WyL.height()+10,B.width(),B.height())

    yMsg=QLabel("y-Label:",self)
    yMsg.setAlignment(Qt.AlignLeft)
    yMsg.setFrameStyle( QFrame.Panel | QFrame.Sunken )
    yMsg.setGeometry(WyL)
    self.yLab = QLineEdit(self)
    self.yLab.setGeometry(WyE)
    self.yLab.setText(self.PPar["y_label"])
    QToolTip.add( self.yLab, "Enter new text for y-Label" )
    self.connect(self.yLab,SIGNAL("returnPressed()"),self.DataChanged )

    WFrame=QRect(2*B.x()+B.width(),B.y(),\
                 WMain.width()-3*B.x()-B.width(),\
	         WMain.height()-2*B.y())

    WPlot=QRect(0,0,WFrame.width()-5,self.h*3/4)
    self.PPar["rect"]=WPlot
    	 
    PlotF = QFrame(self,"PlotFrame")
    PlotF.setGeometry(WFrame)
    PlotF.setFrameStyle( QFrame.Panel | QFrame.Sunken)
    PlotF.setPalette(QPalette(QColor(255,255,255)))
    PlotF.setLineWidth(3)

    self.PlotA=MultiPlot(df,self.PPar,self,name="PlotArea");
    self.PlotA.setGeometry(WFrame.x()+3,WFrame.y()+3,WFrame.width()-3,WFrame.height()-3)
    self.PlotA.setPalette(QPalette(QColor(255,255,255)))
    self.PlotA.show()
    self.show()

#   Push Buttons 
#   Print error bars
    WErrorB=QRect(B.x(), WMain.height()-4*B.y()-4*B.height(),B.width(),B.height())
    self.ErrorB = QPushButton("error bars",self)
    self.ErrorB.setFont(QFont("Helvetica",10))
    if self.PPar['symbol']['type'][0]=='ERROR_BAR':
       text="Symbols"
    else: text="Error bars"
    self.ErrorB.setText(text)
    self.ErrorB.setGeometry(WErrorB)
    self.ErrorB.setFixedHeight( self.ErrorB.sizeHint().height() )
    self.connect(self.ErrorB, SIGNAL("clicked()"), self.ErrorBars )
    self.ErrorB.show()
    QToolTip.add( self.ErrorB, "Toggle error bars/symbols" )

#   Print in archive format
    if len(self.df)==1:
       WArchiveB=QRect(B.x(),WMain.height()-6*B.y()-6*B.height(),B.width(),B.height())
       self.ArchiveB = QPushButton("Form",self)
       self.ArchiveB.setFont(QFont("Helvetica",10))
       self.ArchiveB.setText("Archive Form")
       self.ArchiveB.setGeometry(WArchiveB)
       self.ArchiveB.setFixedHeight(self.ArchiveB.sizeHint().height() )
       self.connect(self.ArchiveB,SIGNAL("clicked()"), self.ArchiveFormat )
       self.ArchiveB.show()
       QToolTip.add( self.ArchiveB,"Toggle archive/normal format" )
    else: self.ArchiveB=None

#   Print on printer
    WPrintB=QRect(B.x(),WMain.height()-2*B.y()-2*B.height(),B.width(),B.height())
    PrintB = QPushButton("Print",self)
    PrintB.setFont(QFont("Helvetica",10))
    PrintB.setText("Print")
    PrintB.setGeometry(WPrintB)
    PrintB.setFixedHeight( PrintB.sizeHint().height() )
    self.connect(PrintB, SIGNAL("clicked()"), self.PrintIt)
    PrintB.show()
    QToolTip.add(PrintB, "Print graphic" )

#   Close window
    WCloseB=QRect(B.x(),WMain.height()-B.y()-B.height(),B.width(),B.height())
    CloseB = QPushButton(self, "close" )
    CloseB.setFont(QFont("Helvetica",10,QFont.DemiBold))
    CloseB.setText("Close")
    CloseB.setGeometry(WCloseB)
    CloseB.setFixedHeight( CloseB.sizeHint().height() )
    self.connect(CloseB, SIGNAL("clicked()"), self.Close )
    CloseB.show()
    QToolTip.add( CloseB, "Close this window" )
    self.quit=END_CLOSE
    self.Format=NORM_FORM
# -----------------------------------------
 def EndClose(self): self.quit=END_CLOSE
# -----------------------------------------
 def EndExit(self): self.quit=END_QUIT
# -----------------------------------------
 def DataChanged(self):
    self.PPar["x_label"]=self.xLab.text()
    self.PPar["y_label"]=self.yLab.text()
    self.PlotA.repaint()
# -----------------------------------------
 def Close(self):
    if self.quit==END_QUIT: qApp.quit()
    else: self.hide()
# -----------------------------------------
 def ErrorBars(self):
    self.setUpdatesEnabled(1)
    for i in range(0,len(self.PPar['symbol']['type'])):
    	if self.PPar['symbol']['type'][i]=='ERROR_BAR': self.PPar['symbol']['type'][i]=i;e=1
    	else: self.PPar['symbol']['type'][i]='ERROR_BAR';e=0
    if e: self.ErrorB.setText("Error bars")
    else: self.ErrorB.setText("Symbols")
    self.PlotA.update()
    self.repaint()
# -----------------------------------------
 def ArchiveFormat(self):
    if self.PPar['n_plots']!=1:
      print "Archive format only for a single plot possible"
      return
    self.setUpdatesEnabled(1)
    if self.Format==NORM_FORM:
       self.PPar['flag'] &= ~BOT_FILELIST
       self.PPar['flag'] |= BOT_FILEHEAD
       self.ArchiveB.setText("Normal Form")
       self.Format=ARCH_FORM   
       self.PPar['rect'].setHeight(self.h*1/2)
    else:
       self.Format=NORM_FORM   
       self.ArchiveB.setText("Archive Form")
       self.PPar['flag'] |= BOT_FILELIST
       self.PPar['flag'] &= ~BOT_FILEHEAD
       self.PPar['rect'].setHeight(self.h*3/4)

    self.PlotA.update()
    self.repaint()
# -----------------------------------------
 def keyPressEvent(self,e):
   global iP
   if e.key() == Qt.Key_1: iP=0 
   elif e.key() == Qt.Key_2: iP=1
   elif e.key() == Qt.Key_3: iP=2
   elif e.key() == Qt.Key_4: iP=3
   elif e.key() == Qt.Key_5: iP=4
   elif e.key() == Qt.Key_6: iP=5
   elif e.key() == Qt.Key_7: iP=6
   elif e.key() == Qt.Key_8: iP=7
   elif e.key() == Qt.Key_C:
       if iP<len(self.df): 
         if self.PPar['symbol']['type'][iP] != 'CONNECT':
            self.PPar['symbol']['type'][iP]='CONNECT'
         else: self.PPar['symbol']['type'][iP]=iP
         self.PlotA.update()
         self.repaint()
   elif e.key() == Qt.Key_E:
       if iP<len(self.df): 
         if self.PPar['symbol']['type'][iP] != 'ERROR_BAR':
            self.PPar['symbol']['type'][iP]='ERROR_BAR'
         else: self.PPar['symbol']['type'][iP]=iP
         self.PlotA.update()
         self.repaint()
   elif e.key() == Qt.Key_Plus:
       if iP<len(self.df): 
         if self.PPar['symbol']['size'][iP] >=0: self.PPar['symbol']['size'][iP]+=1
#         print "Sybolsize[%d] + : %d\n" % (iP,self.PPar['symbol']['size'][iP])
         self.PlotA.update()
         self.repaint()
   elif e.key() == Qt.Key_Minus:
       if iP<len(self.df): 
         if self.PPar['symbol']['size'][iP] >=0: self.PPar['symbol']['size'][iP]-=1
         self.PlotA.update()
         self.repaint()
   elif e.key() == Qt.Key_Escape: self.Close()
   else: e.ignore()
# -----------------------------------------
 def PrintIt(self):
    Prt=QPrinter()
    if Prt.setup(self.PlotA):
#       pdm=QPaintDeviceMetrics( Prt )
#       print pdm.widthMM(),pdm.heightMM()
       paint=QPainter()
       paint.begin(Prt)
       f=self.PPar['font']
       self.PPar['font']=QFont("sansserif",10)
       DrawData(paint,self.df,self.PPar)
       paint.end()
       self.PPar['font']=f
    else: print "Printer setup failed"
# ------------------------------
# MultiPlot
# ------------------------------
class MultiPlot(QWidget):
 w=0
 h=0
 def __init__(self,df,ppar,parent=None,name=None):
    QWidget.__init__(self,parent,name)
    self.df=df
    self.ppar=ppar
    self.setFocusPolicy(QWidget.ClickFocus )
    self.setFocus()
# ------------------------------
 def paintEvent(self,p):
    p= QPainter()
    p.begin(self)
    DrawData(p,self.df,self.ppar)
    p.end()
# ------------------------------
# ManipData
# ------------------------------
class ManipData(QWidget):
  w=0
  h=0
  def __init__(self,df,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     self.TD = [df]
     if self.TD[0].nNLines == 0:
        raise ValueError("Data File empty")
     self.Button=[{"label":"Help F1",       "tip":"Display Help",                "here":1},
                  {"label":"Error bars F2", "tip":"Toggle error bars",  	 "here":1},
	          {"label":"Meanvalue F3",  "tip":"Meanvalue of neigbour points","here":1},
	          {"label":"Remove F4",     "tip":"Remove actual point",	 "here":1},
	          {"label":"New Start F5",  "tip":"Select a new start point", 	 "here":1},
                  {"label":"New End F6",    "tip":"Select a new end point",	 "here":1},
	          {"label":"Previous F7",   "tip":"Goto previous range",         "here":1},
	          {"label":"Integrate F8",  "tip":"Numerical integral",	         "here":1},
	          {"label":"Save Range F9", "tip":"Save range & exit",		 "here":1},
	          {"label":"Exit F10",      "tip":"Exit now",			 "here":1},] 

     self.w=800
     self.h=600
     self.FontSize = 14
     self.Quit=END_CLOSE
     
     self.WMain=QRect(10,10,self.w,self.h)
     self.InitPlot()
     
     self.setGeometry(self.WMain)
     self.setCaption("Change selected data")
     topLayout = QVBoxLayout( self, 5 )
     
     TopLabel = QHBoxLayout()
     topLayout.addLayout(TopLabel)

     self.InfoLT=QLabel("", self, "TopLabel")
     self.InfoLT.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.InfoLT.setFrameStyle( QFrame.Panel | QFrame.Sunken )
     self.InfoLT.setMinimumSize(self.InfoLT.sizeHint())
     self.InfoLT.setFixedHeight(self.InfoLT.height())
     TopLabel.addWidget(self.InfoLT,100)
 
     Middle = QHBoxLayout()
     topLayout.addLayout(Middle) 

     buttons = QVBoxLayout()
     Middle.addLayout(buttons,1)

     self.Bts=[]
     ct=0
     for i in self.Button:
        if i["here"]:
	   B=QPushButton(self,"B%02d" % ct)
           B.setFont(QFont("Helvetica",self.FontSize,QFont.Light))
           B.setText(i["label"])
	   if ct==0: B0=B
           B.setMinimumSize(B0.sizeHint())
           B.setFixedWidth(B0.sizeHint().width()*2)
           buttons.addWidget(B)
           QToolTip.add(B,i["tip"])
           self.Bts.append(B)
        else:
           B=None
           L=QLabel(self)
           L.setFrameStyle(QFrame.NoFrame)
           L.setMinimumSize(self.Bts[0].sizeHint())
           L.setFixedWidth(self.Bts[0].width())
           buttons.addWidget(L)
           self.Bts.append(B)
        ct+=1

     self.connect(self.Bts[1],SIGNAL("clicked()"),self.ErrorBars)
     self.connect(self.Bts[2],SIGNAL("clicked()"),self.MeanVal)
     self.connect(self.Bts[3],SIGNAL("clicked()"),self.RemovePoint)
     self.connect(self.Bts[4],SIGNAL("clicked()"),self.NewStart)
     self.connect(self.Bts[5],SIGNAL("clicked()"),self.NewEnd)
     self.connect(self.Bts[6],SIGNAL("clicked()"),self.RangeBack)
#     self.connect(self.Bts[7],SIGNAL("clicked()"),self.Integrate)
#     self.connect(self.Bts[8],SIGNAL("clicked()"),self.SaveRange)
     self.connect(self.Bts[9],SIGNAL("clicked()"),self.Close)

     self.PlotF = MouseTrackWid(self,"PlotFrame")
     self.PlotF.setMinimumSize( QSize(self.WPl.width(),self.WPl.height()))
     self.PlotF.setFrameStyle( QFrame.Panel | QFrame.Sunken)
     self.PlotF.setPaletteBackgroundColor(Qt.white)
     self.PlotF.setLineWidth(3)
     self.PPar["rect"]=QRect(0,0,self.PlotF.width(),self.PlotF.height())
     self.connect(self.PlotF,PYSIGNAL("RCoords(QPoint)"),self.ShowCoords)
     self.connect(self.PlotF,PYSIGNAL("LCoords(QPoint)"),self.ClosestPoint)
     Middle.addWidget(self.PlotF,20)
     BotLabel = QHBoxLayout()
     topLayout.addLayout(BotLabel)

     self.InfoLB=QLabel("",self,"BotLabel")
     self.InfoLB.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.InfoLB.setFrameStyle(QFrame.Panel | QFrame.Sunken)
     self.InfoLB.setMinimumSize(self.InfoLB.sizeHint())
     self.InfoLB.setFixedHeight(self.InfoLB.height())
     BotLabel.addWidget(self.InfoLB)
     topLayout.activate()
     self.setFocusPolicy(QWidget.ClickFocus)
     self.show()
# -----------------------------------------
  def InitPlot(self):
     self.iStep=0
     self.First=1
     self.iRCount=0
     self.RSelect=0
     self.POld=QPoint(0,0)
     self.CursPM=QPixmap()
     self.PM0=QPoint(0,0)
     self.CurCurs=QPoint(0,0)
     self.Pl=None
     self.Pd=[self.TD[0]]
     self.iColX=self.TD[0].iCx
     self.iColY=self.TD[0].iCy

     self.WPl=QSize(self.WMain.width()-B.y()-B.width()-20, self.WMain.height()-2*B.height()-40)	 		
     self.T=XYData(self.TD[0].Nums.Export(self.TD[0].iCx),self.TD[0].Nums.Export(self.TD[0].iCy))
     self.InitLen=len(self.T)
     self.FontSize=11
     self.PPar={   "font": QFont("sanserif",self.FontSize),
                   "rect": None,
                "x_label":self.TD[0].ColId[self.iColX],
                "y_label":self.TD[0].ColId[self.iColY],
	           "flag":TOP_FILENAME,
	         "istart":0,
	           "iend":len(self.T),
                   "dmin":{"x":self.T.Xmin()[0],"y":self.T.Ymin()},
                   "dmax":{"x":self.T.Xmax()[0],"y":self.T.Ymax()},
                   "pmin":{"x":0,"y":0},
                   "pmax":{"x":0,"y":0},
	        "n_plots":0,
	         "symbol":{"size":[0],"color":[Qt.black],"type":['CIRCLE']}
	      }
     self.start_end=[]      
#     if self.TD[0].FileType == FT_SXSMulti and MOD_ftypes['sxs']: self.PPar["symbol"]["type"][0]="ERROR_BAR"  
# -----------------------------------------
  def ShowCoords(self,p):
 #    print "R: x=%d y=%d w=%d, h=%d; p: x=%d, y=%d" % (self.PPar["rect"].x(),self.PPar["rect"].y(),self.PPar["rect"].width(),self.PPar["rect"].height(),p.x(),p.y())
     if self.PPar["rect"].contains(p):
       x=BackTrafo((self.PPar["pmin"]["x"],self.PPar["pmax"]["x"]),\
                   (self.PPar["dmin"]["x"],self.PPar["dmax"]["x"]),p.x())
       y=BackTrafo((self.PPar["pmax"]["y"],self.PPar["pmin"]["y"]),\
                   (self.PPar["dmax"]["y"],self.PPar["dmin"]["y"]),p.y())
       a="X: %12.5g Y: %12.5g" %(x,y)
     else: a="Out of plot area"
#     print a
     self.InfoLB.setText(a)
# -----------------------------------------
  def ClosestPoint(self,p):
     if self.PPar["rect"].contains(p):
       x=p.x()
       br=self.Pl[0].boundingRect()
       if x<br.x(): s=self.PPar["istart"]
       elif x>br.x()+br.width(): s=self.PPar["iend"]-1
       else:
         size=5
	 r=QRect(p.x()-size,p.y()-size,2*size,2*size)
	 found=-1
	 d=SqDist(QPoint(self.Pl[0].point(self.PPar["istart"])[0],self.Pl[0].point(self.PPar["istart"])[1]),p)
	 for s in range(self.PPar["istart"],self.PPar["iend"]-1):
	     pp=QPoint(self.Pl[0].point(s)[0],self.Pl[0].point(s)[1])
	     if r.contains(pp):
	       dd=SqDist(p,pp)
	       if dd<d: d=dd;found=s
	 if found!=-1: s=found
	 else:
           d=abs(self.Pl[0].point(self.PPar["istart"])[0]-x)
           for i in range(self.PPar["istart"],self.PPar["iend"]-1):
	      dd=abs(self.Pl[0].point(i)[0]-x)
              if dd<d: d=dd; s=i
       self.iStep=s
       self.MoveCursor(Qt.red)
# -----------------------------------------
  def Close(self):
     if self.Quit==END_CLOSE:
       if self.TD[0].FileType == FT_THECAP or self.TD[0].FileType == FT_SXSMulti:
          self.TD[0].ChgTextPar([('Idn','File from - Disp.py')])
	  if self.InitLen != len(self.T.x):
             self.TD[0].ChgTextPar([('DataPoints','%d' % len(self.T.x))])
       self.TD[0].Write(sys.stdout)
     qApp.quit()
#     else: self.hide()
# -----------------------------------------
  def paintEvent(self,p):
     paint=QPainter()
     paint.begin(self.PlotF)
     self.PPar["font"]=QFont("sanserif",self.FontSize)

     self.Pl=DrawData(paint,self.Pd,self.PPar)
     paint.end()
 
     self.MoveCursor(Qt.red)
# -----------------------------------------
  def MoveCursor(self,c):
     p=QPainter()
     size=10
     self.CurCurs.setX(self.Pl[0].point(self.iStep)[0])
     self.CurCurs.setY(self.Pl[0].point(self.iStep)[1])
#     print "********",self.Pl[0].point(self.iStep)[0],self.Pl[0].point(self.iStep)[1]
     self.PM0.setX(self.CurCurs.x()-size-1)
     self.PM0.setY(self.CurCurs.y()-size-1)
  
     p.begin(self.PlotF)
     if not self.First:
       p.setPen(self.PlotF.paletteBackgroundColor())
       DrawSymbol("CROSS", size, p, (self.POld.x(),self.POld.y()))
       p.drawPixmap(self.POld.x()-size-1,self.POld.y()-size-1,self.CursPM)
     else: self.First=0
#     print self.CurCurs.x(),self.CurCurs.y(),self.POld.x(),self.POld.y()
     self.CursPM=QPixmap.grabWindow(self.PlotF.winId(),self.PM0.x(),self.PM0.y(),2*size+2,2*size+2)
     self.POld.setX(self.CurCurs.x())
     self.POld.setY(self.CurCurs.y())
#     p.drawRect(self.PM0.x(),self.PM0.y(),2*size+2,2*size+2)

     p.setPen(c);
     DrawSymbol("CROSS", size, p, (self.CurCurs.x(),self.CurCurs.y()))
     p.end()
     a = "Point: %d of %d X:%12.5g Y:%12.5g prev. ranges:%d" %\
     	 (self.iStep+1,len(self.T),self.T.x[self.iStep][0],\
     	 self.T.y[self.iStep],self.iRCount)
     self.InfoLT.setText(a)
# -----------------------------------------
  def ErrorBars(self):
     self.PlotF.setUpdatesEnabled(1)
     if self.PPar["symbol"]["type"][0] != "ERROR_BAR":
     	self.PPar["symbol"]["type"][0] = "ERROR_BAR"
        self.Bts[1].setText("Symbols F2")
     else: 
     	self.PPar["symbol"]["type"][0] = "CIRCLE"
        self.Bts[1].setText(self.Button[1]["label"])
     self.PlotF.repaint()
     self.repaint()
# -----------------------------------------
  def MeanVal(self):
     self.PlotF.setUpdatesEnabled(1)
     il=self.iStep-1
     if self.iStep==0: il=1
     if self.iStep==len(self.T): il=len(self.T)-2
 
     ir=self.iStep+1;
     if self.iStep==0: ir=2
     if self.iStep==len(self.T):ir=len(self.T)-1

     DeltaX=self.T.x[ir][0] - self.T.x[il][0]
     DeltaY=self.T.y[ir] - self.T.y[il]

     dx=self.T.x[ir][0]-self.T.x[self.iStep][0]
     dy=DeltaY/DeltaX*dx

     oldY=self.T.y[self.iStep]
     oldX=self.T.x[self.iStep][0]
     self.T.y[self.iStep]=self.T.y[self.iStep+1]-dy
     self.TD[0].Nums.Import(self.iColY,self.T.y)
     self.Pd=[self.TD[0]]
     
#     PaintSymbol(DIAMOND, 3, red, PlotF, D);
     a="#%d Y(%f): %12.5g set to %12.5g" % (self.iStep,\
                     self.T.x[self.iStep][0],oldY,self.T.y[self.iStep])
     self.InfoLB.setText(a)
     self.PlotF.repaint()
     self.repaint()
# ----------------------------------------
  def RemovePoint(self):
     t=self.T.pop(self.iStep)
     self.TD[0].Nums.pop_row(self.iStep)
     self.Pd=[self.TD[0]]
     a="#%d Y(%f): %12.5g removed" % (self.iStep+1,\
                     t[0][0],t[1])
     self.InfoLB.setText(a)
     self.PlotF.repaint()
     self.repaint()
# ----------------------------------------
  def NewStart(self):
     self.PlotF.setUpdatesEnabled(1)
     if self.iStep<self.PPar["iend"]:
       self.start_end.append((self.PPar["istart"],self.PPar["iend"]))
       self.PPar["istart"]=self.iStep
       self.iRCount+=1
       a="New start (%d,%d)" % (self.PPar["istart"],self.PPar["iend"])
     else: 
       a="Bad start %d > %d" % (self.PPar["istart"],self.PPar["iend"])
     self.InfoLB.setText(a)
     self.PlotF.repaint()
     self.repaint()
# -----------------------------------------
  def NewEnd(self):
     self.PlotF.setUpdatesEnabled(1)
     if self.iStep>self.PPar["istart"]+3:
       self.start_end.append((self.PPar["istart"],self.PPar["iend"]))
       self.PPar["iend"]=self.iStep
       self.iRCount+=1
       a="New end (%d,%d)" % (self.PPar["istart"],self.PPar["iend"])
     else: 
       a="Bad end %d < %d + 3" % (self.PPar["iend"],self.PPar["istart"])
     self.InfoLB.setText(a)
     self.PlotF.repaint()
     self.repaint()
# -----------------------------------------
  def RangeBack(self):
     self.PlotF.setUpdatesEnabled(1)
     if len(self.start_end)>0:
       t=self.start_end.pop()
       self.PPar["istart"]=t[0]
       self.PPar["iend"]=t[1]
       self.iRCount-=1
       a="Back (%d,%d)" % (self.PPar["istart"],self.PPar["iend"])
     else: 
       a="Last range reached (%d,%d)" % (self.PPar["iend"],self.PPar["istart"])
     self.InfoLB.setText(a)
     self.PlotF.repaint()
     self.repaint()
# -----------------------------------------
  def keyPressEvent(self,e):
     if e.key() == Qt.Key_Escape:
       self.Quit=END_QUIT 
       self.Close()
     elif e.key() == Qt.Key_Right:
     	 if self.iStep<self.PPar["iend"]:self.iStep+=1
     	 if self.iStep==self.PPar["iend"]:self.iStep=self.PPar["istart"]
         self.MoveCursor(Qt.red)
     elif e.key() == Qt.Key_Left:
     	 if self.iStep>self.PPar["istart"]:self.iStep-=1
     	 elif self.iStep==self.PPar["istart"]:self.iStep=self.PPar["iend"]-1
         self.MoveCursor(Qt.red)
     elif e.key() == Qt.Key_Home:
     	 self.iStep=self.PPar["istart"]
         self.MoveCursor(Qt.red)
     elif e.key() == Qt.Key_End:
     	 self.iStep=self.PPar["iend"]-1
         self.MoveCursor(Qt.red)
     elif e.key() == Qt.Key_C:
     	 if self.PPar["symbol"]["type"][0] != "CONNECT":
     	   self.Pd=[self.TD[0],self.TD[0]]
           self.PPar["symbol"]={"size":[0,0],"color":[Qt.black,Qt.red],"type":["CONNECT",'CIRCLE']}
     	 else: 
            self.PPar["symbol"]={"size":[0],"color":[Qt.black],"type":['CIRCLE']}
     	    self.Pd=[self.TD[0]]
     	 self.PlotF.setUpdatesEnabled(1)
     	 self.PlotF.repaint()
     	 self.repaint()
     elif e.key() == Qt.Key_E:self.ErrorBars()
     elif e.key() == Qt.Key_Plus:
     	 if self.PPar["symbol"]["size"][0] >= 0:
     	    self.PPar["symbol"]["size"][0]+=1
     	 self.PlotF.setUpdatesEnabled(1)
     	 self.PlotF.repaint()
     	 self.repaint()
     elif e.key() == Qt.Key_Minus:
     	 if self.PPar["symbol"]["size"][0] > 0:
     	    self.PPar["symbol"]["size"][0]-=1
     	 self.PlotF.setUpdatesEnabled(1)
     	 self.PlotF.repaint()
     	 self.repaint()
     elif e.key() == Qt.Key_0:
     	 if self.PPar['flag'] & ZERO_X_LINE: self.PPar['flag'] &= ~ZERO_X_LINE
     	 else: self.PPar['flag'] |= ZERO_X_LINE
     	 self.repaint()
     elif e.key() == Qt.Key_Up:
     	 self.FontSize+=1
     	 self.PlotF.setUpdatesEnabled(1)
     	 self.PlotF.repaint()
         self.repaint()
     elif e.key() == Qt.Key_Down: 
     	 if self.FontSize>4: self.FontSize-=1
     	 self.PlotF.setUpdatesEnabled(1)
     	 self.PlotF.repaint()
         self.repaint()
     elif e.key() == Qt.Key_F2: self.ErrorBars()
     elif e.key() == Qt.Key_F3: self.MeanVal()
     elif e.key() == Qt.Key_F4: self.NewStart()
     elif e.key() == Qt.Key_F5: self.NewEnd()
     elif e.key() == Qt.Key_F6: self.RangeBack()
     else: e.ignore()

# ------------------------------
# MouseTrackWid
# ------------------------------
class MouseTrackWid(QFrame):
  def __init__(self,parent=None,name=None):
     QFrame.__init__(self,parent,name)
     self.RDown=0
     self.setCursor(Qt.crossCursor)

  def mousePressEvent(self,e):
     if e.button()==Qt.RightButton: self.RDown=1

  def mouseReleaseEvent(self,e):
     if e.button()==Qt.RightButton: self.RDown=0
     if e.button()==Qt.LeftButton: self.emit(PYSIGNAL("LCoords(QPoint)"), (QPoint(e.pos()), ))  # emit LCoords(QPoint(e.pos()))
  
  def mouseMoveEvent(self,e):
     if self.RDown: self.emit(PYSIGNAL("RCoords(QPoint)"), (QPoint(e.pos()), ))  # emit RCoords(QPoint(e.pos()))
# -------------------------------
# RescaleWid
# -------------------------------
class RescaleWid(QWidget):
  w=0
  h=0
  def __init__(self,filename,parent=None,name=None):
     QWidget.__init__(self,parent,name)
      
     self.Button=[{"label":"Help F1",   "tip":"Display Help", "here":1},
                  {"label":"Close F10", "tip":"Close",	      "here":1},
                ] 

     self.RButton=[{"label":"Set Xmin", "tip":"XMin of scale can be set"},
                   {"label":"Set Xmax", "tip":"XMax of scale can be set"},
                   {"label":"Set Ymin", "tip":"YMin of scale can be set"},
                   {"label":"Set Ymax", "tip":"YMax of scale can be set"},
                   {"label":"Pick up",  "tip":"Pick up data points"},
                  ]

     self.w=SCREEN.x()*0.75
     self.h=SCREEN.y()*0.75
     self.FontSize = 10
     self.Filename=filename
     self.InfoTxt=QPixmap.imageFormat(self.Filename)
     if self.InfoTxt==None:
        raise TypeError("%s not a valid grafic file (png,jpg,gif)" % filename)

     self.WMain=QRect(10,10,self.w,self.h)
     
     self.setGeometry(self.WMain)
     self.setCaption("Rescale Plot")
     topLayout = QVBoxLayout( self, 5 )
     
     TopLabel = QHBoxLayout()
     topLayout.addLayout(TopLabel)

     self.InfoLT=QLabel("", self, "TopLabel")
     self.InfoLT.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.InfoLT.setFrameStyle( QFrame.Panel | QFrame.Sunken )
     self.InfoLT.setMinimumSize(self.InfoLT.sizeHint())
     self.InfoLT.setFixedHeight(self.InfoLT.height())
     self.InfoLT.setFont(QFont("sansserif",10))
     TopLabel.addWidget(self.InfoLT,100)
     
     Middle = QHBoxLayout()
     topLayout.addLayout(Middle) 

     buttons = QVBoxLayout()
     Middle.addLayout(buttons,1)

     self.Bts=[]
     ct=0
     for i in self.Button:
        if i["here"]:
	   B=QPushButton(self,"B%02d" % ct)
           B.setFont(QFont("Helvetica",self.FontSize,QFont.Light))
           B.setText(i["label"])
	   if ct==0: B0=B
           B.setMinimumSize(B0.sizeHint())
           B.setFixedWidth(B0.sizeHint().width()*2)
           buttons.addWidget(B)
           QToolTip.add(B,i["tip"])
           self.Bts.append(B)
        else:
           B=None
           L=QLabel(self)
           L.setFrameStyle(QFrame.NoFrame)
           L.setMinimumSize(self.Bts[0].sizeHint())
           L.setFixedWidth(self.Bts[0].width())
           buttons.addWidget(L)
           self.Bts.append(B)
        ct+=1

     self.connect(self.Bts[1],SIGNAL("clicked()"),self.Close)

     self.bg = QButtonGroup( self, "radioGroup" )
#     self.bg.setTitle( "Radio buttons" )

     buttons.addWidget( self.bg)
     self.vbox = QVBoxLayout( self.bg, 10 )
     self.vbox.addSpacing( self.bg.fontMetrics().height() )

     self.Edit=[]
     ct=0
     for i in self.RButton:
         self.rb = QRadioButton(self.bg)
         self.rb.setText(i["label"])
         if ct==0: self.rb.setChecked(1)
         self.vbox.addWidget( self.rb )
         self.rb.setMinimumSize( self.rb.sizeHint() )
         QToolTip.add( self.rb, i["tip"] )
         le=QLineEdit( self.bg, "lineEdit%d" % ct)
         self.Edit.append(le)
         self.vbox.addWidget( le )
         le.setFixedHeight( le.sizeHint().height() )
	 if ct==4: QToolTip.add(le, "Enter filename to save data" )
         else: QToolTip.add(le, "Enter data value and click left" )
         
	 ct+=1

     self.vbox.activate()

     self.RB=0
     self.connect( self.bg, SIGNAL("clicked(int)"), self.radioButtonClicked )
     self.connect( self.Edit[4], SIGNAL("returnPressed()"), self.NewFilename )


     self.WPl=QSize(self.WMain.width()-B.y()-B.width()-20, self.WMain.height()-2*B.height()-40)	 		

     self.PlotF = PickUpWid(self,"PlotLabel")
     self.PlotF.setMinimumSize( QSize(self.WPl.width(),self.WPl.height()))
     self.PlotF.setFrameStyle( QFrame.Panel | QFrame.Sunken)
#     self.PlotF.setPaletteBackgroundColor(Qt.white)
#     self.PlotF.setLineWidth(3)
     self.connect(self.PlotF,PYSIGNAL("MCoords(QPoint)"),self.ShowCoords)
     self.connect(self.PlotF,PYSIGNAL("LCoords(QPoint)"),self.GetMinMax)
     Middle.addWidget(self.PlotF,20)
     BotLabel = QHBoxLayout()
     topLayout.addLayout(BotLabel)

     self.InfoLB=QLabel("",self,"BotLabel")
     self.InfoLB.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.InfoLB.setFrameStyle(QFrame.Panel | QFrame.Sunken)
     self.InfoLB.setMinimumSize(self.InfoLB.sizeHint())
     self.InfoLB.setFixedHeight(self.InfoLB.height())
     BotLabel.addWidget(self.InfoLB)
     topLayout.activate()

     self.pm = QPixmap(self.Filename)
     img=self.pm.convertToImage()
     if img.width()>self.PlotF.width() or img.height()>self.PlotF.height():
        nimg=img.smoothScale(self.PlotF.width(),self.PlotF.height(),QImage.ScaleMin)
     else: nimg=img
     
     self.Size=(self.pm.width(),self.pm.height())
     self.InfoTxt="%s %s %sx%s -> %sx%s" % \
                (self.Filename,self.InfoTxt,self.Size[0],self.Size[1],self.PlotF.width(),self.PlotF.height())
     self.InfoLT.setText(self.InfoTxt) 
     self.pm.convertFromImage(nimg)
     self.PlotF.setPixmap(self.pm)
     self.PlotF.setFixedSize(self.PlotF.size())

     self.setFocusPolicy(QWidget.ClickFocus)

     self.PPar={   "dmin":{"x":None,"y":None},
                   "dmax":{"x":None,"y":None},
                   "pmin":{"x":0,"y":0},
                   "pmax":{"x":0,"y":0},
                   "rect":QRect(0,0,self.PlotF.width(),self.PlotF.height())
	       }
    
     self.PickW=None
     self.show()
# -----------------------------------------
  def radioButtonClicked( self, id ):
      self.InfoLB.setText( "#%d: %s" % (id,self.RButton[id]["label"] ))
      self.RB=id
# -----------------------------------------
  def ToFloat(self,text):
      try: r=atof(text)
      except: return None
      return r
# -----------------------------------------
  def NewFilename(self):
     fn=self.Edit[4].text()
     if self.PickW!=None and fn!="": self.PickW.Filename=str(fn)
# -----------------------------------------
  def GetMinMax(self,p):
      Text=str(self.Edit[self.RB].text())
      if self.RB==0:
         self.PPar["pmin"]["x"]=p.x()
         self.PPar["dmin"]["x"]=self.ToFloat(Text)
      elif self.RB==1:
         self.PPar["pmax"]["x"]=p.x()
         self.PPar["dmax"]["x"]=self.ToFloat(Text)
      elif self.RB==2:
         self.PPar["pmax"]["y"]=p.y()
         self.PPar["dmax"]["y"]=self.ToFloat(Text)
      elif self.RB==3:
         self.PPar["pmin"]["y"]=p.y()
         self.PPar["dmin"]["y"]=self.ToFloat(Text)
      elif self.RB==4:
         if self.PPar["dmin"]["x"]==None or self.PPar["dmax"]["x"]==None or \
	    self.PPar["dmax"]["y"]==None or self.PPar["dmin"]["y"]==None:
	    QMessageBox.about(self,"Rescale","Not all x-y min/max are set.\nPlease set min/max first")
         else:
	    if self.PickW==None:
	       self.PickW=EditWid(Text, None, "multiLineEdit" )
	       self.PickW.setWFlags(Qt.WType_TopLevel)
	       self.PickW.show()
	    if self.PickW.isHidden():
	       self.PickW.clear()
	       self.PickW.show()
	    self.PickW.Filename=Text
            X=BackTrafo((self.PPar["pmin"]["x"],self.PPar["pmax"]["x"]),\
                      (self.PPar["dmin"]["x"],self.PPar["dmax"]["x"]),p.x())
            Y=BackTrafo((self.PPar["pmin"]["y"],self.PPar["pmax"]["y"]),\
                      (self.PPar["dmin"]["y"],self.PPar["dmax"]["y"]),p.y())
	    self.PickW.insert("%g  %g\n" %(X,Y))
      lText=self.InfoTxt
      if self.PPar["dmin"]["x"]!=None: lText+=" x(min): %g @ %d" % (self.PPar["dmin"]["x"],self.PPar["pmin"]["x"])	 
      if self.PPar["dmax"]["x"]!=None: lText+=" x(max): %g @ %d" % (self.PPar["dmax"]["x"],self.PPar["pmax"]["x"])	 
      if self.PPar["dmax"]["y"]!=None: lText+=" y(min): %g @ %d" % (self.PPar["dmax"]["y"],self.PPar["pmax"]["y"])	 
      if self.PPar["dmin"]["y"]!=None: lText+=" y(max): %g @ %d" % (self.PPar["dmin"]["y"],self.PPar["pmin"]["y"])	 
      self.InfoLT.setText(lText)
# -----------------------------------------
  def ShowCoords(self,p):
     if self.PPar["rect"].contains(p):
       a="x: %d y: %d" %(p.x(),p.y())
       if not (self.PPar["dmin"]["x"]==None or self.PPar["dmax"]["x"]==None):
          X=BackTrafo((self.PPar["pmin"]["x"],self.PPar["pmax"]["x"]),\
                      (self.PPar["dmin"]["x"],self.PPar["dmax"]["x"]),p.x())
          a=a+" X: %g" % X
       if not (self.PPar["dmin"]["y"]==None or self.PPar["dmax"]["y"]==None):
          Y=BackTrafo((self.PPar["pmin"]["y"],self.PPar["pmax"]["y"]),\
                      (self.PPar["dmin"]["y"],self.PPar["dmax"]["y"]),p.y())
          a=a+" Y: %g" % Y
     else: a="Out of plot area"
#     print a
     self.InfoLB.setText(a)
# -----------------------------------------
  def Close(self):
     qApp.quit()
# -----------------------------------------
# ------------------------------
# PickUpWid
# ------------------------------
class PickUpWid(QLabel):
  def __init__(self,parent=None,name=None):
     QLabel.__init__(self,parent,name)
#     self.LDown=0
     self.setMouseTracking(1)
#     self.RDown=0
     self.setCursor(Qt.crossCursor)

  def mousePressEvent(self,e):
#     if e.button()==Qt.RightButton: self.RDown=1
     if e.button()==Qt.LeftButton:  self.LDown=1

  def mouseReleaseEvent(self,e):
#     if e.button()==Qt.RightButton: self.RDown=0
#     if e.button()==Qt.LeftButton:  self.LDown=0
     if e.button()==Qt.LeftButton: self.emit(PYSIGNAL("LCoords(QPoint)"), (QPoint(e.pos()), ))  # emit LCoords(QPoint(e.pos()))
#  
  def mouseMoveEvent(self,e):
#     if self.LDown: self.emit(PYSIGNAL("MCoords(QPoint)"), (QPoint(e.pos()), ))  # emit RCoords(QPoint(e.pos()))
#     if self.RDown: self.emit(PYSIGNAL("RCoords(QPoint)"), (QPoint(e.pos()), ))  # emit RCoords(QPoint(e.pos()))
     self.emit(PYSIGNAL("MCoords(QPoint)"), (QPoint(e.pos()), ))  # emit RCoords(QPoint(e.pos()))
# ------------------------------
# PickUpWid
# ------------------------------
class EditWid(QMultiLineEdit):
  def __init__(self,filename,parent=None,name=None):
     QMultiLineEdit.__init__(self,parent,name)
     self.Filename=filename
#     self.setWFlags(Qt.WType_TopLevel)
     self.setBaseSize(SCREEN.x()*0.25,SCREEN.y()*0.7)
  def closeEvent(self,e):
     if self.Filename=="" and self.text()!="":
        a=QMessageBox.question(self,"Rescale","No filename specified all data will be lost.\nContinue ?",\
	                       QMessageBox.Yes,QMessageBox.No)
        if a==QMessageBox.No: return
	else: e.accept()
	
     try: fp=open(self.Filename,"r")
     except: a=QMessageBox.Yes
     else:
        a=QMessageBox.question(self,"Rescale","Output file %s already exists!\nOverwrite ?" % self.Filename,\
	                       QMessageBox.Yes,QMessageBox.No)
        fp.close()
     if a==QMessageBox.Yes:
        fp=open(self.Filename,"w")
        fp.write(str(self.text()))
        fp.close()
        e.accept()
     else: 
        QMessageBox.about(self,"Rescale","Please change filename in pick up section")
        return
#if __name__=='__main__':
