import sys,types,os
from qt import *
from qttable import *

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
from grafic import *
from popen2 import *

#$Log: wid.py,v $
#Revision 1.1  2006/01/04 14:41:36  herbie
#Initial revision
#
#
CVS_ID="$Id: wid.py,v 1.1 2006/01/04 14:41:36 herbie Exp herbie $"

lib_path='PAT_LIB_PATH'
if os.name=='posix': PATH_DEL="/"
else: PATH_DEL="\\"
HIST_FILE=os.environ[lib_path]+PATH_DEL+".cmd_history"

# -------------------------------
# ColInfoWid
# -------------------------------
class ColInfoWid(QWidget):
  def __init__(self,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     self.FontSize = 8
     self.Font=QFont("sanserif",10)
     self.td=None
     
     # Info label, x-col info, y-col info
     topLayout = QVBoxLayout(self, 3)
     
     # No of Cols, No of Rows
     Info = QHBoxLayout(2)
     topLayout.addLayout(Info)

     # Label, x-column box
     xBox = QHBoxLayout(2)
     topLayout.addLayout(xBox)
 
     # Label, y-column box
     yBox = QHBoxLayout(2)
     topLayout.addLayout(yBox)
 
     self.RowLabel=QLabel("", self, "RowLabel")
     self.RowLabel.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.RowLabel.setFrameStyle( QFrame.Panel | QFrame.Sunken )
     self.RowLabel.setFont(self.Font)
     self.RowLabel.setText("Rows: -")
     self.RowLabel.setMinimumSize(self.RowLabel.sizeHint())
     self.RowLabel.setFixedHeight(self.RowLabel.height())
     self.RowLabel.show()
     Info.addWidget(self.RowLabel)
     
     self.ColLabel=QLabel("", self, "ColLabel")
     self.ColLabel.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.ColLabel.setFrameStyle( QFrame.Panel | QFrame.Sunken )
     self.ColLabel.setFont(self.Font)
     self.ColLabel.setText("Cols: -")
     self.ColLabel.setMinimumSize(self.ColLabel.sizeHint())
     self.ColLabel.setFixedHeight(self.ColLabel.height())
     self.RowLabel.show()
     Info.addWidget(self.ColLabel)

     self.xCl=QLabel("x:", self, "xcl")
     self.xCl.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
     self.xCl.setFont(self.Font)
     self.xCl.setMinimumSize(self.xCl.sizeHint())
     self.xCl.setFixedHeight(self.xCl.height())
     self.xCl.show()
     xBox.addWidget(self.xCl)
     self.xCombo = QComboBox( 0, self, "xComboBox" )
     self.xCombo.setFixedHeight( self.xCombo.sizeHint().height() )
     self.xCombo.setFont(self.Font)
     self.xCombo.show()
     xBox.addWidget(self.xCombo)

     self.yCl=QLabel("y:", self, "ycl")
     self.yCl.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
     self.yCl.setFont(self.Font)
     self.yCl.setMinimumSize(self.yCl.sizeHint())
     self.yCl.setFixedHeight(self.yCl.height())
     self.yCl.show()
     yBox.addWidget(self.yCl)
     self.yCombo = QComboBox( 0, self, "yComboBox" )
     self.yCombo.setFixedHeight( self.yCombo.sizeHint().height() )
     self.yCombo.setFont(self.Font)
     self.yCombo.show()
     yBox.addWidget(self.yCombo)
     topLayout.activate()
     self.show()
# -------------------------------
  def SetValues(self,td):
     self.td=td
     self.RowLabel.setText("Rows: %d" % len(td.Nums[0]))
     self.ColLabel.setText("Cols: %d" % len(td.Nums))
     self.x=td.iCx
     self.y=td.iCy

     ix=self.xCombo.count()
     iy=self.yCombo.count()
#     print ix,iy
     for i in range(len(td.ColId)):
        t="%1d: %s" % (i+1, td.ColId[i])
        if i>=ix: self.xCombo.insertItem(t)
	else: self.xCombo.changeItem(t,i)
        if i>=iy: self.yCombo.insertItem(t)
	else: self.yCombo.changeItem(t,i)

     self.connect( self.xCombo, SIGNAL("activated(int)"),self.xColChanged)
     self.connect( self.yCombo, SIGNAL("activated(int)"),self.yColChanged)
#     print td.iCx,td.iCy
     self.xCombo.setCurrentItem(td.iCx)
     self.yCombo.setCurrentItem(td.iCy)
# -------------------------------
  def xColChanged(self,x):
     self.td.SetXYCol(x,self.y)
     self.x=x
#     print "x",x
     self.emit(PYSIGNAL("NewXCol(int)"), (x, ))
# -------------------------------
  def yColChanged(self,y):
     self.td.SetXYCol(self.x,y)
     self.y=y
#     print "y",y
     self.emit(PYSIGNAL("NewYCol(int)"), (y, ))
# -------------------------------
# DrawWid
# -------------------------------
class DataWid(QFrame):
  def __init__(self,parent=None,name=None):
     QFrame.__init__(self,parent,name)
     self.iColX = -1
     self.iColY = -1
     self.FontSize = self.height()/24
     self.Table=None
     self.Header=None
     self.TD = [None]
# -------------------------------
  def InitPlot(self,df):
     self.TD[0]=df
     if self.TD[0].nNLines == 0:
       raise ValueError("Data File empty")

     self.iColX=self.TD[0].iCx
     self.iColY=self.TD[0].iCy

#   Change here if font size of Show.py does not fit
#   Fontsize in parts of the window height
     self.FontSize = self.height()/24
     
     self.PPar={   "font": QFont("Helvetica",10),
                   "rect": QRect(),
                "x_label":self.TD[0].ColId[self.iColX],
                "y_label":self.TD[0].ColId[self.iColY],
	           "flag":TOP_FILENAME,
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

#     self.PPar["rect"].setRect(3,3,self.width()-3,self.height()-3)
#     pprint(self.PPar)
     self.setFrameStyle( QFrame.Panel | QFrame.Sunken)
     try: self.setPaletteBackgroundColor(Qt.white)
     except AttributeError: self.setBackgroundColor(Qt.white)
     self.show()
#     print self.TD[0].iCx,self.TD[0].iCy
# -------------------------------
  def keyPressEvent(self,e):
      if e.key() == Qt.Key_C:
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
# -------------------------------
  def paintEvent(self,p):
      paint = QPainter()
      if not paint.begin(self):
         print "Painter: begin failed"
      if self.TD[0]==None:
         paint.setFont(QFont("sanserif",8))
         paint.drawText(3,self.height()/2,"Click button to load file")
      else:
         self.PPar["font"]=QFont("sansserif",self.FontSize);
#         print "Fontsize:",iF
         self.PPar["rect"].setRect(3,3,self.width()-3,self.height()-3)
         DrawData(paint,self.TD,self.PPar)
      paint.end()
# -------------------------------
  def mouseDoubleClickEvent(self,e):
      #print "Double clicked",str(e)
      if e.button() == Qt.LeftButton:
         if self.Table == None or self.Table.isHidden(): self.Table=DataTable(self.TD[0])
         else: 	 QMessageBox.about(self,"Pat","Table already open") 
      if e.button() == Qt.RightButton:
         if self.Header == None or self.Header.isHidden(): self.Table=DataHeader(self.TD[0])
         else: 	 QMessageBox.about(self,"Pat","Table already open") 

# -------------------------------
# DataTable
# -------------------------------
class DataTable(QWidget):
  def __init__(self,td,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     if td == None:
	 QMessageBox.about(self,"Pat","No data to display.\nOpen first.") 
	 return
     self.nCols=len(td.Nums)
     self.nRows=td.nNLines
     self.Font=QFont("mono",10)
#     print self.nCols,self.nRows,td.Nums[0]
     self.table = QTable(self.nRows,self.nCols)
     self.table.setCaption("Data Columns")
     self.table.setFont(self.Font)
     self.table.setReadOnly(True)
     self.header = self.table.horizontalHeader()
     for i in range(len(td.ColId)):
        self.header.setLabel(i, "%1d: %s" % (i+1, td.ColId[i]))
     
     for i in range(self.nRows):
         for j in range(self.nCols):
	     self.table.setText(i,j,td.Nums[j][i])
     self.table.show()  

# -------------------------------
# Dataheader
# -------------------------------
class DataHeader(QWidget):
  def __init__(self,td,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     if td == None:
	 QMessageBox.about(self,"Pat","No data to display.\nOpen first.") 
	 return
     self.Font=QFont("mono",10)
     self.head=QTextEdit()
     self.head.setFont(self.Font)
     self.head.setReadOnly(True)
     self.head.setWordWrap(QTextEdit.NoWrap)
     self.head.setText(join(td.FText,"\n"))
     self.head.setCaption("Header lines")
     self.head.setGeometry(30,30,SCREEN.x()/3,len(td.FText)*10)
     self.head.show()  

# -------------------------------
# PatWid
# -------------------------------
class PatWid(QWidget):
  def __init__(self,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     self.DFrames =[{"label":"First data file", "tip":"First data file, double click: right: show data; left: show text ", "name":"DFile1"},
                    {"label":"Second data file","tip":"Second data file","name":"DFile2"}
#                    {"label":"Result",          "tip":"Result",          "name":"Result"}
                  ]
     self.DBts =  [[{"label":"Open 1st","name":"BFile1"},{"label":"Save","name":"BSFile1"}],
                   [{"label":"Open 2nd","name":"BFile2"}]
#                   {"label":"Save result",       "name":"BFile3"}
                  ]
# type: text, int, list, box
# max: nlen#: number of data points, #: 0,1: number of file
#      tlen#: number of text lines
#      clen#: number of columns
#xy:   0: -x, -y not passed to command
#      1: -x, -y are passed to command
#      2: -x, -y of two files are passed to command
     self.Cmds =  [ {"cmd":"Calc.py",   "tip": "Perform a calculation with one column", "xy":0, "files":1, 
                           "options":{"-f": {"type":["text"], "name":["Formula"], "min":[None], "max":[None]}}          },
   	   	    {"cmd":"Deriv.py",  "tip": "Calculate the numeric derivative", "xy":1, "files":1,
                           "options":{"-n": {"type":["int"], "name":["Neighbour points"], "min":[1], "max":["nlen0"]}}          },
   	   	    {"cmd":"Exchg.py",  "tip": "Swap columns", "xy":0, "files":1,
                            "options":{"-c": {"type":["int","int"], "name":["Column","Column"], "min":[1,1], "max":["clen0","clen0"],}}  },
  	   	    {"cmd":"Fileop.py", "tip": "Add, subtract, multiply or divide two data files", "xy":2, "files":2,
                            "options":{"-f": {"type":["list"], "name":["Operation"], "min":[],"max":["+,-,'*',/"]}}  },
   	   	    {"cmd":"Integ.py",  "tip": "Numerical integration", "xy":1, "files":1,
                            "options":{}  },
#   	   	    {"cmd":"Intpol.py", "tip": "Interpolate one value from a table"},
   	   	    {"cmd":"Look.py",   "tip": "Interpolate all values of data file from a table file", "xy":2, "files":2,
                            "options":{"-a": {"type":["box"], "name":["Column"], "min":["append"],"max":[]}}  },
   	   	    {"cmd":"Mean.py",   "tip": "Calculate the mean value in bins with size of delta x.", "xy":1, "files":1,
                            "options":{"-f": {"type":["text"], "name":["Delta X"], "min":[None],"max":[None]}}  },
   	   	    {"cmd":"Merge.py",  "tip": "Merge two data files", "xy":0, "files":2,
                            "options":{"-c": {"type":["box"], "name":["Merge by"], "min":["column"],"max":[]},
                                       "-r": {"type":["box"], "name":["Merge by"], "min":["row"],"max":[]}}  },
   	   	    {"cmd":"Select.py", "tip": "Obtain a part from a data file", "xy":0, "files":1,
                            "options":{# "-c": {"type":["int","int"], "name":["From column","To column"], "min":[1,1], "max":["clen0","clen0"]},
			               #"-r": {"type":["int","int","int"], "name":["Start row","End row","Row step"], "min":[1,1,1], "max":["nlen0","nlen0","nlen0"]},
			               "-f": {"type":["int","text","text"], "name":["Column","Start value","End value"], "min":[1,None,None], "max":["clen0",None,None]}}  },
   	   	    {"cmd":"Sort.py",   "tip": "Sort columns", "xy":0, "files":1,
                            "options":{"-c": {"type":["int"], "name":["Sort column"], "min":[1], "max":["clen0"]},
                                       "-r": {"type":["box"], "name":["Sort"], "min":["reverse"],"max":[]}}  },
#   	   	    {"cmd":"Stext.py",  "tip": "Change one text line"},
   	   	    {"cmd":"Show.py",   "tip": "Display a graphic of one data file on the screen", "xy":1, "files":1,
                           "options":{}  },
#   	   	    {"cmd":"Plot.py",   "tip": "Display and print up to 8 data files in one plot"},
   	   	    {"cmd":"Disp.py",   "tip": "Detailed graphic display with the possibility to change data", "xy":1, "files":1,
                           "options":{}  }
   	   	 ]
     self.DF=[None,None]

     self.w_show=SCREEN.x()/4
     self.h_show=SCREEN.y()/4
     self.w=self.w_show+10
     self.h=3*(self.h_show+6)
     self.FontSize = 10
     self.WMain=QRect(50,50,self.w,self.h)
 
     self.setCaption(str("Pat"))
     self.setGeometry(self.WMain)
     self.CmdW=None
     
     #Command box, column info, show frames
     widLayout = QVBoxLayout( self, 2 )
     topLayout = QHBoxLayout( 3 )
     
     widLayout.addLayout(topLayout)

     self.CmdBox = QListBox( self, "cmdBox" )
     for i in self.Cmds:
        self.CmdBox.insertItem(i["cmd"])
     topLayout.addWidget(self.CmdBox)
     self.connect( self.CmdBox, SIGNAL("selected(int)"), self.CmdSelected )
     self.connect( self.CmdBox, SIGNAL("highlighted(int)"), self.CmdHighlighted )
 
	
     colInfo = QVBoxLayout(2)
     topLayout.addLayout(colInfo)
     self.Bts=[]
     self.Cols=[]
     ct=0
     for i in self.DBts:
        C=ColInfoWid(self,"C%1d" % ct)
        self.Cols.append(C)
        colInfo.addWidget(C)
	buttons=QHBoxLayout(2)
	colInfo.addLayout(buttons)
        for j in i:
 	   B=QPushButton(self,j['name'])
           B.setText(j["label"])
           buttons.addWidget(B)
           self.Bts.append(B)
	   ct+=1
     self.connect(self.Bts[0],SIGNAL("clicked()"),self.ReadFirst)
     self.connect(self.Bts[1],SIGNAL("clicked()"),self.SaveFirst)
     self.connect(self.Bts[2],SIGNAL("clicked()"),self.ReadSecond)

     self.ShowFrames = QVBoxLayout(3)
     topLayout.addLayout(self.ShowFrames)

     self.DataFrames=[]
     for i in self.DFrames:
	 F=DataWid(self,i['name'])
         F.setFixedHeight(self.h_show)
         F.setFixedWidth(self.w_show)
         F.setFrameStyle( QFrame.Panel | QFrame.Sunken )
         self.ShowFrames.addWidget(F)
         QToolTip.add(F,i["tip"])
         F.setFocusPolicy(QWidget.ClickFocus )
         self.DataFrames.append(F)

     #Tip Label
     self.TipL=QLabel("", self, "TipLabel")
     self.TipL.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.TipL.setFrameStyle( QFrame.Panel | QFrame.Sunken )
     self.TipL.setFont(QFont("sanserif",self.FontSize))
     self.TipL.setMinimumSize(self.TipL.sizeHint())
     self.TipL.setFixedHeight(self.TipL.height())
     widLayout.addWidget(self.TipL)
     
     self.setUpdatesEnabled(1)
     widLayout.activate()
     self.connect(self.Cols[0],PYSIGNAL("NewXCol(int)"),self.NewColFirst)
     self.connect(self.Cols[0],PYSIGNAL("NewYCol(int)"),self.NewColFirst)
     self.show()
# -------------------------------------------
  def ReadFirst(self):
     self.Read(0)
# -------------------------------------------
  def ReadSecond(self):
     self.Read(1)
# -------------------------------------------
  def Read(self,n):
      fn = QFileDialog.getOpenFileName( QString.null, QString.null, self )
      if not fn.isEmpty():
         Filename=str(fn)
         fp=open(Filename,"r")
         Type=CheckFileType(fp)
         if Type[0] == None:
	    QMessageBox.about(self,"Pat","Not a valid data file")
	    return
         if Type[0] == FT_THECAP and MOD_ftypes['the']: self.DF[n]=TheFile(Filename)
         elif Type[0] == FT_SXSMulti and MOD_ftypes['sxs']: self.DF[n]=SxSFile(Filename)
         else: self.DF[n]=AsciiFile(Filename)
	 try: self.DF[n].Read(1,fp)
	 except:
	    QMessageBox.about(self,"Pat","%s\nNot a valid data file" % Filename) 
	    return
	 fp.close()
	 self.DataFrames[n].InitPlot(self.DF[n])
	 self.Cols[n].SetValues(self.DF[n])
	 self.DataFrames[n].repaint()
#	 self.repaint()
# -------------------------------------------
  def SaveFirst(self):
      if self.DF[0] == None:
	 QMessageBox.about(self,"Pat","No data to save.\n") 
	 return
         
      fn = QFileDialog.getSaveFileName( QString.null, QString.null, self )
      if not fn.isEmpty():
         Filename=str(fn)
	 try: fp=open(Filename,"r")
	 except IOError, (errno, strerror):
#	    print errno
	    if errno==2: #File does not exist; try to save
              self.SaveFile(Filename)
	      self.DataFrames[0].repaint()
	      return
	    else:
	      QMessageBox.about(self,"Pat","%s\nCannot save file\n%d %s" % (Filename,errno, strerror)) 
	      return
	 else:
	    a=QMessageBox.question(self,"Pat","%s\n already exists!\nOverwrite ?" % (Filename),QMessageBox.Yes,QMessageBox.No)
            fp.close()
            if a==QMessageBox.Yes:
              self.SaveFile(Filename)
	      self.DataFrames[0].repaint()
	      return
	    else: print "File %s not saved" % Filename; return
	 print "Hmmm: Save File"
#	 self.repaint()
# -------------------------------------------
  def SaveFile(self,Filename):
     try: fo=open(Filename,"w")
     except IOError, (errno, strerror):
     	QMessageBox.about(self,"Pat","%s\nCannot save file\n%d %s" % (Filename,errno, strerror)) 
     	return 0
     self.DF[0].Write(fo)
     self.DF[0].FileName=Filename
     fo.close()
     return 1
# -------------------------------------------
  def NewColFirst(self,c):
     self.DataFrames[0].InitPlot(self.DF[0])
     self.DataFrames[0].repaint()
# -------------------------------------------
  def CmdSelected(self,cmd):
#     print cmd,self.Cmds[cmd]
#     for k,v in self.Cmds[cmd]["options"].iteritems():
#         print k,v
     try: self.Cmds[cmd]["options"]
     except:
       QMessageBox.information(self,"Pat-Warning","%s\nNot implemented" % self.Cmds[cmd]["cmd"])
       return
     if self.DF[0]==None:
       QMessageBox.information(self,"Pat-Error","No 1st file opened")
       return

     if self.Cmds[cmd]["files"]==2 and self.DF[1]==None:
       QMessageBox.information(self,"Pat-Error","No 2nd file opened")
       return
    
     if self.CmdW == None or self.CmdW.isHidden():
        self.CmdW = CmdWid(self.Cmds[cmd], self.DF, None, "actCmd" )
	self.CmdW.setWFlags(Qt.WType_TopLevel)
        self.connect(self.CmdW,PYSIGNAL("NewFile(char *)"),self.NewResult)
	self.CmdW.show()
     else:
        QMessageBox.about(self,"Pat","Another command window is already open\nClose it first")
# -------------------------------------------
  def CmdHighlighted(self,cmd):
      self.TipL.setText("%s: %s" % (self.Cmds[cmd]["cmd"],self.Cmds[cmd]["tip"]))
# -------------------------------------------
  def NewResult(self,Result):
      Filename=Result[0]
#      print Result[1]
#      print "Result:",Filename
      self.TipL.setText("%s finished " % Result[1])
      if Filename == None:
         self.CmdW.close()
	 return

      x=self.DF[0].iCx
      y=self.DF[0].iCy
      fp=open(Filename,"r")
      try: Type=CheckFileType(fp)
      except: Type=(None,None)
      if Type[0] == None:
         QMessageBox.about(self,"Pat","Not a valid data file")
         return
      if Type[0] == FT_THECAP and MOD_ftypes['the']: self.DF[0]=TheFile(Filename)
      elif Type[0] == FT_SXSMulti and MOD_ftypes['sxs']: self.DF[0]=SxSFile(Filename)
      else: self.DF[0]=AsciiFile(Filename)
      self.DF[0].Read(1,fp)
      self.DF[0].SetXYCol(x,y)
      self.DataFrames[0].InitPlot(self.DF[0])
      self.Cols[0].SetValues(self.DF[0])
      self.DataFrames[0].repaint()
#      os.remove(Filename)
#      print self.DF[0].iCx, self.DF[0].iCy
      self.CmdW.close()
# -------------------------------
# CmdWid
# -------------------------------
class CmdWid(QWidget):
  def __init__(self,cmd,df_list,parent=None,name=None):
     QWidget.__init__(self,parent,name)
#     print cmd
     self.Cmd=cmd
     self.DF=df_list
     self.setCaption(self.Cmd["cmd"])
     self.TmpFile=None

     #Label, command parameters, go & help button, (help)
     self.topLayout = QVBoxLayout( self, 4 )
     self.Font=QFont("sanserif",8)

     #Label
     self.Info=QLabel("", self, "InfoLabel")
     self.Info.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
     self.Info.setFrameStyle( QFrame.Panel | QFrame.Sunken )
     self.Info.setFont(self.Font)
     self.Info.setText(cmd["tip"])
     self.Info.setMinimumSize(self.Info.sizeHint())
     self.Info.setFixedHeight(self.Info.height())
     self.Info.show()
     self.topLayout.addWidget(self.Info)

     cmd_pars = QVBoxLayout()
     self.topLayout.addLayout(cmd_pars)

     self.Opts={}
     for option,v in self.Cmd["options"].iteritems():		    
#     	print option,":",v					    
     	param = QHBoxLayout(2)  				    
     	W = None
        # type, name, min, max: are lists ---> -f 3,2.5,4.5
        self.Opts[option]=[]
	for i in range(len(v["name"])):						    
     	   L=QLabel("", self, option+"_l") 			    
     	   L.setAlignment(Qt.AlignRight|Qt.AlignVCenter)		    
#    	    L.setFrameStyle( QFrame.Panel | QFrame.Sunken )	    
     	   L.setFont(self.Font)					    
     	   L.setText(v["name"][i])					    
     	   L.setMinimumSize(L.sizeHint())  			    
     	   L.setFixedHeight(L.height())				    
     	   param.addWidget(L)					    
     	   if v["type"][i] == "box":
	      W = QCheckBox( self, name )
	      W.setFont(self.Font)
	      W.setText(v["min"][i])
	      rfunc=W.isChecked
	      otype="novalue"
     	   if v["type"][i] == "list":
	      W = QComboBox( self, name )
	      rfunc=W.currentText
	      otype="value"
	      for t in v["max"][i].split(","):W.insertItem(t)             
     	   if v["type"][i] == "text":
	      W = QLineEdit( self, name )	    
              if v["min"][i]:
	         W.setMaxLength(v["min"][i])
     	         W.setFixedWidth(v["min"][i]*30)				    
              rfunc=W.text
	      otype="value"
     	   if v["type"][i] == "int":
	      if v["max"][i][:-1] == "nlen":
	         fno=atoi(v["max"][i][-1])
	         ma=len(self.DF[fno].Nums[0])
	      elif v["max"][i][:-1] == "clen":
	         fno=atoi(v["max"][i][-1])
	         ma=len(self.DF[fno].Nums)
	      else: ma=v["max"][i]
              W = QSpinBox(v["min"][i], ma, 1, self, name )
              rfunc=W.text
	      otype="value"
	   #self.Opts[option]=W
     	   self.Opts[option].append({"obj":W, "func":rfunc, "type": otype})					    
     	   param.addWidget(W)					    
     	cmd_pars.addLayout(param)				    
    							    
     #Buttons
     buttons = QHBoxLayout(2)
     self.topLayout.addLayout(buttons)
     
     if self.DF[0]!=None: 
       self.Go=QPushButton(self,"goButton")
       self.Go.setText("Go")
       self.Go.setMinimumSize(self.Go.sizeHint())
       self.Go.setFixedHeight(self.Go.height())
       buttons.addWidget(self.Go)
       self.connect(self.Go,SIGNAL("clicked()"),self.RunCmd)
     else:
       L=QLabel("", self, "noFile")			    
       L.setAlignment(Qt.AlignCenter|Qt.AlignVCenter)		 
       L.setFrameStyle( QFrame.Panel | QFrame.Sunken ) 	 
       L.setFont(self.Font)					 
       L.setText("No file selected")					 
       L.setMinimumSize(L.sizeHint())				 
       L.setFixedHeight(L.height())				 
       buttons.addWidget(L)					 

     self.Help=QPushButton(self,"helpButton")
     self.Help.setText("Help")
     self.Help.setMinimumSize(self.Help.sizeHint())
     self.Help.setFixedHeight(self.Help.height())
     buttons.addWidget(self.Help)
     self.connect(self.Help,SIGNAL("clicked()"),self.DispHelp)
# -------------------------------
  def DispHelp(self):
      fp=popen3(self.Cmd["cmd"]+" -h")
      h_mesg=fp[0].read()
      e_mesg=fp[2].read()
#      if h_mesg != "": QMessageBox.information(self,"Pat-Help",h_mesg)
      Info=QLabel("", self, "InfoLabel")
      Info.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
      Info.setFrameStyle( QFrame.Panel | QFrame.Sunken )
      Info.setFont(QFont("mono",8))
      Info.setText(h_mesg)
      Info.setMinimumSize(Info.sizeHint())
      Info.setFixedHeight(Info.height())
      Info.show()
      self.topLayout.addWidget(Info)
      if e_mesg != "": QMessageBox.about(self,"Pat-Error",e_mesg)
# -------------------------------
  def RunCmd(self):
#      print self.Opts
      col=""
#      if self.Cmd["xy"]: col=" -x %1d -y %1d " % (self.DF[0].iCx+1, self.DF[0].iCy+1)
#      else: col=" "
      xopt=[]
      yopt=[]
      for i in range(self.Cmd["xy"]):
         xopt.append(str(self.DF[i].iCx+1))
         yopt.append(str(self.DF[i].iCy+1))
      
      if xopt!=[]:  col = " -x " +  join(xopt,",")
      if yopt!=[]: col += " -y " +  join(yopt,",")
      col += " "
      
      cmd = self.Cmd["cmd"]+col
      for option,v in self.Opts.iteritems():
          #print option,v
          ol=[]
          for i in v:
	     #print i["obj"],i["func"]
	     ol.append(str(i["func"]()).strip())
	  o=join(ol,",")
          if i["type"] == "novalue":
	     #print o
	     if o == "True": cmd += "%s " % option
             else: pass
	  else:
	     cmd += "%s" % option +  " %s " % o 
      cmd += ' "%s"' % self.DF[0].FileName
      if self.Cmd["files"] == 2: cmd += " %s" % self.DF[1].FileName
      fp=popen3(cmd)

      try: fh=open(HIST_FILE,"a")
      except: print "Hmmm: Cannot write history: %s", HIST_FILE
      else: fh.write(cmd); fh.close()
      #print cmd
      out = fp[0].read()
      err = fp[2].read()
      if self.TmpFile:
         try: os.remove(self.TmpFile[1])
	 except: print "Error removing tmp-file:%s" % self.TmpFile[1]
      if out !="":	 
         self.TmpFile = mkstemp(suffix="tmp",prefix=self.Cmd["cmd"].split(".")[0], text=True)
         #print self.TmpFile[1]
         os.write(self.TmpFile[0],out)
         #print out
      else: self.TmpFile=(None,None)
       	 
      if err != "":
         QMessageBox.about(self,"Pat-Error",err)
	 print err
	 return
      self.emit(PYSIGNAL("NewFile(char*)"), ((self.TmpFile[1],cmd), ))
