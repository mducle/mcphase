# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'patwid.ui'
#
# Created: Fri Jan 6 14:40:52 2006
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
#

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

#$Log$
#
CVS_ID="$Id$"

lib_path='PAT_LIB_PATH'
if os.name=='posix': PATH_DEL="/"
else: PATH_DEL="\\"
HIST_FILE=os.environ[lib_path]+PATH_DEL+".cmd_history"
ICON_PATH=os.environ[lib_path]+PATH_DEL+"icons/"

# -----------------------------------------------------
# PatWid
# -----------------------------------------------------
class PatWid(QWidget):
  def __init__(self,parent = None,name = None,fl = 0):
     QWidget.__init__(self,parent,name,fl)
     if not name:
     	 self.setName("PatWid")
     # ------------------------------------------
     self.DF=[None,None]
     self.Table=[None,None]
     self.Header=[None,None]
     self.FontSize = 10
     self.CmdW=None

     # type: text, int, list, box
     # max: nlen#: number of data points, #: 0,1: number of file
     #      tlen#: number of text lines
     #      clen#: number of columns
     #xy:   0: -x, -y not passed to command
     #      1: -x, -y are passed to command
     #      2: -x, -y of two files are passed to command
     if os.name=='posix': file_ops=["+,-,'*',/"]
     else: file_ops=["+,-,*,/"]
     
     self.Cmds =  [ {"cmd":"Calc.py",   "tip": "Perform a calculation with one column", "xy":0, "files":1, 
                           "options":{"-f": {"type":["text"], "name":["Formula"], "min":[None], "max":[None]}}          },
   	   	    {"cmd":"Deriv.py",  "tip": "Calculate the numeric derivative", "xy":1, "files":1,
                           "options":{"-n": {"type":["int"], "name":["Neighbour points"], "min":[1], "max":["nlen0"]}}          },
   	   	    {"cmd":"Exchg.py",  "tip": "Swap columns", "xy":0, "files":1,
                            "options":{"-c": {"type":["int","int"], "name":["Column","Column"], "min":[1,1], "max":["clen0","clen0"],}}  },
  	   	    {"cmd":"Fileop.py", "tip": "Add, subtract, multiply or divide two data files", "xy":2, "files":2,
                            "options":{"-f": {"type":["list"], "name":["Operation"], "min":[],"max":file_ops}}  },
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
     # ------------------------------------------
     self.setCaption(str("Pat"))

     self.imageSave = QPixmap(ICON_PATH+"filesave.png")
     self.imageOpen = QPixmap(ICON_PATH+"fileopen.png")
     self.imageSpread = QPixmap(ICON_PATH+"spread.png")
     self.imageHead = QPixmap(ICON_PATH+"head.png")

     LayoutWidget = QWidget(self,"WidLayout")
     LayoutWidget.setGeometry(QRect(12,0,644,534))
     WidLayout = QVBoxLayout(LayoutWidget,11,6,"WidLayout")

     TopLayout = QHBoxLayout(None,0,6,"TopLayout")

     Cmd = QVBoxLayout(None,0,6,"Cmd")

     self.CmdBox = QListBox(LayoutWidget,"CmdBox")
     try: self.CmdBox.setSizePolicy(QSizePolicy(QSizePolicy.Minimum,QSizePolicy.Expanding,0,0,self.CmdBox.sizePolicy().hasHeightForWidth()))
     except: pass
     self.CmdBox.setMinimumSize(QSize(110,320))
     CmdBox_font = QFont(self.CmdBox.font())
     CmdBox_font.setPointSize(11)
     self.CmdBox.setFont(CmdBox_font)
     Cmd.addWidget(self.CmdBox)

     spacer1 = QSpacerItem(110,30,QSizePolicy.Minimum,QSizePolicy.Expanding)
     Cmd.addItem(spacer1)
     TopLayout.addLayout(Cmd)

     View = QVBoxLayout(None,0,6,"View")

     Data1 = QHBoxLayout(None,0,6,"Data1")

     Handling1 = QVBoxLayout(None,0,6,"Handling1")

     InfoLayout1 = QVBoxLayout(None,0,6,"InfoLayout1")

     Info = QHBoxLayout(None,0,6,"Info")

     self.RowLabel = [None,None]
     self.RowLabel[0] = QLabel(LayoutWidget,"RowLabel")
     RowLabel_font = QFont(self.RowLabel[0].font())
     RowLabel_font.setPointSize(12)
     self.RowLabel[0].setFont(RowLabel_font)
     self.RowLabel[0].setFrameShape(QLabel.Panel)
     self.RowLabel[0].setFrameShadow(QLabel.Sunken)
     self.RowLabel[0].setAlignment(QLabel.AlignVCenter | QLabel.AlignLeft)
     Info.addWidget(self.RowLabel[0])

     self.ColLabel = [None,None]
     self.ColLabel[0] = QLabel(LayoutWidget,"ColLabel")
     ColLabel_font = QFont(self.ColLabel[0].font())
     ColLabel_font.setPointSize(12)
     self.ColLabel[0].setFont(ColLabel_font)
     self.ColLabel[0].setFrameShape(QLabel.Panel)
     self.ColLabel[0].setFrameShadow(QLabel.Sunken)
     self.ColLabel[0].setAlignment(QLabel.AlignVCenter | QLabel.AlignLeft)
     Info.addWidget(self.ColLabel[0])
     InfoLayout1.addLayout(Info)

     xInfo = QHBoxLayout(None,0,6,"xInfo")

     self.xCol = [None,None]
     self.xCol[0] = QLabel(LayoutWidget,"xCol")
     try: self.xCol[0].setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.MinimumExpanding,0,0,self.xCol[0].sizePolicy().hasHeightForWidth()))
     except: pass
     self.xCol[0].setMaximumSize(QSize(30,20))
     xCol_font = QFont(self.xCol[0].font())
     xCol_font.setPointSize(12)
     self.xCol[0].setFont(xCol_font)
     self.xCol[0].setFrameShadow(QLabel.Sunken)
     self.xCol[0].setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
     xInfo.addWidget(self.xCol[0])

     self.xComboBox = [None,None]
     self.xComboBox[0] = QComboBox(0,LayoutWidget,"xComboBox")
     xComboBox_font = QFont(self.xComboBox[0].font())
     xComboBox_font.setPointSize(12)
     self.xComboBox[0].setFont(xComboBox_font)
     self.connect( self.xComboBox[0], SIGNAL("activated(int)"),self.x0ColChanged)
     xInfo.addWidget(self.xComboBox[0])
     InfoLayout1.addLayout(xInfo)

     yInfo = QHBoxLayout(None,0,6,"yInfo")

     self.yCol = [None,None]
     self.yCol[0] = QLabel(LayoutWidget,"yCol")
     try:self.yCol[0].setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.MinimumExpanding,0,0,self.yCol[0].sizePolicy().hasHeightForWidth()))
     except: pass
     self.yCol[0].setMaximumSize(QSize(30,20))
     yCol_font = QFont(self.yCol[0].font())
     yCol_font.setPointSize(12)
     self.yCol[0].setFont(yCol_font)
     self.yCol[0].setFrameShadow(QLabel.Sunken)
     self.yCol[0].setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
     yInfo.addWidget(self.yCol[0])

     self.yComboBox = [None,None]
     self.yComboBox[0] = QComboBox(0,LayoutWidget,"yComboBox")
     yComboBox_font = QFont(self.yComboBox[0].font())
     yComboBox_font.setPointSize(12)
     self.yComboBox[0].setFont(yComboBox_font)
     self.connect( self.yComboBox[0], SIGNAL("activated(int)"),self.y0ColChanged)
     yInfo.addWidget(self.yComboBox[0])
     InfoLayout1.addLayout(yInfo)
     Handling1.addLayout(InfoLayout1)

     self.Actions = [None,None]
     self.Actions[0] = QButtonGroup(LayoutWidget,"Actions")
     try: self.Actions[0].setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed,0,0,self.Actions[0].sizePolicy().hasHeightForWidth()))
     except: pass
     self.Actions[0].setMinimumSize(QSize(210,0))

     self.SaveFile = QToolButton(self.Actions[0],"SaveFile")
     self.SaveFile.setGeometry(QRect(60,20,40,40))
     self.SaveFile.setIconSet(QIconSet(self.imageSave))
     self.SaveFile.setUsesBigPixmap(0)
     self.Actions[0].insert( self.SaveFile,1)

     self.OpenFile = QToolButton(self.Actions[0],"OpenFile")
     self.OpenFile.setGeometry(QRect(10,20,40,40))
     self.OpenFile.setIconSet(QIconSet(self.imageOpen))
     self.OpenFile.setUsesBigPixmap(0)
     self.Actions[0].insert( self.OpenFile,0)
     self.connect(self.Actions[0],SIGNAL("clicked(int)"),self.Action1)

     self.ViewCols = QToolButton(self.Actions[0],"ViewCols")
     self.ViewCols.setGeometry(QRect(110,20,40,40))
     self.ViewCols.setIconSet(QIconSet(self.imageSpread))
     self.ViewCols.setUsesBigPixmap(0)
     self.Actions[0].insert( self.ViewCols,2)

     self.ShowHeader = QToolButton(self.Actions[0],"ShowHeader")
     self.ShowHeader.setGeometry(QRect(160,20,40,40))
     self.ShowHeader.setIconSet(QIconSet(self.imageHead))
     self.ShowHeader.setUsesBigPixmap(0)
     Handling1.addWidget(self.Actions[0])
     Data1.addLayout(Handling1)
     self.Actions[0].find(1).setEnabled(0)
     self.Actions[0].find(2).setEnabled(0)
     self.Actions[0].find(3).setEnabled(0)
     
     self.DrawWid = [None,None]
     self.DrawWid[0] = DataWid(LayoutWidget,"DrawWid1")
     try: self.DrawWid[0].setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding,0,0,self.DrawWid[0].sizePolicy().hasHeightForWidth()))
     except: pass
     self.DrawWid[0].setMinimumSize(QSize(300,233))
     self.DrawWid[0].setFrameShape(QFrame.StyledPanel)
     self.DrawWid[0].setFrameShadow(QFrame.Sunken)
     #self.DrawWid[0].setPaletteBackgroundColor(QColor(255,255,255))
     self.DrawWid[0].setFocusPolicy(QFrame.ClickFocus)
     Data1.addWidget(self.DrawWid[0])
     View.addLayout(Data1)

     Data2 = QHBoxLayout(None,0,6,"Data2")

     Handling2 = QVBoxLayout(None,0,6,"Handling2")

     InfoLayout2 = QVBoxLayout(None,0,6,"InfoLayout2")

     Info_2 = QHBoxLayout(None,0,6,"Info_2")

     self.RowLabel[1] = QLabel(LayoutWidget,"RowLabel_2")
     RowLabel_2_font = QFont(self.RowLabel[1].font())
     RowLabel_2_font.setPointSize(12)
     self.RowLabel[1].setFont(RowLabel_2_font)
     self.RowLabel[1].setFrameShape(QLabel.Panel)
     self.RowLabel[1].setFrameShadow(QLabel.Sunken)
     self.RowLabel[1].setAlignment(QLabel.AlignVCenter | QLabel.AlignLeft)
     Info_2.addWidget(self.RowLabel[1])

     self.ColLabel[1] = QLabel(LayoutWidget,"ColLabel_2")
     ColLabel_2_font = QFont(self.ColLabel[1].font())
     ColLabel_2_font.setPointSize(12)
     self.ColLabel[1].setFont(ColLabel_2_font)
     self.ColLabel[1].setFrameShape(QLabel.Panel)
     self.ColLabel[1].setFrameShadow(QLabel.Sunken)
     self.ColLabel[1].setAlignment(QLabel.AlignVCenter | QLabel.AlignLeft)
     Info_2.addWidget(self.ColLabel[1])
     InfoLayout2.addLayout(Info_2)

     xInfo_2 = QHBoxLayout(None,0,6,"xInfo_2")

     self.xCol[1] = QLabel(LayoutWidget,"xCol_2")
     try: self.xCol[1].setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.MinimumExpanding,0,0,self.xCol[1].sizePolicy().hasHeightForWidth()))
     except: pass
     self.xCol[1].setMaximumSize(QSize(30,20))
     xCol_2_font = QFont(self.xCol[1].font())
     xCol_2_font.setPointSize(12)
     self.xCol[1].setFont(xCol_2_font)
     self.xCol[1].setFrameShadow(QLabel.Sunken)
     self.xCol[1].setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
     xInfo_2.addWidget(self.xCol[1])

     self.xComboBox[1] = QComboBox(0,LayoutWidget,"xComboBox_2")
     xComboBox_2_font = QFont(self.xComboBox[1].font())
     xComboBox_2_font.setPointSize(12)
     self.xComboBox[1].setFont(xComboBox_2_font)
     self.connect( self.xComboBox[1], SIGNAL("activated(int)"),self.x1ColChanged)
     xInfo_2.addWidget(self.xComboBox[1])
     InfoLayout2.addLayout(xInfo_2)

     yInfo_2 = QHBoxLayout(None,0,6,"yInfo_2")

     self.yCol[1] = QLabel(LayoutWidget,"yCol_2")
     try: self.yCol[1].setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.MinimumExpanding,0,0,self.yCol[1].sizePolicy().hasHeightForWidth()))
     except: pass
     self.yCol[1].setMaximumSize(QSize(30,20))
     yCol_2_font = QFont(self.yCol[1].font())
     yCol_2_font.setPointSize(12)
     self.yCol[1].setFont(yCol_2_font)
     self.yCol[1].setFrameShadow(QLabel.Sunken)
     self.yCol[1].setAlignment(QLabel.AlignVCenter | QLabel.AlignRight)
     yInfo_2.addWidget(self.yCol[1])

     self.yComboBox[1] = QComboBox(0,LayoutWidget,"yComboBox_2")
     yComboBox_2_font = QFont(self.yComboBox[1].font())
     yComboBox_2_font.setPointSize(12)
     self.yComboBox[1].setFont(yComboBox_2_font)
     self.connect( self.yComboBox[1], SIGNAL("activated(int)"),self.y1ColChanged)
     yInfo_2.addWidget(self.yComboBox[1])
     InfoLayout2.addLayout(yInfo_2)
     Handling2.addLayout(InfoLayout2)

     self.Actions[1] = QButtonGroup(LayoutWidget,"Actions_2")
     try: self.Actions[1].setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed,0,0,self.Actions[1].sizePolicy().hasHeightForWidth()))
     except: pass
     self.Actions[1].setMinimumSize(QSize(210,0))

     self.OpenFile_3 = QToolButton(self.Actions[1],"OpenFile_3")
     self.OpenFile_3.setGeometry(QRect(10,20,40,40))
     self.OpenFile_3.setIconSet(QIconSet(self.imageOpen))
     self.OpenFile_3.setUsesBigPixmap(0)
     self.Actions[1].insert( self.OpenFile_3,0)

     self.SaveFile_3 = QToolButton(self.Actions[1],"SaveFile_3")
     self.SaveFile_3.setGeometry(QRect(60,20,40,40))
     self.SaveFile_3.setIconSet(QIconSet(self.imageSave))
     self.SaveFile_3.setUsesBigPixmap(0)
     self.Actions[1].insert( self.SaveFile_3,1)

     self.ViewCols_3 = QToolButton(self.Actions[1],"ViewCols_3")
     self.ViewCols_3.setGeometry(QRect(110,20,40,40))
     self.ViewCols_3.setIconSet(QIconSet(self.imageSpread))
     self.ViewCols_3.setUsesBigPixmap(0)
     self.Actions[1].insert( self.ViewCols_3,2)

     self.ShowHeader_3 = QToolButton(self.Actions[1],"ShowHeader_3")
     self.ShowHeader_3.setGeometry(QRect(160,20,40,40))
     self.ShowHeader_3.setIconSet(QIconSet(self.imageHead))
     self.ShowHeader_3.setUsesBigPixmap(0)
     self.connect(self.Actions[1],SIGNAL("clicked(int)"),self.Action2)

     Handling2.addWidget(self.Actions[1])
     Data2.addLayout(Handling2)
     self.Actions[1].find(1).setEnabled(0)
     self.Actions[1].find(2).setEnabled(0)
     self.Actions[1].find(3).setEnabled(0)

     self.DrawWid[1] = DataWid(LayoutWidget,"DrawWid1")
     try: self.DrawWid[1].setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding,0,0,self.DrawWid[1].sizePolicy().hasHeightForWidth()))
     except: pass
     self.DrawWid[1].setMinimumSize(QSize(300,233))
     self.DrawWid[1].setFrameShape(QFrame.StyledPanel)
     self.DrawWid[1].setFrameShadow(QFrame.Sunken)
     #self.DrawWid[1].setPaletteBackgroundColor(QColor(255,255,255))
     self.DrawWid[1].setFocusPolicy(QFrame.ClickFocus)
     Data2.addWidget(self.DrawWid[1])
     View.addLayout(Data2)
     TopLayout.addLayout(View)
     WidLayout.addLayout(TopLayout)

     self.CmdInfo = QLabel(LayoutWidget,"CmdInfo")
     try: self.CmdInfo.setSizePolicy(QSizePolicy(QSizePolicy.Expanding,QSizePolicy.Fixed,0,0,self.CmdInfo.sizePolicy().hasHeightForWidth()))
     except: pass
     self.CmdInfo.setMinimumSize(QSize(0,20))
     self.CmdInfo.setMaximumSize(QSize(32767,20))
     CmdInfo_font = QFont(self.CmdInfo.font())
     CmdInfo_font.setPointSize(9)
     self.CmdInfo.setFont(CmdInfo_font)
     self.CmdInfo.setFrameShape(QLabel.Panel)
     self.CmdInfo.setFrameShadow(QLabel.Sunken)
     WidLayout.addWidget(self.CmdInfo)

     self.languageChange()
     for i in self.Cmds:
        #print i["cmd"]
        #self.CmdBox.insertItem(self.__tr("New Item"))
        self.CmdBox.insertItem(i["cmd"])
	
     self.connect( self.CmdBox, SIGNAL("selected(int)"), self.CmdSelected )
     self.connect( self.CmdBox, SIGNAL("highlighted(int)"), self.CmdHighlighted )

     self.resize(QSize(661,529).expandedTo(self.minimumSizeHint()))
     try: self.clearWState(Qt.WState_Polished)
     except: pass
# -----------------------------------------------------
  def Action1(self,button):
     #print "Pressed 1",button
     self.Action(0,button)
# -----------------------------------------------------
  def Action2(self,button):
     self.Action(1,button)
     #print "Pressed 2",button
# -----------------------------------------------------
  def Action(self,n,button):
     #print n,button
     if button == 0:  # OpenFile
        self.Read(n)
	return 
     elif button == 1:  # SaveFile
        self.GetFilename(n)
	return 
     elif button == 2:  # TableReq
        self.TableReq(n)
	return 
     elif button == 3:  # HeadReq
        self.HeaderReq(n)
	return 
# -----------------------------------------------------
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
	 #self.DrawWid[n].InitPlot(self.DF[n])
	 self.Redraw(n)
	 self.SetColValues(n)
         self.Actions[n].find(1).setEnabled(1)
         self.Actions[n].find(2).setEnabled(1)
	 self
         self.Actions[n].find(3).setEnabled(1)
	 #self.DataFrames[n].repaint()
         #self.repaint()
# -------------------------------------------
  def GetFilename(self,n):
      if self.DF[n] == None:
	 QMessageBox.about(self,"Pat","No data to save.\n") 
	 return
         
      fn = QFileDialog.getSaveFileName( QString.null, QString.null, self )
      if not fn.isEmpty():
         Filename=str(fn)
	 try: fp=open(Filename,"r")
	 except IOError, (errno, strerror):
#	    print errno
	    if errno==2: #File does not exist; try to save
              #print n,Filename
              self.SaveIt(n,Filename)
	      #self.DataFrames[0].repaint()
	      return
	    else:
	      QMessageBox.about(self,"Pat","%s\nCannot save file\n%d %s" % (Filename,errno, strerror)) 
	      return
	 else:
	    a=QMessageBox.question(self,"Pat","%s\n already exists!\nOverwrite ?" % (Filename),QMessageBox.Yes,QMessageBox.No)
            fp.close()
            if a==QMessageBox.Yes:
              self.SaveIt(n,Filename)
	      #self.DataFrames[0].repaint()
	      return
	    else: print "File %s not saved" % Filename; return
	 print "Hmmm: Save File"
#	 self.repaint()
# -------------------------------------------
  def SaveIt(self,n,Filename):
     try: fo=open(Filename,"w")
     except IOError, (errno, strerror):
     	QMessageBox.about(self,"Pat","%s\nCannot save file\n%d %s" % (Filename,errno, strerror)) 
     	return 0
     self.DF[n].Write(fo)
     self.DF[n].FileName=Filename
     fo.close()
     return 1
# -------------------------------------------
  def x0ColChanged(self,i):self.xColChanged(0,i)
# ------------------------------------------
  def x1ColChanged(self,i):self.xColChanged(1,i)
# ------------------------------------------
  def y0ColChanged(self,i):self.yColChanged(0,i)
# ------------------------------------------
  def y1ColChanged(self,i):self.yColChanged(1,i)
# -------------------------------------------
  def xColChanged(self,n,x):
     self.DF[n].SetXYCol(x,self.y)
     self.x=x
     #print "x",x
     self.Redraw(n)
# -------------------------------------------
  def yColChanged(self,n,y):
     self.DF[n].SetXYCol(self.x,y)
     self.y=y
     #print "y",y
     self.Redraw(n)
# -------------------------------------------
  def SetColValues(self,n):
     self.RowLabel[n].setText("Rows: %d" % len(self.DF[n].Nums[0]))
     self.ColLabel[n].setText("Cols: %d" % len(self.DF[n].Nums))
     self.x=self.DF[n].iCx
     self.y=self.DF[n].iCy

     ix=self.xComboBox[n].count()
     iy=self.yComboBox[n].count()
#     print ix,iy
     for i in range(len(self.DF[n].ColId)):
        t="%1d: %s" % (i+1, self.DF[n].ColId[i])
        if i>=ix: self.xComboBox[n].insertItem(t)
	else: self.xComboBox[n].changeItem(t,i)
        if i>=iy: self.yComboBox[n].insertItem(t)
	else: self.yComboBox[n].changeItem(t,i)

#     print td.iCx,td.iCy
     self.xComboBox[n].setCurrentItem(self.DF[n].iCx)
     self.yComboBox[n].setCurrentItem(self.DF[n].iCy)
# -------------------------------------------
  def Redraw(self,n):
     self.DrawWid[n].InitPlot(self.DF[n])
     self.DrawWid[n].repaint()
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
      self.CmdInfo.setText("%s: %s" % (self.Cmds[cmd]["cmd"],self.Cmds[cmd]["tip"]))
# -------------------------------------------
  def NewResult(self,Result):
      Filename=Result[0]
#      print Result[1]
#      print "Result:",Filename
      self.CmdInfo.setText("%s finished " % Result[1])
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
      self.DrawWid[0].InitPlot(self.DF[0])
      self.SetColValues(0)
      self.DrawWid[0].repaint()
#      os.remove(Filename)
#      print self.DF[0].iCx, self.DF[0].iCy
      self.CmdW.close()
# -------------------------------------------
  def TableReq(self,n):
     if self.Table[n] == None or self.Table[n].isHidden(): self.Table[n]=DataTable(self.DF[n])
     else:   QMessageBox.about(self,"Pat","Table already open") 
# -------------------------------------------
  def HeaderReq(self,n):
     if self.Header[n] == None or self.Header[n].isHidden(): self.Header[n]=DataHeader(self.DF[n])
     else:   QMessageBox.about(self,"Pat","Header already open") 
# -------------------------------------------
  def languageChange(self):
     self.setCaption(self.__tr("Pat"))
     self.CmdBox.clear()
#     self.CmdBox.insertItem(self.__tr("New Item"))
     self.RowLabel[0].setText(self.__tr("Rows: ----"))
     self.ColLabel[0].setText(self.__tr("Cols: --"))
     self.xCol[0].setText(self.__tr("x:"))
     QToolTip.add(self.xComboBox[0],self.__tr("Select x column"))
     self.yCol[0].setText(self.__tr("y:"))
     QToolTip.add(self.yComboBox[0],self.__tr("Select y column"))
     self.Actions[0].setTitle(self.__tr("File 1"))
     self.SaveFile.setText(QString.null)
     QToolTip.add(self.SaveFile,self.__tr("Save File"))
     self.OpenFile.setText(QString.null)
     QToolTip.add(self.OpenFile,self.__tr("Open File"))
     self.ViewCols.setText(QString.null)
     QToolTip.add(self.ViewCols,self.__tr("View Columns"))
     self.ShowHeader.setText(QString.null)
     QToolTip.add(self.ShowHeader,self.__tr("Show Header"))
     self.RowLabel[1].setText(self.__tr("Rows: ----"))
     self.ColLabel[1].setText(self.__tr("Cols: --"))
     self.xCol[1].setText(self.__tr("x:"))
     QToolTip.add(self.xComboBox[1],self.__tr("Select x column"))
     self.yCol[1].setText(self.__tr("y:"))
     QToolTip.add(self.yComboBox[1],self.__tr("Select y column"))
     self.Actions[1].setTitle(self.__tr("File 2"))
     self.OpenFile_3.setText(QString.null)
     QToolTip.add(self.OpenFile_3,self.__tr("Open File"))
     self.SaveFile_3.setText(QString.null)
     QToolTip.add(self.SaveFile_3,self.__tr("Save File"))
     self.ViewCols_3.setText(QString.null)
     QToolTip.add(self.ViewCols_3,self.__tr("View Columns"))
     self.ShowHeader_3.setText(QString.null)
     QToolTip.add(self.ShowHeader_3,self.__tr("Show Header"))
     self.CmdInfo.setText(self.__tr("textLabel1"))
# -----------------------------------------------------
  def __tr(self,s,c = None):
     return qApp.translate("PatWid",s,c)
# -----------------------------------------------------
# DataTable
# -----------------------------------------------------
class DataTable(QWidget):
  def __init__(self,td,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     if td == None:
	 QMessageBox.about(self,"Pat","No data to display.\nOpen first.") 
	 return
     self.nCols=len(td.Nums)
     self.nRows=td.nNLines
     self.Font=QFont("mono",10)
     #print self.nCols,self.nRows,td.Nums[0]
     self.table = QTable(self.nRows,self.nCols,self,"table")
     self.setCaption("Data Columns")
     self.table.setFont(self.Font)
     #self.table.setReadOnly(True)
     self.header = self.table.horizontalHeader()
     for i in range(len(td.ColId)):
        self.header.setLabel(i, "%1d: %s" % (i+1, td.ColId[i]))
     
     for i in range(self.nRows):
         for j in range(self.nCols):
	     self.table.setText(i,j,td.Nums[j][i])
     self.setFocusPolicy(QWidget.ClickFocus )
     self.connect(self.table,SIGNAL("contextMenuRequested(int,int,const QPoint &)"),self.RightMouse)
     self.pm= QPopupMenu ( self, "copy")
     self.pm.insertItem("copy")
     self.connect(self.pm,SIGNAL("activated(int)"),self.PopupAction)
     self.show()
# -------------------------------------------
  def RightMouse(self,row,col,pos):
      self.pm.popup(pos)
      print "context",row,col,pos.x(),pos.y()
# -------------------------------------------
  def PopupAction(self,item):
      #print "Popup",item
      if item==-2:
         cs=self.table.currentSelection()
	 s=self.table.selection(cs)
	 #print "selection:",cs,self.table.numSelections()
	 if not s.isEmpty():
	    t=""
	    if s.leftCol()==s.rightCol(): sep=''
	    else :sep='\t'
	    #print s.topRow(),s.bottomRow(),s.leftCol(),s.rightCol()
	    for r in range(s.topRow(),s.bottomRow()+1):
	        for c in range(s.leftCol(),s.rightCol()+1):
		    t+="%s%s" % (str(self.table.text(r,c)),sep)
		t+='\n'
	    #print t    
            cb = QApplication.clipboard()
            #print cb.supportsSelection ()
	    cb.setText(t,QClipboard.Clipboard)
	 else: print "Mhhh: Popupaction: empty selection"
      else: print "Mhhh: Popupaction: unexpected item"
# -----------------------------------------------------
# DataHeader
# -----------------------------------------------------
class DataHeader(QWidget):
  def __init__(self,td,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     if td == None:
	 QMessageBox.about(self,"Pat","No data to display.\nOpen first.") 
	 return
     self.Font=QFont("mono",10)

     if os.name=='posix':
        self.head=QTextEdit()
        self.head.setWordWrap(QTextEdit.NoWrap)
     else:
        self.head=QMultiLineEdit()

     self.head.setFont(self.Font)
     self.head.setReadOnly(True)
     self.head.setText(join(td.FText,"\n"))
     self.head.setCaption("Header lines")
     self.head.setGeometry(30,30,SCREEN.x()/3,len(td.FText)*10)
     self.head.show()  
# -------------------------------
# DrawWid
# -------------------------------
class DataWid(QFrame):
  def __init__(self,parent=None,name=None):
     QFrame.__init__(self,parent,name)
     self.iColX = -1
     self.iColY = -1
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
     try: self.setPaletteBackgroundColor(Qt.white)
     except AttributeError: self.setBackgroundColor(Qt.white)
     #self.setFrameShape(QFrame.StyledPanel)
     #self.setFrameShadow(QFrame.Sunken)
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
         #print "Fontsize:",iF
	 #print self.width(),self.height()
         self.PPar["rect"].setRect(5,5,self.width()-5,self.height()-5)
         DrawData(paint,self.TD,self.PPar)
      paint.end()
# -------------------------------
# -------------------------------
# CmdWid
# -------------------------------
class CmdWid(QWidget):
  def __init__(self,cmd,df_list,parent=None,name=None):
     QWidget.__init__(self,parent,name)
     self.imageGo = QPixmap(ICON_PATH+"go_b.png")
     self.imageHelp = QPixmap(ICON_PATH+"help.png")

#     print cmd
     self.Cmd=cmd
     self.DF=df_list
     self.setCaption(self.Cmd["cmd"])
     self.TmpFile=None

     #Label, command parameters, go & help button, (help)
     self.topLayout = QVBoxLayout( self, 4 )
     self.Font=QFont("freesans",11)

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
     
     if self.DF[0] != None: 
       self.Go = QToolButton(self,"goButton")
       self.Go.setGeometry(QRect(0,0,32,32))
       self.Go.setIconSet(QIconSet(self.imageGo))
       self.Go.setUsesBigPixmap(0)

       self.Go.setTextLabel("Go")
       try: self.Go.setTextPosition(QToolButton.BelowIcon)
       except: pass
       #self.Go.setMinimumSize(self.Go.sizeHint())
       self.Go.setFixedHeight(self.Go.height())
       try: self.Go.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.MinimumExpanding,0,0,self.Go.sizePolicy().hasHeightForWidth()))
       except: pass
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

     self.Help = QToolButton(self,"helpButton")
     self.Help.setGeometry(QRect(40,0,32,32))
     self.Help.setIconSet(QIconSet(self.imageHelp))
     self.Help.setUsesBigPixmap(0)
     try: self.Help.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.MinimumExpanding,0,0,self.Help.sizePolicy().hasHeightForWidth()))
     except: pass
     #self.Help.setText("Help")
     #self.Help.setMinimumSize(self.Help.sizeHint())
     self.Help.setFixedHeight(self.Help.height())
     self.Help.setTextLabel("Help")
     try: self.Help.setTextPosition(QToolButton.BelowIcon)
     except: pass
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
      if self.Cmd["files"] == 2: cmd += ' "%s"' % self.DF[1].FileName
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


