#! /usr/bin/env /usr/bin/python
"""Usage: Pat [-v #] [-h]
          -V: Print Version number
          -h: Print this help message
  RESULT:
  Shows a graphic window where all pat commands are available
  Keys changing the layout:
        +-: Change symbol size
  up, down: Change font size
         C: Connect points     
         E: Plot error bars for count rate
	 P: Print graph
	 0: Plot zero line
    Q, Esc: Quit 
"""
import sys,string,os 
from getopt import *

lib_path='PAT_LIB_PATH'
try: path=os.environ[lib_path]
except KeyError:
  sys.stderr.write('Set environment variable %s="path to lib-files of pat-package"\n' % lib_path)
  sys.exit(0)
sys.path.append(path)
#print sys.path
from SysLog import *
from stdfunc import *
from asciifile import *
try: from thefile import *
except: pass
from xydata import *
from grafic import *
from wid_ui import *

#$Log: Pat.py,v $
#Revision 1.1  2006/01/04 14:45:31  herbie
#Initial revision
#
#
CVS_ID="$Id: Pat.py,v 1.1 2006/01/04 14:45:31 herbie Exp herbie $"

iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
if len(sys.argv)<1: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE)

opt=getopt (sys.argv[1:], "v:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)

QApplication.setColorSpec( QApplication.CustomColor )
a=QApplication(sys.argv)
QApplication.setFont( QFont("sanserif") )
  
d=QApplication.desktop()

dw=d.width()
if dw > 1280: Swidth=1280
else: Swidth=dw

SCREEN.setX(Swidth)
SCREEN.setY(d.height())

#print SCREEN.x(),SCREEN.y()

sw=PatWid()
#sw.SetEndExit()
a.setMainWidget( sw )
sw.show()
a.exec_loop()
