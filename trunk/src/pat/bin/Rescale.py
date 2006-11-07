#! /usr/bin/env /usr/bin/python
"""Usage: Rescale [-h] [-v] [-V] GraficFile
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
  Plots a grafic file (jpg,gif,png) on the screen.
  Performs a scaling procedure to obtain original data
  No output is written to stdout.
  Keys changing the layout:
        +-: Change picture size
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

from SysLog import *
from stdfunc import *
from grafic import *

#$Log:$
CVS_ID="$Id:$"

iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
Filename=None
fp=None
if len(sys.argv)<1: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE)

opt=getopt (sys.argv[1:], "v:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)

Filename=opt[-1][0]
fp=open(Filename,"r")
fp.close()

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

sw=RescaleWid(Filename)
a.setMainWidget( sw )
sw.show()
a.exec_loop()
