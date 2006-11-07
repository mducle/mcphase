#! /usr/bin/env /usr/bin/python
"""Usage: Disp [-x #] [-y #] [-s #] [-v #] [-h] InputFile
        -x #:
        -y #: (int) number of x-, y-columns used for plot.
              default: -x 1 -y 2
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
  x versus y column of the data file are displayed on screen.
  Data can be changed "on line" by several functions.
  X-Windows is needed.
  The following keys can be used to control the layout of the display:
  Function keys (F1- F10) allow to manipulate data
    + Key: increase symbolsize
    - Key: decrease symbolsize to default value
    E Key: toggle error bars/symbols
    C Key: connect points with line
  ErrorBar: ONLY for spectroscopic data !! (Y: count rate)
            displays error bars: sqrt(Y)
  Changed data is written to stdout. 
  InputFile can be piped (>).
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
from datafile import *
try: from thefile import *
except: pass
from thefile import *
from xydata import *
from grafic import *

#$Log:$
CVS_ID="$Id:$"

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
Filename=None
fp=None
ix=1
iy=2
if len(sys.argv)<1: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE)

opt=getopt (sys.argv[1:], "s:x:y:v:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-x': ix=string.atoi(o[1])
  if o[0]=='-y': iy=string.atoi(o[1])
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)

if len(opt[-1]): Filename=opt[-1][0]

if Filename == None:
   fp=StdinToFile(sys.argv[0])
else:
   fp=open(Filename,"r")

Type=CheckFileType(fp)
S.Log(16,str(Type))

if Type[0] == None:
   raise ValueError('%s: Bad input file type' % Filename)

if Type[0] == FT_THECAP and MOD_ftypes['the']: DF=TheFile(Filename)
elif Type[0] == FT_SXSMulti and MOD_ftypes['sxs']: DF=SxSFile(Filename)
else: DF=AsciiFile(Filename)

S.Log(16,str(MOD_ftypes))

DF.Read(iSet,fp)

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

DF.SetXYCol(ix-1,iy-1)

mp=ManipData(DF)
#sw.SetEndExit()
a.setMainWidget( mp )
mp.show()
a.exec_loop()
