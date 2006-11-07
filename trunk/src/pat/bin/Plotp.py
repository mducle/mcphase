#! /usr/bin/env /usr/bin/python
"""Usage: Plotp [-x #[,...]] [-y #[,...]] [-v] [-h]  (File1 or -) File2 ...
        -x #,...:             
        -y #,...#: (int) number of x-, y-columns used for the plot
	        Numbers separated by commas can be specified, refering 
		to first file, second file, ...  respectively   
	        If ommited the default values x:1, y:2 are assumed.
            -v: verify -> print header before and after operation (stderr)          
            -V: Print Version number
            -h: Print this help message 
     InputFile: Input data file
RESULT:
The x-, y-columns of te data files are plotted on screen.
A plot on a printer can be requested.
X-Windows is needed.
If a plot is displayed the following keys can be used to control
the layout:
    # + Key: increase symbolsize of curve #
    # - Key: decrease symbolsize of curve #
    # E Key: set/reset error bars for count rate of curve #
    # C Key: connect points with line
      Q,Esc: Quit
If only one file, the file name can be ommited to be piped (|) from stdin.
If more files use "-" for the first file name then it can be piped from stdin.
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
try: from sxsfile import *
except: pass
from xydata import *
from grafic import *

#$Log:$
CVS_ID="$Id:$"

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
Filename=None
fp=None
ix=[]
iy=[]
s=[]
if len(sys.argv)<1: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE)

opt=getopt (sys.argv[1:], "s:x:y:v:Vh")
for o in opt[0]:
  if o[0]=='-s': s=o[1].split(',')
  if o[0]=='-x': ix=o[1].split(',')
  if o[0]=='-y': iy=o[1].split(',')
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str(opt))

if len(opt[-1]): Filename=list(opt[-1])
if len(opt[-1]) > 8:
   raise ValueError('Plot list to long (<=8)')
   
fp=[]
for i in Filename:
    S.Log(16,str(i))
    if i == '-':
       fp.append(StdinToFile(sys.argv[0]))
    else :
       fp.append(open(i,"r"))
   
DF=[]
for i in range(0,len(fp)):
    r=CheckFileType(fp[i])
    if r == None:
       raise ValueError('%s: Bad input file type' % Filename[i])
    S.Log(16,str(r))
    if r[0] == FT_THECAP and MOD_ftypes['the']: DF.append(TheFile(Filename[i]))
    elif r[0] == FT_SXSMulti and MOD_ftypes['sxs']: DF.append(SxSFile(Filename[i]))
    else: DF.append(AsciiFile(Filename[i]))
xl=[]
yl=[]
sl=[]
for i in range(0,len(Filename)):
    xl.append(1)
    yl.append(2)
    sl.append(1)
for i in range(0,len(ix)): xl[i]=string.atoi(ix[i])
for i in range(0,len(iy)): yl[i]=string.atoi(iy[i])
for i in range(0,len(s)):  sl[i]=string.atoi(s[i])
S.Log(16,str(xl))
S.Log(16,str(yl))
S.Log(16,str(sl))
S.Log(16,str(MOD_ftypes))
for i in range(0,len(Filename)):
    DF[i].Read(sl[i],fp[i])
    DF[i].SetXYCol(xl[i]-1,yl[i]-1)
QApplication.setColorSpec( QApplication.CustomColor )
a=QApplication(sys.argv)
QApplication.setFont( QFont("sansserif") )
  
d=QApplication.desktop()

dw=d.width()
if dw > 1280: Swidth=1280
else: Swidth=dw

SCREEN.setX(Swidth)
SCREEN.setY(d.height())

#print SCREEN.x(),SCREEN.y()

pw=PlotData(DF)
pw.EndExit()
a.setMainWidget( pw )
pw.show()
a.exec_loop()
