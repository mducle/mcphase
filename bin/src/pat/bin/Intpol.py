#! /usr/bin/env /usr/bin/python
"""Usage: Intpol -l -c [-f string] [-x #] [-y #] [-s #] [[-o] or [-u]] [-v #] [-h] TableFile
          -c: perform a cubic spline interpolation
          -l: perform a linear interpolation
        -x #:             
        -y #: (int) number of x-, y-columns used for derivation.
              default: -x 1 -y 2
   -f string: C-format string to print x and the result
              default '%f %f'
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
   TableFile: Interpolation table data file
  RESULT:
  Reads a x value from stdin. Prints x and the interpolated y value for x to stdout.
  The interpolation table is obtained from TableFile
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
from asciifile import *
from xydata import *
from spline import *

#$Log: Intpol.py,v $
#Revision 1.1  2005/12/15 09:16:48  herbie
#Initial revision
#
#
CVS_ID="$Id: Intpol.py,v 1.1 2005/12/15 09:16:48 herbie Exp $"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=None
SplineMeth=1
Format=None
ix=1
iy=2
if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:clf:x:y:ouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-l': SplineMeth*=2
  if o[0]=='-c': SplineMeth*=4
  if o[0]=='-f': Format=o[1]
  if o[0]=='-x': ix=string.atoi(o[1])
  if o[0]=='-y': iy=string.atoi(o[1])
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

if SplineMeth > 4 or SplineMeth < 2:
  raise GetoptError("You must specify exact one spline method (-l or -c)")
  sys.exit(EXIT_FAILURE)

Filename=opt[-1][0]

if sys.stdin.isatty():
   raise IOError("No input data on stdin")

xVal=string.atof(sys.stdin.readline())
S.Log(16,str([iSet,xVal,SplineMeth,ix,iy,LineT,Filename,iDbg]))

DF=AsciiFile(Filename)
		    
DF.Read(iSet)
XY=XYData(DF.Nums.Export(ix-1),DF.Nums.Export(iy-1))

if SplineMeth == 2:
  yVal=XY.LinIntpol(xVal)
elif SplineMeth == 4:
  S=SplineTable(XY.X(),XY.y)
  S.Spline()
  yVal=S.Intpol(xVal)
else:
  raise ValuError('Hmmm, allowed to be 2 or 4 but SplineMeth = %d' % SplineMeth)

if Format == None: print  xVal,yVal
else: print Format % (xVal,yVal)

sys.exit(EXIT_SUCCESS)
