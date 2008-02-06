#! /usr/bin/env /usr/bin/python
"""Usage: Deriv -n # [-x #] [-y #] [-s #] [[-o] or [-u]] [-v #] [-h] InputFile
        -n #: Number of neigbourpoints to be concerned for derivation
        -x #:             
        -y #: (int) number of x-, y-columns used for derivation.
              default: -x 1 -y 2
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
  dy/dx is build by an algorithm using n neigbour points.
  Values will be sorted by x-column. Other columns than x/y are not affectd
  InputFile can be piped (|). Output is written to stdout.
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

#$Log: Deriv.py,v $
#Revision 1.1  2005/12/15 09:16:48  herbie
#Initial revision
#
#
CVS_ID="$Id: Deriv.py,v 1.1 2005/12/15 09:16:48 herbie Exp $"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=None
fp=None
nP=None
ix=1
iy=2
if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:n:x:y:ouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-n': nP=string.atoi(o[1])
  if o[0]=='-x': ix=string.atoi(o[1])
  if o[0]=='-y': iy=string.atoi(o[1])
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,nP,ix,iy,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

if nP == None:
  raise GetoptError('Missing option -n #')

if len(opt[-1]): Filename=opt[-1][0]

if Filename == None:
   fp=StdinToFile(sys.argv[0])
else:
   fp=open(Filename,"r")

Type=CheckFileType(fp)
S.Log(16,str(Type))

if Type[0] == None:
   raise ValueError('%s: Bad input file type' % Filename)

DF=AsciiFile(Filename)
		    
DF.Read(iSet,fp)
XY=XYData(DF.Nums.Export(ix-1),DF.Nums.Export(iy-1))
XY.Deriv(nP)

for i in range(0,nP): DF.Nums.pop_row()

DF.Nums.Import(ix-1,XY.X())
DF.Nums.Import(iy-1,XY.y)

if Type[0] == FT_THECAP or Type[0] == FT_SXSMulti:
   DF.ChgTextPar([('Idn','File from %s - Deriv x:%d, y:%d, n:%d' % (ID,ix,iy,nP))])

if LineT !=1: DF.SetDelim(LineT/2)

DF.Write(sys.stdout)

sys.exit(EXIT_SUCCESS)
