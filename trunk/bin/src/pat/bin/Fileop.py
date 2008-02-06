#! /usr/bin/env /usr/bin/python
"""Usage: Fileop f # [-x #] [-y #] [-s #] [[-o] or [-u]] [-v #] [-h] File1 File2
       -f op: A single character defining operation one of <+-*/>
    -x #[,#]:             
    -y #[,#]: (int) number of x-, y-columns used for calculation.
              Two numbers separated by commas can be specified, refering 
              to File1 or File2, respectively   
              default: -x 1,1 -y 2,2
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
    The operation defined by the character is performed:
    y1(x1) = y1(x1) <op> y2(x2). 
    x(File1) and x(File2) need not be same, linear interpolation is used to get
    the value y2(x1) from y2(x2).
    Columns not involved in the operation are unchanged.
    File1 can be piped (|). Output is written to stdout.
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

#$Log: Fileop.py,v $
#Revision 1.2  2005/12/19 09:14:57  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 09:16:48  herbie
#Initial revision
#
#
CVS_ID="$Id: Fileop.py,v 1.2 2005/12/19 09:14:57 herbie Exp $"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=[None,None]
fp=[None,None]
op=None
ix=[1,1]
iy=[2,2]
Type=[None,None]
DF=[None,None]
XY=[None,None]

if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:f:x:y:ouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-f': op=o[1]
  if o[0]=='-x': 
     t=o[1].split(',')
     for i in [0,1]:
         if t[i]!='': ix[i]=string.atoi(t[i])
  if o[0]=='-y':
     t=o[1].split(',')
     for i in [0,1]:
        if t[i]!='': iy[i]=string.atoi(t[i])
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,op,ix,iy,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError('Missing option -f "formula"')

if len(opt[-1]) > 2 : print __doc__; sys.exit(EXIT_SUCCESS)
if len(opt[-1]) == 0 : raise GetoptError('At least one file name necessary')
if len(opt[-1]) == 1 : Filename[1]=opt[-1][0]
if len(opt[-1]) == 2 : Filename[0]=opt[-1][0];Filename[1]=opt[-1][1];

S.Log(16,str([iSet,op,ix,iy,LineT,Filename,iDbg]))

if Filename[0] == None:
   fp[0]=StdinToFile(sys.argv[0])
else:
   fp[0]=open(Filename[0],"r")

fp[1]=open(Filename[1],"r")

for i in [0,1]:
  Type[i]=CheckFileType(fp[i])
  DF[i]=AsciiFile(Filename[i])
		    
  DF[i].Read(iSet,fp[i])
  XY[i]=XYData(DF[i].Nums.Export(ix[i]-1),DF[i].Nums.Export(iy[i]-1))


XY[0].MathOper(XY[1],op)

w=XY[0].Weight()
DF[0].Nums.Import(ix[0]-1,XY[0].X())
DF[0].Nums.Import(iy[0]-1,XY[0].y)

if Type[0][0] == FT_THECAP or Type[0][0] == FT_SXSMulti:
   Id='File from %s - Fileop "%s" x:%s y:%s File:%s' % \
           (ID,op,str(ix),str(iy),Filename[1])
   DF[0].ChgTextPar([('Idn',Id),('DataPoints','%d' % w[1])])

if LineT !=1: DF.SetDelim(LineT/2)

DF[0].Write(sys.stdout,Sort=w[0])

sys.exit(EXIT_SUCCESS)
