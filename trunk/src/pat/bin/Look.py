#! /usr/bin/env /usr/bin/python
"""Usage: Look [-x #,#] [-y #,#] [-s #,#] [-a] [[-o] or [-u]] [-v #] [-h] (Data) Table
      -x #,#:             
      -y #,#: (int) number of x-, y-columns used for table lookup.
              Two numbers separated by commas can be specified, refering 
              to Data file or Table file, respectively   
              default: -x 1,1 -y 2,2
          -a: append interpolated values as last column to Data file
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
      -s #,#: (int) number of data sets to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
    A table lookup is performed:
    For all x-values of the Data file a table lookup (cubic spline)
    in the x-values are performed and the corresponding y-values of the
    Data file are replaced by the interpolated y-values of the Table file
    Columns not involved in the operation are unchanged.
    1st x: x-column of DataFile (data values)
    2nd x: x-column of TableFile (data values)
    1st y: y-column of DataFile (interpolated values) will be replaced or inserted
    2nd y: y-column of  TableFile (table values)
       -a: Append; interpolated values are appended as last column,
           else replaced
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

#$Log: Look.py,v $
#Revision 1.1  2005/12/15 09:16:48  herbie
#Initial revision
#
#
CVS_ID="$Id: Look.py,v 1.1 2005/12/15 09:16:48 herbie Exp $"
ID=string.join(CVS_ID.split()[1:4])

iSet=[1,1]
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=[None,None]
fp=[None,None]
ix=[1,1]
iy=[2,2]
Type=[None,None]
DF=[None,None]
XY=[None,None]
Append=0

if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:x:y:aouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s':
     t=o[1].split(',')
     for i in [0,1]:
         if t[i]!='': iSet[i]=string.atoi(t[i])
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
  if o[0]=='-a': Append=1
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,Append,ix,iy,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

if len(opt[-1]) > 2 : print __doc__; sys.exit(EXIT_SUCCESS)
if len(opt[-1]) == 0 : raise GetoptError('At least one file name necessary')
if len(opt[-1]) == 1 : Filename[1]=opt[-1][0]
if len(opt[-1]) == 2 : Filename[0]=opt[-1][0];Filename[1]=opt[-1][1];

S.Log(16,str([iSet,Append,ix,iy,LineT,Filename,iDbg]))

if Filename[0] == None:
   fp[0]=StdinToFile(sys.argv[0])
else:
   fp[0]=open(Filename[0],"r")

fp[1]=open(Filename[1],"r")

for i in [0,1]:
  Type[i]=CheckFileType(fp[i])
  DF[i]=AsciiFile(Filename[i])
		    
  DF[i].Read(iSet[i],fp[i])
  XY[i]=XYData(DF[i].Nums.Export(ix[i]-1),DF[i].Nums.Export(iy[i]-1))


XY[0].LookTable(XY[1].X(),XY[1].y)

w=XY[0].Weight()
DF[0].Nums.Import(ix[0]-1,XY[0].X())
if Append: DF[0].Nums.Import(-1,XY[0].y)
else: DF[0].Nums.Import(iy[0]-1,XY[0].y)

if Type[0][0] == FT_THECAP or Type[0][0] == FT_SXSMulti:
   Id='File from %s - Look x:%s y:%s File:%s' % \
           (ID,str(ix),str(iy),Filename[1])
   DF[0].ChgTextPar([('Idn',Id),('DataPoints','%d' % w[1])])

if LineT !=1: DF.SetDelim(LineT/2)

n=DF[0].Write(sys.stdout,Sort=w[0])
d=len(DF[0].Nums[0]) - n[0]
if d != 0:
   S.Log(0,"%d values not in table; skipped" % d)
   
sys.exit(EXIT_SUCCESS)
