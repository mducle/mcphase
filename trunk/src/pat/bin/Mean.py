#! /usr/bin/env /usr/bin/python
"""Usage: Mean -f # [-x #] [-y #] [-s #] [[-o] or [-u]] [-v #] [-h] InputFile
        -f #: Delta X to be considered for building mean values
        -x #:             
        -y #: (int) number of x-, y-columns used for mean values.
              default: -x 1 -y 2
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
  Put all values within delta X in a bin.
  Caluclate mean values and deviations for x and y within this bin .
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
from thefile import *
from xydata import *

#$Log: Mean.py,v $
#Revision 1.2  2006/01/04 14:45:31  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/19 09:14:57  herbie
#Initial revision
#
#
CVS_ID="$Id: Mean.py,v 1.2 2006/01/04 14:45:31 herbie Exp herbie $"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=None
fp=None
delta=None
ix=1
iy=2
if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:f:x:y:ouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-f': delta=string.atof(o[1])
  if o[0]=='-x': ix=string.atoi(o[1])
  if o[0]=='-y': iy=string.atoi(o[1])
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,delta,ix,iy,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

if delta == None:
  raise GetoptError('Missing option -f #')

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
else: DF=AsciiFile(Filename)
		    
DF.Read(iSet,fp)
XY=XYData(DF.Nums.Export(ix-1),DF.Nums.Export(iy-1))
out=XY.MeanBin(delta)

if Type[0] == FT_THECAP and MOD_ftypes['the']:
   DF.ChgTextPar([('Idn','File from %s - Mean x:%d, y:%d, delta:%g' % (ID,ix,iy,delta))])
   DF.ChgTextPar([('DataPoints','%d' % len(out))])
   DF.ChgTextPar([('NoofCols','4')])
   DF.ChgTextPar([('ColID1',DF.ColId[ix])])
   DF.ChgTextPar([('ColID2',DF.ColId[iy])])
   DF.ChgTextPar([('ColID3','Sigma(Col1)')])
   DF.ChgTextPar([('ColID4','Sigma(Col2)')])
   DF.ChgTextPar([('ColID5','# of values in bin')])

if LineT !=1: DF.SetDelim(LineT/2)
for i in DF.FText: sys.stdout.write(i+DF.DeLim)
f1=DF.Nums.GetColumnFormat(ix-1,2)
f2=DF.Nums.GetColumnFormat(iy-1,2)
print DF.Nums.ColDesc
#print f1,f2
f=("%s  %s  %s  %s  %%d%%s" %(f1,f2,f1,f2))
for i in out: sys.stdout.write(f %(i[0],i[1],i[2],i[3],i[4],DF.DeLim))

sys.exit(EXIT_SUCCESS)
