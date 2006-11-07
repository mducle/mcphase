#!/usr/bin/env python
"""Usage: Sort  [-s #] -c # [-t] [-v] [-h] InputFile
        -c #: Number of column used for sort
	  -r: sort reverse
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
    RESULT:
      All columns are sorted with respect to column number specified by -c.
      InputFile can be piped (|). Output is written to stdout.
"""
import sys,string,types,os
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
from xydata import *

def rev_sort(x,y):
  if x>y: return -1
  else: return 0

#$Log$
CVS_ID="$Id$"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=None
fp=None
ic=None
DF=None
sort_f=None
if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:c:rouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-c': ic=string.atoi(o[1])
  if o[0]=='-r': sort_f=rev_sort
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,ic,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

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

l=DF.Nums.Export(ic-1,types.TupleType)
l.sort(sort_f)
ll=[]
for i in l: ll.append(i[1])

if Type[0] == FT_THECAP or Type[0] == FT_SXSMulti:
   Id='File from %s - Sort c:%s File:%s' % \
           (ID,str(ic),Filename)
   DF.ChgTextPar([('Idn',Id)])

if LineT !=1: DF.SetDelim(LineT/2)

for i in DF.FText:  sys.stdout.write("%s%s" % (i,DF.DeLim))

DF.Nums.Write(fp=sys.stdout,Sort=ll)
sys.exit(EXIT_SUCCESS)
