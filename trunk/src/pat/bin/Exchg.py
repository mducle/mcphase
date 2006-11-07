#! /usr/bin/env /usr/bin/python
"""Usage: Exchg -c #,# [-s #] [[-o] or [-u]] [-v #] [-h] InputFile
        -c #,#: (ints) number of columns to be exchanged
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
  Swaps two columns
  Other columns are not affectd
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
from datafile import *
from xydata import *

#$Log: Exchg.py,v $
#Revision 1.1  2004/08/04 11:28:31  herbie
#*** empty log message ***
#
CVS_ID="$Id: Exchg.py,v 1.1 2004/08/04 11:28:31 herbie Exp $"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=None
fp=None
l0=None
c=[]
if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:c:ouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-c':
    t=o[1].split(',')
    c=[string.atoi(x) for x in t] 
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

if len(c) < 2:
  raise GetoptError("Option -c #,# missing or incomplete")

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

S.Log(16,str([iSet,c[0],c[1],LineT,Filename,iDbg]))

DF.Nums.swap_col(c[0]-1,c[1]-1)

if Type[0] == FT_THECAP or Type[0] == FT_SXSMulti:
   t=DF.ColId[c[0]-1]
   DF.ColId[c[0]-1]=DF.ColId[c[1]-1]
   DF.ColId[c[1]-1]=t
   Id='File from %s - Exchg c%d <-> c%d' % (ID,c[0],c[1])
   DF.ChgTextPar([('Idn',Id)])

if LineT !=1: DF.SetDelim(LineT/2)

DF.Write(sys.stdout)

sys.exit(EXIT_SUCCESS)
