#! /usr/bin/env /usr/bin/python
"""Usage: Calc -f string [-s #] [[-o] or [-u]] [-v #] [-h] InputFile
   -f string: Formula defining calculation
              Allowed formulas  (#: Column number : 1 - n)
                c# = 1.23 <+-*/^> c#
                c# = c# <+-*/^> 1.234
                c# = c# <+-*/^> c#
                c# = FUN(c#)
                c# = FUN(1.234)
                FUN: ABS, COS, EXP, LOG, SIN,  TAN, SEQ,
                     ACOS, ASIN, ATAN, SQRT
                Not case sensitive; blancs are ignored
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
  The operation defined by the formula is performed.
  Columns not involved in the formula are unchanged.
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

#$Log: Calc.py,v $
#Revision 1.1  2004/08/04 11:26:03  herbie
#*** empty log message ***
#
CVS_ID="$Id: Calc.py,v 1.1 2004/08/04 11:26:03 herbie Exp $"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Formula=None
Filename=None
fp=None

if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:f:ouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-f': Formula=o[1]
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,Formula,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

if Formula == None:
  raise GetoptError('Missing option -f "formula"')

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
DF.Nums.Calculate(Formula)

if Type[0] == FT_THECAP or Type[0] == FT_SXSMulti:
   DF.ChgTextPar([('Idn','File from %s - Calc Formula: %s' % (ID,Formula))])
 
if LineT !=1: DF.SetDelim(LineT/2)

DF.Write(sys.stdout)

S.Log(8,str([iSet,Formula,LineT,Filename,iDbg]))

sys.exit(EXIT_SUCCESS)
