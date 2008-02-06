#! /usr/bin/env /usr/bin/python
"""Usage: SetText -l # -t string [[-o] or [-u]] [-v #] [-h]  InputFile
   -t string: New text
        -l #: Number of text line
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
  The line specified by -l # will be changed.
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

#$Log: Stext.py,v $
#Revision 1.1  2005/12/15 09:16:48  herbie
#Initial revision
#
#
CVS_ID="$Id: Stext.py,v 1.1 2005/12/15 09:16:48 herbie Exp $"
ID=string.join(CVS_ID.split()[1:4])

iLine=None
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Text=None
Filename=None
fp=None

if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "t:l:uv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-l': iLine=string.atoi(o[1])
  if o[0]=='-t': Text=o[1]
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iLine,Text,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")

if Text == None:
  raise GetoptError('Missing option -t "text"')

if iLine == None:
  raise GetoptError('Missing option -l #')
  sys.exit(EXIT_FAILURE)

if len(opt[-1]): Filename=opt[-1][0]

if Filename == None:
   fp=StdinToFile(sys.argv[0])
else:
   fp=open(Filename,"r")

DF=AsciiFile(Filename); FErr=0
 		    
DF.Read(1,fp)

try: DF.FText[iLine-1]=Text
except IndexError:
  raise IndexError("Bad line number %d" % iLine)
  sys.exit(EXIT_FAILURE)

if LineT !=1: DF.SetDelim(LineT/2)

DF.Write(sys.stdout)

sys.exit(EXIT_SUCCESS)
