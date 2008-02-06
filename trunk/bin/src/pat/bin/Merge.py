#! /usr/bin/env /usr/bin/python
"""Usage: Merge [[-c] or [-r]] [-s #] [-a] [[-o] or [-u]] [-v #] [-h] (File1) File2
          -c: Merge files by columns
          -r: Merge files by rows
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
      -s #,#: (int) number of data sets to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
  RESULT:
    The two file are merged to one data file. The files can have different 
    file types. The output file inherits the type of File1.
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

#$Log: Merge.py,v $
#Revision 1.1  2005/12/15 09:16:48  herbie
#Initial revision
#
#
CVS_ID="$Id: Merge.py,v 1.1 2005/12/15 09:16:48 herbie Exp $"
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
Merge=1
if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:crouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s':
     t=o[1].split(',')
     for i in [0,1]:
         if t[i]!='': iSet[i]=string.atoi(t[i])
  if o[0]=='-c': Merge=Merge*2
  if o[0]=='-r': Merge=Merge*4
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,Merge,LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")
if Merge > 4 :
  raise  GetoptError("-c or -r: only one option possible")

if len(opt[-1]) > 2 : print __doc__; sys.exit(EXIT_SUCCESS)
if len(opt[-1]) == 0 : raise GetoptError('At least one file name necessary')
if len(opt[-1]) == 1 : Filename[1]=opt[-1][0]
if len(opt[-1]) == 2 : Filename[0]=opt[-1][0];Filename[1]=opt[-1][1];

S.Log(16,str([iSet,Merge,LineT,Filename,iDbg]))

if Filename[0] == None:
   fp[0]=StdinToFile(sys.argv[0])
else:
   fp[0]=open(Filename[0],"r")

fp[1]=open(Filename[1],"r")

for i in [0,1]:
  Type[i]=CheckFileType(fp[i])
  DF[i]=AsciiFile(Filename[i])
		    
  DF[i].Read(iSet[i],fp[i])

n1=len(DF[0].Nums[0])
n2=len(DF[1].Nums[0])
c1=len(DF[0].Nums)
c2=len(DF[1].Nums)

if Merge == 2:  # append columns
  if n1 != n2:
     raise ValueError('Number of data lines (%d,%d) differ' % (n1,n2))
  for i in  DF[1].Nums:
     DF[0].Nums.append(i)

  if Type[0][0] == FT_THECAP or Type[0][0] == FT_SXSMulti:
        nc=[]
        for i in range(0,c2):
            nc.append("ColID%d=%s" %((c1+i+1),DF[1].ColId[i]))
        for i in range(0,c2): DF[0].InsTextPar('ColID%d' % (c1+i),nc[i])
elif Merge == 4: # append rows
  if c1 != c2:
     raise ValueError('Number of data columns (%d,%d) differ' % (c1,c2))
  DF[0].Nums.append_rows(DF[1].Nums)

  if Type[0][0] == FT_THECAP or Type[0][0] == FT_SXSMulti:
     DF[0].ChgTextPar([('DataPoints',n1+n2)])
     
else:
  raise ValueError('Hmmm, allowed to be 2 or 4 but Merge = %d' % Merge)

if Type[0][0] == FT_THECAP or Type[0][0] == FT_SXSMulti:
   Id='File from %s - Merge %d File:%s' % \
           (ID,Merge,Filename[1])
   DF[0].ChgTextPar([('Idn',Id)])

if LineT !=1: DF.SetDelim(LineT/2)

DF[0].Write(sys.stdout)

sys.exit(EXIT_SUCCESS)
