#!/usr/bin/env python
"""Usage: Select  [-s #] [-c #,#] [-r #,#[,#]] [-f #,#,#] [-t] [-v] [-h] InputFile
      -c #,#: Range of columns (start,end) to be selected.
              If -c ommited all columns are selected
    -r #,#,#: Range of rows to be selected.
              If -r ommited all rows are selected
              Range: start,end,step, if step ommited each row is selected
    -f #,#,#: Range of values from column (column,min,max) to be selected.
              Range: column,min,max
          -o: OutFile in DOS <lf><cr> format
          -u: OutFile in UNIX <lf> format
        -s #: (int) number of data set to perform operation  
              (Only for file types supporting multiple data sets)
        -v #: Print some debug info (stderr)
          -V: Print Version number
          -h: Print this help message
    RESULT:
    Selects the specified rows/columns of a file.
    EXAMPLE: 
    Select -c 2,3 -r 1,50 Data.cnt
      Writes column 2 and 3 and the rows from 1 - 50 of the file Data.cnt
      to stdout

    Select -f 1,2.,30.5 Data.cnt
      Writes all rows with values between 2 and 3.5 in coumn 1 to stdout
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
try: from sxsfile import *
except: pass
from xydata import *

#$Log: Select.py,v $
#Revision 1.2  2006/01/04 14:45:31  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 09:16:48  herbie
#Initial revision
#
CVS_ID="$Id: Select.py,v 1.2 2006/01/04 14:45:31 herbie Exp herbie $"
ID=string.join(CVS_ID.split()[1:4])

iSet=1
iDbg=0
S=SysLog(iDbg,sys.argv[0],2)
LineT=1
Filename=None
fp=None
c=[-1,-1]
r=[-1,-1,-1]
v=[-1,None,None]
DF=None

if len(sys.argv)<1: sys.exit(EXIT_FAILURE);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(EXIT_FAILURE);

opt=getopt (sys.argv[1:], "s:c:r:f:ouv:Vh")
S.Log(16,str(opt))
for o in opt[0]:
  if o[0]=='-s': iSet=string.atoi(o[1])
  if o[0]=='-f': 
    t=o[1].split(',')
    if t[0]!='': v[0]=string.atoi(t[0])-1
    if t[1]!='': v[1]=string.atof(t[1])
    if t[2]!='': v[2]=string.atof(t[2])
    
  if o[0]=='-c':
    t=o[1].split(',')
    if len(t)==1:
       x=string.atoi(t[0])-1
       c=[x,x]
    elif len(t)==2: 
       for i, x in enumerate(t): 
           if x!='': c[i]=string.atoi(x)-1
    else:
       raise ValueError("Bad range specifier for -c")
  if o[0]=='-r':
    t=o[1].split(',')
    if len(t)==1: r[2]=string.atoi(t[0])
    elif len(t)==2:
       r[0]=string.atoi(t[0])-1
       r[1]=string.atoi(t[1])-1
    elif len(t)==3:
       if t[0]!='': r[0]=string.atoi(t[0])-1
       if t[1]!='': r[1]=string.atoi(t[1])-1
       if t[2]!='': r[2]=string.atoi(t[2])
    else:
       raise ValueError("Bad range specifier for -r")
       sys.exit(EXIT_SUCCESS)
  if o[0]=='-o': LineT=LineT*4
  if o[0]=='-u': LineT=LineT*2
  if o[0]=='-v': iDbg=string.atoi(o[1])
  if o[0]=='-V': print CVS_ID; sys.exit(EXIT_SUCCESS)
  if o[0]=='-h': print __doc__; sys.exit(EXIT_SUCCESS)
S.SetLevel(iDbg)
S.Log(16,str([iSet,str(c),str(r),str(v),LineT,Filename,iDbg]))

if LineT > 4 :
  raise GetoptError("Options -o and -u are mutual exclusive")
  sys.exit(EXIT_FAILURE)

if len(opt[-1]): Filename=opt[-1][0]
if Filename == None:
   fp=StdinToFile(sys.argv[0])
else:
   fp=open(Filename,"r")

Type=CheckFileType(fp)
S.Log(16,str(Type))
if Type[0] == None:
   raise ValueError('%s: Bad input file type' % Filename)

if Type[0] == FT_SXSMulti and MOD_ftypes['sxs']: DF=SxSFile(Filename)
else: 
  DF=AsciiFile(Filename)

S.Log(16,str(MOD_ftypes))
		    
DF.Read(iSet,fp)
w=DF.Nums.Write(sys.stdout,0,r,c,v,quiet=1)

if Type[0] == FT_THECAP or Type[0] == FT_SXSMulti:
   Id='File from %s - Select c:%s r:%s f:%s File:%s' % \
           (ID,str(c),str(r),str(v),Filename)
   n_ci=[]
   for i,x in enumerate(w[1]): n_ci.append(("ColID%1d" % (i+1),DF.ColId[x])) 
   for i in range(len(n_ci,),len(DF.ColId)):n_ci.append(("ColID%1d" % (i+1),"-"))
   n_ci.append(("NoofCols","%d" % len(w[1])))
   DF.ChgTextPar([('Idn',Id),('DataPoints','%d' % w[0])]+n_ci)
   

if LineT !=1: DF.SetDelim(LineT/2)

for i in DF.FText:  sys.stdout.write("%s%s" % (i,DF.DeLim))

DF.Nums.Write(sys.stdout,0,r,c,v)
sys.exit(EXIT_SUCCESS)
