#! /usr/bin/env /usr/bin/python
"""Usage: CTime UnixTimestamp
          -l: Long (Numeric format)
          -V: Print Version number
          -h: Print this help message
  RESULT:
     Print UnixTimestamp im human readable format
"""

import sys,string,os,time
from getopt import *

#$Log:$
#
CVS_ID="$Id:$"
ID=string.join(CVS_ID.split()[1:4])
if len(sys.argv)<1: sys.exit(1);
if len(sys.argv)<2: sys.stderr.write(__doc__);sys.exit(1);

out=1
opt=getopt (sys.argv[1:], "lVh")
for o in opt[0]:
  if o[0]=='-l':
     for i in time.localtime(float(opt[-1][0])):
       print "%02d" % i,
     print
     sys.exit(0)
  if o[0]=='-V': print CVS_ID; sys.exit(0)
  if o[0]=='-h': print __doc__; sys.exit(0)

print time.ctime(float(opt[-1][0]))
