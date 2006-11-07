import sys,string
from tempfile import *
from types import *

import smtplib
import posixpath
from email.MIMEText import MIMEText

#$Log: stdfunc.py,v $
#Revision 1.2  2004/12/22 16:19:56  herbie
#*** empty log message ***
#
#Revision 1.1  2004/08/04 11:17:18  herbie
#*** empty log message ***
#
CVS_ID="$Id: stdfunc.py,v 1.2 2004/12/22 16:19:56 herbie Exp $"

_HUGE_VAL = 1.79769313486231e+308

EPSILON0=0.8854188E-11 # A sec/V m Dielektrizitaetskonstante

LT_UNIX = 1
LT_DOS  = 2
LINE_T=['','DOS','UNIX']

DOS_END='\r\n'
UNIX_END='\n'

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

def sign(x):
   if not ( type(x) == FloatType or type(x) == IntType):
      raise ValueError('Bad input parameter x')
      return x
   if x==0: return 0
   if x<0:  return -1
   if x>0:  return 1


def StdinToFile(name):
   if sys.stdin.isatty():
     print "No input file (%s)\n" % name
     sys.exit(EXIT_FAILURE)

   suff=posixpath.basename(sys.argv[0].split(".")[0])
   fp=TemporaryFile('w+',-1,suff)
   b=sys.stdin.readlines()
   fp.writelines(b)
   fp.seek(0,0)
   return fp

def swap(a,b): t=a; a=b; b=t

def ReadSockWErr(sock,msg,log_msg,SL,log_pr=3,exit=0):
  data=sock.read(msg)
  if string.find(data,"-ERR") != -1:
    if log_pr & 1: SL.Log(0, log_msg+"%s" % data)
    if log_pr & 2: print  log_msg+"%s" % data
    if exit: sys.exit(1)
  return data

def mail(sender,recip,subj,msg):
  if type(recip)==StringType: r=[recip]
  else: r=list(recip)
  m = MIMEText(msg)

  m['Subject'] = subj
  m['From'] = sender
  m['To'] = r[0]
  s = smtplib.SMTP()
  s.connect()
  s.sendmail(sender, r , m.as_string())
 
def LinTrafo(pmin_max,dmin_max, d):
    diffz= pmin_max[1]-pmin_max[0]
    diffn= dmin_max[1]-dmin_max[0]
    return  (diffz/diffn*d + (pmin_max[1]-diffz/diffn*dmin_max[1]))

def BackTrafo(pmin_max,bmin_max, b):
    diffz= pmin_max[1]-pmin_max[0]
    diffn= bmin_max[1]-bmin_max[0]
    return  ((b-(pmin_max[1]-(diffz/diffn)*bmin_max[1]))*diffn/diffz)
