import sys,string,os,popen2,time
from tempfile import *
from types import *

import smtplib
import posixpath
from email.MIMEText import MIMEText

#$Log: stdfunc.py,v $
#Revision 1.4  2006/07/11 12:10:08  herbie
#*** empty log message ***
#
#Revision 1.3  2006/01/04 14:41:36  herbie
#*** empty log message ***
#
#Revision 1.2  2005/12/19 09:13:05  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 08:57:36  herbie
#Initial revision
#
CVS_ID="$Id: stdfunc.py,v 1.4 2006/07/11 12:10:08 herbie Exp herbie $"

_HUGE_VAL = 1.79769313486231e+308

EPSILON0=0.8854188E-11 # A sec/V m Dielektrizitaetskonstante

LT_UNIX = 1
LT_DOS  = 2
LINE_T=['','DOS','UNIX']

DOS_END='\r\n'
UNIX_END='\n'

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

# -------------------------------------
def sign(x):
   if not ( type(x) == FloatType or type(x) == IntType):
      raise ValueError('Bad input parameter x')
      return x
   if x==0: return 0
   if x<0:  return -1
   if x>0:  return 1
# -------------------------------------
def StdinToFile(name):
   if sys.stdin.isatty():
      raise IOError('No input file %s' % name)
   
#     print "No input file (%s)\n" % name
#     sys.exit(EXIT_FAILURE)

   suff=posixpath.basename(sys.argv[0].split(".")[0])
   fp=TemporaryFile('w+',-1,suff)
   b=sys.stdin.readlines()
   fp.writelines(b)
   fp.seek(0,0)
   return fp
# -------------------------------------
def swap(a,b): t=a; a=b; b=t
# -------------------------------------
def ReadSockWErr(sock,msg,log_msg,SL,log_pr=3,exit=0):
  data=sock.read(msg)
  if string.find(data,"-ERR") != -1:
    if log_pr & 1: SL.Log(0, log_msg+"%s" % data)
    if log_pr & 2: print  log_msg+"%s" % data
    if exit: sys.exit(1)
  return data
# -------------------------------------
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
# -------------------------------------
def ReadPipe(cmd,end_string=None):
  """Reads stdin from a command and waits until end_str appears
     or the pipe is closed.
  """
  fp=os.popen(cmd)

  r=''
  while 1:
    a=fp.readline(-1)
    #print a
    r+=a
    if end_string:
      if a.find(end_string)!=-1:break
    if a==None: break
    if a=="": break
  fp.close()
  return r
# -----------------------------------------
def ReadPipe2(cmd,end_string=None,tmo=None):
  """Reads stdin & stderr from a command and waits until end_str appears
     or the pipe is closed.
  """
  fp=popen2.popen3(cmd)
  b=''

  s=select.select([fp[0],fp[2]],[],[],tmo)
  r=['','']
  #print s
  for i in range(len(s[0])):
     while 1:
       a=s[0][i].readline()
       #print "*",a
       if s[0][i]==fp[0]:r[0]+=a
       if s[0][i]==fp[2]:r[1]+=a
       if end_string:
         if a.find(end_string)!=-1: break #return list(r)
       if a==None: break
       if a=="": break
     #print i,s[0][i],'finished'
     #print '***',i,
     s[0][i].close()
  return list(r)
# -----------------------------------------
def ReadFile(name,end_string=None,tmo=None):
  """Reads fa file and waits until end_str appears
     or timout appears.
  """
  try: fp=open(name,'r+')
  except: return (1,"Error opening %s" % name)
  r=''
  t=None
  while 1:
     a=fp.readline()
     if a == '' or a==None: 
       if tmo:
         if t==None: t=time.time()
         if time.time()-t > tmo:
            fp.close()
            return (2,r)
	    break  
     else:
       t=None
       r+=a
       if end_string:
          if type(end_string)==StringType:
	     if a.find(end_string)!=-1:break
          elif type(end_string)==ListType:
	     for i in end_string:
	        #print a,i
	        if a.find(i)!=-1:fp.close();return(0,r)
	  else: 
	     return (3,"end_string must be a list or string")
  fp.close()
  return (0,r)
# -----------------------------------------
def LinTrafo(pmin_max,dmin_max, d):
    diffz= pmin_max[1]-pmin_max[0]
    diffn= dmin_max[1]-dmin_max[0]
    return  (diffz/diffn*d + (pmin_max[1]-diffz/diffn*dmin_max[1]))
# -------------------------------------
def BackTrafo(pmin_max,bmin_max, b):
    diffz= pmin_max[1]-pmin_max[0]
    diffn= bmin_max[1]-bmin_max[0]
    return  ((b-(pmin_max[1]-(diffz/diffn)*bmin_max[1]))*diffn/diffz)
# -------------------------------------
def Horner(x,a):
    """ Calculates sum(a[i]*x**i) using Horners algorithm
        [a0,a1,a2,...]
	-5*x*x+2*x+1=(-5*x+2)x+1
    """
    s=a[-1]
    for i in range(len(a)-2,-1,-1):
     s=s*x+a[i]
    return s
# ------------------------------------
def LorentzFn(x,p):
    """Calculates Lorentian function:
       
       Sum of Lorentz functions: BG + A / { 1 + [ 2* (x - x0) / HBW]**2 } 
         A: Amplitude
        x0: Position of peak (maximum)
       HBW: Half band width
        BG: Background
        nP: number of peaks
	
      p[0]: Amplitude[0] -> [n*3]
      p[1]: x0[0]
      p[2]: HBW[0]
      p[3]: Amlidude[1]
      p[4]: x0[1]
      p[5]: HBW[1]
         ... 
      p[3*nP] Background
    """
    k=len(p)-1
    if k%3:
       raise ValueError, "Number of parameters (%d) must be 3*n+1" % len(p)
       
    s=0
    for i in range(k/3):
        j=3*i
        s+=( p[j] / (1. + (2. * ( x - p[j+1])/p[j+2])**2 ) )
    return s+p[-1]
# -------------------------------------

