from types import *
from array import *
from stdfunc import *
from spline import *
import string,sys
from math import *
from decimal import *
#$Log: xydata.py,v $
#Revision 1.2  2005/12/19 09:13:05  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 08:57:36  herbie
#Initial revision
#
CVS_ID="$Id: xydata.py,v 1.2 2005/12/19 09:13:05 herbie Exp $"

class XYData:
 """ class XYData

     Represents x/y data pairs by two lists (columns) with the same count of float numbers.
     Each data pair can have a 'weight'; this is an integer number.
     This 'weight' can be used to control operations. (Default: Running index; first=1)
     E.g. If you sort the list, the weight is not sorted. so the original order is not lost.
     The data pairs can be interpreted as y=f(x) and some standard operations can be
     performed.

 """
 def __init__(self,x,y):
    """ Initializes x and y data.
        Parameters: x: x values (list of float)
	            y: y values (list of float); same size as x
        
    """
    if not (type(x) == ListType and type(y) == ListType and len(x) == len(y)):
       raise TypeError('x and y must be lists with the same size')
#    self.length=len(x)
    
    self.x=[]
    self.y=[]
    for i in range(0,len(x)):
	if not (type(x[i]) == FloatType and type(y[i])==FloatType):
	   raise TypeError(i,x[i],y[i],'x and y must be lists of floats')
	self.x.append((x[i],i+1))
	self.y.append(y[i])
#    for i, v in x.enumerate():
#	if not (type(v) == FloatType and type(y[i])==FloatType):
#	   raise TypeError(i,v,y[i],'x and y must be lists of floats')
#	self.x.append((v,i))
#	self.y.append(y[i])
# --------------------------------------------------
 def __len__(self):
     """ Gives the length of the data pair list"""
     return len(self.x)
# --------------------------------------------------
 def X(self):
     """ Returns a list with the x-values"""
     r=[]
     for i in self.x: r.append(i[0])
     return r
# --------------------------------------------------
 def Weight(self):
     """ Returns a tuple:
         First: Weight of the data pairs as a list
	 Second: Count where weight > 0
     """
     r=[]
     c=0
     for i in self.x:
       r.append(i[1])
       if i[1] :c+=1
     return (r,c) 
# --------------------------------------------------
 def Xmax(self,start=0,end=None):
     """ Gives the maximum of the x-values between start and end.
         If start omitted starting from first value.
	 If end omitted end with lats value.
     """
     if end==None: end=len(self.x)
     return max(self.x[start:end])
# --------------------------------------------------
 def Ymax(self,start=0,end=None):
     """ Gives the maximum of the y-values between start and end.
         If start omitted starting from first value.
	 If end omitted end with last value.
     """
     if end==None: end=len(self.y)
     return max(self.y[start:end])
# --------------------------------------------------
 def Xmin(self,start=0,end=None):
     """ Gives the minimum of the x-values between start and end.
         If start omitted starting from first value.
	 If end omitted end with last value.
     """
     if end==None: end=len(self.x)
     return min(self.x[start:end])
# --------------------------------------------------
 def Ymin(self,start=0,end=None):
     """ Gives the minimum of the y-values between start and end.
         If start omitted starting from first value.
	 If end omitted end with last value.
     """
     if end==None: end=len(self.y)
     return min(self.y[start:end])
# --------------------------------------------------
 def pop(self,i):
     """remove and return xy-data-pair at index i
        return: [x[i],y[i]]
     """
     r1=self.x.pop(i) 
     r2=self.y.pop(i)
#     self.length-=1
     return [r1,r2]
# --------------------------------------------------
 def Write(self,fp=sys.stdout,ctrl=-1):
    """Write data pairs to file.
       Parameters:
	     fp: File to write data to. default: stdout
	   ctrl: Integer to control output. Values must be or-ed.
	         1: Write running index.
		 2: Write x and y
		 4: Write weight
		 default: Write all
    """
    for i in range(0,len(self.x)):
        if ctrl & 1: fp.write("%d:" % i)
        if ctrl & 4: fp.write(" %d" % self.x[i][1])
        if ctrl & 2: fp.write(" %14.7g" % self.x[i][0])
        if ctrl & 2: fp.write(" %14.7g" % self.y[i])
        fp.write('\n')
# --------------------------------------------------
 def Sort(self):
    """ Sort data pairs by ascending x values.
        Weight is not changed.
    """
    self.x.sort()
    t=list(self.y)
    for i in range(0,len(self.x)):
          self.y[i]=t[self.x[i][1]-1]
# --------------------------------------------------           
 def LinIntpol_sorted(self,x):
    """ Calculate the corresponding y value using linear interpolation
        between the neigbour ppints.
	Asume that the data pairs are sorted by ascending x-values.
	Return: y=f(x)
    """
    if len(self.x) < 2:
      raise IndexError('More than 2 points needed')
    iLo = 0
    iHi = len(self.x)-1
    Outside = 0;
    if x > self.x[len(self.x)-1][0]:
      iHi=len(self.x)-1
      iLo=len(self.x)-2
      Outside=1
    if x < self.x[0][0]:
      iHi=1
      iLo=0
      Outside=1
    if not Outside:
      while iHi-iLo>1:
        k=(iHi+iLo) >> 1
        if  self.x[k][0] > x: iHi=k
        else: iLo=k
  
    if(iHi==iLo): raise ValueError('iHi == iLo')

    yk =( self.y[iHi] - self.y[iLo] ) /\
	   ( self.x[iHi][0] - self.x[iLo][0] )*(x - self.x[iLo][0])
    if Outside:
      raise UserWarning('x: %f outside table' % x)
    return (self.y[iLo]+yk)
# --------------------------------------------------
 def LinIntpol(self,x):
    """ Calculate the corresponding y value using linear interpolation
        between the neigbour ppints.
	Asume that the data pairs are NOT sorted by ascending x-values.
	(Slower than LinIntpol_sorted)
	Return: y=f(x)
    """
    if len(self.x) < 2:
      raise IndexError('More than 2 points needed')
    iLo = -1
    iHi = -1
    for k in range(0,len(self.x)):
       d=x-self.x[k][0]
       if d==0: return self.y[k]
       if d < 0:
	 if iHi==-1 or fabs(d) < fabs(dminm): dminm=d; iHi=k; continue

       if d > 0:
	 if iLo==-1 or fabs(d) < fabs(dminp): dminp=d; iLo=k; continue

    if iHi == iLo: raise ValueError('Value %14.7g not in table' % x)
    if iHi==-1 or iLo==-1: raise ValueError('Value %14.7g not in table' % x)

    yk =( self.y[iHi] - self.y[iLo] ) /\
	   ( self.x[iHi][0] - self.x[iLo][0] )*(x - self.x[iLo][0])

    return (self.y[iLo]+yk)
# --------------------------------------------------
 def MathOper(self,S,Sy):
    """ Perform a basic mathematical operation y1(x1) = y1(x1) <+-*/> y2(x2).
	x1 and x2 need not to be the same values, linear interpolation is used to get
        the value y2(x1) from y2(x2)
	Parameters:
	     S: XYData object representing y2(x2)
	    Sy: Operation (one of +-*/) (one character string)
	Result: self=self <+-*/> S
	        values of x2 which are not within x1 are skipped.
	Return: Number of skipped points.
    """
    if not (Sy=='+' or Sy=='-' or Sy=='*' or Sy=='/'):
       raise ValueError('Illegal operation char <%s>; must be "+*-/"'%Sy)

    if len(S)<3:
       raise ValueError('More than 3 points needed')
    
    if S.Xmin()[0] > self.Xmax()[0] or S.Xmax()[0] < self.Xmin()[0]:
       raise ValueError("x ranges do not overlap")

    sc=0
    nx=[]
    ny=[]
    for i in range(0,len(self.x)):
        try: yn=S.LinIntpol(self.x[i][0])
        except ValueError:
          sc+=1
          nx.append((self.x[i][0],0))
          ny.append(self.y[i])
          continue
        nx.append(self.x[i])
        if Sy == '+': ny.append(self.y[i]+yn); continue
        if Sy == '-': ny.append(self.y[i]-yn); continue
        if Sy == '*': ny.append(self.y[i]*yn); continue
        if Sy == '/': 
           if yn==0: ny.append(_HUGE_VAL)
           else: ny.append(self.y[i]/yn)
#           continue
    self.x=list(nx)
    self.y=list(ny)
#    self.length=len(nx)
    return sc
# --------------------------------------------------
 def Deriv(self,nP):
    """  dy/dx is build by an algorithm using n neigbour points.
         Values will be sorted ascending by x.
    """
    if nP>len(self.x)-3 or nP<1:
     raise ValueError("Invalid number of neighbours %d (1 ... %d)" %\
                       (nP,len(self.x)-3))
    self.Sort()
    nErrs=0
    nx=[]
    ny=[]
    for i in range(0,len(self.x)-nP):
       Dx=0.;Dy=0.
       for j in range (i+nP,i,-1):
          Dx+=self.x[j][0] - self.x[i][0]
          Dy+=self.y[j] - self.y[i]
#       print Dx,Dy
       Dx/=nP
       Dy/=nP
       if Dx==0: nErrs+=1; continue
       else: nx.append(self.x[i]); ny.append(Dy/Dx)
    self.x=list(nx)
    self.y=list(ny)
#    self.length=len(nx)
    return nErrs
# --------------------------------------------------
 def Integrate(self):
    """  Intgral y.dx is build by 0.5('Obersumme' + 'Untersumme').
         Values will be sorted ascending by x.
    """
    self.Sort()
    ny=[]
    SumI = 0
    for i in range(1,len(self.x)):
      Dx = self.x[i][0]-self.x[i-1][0]
      y = 0.5*(self.y[i]+self.y[i-1])*Dx
      SumI+=y
      ny.append(SumI)
    self.pop(0)
    self.y=list(ny)
# --------------------------------------------------
 def ZeroShift(self,x0):
    """  Shift the zero point of y(x) to x0 by substracting y(x0) from all y values.
    """
    if x0 < self.Xmin()[0] or x0 > self.Xmax()[0]:
       raise ValueError("x value %14.7g out of range" % x0)
    y0=self.LinIntpol(x0)
    for i in range(0,len(self.x)): self.y[i]-=y0
# --------------------------------------------------
 def DeltaL_l(self,l0):
   """  Calculate Delta L/l0 -> (y-l0)/l0.
   """
   if(l0 == 0.0): raise ZeroDivisionError
   lmin=self.y[self.x.index(self.Xmin())]
   for i in range(0,len(self.x)):
       self.y[i]=(self.y[i]-lmin)/l0
# --------------------------------------------------
 def CalcGap(self,CDiam,GGap=None):
   """Calculates gap of a parallel plate disk condensor
      y-colummn is assumed to be the capacity in pF at input,
      gap in [mm]!! at output.
      x-column not relevant/changed
      Paramters: 
        CDiam:Cell (condesor plate) diameter in [cm]
         GGap: Gap between capacitor plate and ground in [cm]
               None if not not used
   """
   CKonst = 10.*0.0885*0.25*CDiam*CDiam*pi*1e-3
   if CDiam<=0:
      raise ValueError('Diameter must be greater than zero')
   if GGap == None or GGap < 0:
     for i in range(0,len(self.y)): self.y[i]=CKonst/self.y[i];
   else:
     r = 0.5*CDiam
     for i in range(0,len(self.y)):
       a=self.y[i]/0.0885
       b=0.22/0.0885*self.y[i]*GGap - r*r*pi - r*GGap*pi*(1.+GGap/CDiam)
       c=-0.22*r*r*pi*GGap
       d=b*b-4.*a*c
       l1=(-b+sqrt(d))/2/a
       l2=(-b-sqrt(d))/2/a
       if l1<0 and l2<0: self.y[i]=0.
       elif l1<0: self.y[i]=l2*10.
       elif l2<0: self.y[i]=l1*10.
# --------------------------------------------------
 def LookTable(self,t,t1):
   """Looks for all x-values in spline table t, t1 and replaces
      the y-values with lookup values
       t: List of x-values for lookup
      t1: List of y-values for lookup
   """  
   S=SplineTable(t,t1)
   S.Spline()
   for i in range(0,len(self.x)):
       try: ny=S.Intpol(self.x[i][0])
       except ValueError: 
         self.x[i]=(self.x[i][0],0)
       else: self.y[i]=ny
# --------------------------------------------------
 def MeanBin(self,delta_x):
   """Put all values within delta_x in a bin.
      Caluclate mean values and deviations for x and y within this bin
      List will be sorted.
      Return: List containing ( mean(x(bin)), mean(y(bin)), sigma(x(bin)), sigma(y(bin), len(bin)) )
   """  
   self.Sort()
   xmin=self.Xmin()[0]
   xmax=self.Xmax()[0]
   dd=Decimal(("%g" % (delta_x/2.)).strip())
   d=0.5*delta_x
   imin=int(Decimal(("%g" % xmin).strip())/dd)
   imax=int(Decimal(("%g" % (xmax)).strip())/dd)
   xnew=[(i*d,d*(i-1),d*(i+1)) for i in range(imin,imax+1,2)]
#   print xnew
   jstart=0
   bin=[]
   out=[]
   for i in xnew:
       for j in range(jstart,len(self.x)):
           if  self.x[j][0]<i[2]:
#	      print i[0],self.x[j][0]
	      bin.append( (self.x[j][0],self.y[j],))
	   else:
	      mv=[0.,0.]
	      for k in bin:mv[0]+=k[0];mv[1]+=k[1]
	      l=len(bin)
              if l!=0:
	         mv[0]/=l
	         mv[1]/=l
	         sigma=[0.,0.]
		 if l>1:
	            for k in bin:sigma[0]+=(mv[0]-k[0])**2;sigma[1]+=(mv[1]-k[1])**2
	            sigma[0]=sqrt(1./(l-1)*sigma[0])
	            sigma[1]=sqrt(1./(l-1)*sigma[1])
              if l>0: out.append( (mv[0],mv[1],sigma[0],sigma[1],l) )   
	      bin=[]
	      jstart=j 
	      break
   return out
# --------------------------------------------------
if __name__=='__main__':
 X=XYData([5.,2.,6.,4.,1.,3.],[25.,4.,36.,16.,1.,9.])
 Y=XYData([1.5,2.5,3.5,4.5],[2.25, 6.25, 12.25, 20.25])
 Z=XYData([1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,6.],[1.,2.25,4.,6.25,9.,12.25,16., 20.25,25.,36.])
 A=XYData([1.,1.1,1.5,1.6,1.8,2.1,2.5,3.,3.3,3.35,3.5,4.,4.5,4.55,4.66],[1.,1.1,1.5,1.6,1.8,2.1,2.5,3.,3.3,3.35,3.5,4.,4.5,4.55,4.66])
# X.Sort()
# Z.Write()
# print X.MathOper(Y,'/')
# print Z.Deriv(8)
# Z.DeltaL_l(1)
# print
# Z.Write()
#---
# print Z.Xmin(),Z.Xmax()
# print Z.Xmin(0,10),Z.Xmax(0,10)
# print Z.x[9]
#---
 print A.MeanBin(1)
