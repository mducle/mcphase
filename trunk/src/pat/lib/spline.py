import sys,types,string
from pprint import *

#$Log: spline.py,v $
#Revision 1.1  2004/02/24 14:19:07  herbie
#Initial version
#
CVS_ID="$Id: spline.py,v 1.1 2004/02/24 14:19:07 herbie Exp $"

class SplineTable:
  """ Class for cubic spline interpolation
      Reads a table (file or list of x,y-tuples)
      Layout of the file:
        1st line: SPLINETABLE any text
        2nd line: number of following x,y-tuples
      next lines: x/y-tuples (at least 3 tuples)

      Methods:
      __init__(self,table)
     ReadFile(self,filename)
     Spline(self,D0=None, Dn=None)
     Intpol(self,x)
  """
  def __init__(self,table,t2=None):
     """ Initalizes the lookup
         Reads the table, calculates 2nd derivative
         Parameters:
           table:  Either a filename, a list of tuples or a list representing the table
              t2:  If table is a list it gives the x-values, then
                   t2 is also a list representing the y-values
         Errors:
           TypeError: If the table is neither a filename or a list of tuples
          ValueError: If the file is no a spline file according to the specs above
     """
     self.table=[]
     tp1=type(table)
     tp2=type(t2)
     if type(table) == types.StringType and t2 == None: # Its a filename
         self.ReadFile(table)
         self.Spline()
         return
     if type(table) == types.ListType and t2 == None: # Its a list
        if type(table[0]) == types.TupleType: # Its a list of tuples
           self.table=table
           self.table.sort()
           self.Spline()
           return
     if type(table) == types.ListType and type(t2) == types.ListType: # There are two lists
       if len(table) == len(t2): # There are two lists
           for i in range(0,len(table)):
             self.table.append((table[i],t2[i]))
           self.table.sort()
           self.Spline()
           return
#     print len(table),len(t2)
     raise TypeError('Input: filename or list of tuples or two lists (%s,%s)' % (tp1,tp2))

  def ReadFile(self,filename):
      """ ReadFile(self,filename)
          Reads the table from a file
          for internal use only!
          Errors:
          ValueError: If the file is no a spline file according to the specs above
      """
      self.buffer=None
      inp=open(filename,'r')
      self.buffer=inp.readlines()
      if string.find(self.buffer[0],'SPLINETABLE',0)==-1:
         raise ValueError(' %s is not a valid spline table' % filename)
      self.header=self.buffer[0]
      try:
        self.length=string.atof(self.buffer[1])
      except NameError:
        raise ValueError(' %s is not a valid spline table' % filename)
      self.buffer[0:2]=[]
      if len(self.buffer) < 3:
         raise ValueError(' %s spline table to small (>3)' % filename)
 
      for i in self.buffer:
          t=string.split(string.strip(i))
          if len(t) !=2 :
             raise ValueError(' %s is not a valid spline table' % filename)
          self.table.append((string.atof(t[0]) ,string.atof(t[1])))
          self.table.sort()

  def Spline(self,D0=None, Dn=None):
      """ Spline(self,D0=None, Dn=None)
          Calculates the 2nd derivative
          for internal use only!
          Parameters:
           D0,D1: curvature at both ends of the table
          Errors:
          ValueError: If D0,Dn are not floats or None
      """
      if not ( (type(D0) == types.FloatType or type(D0) == types.NoneType) and \
              (type(Dn) == types.FloatType or type(Dn) == types.NoneType)):
          raise ValueError('Bad input parameter D0 or Dn')
      u=[0]*len(self.table)
      self.d2=[0]*len(self.table)
      if D0 == None:
         self.d2[0]=0
         u[0]=0
      else:
         self.d2[0]=-0.5
         u[0]=(3.0/(self.table[1][0]-self.table[0][0])) *\
              ((self.table[1][1]-self.table[0][1])/(self.table[1][0]-self.table[0][0])-D0)
      for i in range(1,len(self.table)-1):
          sig=(self.table[i][0]-self.table[i-1][0])/(self.table[i+1][0]-self.table[i-1][0])
          p=sig*self.d2[i-1]+2.0;
          self.d2[i]=(sig-1.0)/p;
          u[i]=(self.table[i+1][1]-self.table[i][1])/(self.table[i+1][0]-self.table[i][0]) -\
               (self.table[i][1]-self.table[i-1][1])/(self.table[i][0]-self.table[i-1][0])
          u[i]=(6.0*u[i]/(self.table[i+1][0]-self.table[i-1][0])-sig*u[i-1])/p

      qn=un=0.0
      if not Dn == None:
         qn=0.5
         un=(3.0/(self.table[-1][0]-self.table[-2][0])) *\
            (Dn-self.table[-1][1]-self.table[-2][1])/(self.table[-1][0]-self.table[-2][0])

      self.d2[-1]=(un-qn*u[-2])/(qn*self.d2[-2]+1.0)

      for k in range(len(self.table)-2,0,-1):
         self.d2[k]=self.d2[k]*self.d2[k+1]+u[k]

  def __str__(self): return self.table.__str__()

  def Intpol(self,x):
     """ Performs the interpolation
          Parameters:
            x: x-value (1st column) for that the corresponding y-value (2nd column) is calculated
          Return:
            corresponding y-value
          Errors:
          ValueError: If x is less than first table value or greater than last table value
     """
     if x>self.table[-1][0]:
        raise ValueError('x value:%f greater than last in table:%f' % (x,self.table[-1][0]))
        return self.table[-1][1]
     if x<self.table[0][0]:
        raise ValueError('x value:%f less than first in table:%f' % (x,self.table[0][0]))
        return self.table[0][1]
    
     Lo=0
     Hi=len(self.table)-1

     while(Hi-Lo>1):
          k=(Hi+Lo) >> 1;
          if self.table[k][0] > x: Hi=k
          else: Lo=k
 
     h=self.table[Hi][0]-self.table[Lo][0]
     if h == 0.0: 
        raise ValueError('Delta x = 0')
        return None;
     a=(self.table[Hi][0]-x)/h;
     b=(x-self.table[Lo][0])/h;
     return a*self.table[Lo][1]+b*self.table[Hi][1]+\
         ((a*a*a-a)*self.d2[Lo]+(b*b*b-b)*self.d2[Hi])*(h*h)/6.0

if __name__=='__main__':

 print 'Usage:', sys.argv[0],'file value'
 S=SplineTable(sys.argv[1])
# print sys.argv[1],sys.argv[2] 
 print S.Intpol(string.atof(sys.argv[2]))
####
# L=[(0,0),(1,1),(2,4),(3,9),(4,16),(5,25),(10,100)]
# print 'Usage:', sys.argv[0],'value'
# S=SplineTable(L)
# print S.Intpol(string.atof(sys.argv[1]))

# L1=[0,1,2,3,4,5,10]
# L2=[0,1,4,9,16,25,100]
# S1=SplineTable(L1,L2)
# print S1
# print S1.Intpol(string.atof(sys.argv[1]))
