from time import *
from math import *
import types,sys

#$Log: RoundRobin.py,v $
#Revision 1.3  2006/07/11 11:57:52  herbie
#Print ctrl added
#
#Revision 1.2  2005/12/19 09:13:05  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 08:57:36  herbie
#Initial revision
#
CVS_ID="$Id: RoundRobin.py,v 1.3 2006/07/11 11:57:52 herbie Exp herbie $"

__all__=['RoundRobin']

class RoundRobin(list):
   """ Defines a simple round robin data archive (cylic list)
       holding a time value and corresponding a data value
       Some simple basic operations can be performed:
       __init__(self,max_size=3)
       append(self,item,t=None)
       *Insert(self,item,t=None) -> same as append(...)
       *GetSize(self) -> same as len(...)
       GetRate(self,last_n=None)
       When(self,value,last_n=None)

   """
   def __init__(self,max_size=3):
     """ __init__(self,max_size)
        Initializes the data grave
        Parameters:
           max_size: Number of data pairs after that cycling starts
                     must be greater than 3
        Errors:
           ValueError: raised if max_size <= 3
     """
     if max_size < 3: raise ValueError('max_size > 3 required')
     list.__init__(self)
     self.max_size=max_size

   def append(self,item,t=None):
     """ append(self,item,t=None)
         Insert a new data pair in the archive
         Parameters:
	        item: either a data value
                      or a tuple: (time,value)
                   t: corresponding time value
                      if omitted UNIX time stamp is used
                      Consecutive time values must be in growing order
              Error:
         ValueError: raised if t(n+1) <= t(n) (n: order of insert)
     """
     if not type(item)==types.TupleType:
      if t == None: t=time()
      item=(t,item)
     if len(self) and item[0]<self[-1][0]:
       raise ValueError('Time must grow %f < %f' % (item[0],self[-1][0]))
     list.append(self,item)
     if len(self) > self.max_size: self.pop(0)

   Insert=append
   
   GetSize=list.__len__

   def Print(self,last_n=None,fh=None,ctrl=None):
     """Print(self,last_n=None,ctrl=None)
        Prints time/values data pairs
           last_n: The number of last data points to be printed (None: all)
               fh: File handle for printing (None: stdout) or string
	     ctrl: Controls output
	           None: First point has time zero
		      1: Time in UNIX timestamp (secs since 1.1.1970)
		      2: Last value has time zero; others negative
     """
     if last_n == None: last_n=len(self)#-1
     if last_n < 1 or last_n > len(self):
        raise IndexError, 'last_n must be 1 - %d (%d)' % (len(self)-1,last_n)
        return 0
     if fh==None: fh=sys.stdout
     if type(fh) == types.StringType: r=''
     else: r=None
     for i in range(len(self)-last_n,len(self)):
       if not ctrl:
           a="%.2f  %f\n" % (self[i][0]-self[0][0],self[i][1])
       elif ctrl == 1:
           a="%.2f  %f\n" % (self[i][0],self[i][1])
       else: 
           a="%.2f  %f\n" % (self[i][0]-self[-1][0],self[i][1])
           
       if type(fh) == types.StringType: r+=a
       else: fh.write(a)
     return r

   def GetRate(self,last_n=None):
     """ GetRate(self,last_n=None)
         Returns a mean value of the time-derivative of th data values (dv/dt)
         Parameters:
           last_n: The number of last data points over that
                   the mean value of Delta v/Delta t is calculated
                   If omitted the whole archive size is used
             Return: Rate (derivative) in 1/seconds
             Errors:
               IndexError: raised if last_n not between 1 and self.size-1 
     """
     if last_n == None: last_n=len(self)-1
     if last_n < 1 or last_n >= len(self):
        raise IndexError, 'last_n must be 1 - %d' % (len(self)-1)
        return 0
     sum=0
     for i in range(len(self)-last_n,len(self)):
        delta_t=self[i][0]-self[i-1][0]  
        sum=sum+(self[i][1]-self[i-1][1])/delta_t  
     return sum/last_n

   def When(self,value,last_n=None):
     """ When(self,value,last_n=None)
	 Calculates the time when a data value will be reached
	 (linear extra- / interpolation)
	 Parameters:
	   value: Value to be reached
	  last_n: The number of last data points over that
		  will be extra- / interpolated
		  If omitted the whole archive size is used
	 Return: A tuple; first: Starting time of archive
			 second: Calculated time when value will
				 (has) be(en) reached
     """
     if last_n == None: last_n=len(self)-1
     k=self.GetRate(last_n)
     if k == 0:
       raise ZeroDivisionError('rate is zero')
       return (0,0)
     #y=k.x + d; d=y-k.x  
     d=self[len(self)-1][1] - k*self[len(self)-1][0]
     #x=(y-d)/k
     return (self[0][0],(value-d)/k)

   def At(self,tim):
     """ At(self,tim)
	 Calculates the value at the time specified
	 (linear extra- / interpolation)
	 Parameters:
	     tim: time for that the value will be calculated
	  last_n: The number of last data points over that
		  will be extra- / interpolated
		  If omitted the whole archive size is used
	 Return: Tuple: Interpolated value, slope of interpolation
     """
     last_n = len(self)
     if last_n < 2:
        raise ValueError, 'Less than 2 values in archive'
        return (0,0)
     if tim < self[0][0] or tim > self[len(self)-1][0]:
        raise ValueError, 'Value %d not in archive (%d ... %d)' % (tim,self[0][0],self[len(self)-1][0])
        return (0,0)
     
     for i in self:
        if i[0] == tim: return (i[1],0)
        if i[0] > tim:
           k=(i[1]-v[1])/(i[0]-v[0])
           d=i[1]-k*i[0]
           return (k*tim+d,k)
        v=i
     raise ValueError, '%d greater than last in archive' % tim
     return (0,0)

   def Mean(self,last_n=None):
     """ Mean(self,last_n=None)
	 Calculates mean value and mean deviation (sigma)
	 Parameters:
	  last_n: The number of last data points over that
		  the mean value is calculated
		  If omitted the whole archive size is used
	 Return: A tuple; first: mean value
			 second: mean deviation
     """
     if last_n == None: last_n=len(self)
     if last_n < 2 or last_n > len(self):
        raise IndexError, 'last_n must be 2 - %d' % (len(self))
        return (0,0)
     sum=0
     for i in range(len(self)-last_n,len(self)):
        sum=sum+self[i][1]
     m=float(sum)/last_n
     sum=0
     for i in range(len(self)-last_n,len(self)):
        sum=sum+(m-self[i][1])**2
     return (m,sqrt(sum/float(last_n-1)))

   def MaxDev(self,last_n=None):
     """ MaxDev(self,last_n=None)
	 Calculates max deviation
	 Parameters:
	  last_n: The number of last data points over that
		  the maximal deviation is calculated
		  If omitted the whole archive size is used
	 Return: A tuple; first: max-min
			 second: max
			  third: min
     """
     if last_n == None: last_n=len(self)
     if last_n < 2 or last_n > len(self):
        raise IndexError, 'last_n must be 2 - %d' % (len(self))
        return (0,0)
     maxv=self[len(self)-last_n][1]
     minv=maxv
     for i in range(len(self)-last_n,len(self)):
      if self[i][1] > maxv : maxv=self[i][1]
      if self[i][1] < minv : minv=self[i][1]
     return (maxv-minv,maxv,minv)
     
if __name__=='__main__':

 from whrandom import *
#
# R=RoundRobin(10)
# for i in range(10):
#   R.Insert(10*i+random(),i)
# print R.GetRate(5)
# print R.When(0)
# print R
# R=RoundRobin(10000)
# for i in range(10000):
#   R.append(10+uniform(-1.,1.))
# print R[-2:]
# print R.Mean()
 R=RoundRobin(10000)
 for i in range(100):
   R.append(i/10.,i)
# print R.MaxDev()
# print R.MaxDev(100)
 print R.At(98)
 print R.At(98)
 print R.At(100)
# R=RoundRobin(10)
# for i in range(10):
#   R.Insert(200-10*i+random())
#   sleep(0.2)
# print R.GetRate(5)
# e=R.When(0)
# print ctime(e[0]),ctime(e[1])
# print R 
