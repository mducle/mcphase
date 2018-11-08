from math import *
from decimal import *

#$Log: axis.py,v $
#Revision 1.2  2006/04/26 12:25:35  herbie
##.#e+-## representation improved
#
#Revision 1.1  2005/12/15 08:56:43  herbie
#Initial revision
#
CVS_ID="$Id: axis.py,v 1.2 2006/04/26 12:25:35 herbie Exp herbie $"

class axis:
   """ Scales an axis
   """
   def __init__(self,min_max,nsig,len_chr,typ):
      """Automatic scale of axis in decades 1-2-5
         Parameters:
            min_max: touple (min,max) of values (float)
               nsig: number significant digits
            len_chr: length of axis in characters (int)
        	typ: type of axis: 'x' or 'y'
      """
      self.delta=(min_max[1]-min_max[0])
      if self.delta <=0:
        raise ValueError('Min (%f) >= Max (%f)' % (min_max[0],min_max[1]))

      if nsig <= 0: raise ValueError("nsig must be > 0")
      self.nsig=nsig
      if len_chr <= 0: raise ValueError("len_chr must be > 0")
      self.len_chr=len_chr

      if not (typ=='y' or typ == 'x'):
        raise ValueError('Typ (%s) must be "x" or "y"' % typ)
	
      self.exp=0
      self.scale=None
      lgDelta=[None,None]
   

      if typ == 'x':
         min_dig = min(self.left_digits(min_max[1])[0],self.left_digits(min_max[0])[0])
         ntic_max=len_chr/(min_dig+3)
#         print "min_dig: %f" % min_dig
#         print "ntic_max: %f" % ntic_max

         delta_min=float(self.delta)/ntic_max
#         print "delta_min: %f" % delta_min
         lgDelta[1]=log10(self.delta)
         lgDelta[0]=log10(delta_min)

         s=[1,2,5]
         pDelta=[]    # possible deltas
#         print lgDelta
         for i in range(int(floor(lgDelta[0])),int(floor(lgDelta[1]))+1):
#            print i
            for j in s:
        	t=j*10**i
        	if t > delta_min and t<self.delta: pDelta.append(t)
#         print pDelta
  
         pMM=[]    # possible [delta, min, max]
         for i in pDelta:
#         print i,float(min_max[0])/i,float(min_max[1])/i
            pMM.append((i,floor(float(min_max[0])/i)*i,ceil(float(min_max[1])/i)*i))

         for i in pMM:
             steps=int((i[2]-i[1])/i[0])
#             print i,steps
             all_d=0
             for j in range(steps):
        	 v=i[1]+j*i[0]
#       	  print self.all_digits(v,i[0],nsig)
        	 all_d+=self.all_digits(v)[0]+2
#             print all_d,len_chr
             if all_d <= len_chr: self.scale=self.scale_list(i); return
         raise ValueError('No scale: x-axis')

      if typ == 'y':
         ntic_max=len_chr/2
#         print "ntic_max: %f, len_chr:%f" % (ntic_max,len_chr)
         delta_min=float(self.delta)/ntic_max
#         print "delta_min: %f, delta: %f" % (delta_min,self.delta)
#         print min_max
         lgDelta[1]=log10(self.delta)
         lgDelta[0]=log10(delta_min)

         s=[1,2,5]
#         pDelta=[]    # possible deltas
#        print lgDelta
         Delta=self.delta
         brk=0
#         print "Range:",int(floor(lgDelta[0])),int(floor(lgDelta[1]))+1
         for i in range(int(floor(lgDelta[0])),int(floor(lgDelta[1]))+1):
#            print i
            for j in s:
        	t=j*10**i
#                print t
        	if t > delta_min and t<self.delta: Delta=t; brk=1; break
            if brk: break

         pMM=(Delta,floor(float(min_max[0])/Delta)*Delta,ceil(float(min_max[1])/Delta)*Delta)

	 if pMM[2]>0:
	    if pMM[2] <= 1.e-3 or pMM[2] >= 1.e4: self.exp=1

	 if pMM[1]<0:
	    if fabs(pMM[1]) <= 1.e-3 or fabs(pMM[1]) >= 1.e4: self.exp=1
         #print pMM,self.exp
         self.scale=self.scale_list(pMM)
#         return 'No scale y'


   def left_digits(self,v):
       """ Calculates max. number of digits left from decimal point
           including sign.
	   parameters:
	    v: value (float)
	   return: tuple (number of left digits, 'f':float or 'e':expontial)
       """
       if v == 0: return (1,'f')
       l=fabs(log10(fabs(v)))
       f='f' 
       if l > self.nsig or self.exp: r=1;f='e'
       else: r=int(l+1)
       if v<0:r+=1
       return (r,f)

   def right_digits(self,v,exp=0):
       """ Calculates max. number of digits right from decimal point
	   parameters:
	    v: value (float)
	   return: int: number of right digits
       """
       if v == 0: return 0
       l=fabs(log10(fabs(v))) 
       #if l > self.nsig: return self.nsig-1
       t=log10(fabs(self.delta))
       r=0
       if l > self.nsig or exp:
         d=self.delta*10**-floor(t)
         t=log10(fabs(d))
         #print v,d,t,self.delta
	 r=int(log10((self.max-self.min)/self.delta)+.3)
         #print v,d,t,(self.max-self.min)/self.delta,r
       #print v,t
       if t<0: return int(ceil(-t))
       else: return r

   def all_digits(self,v,exp=0):
       """ Calculates max. number of digits including decimal point
	   parameters:
	    v: value (float)
	   return: tuple (number of all digits including dp, 'f' or 'e' )
       """
       r=self.right_digits(v,exp)
       if r: r+=1  # decimal point
       rr=self.left_digits(v)
       r+=rr[0]
       if v:
         if fabs(log10(fabs(v))) > self.nsig or exp: r+=4 # exponential repr 1.234E-05
       return (r,rr[1])

   def scale_list(self,a):
     """ Makes a list containing the scale
         parameters:
	    a: touple containing (delta,min,max) of the scale
	 return:
	    list containig tuples for for each scale tic
	    (number (string), format (%5.2d), value (float))
     """
     self.max=a[2]
     self.min=a[1]
     self.delta=a[0]
#     steps=int((a[2]-a[1])/a[0])
     steps=int( (Decimal(self.make_number(a[2]))-Decimal(self.make_number(a[1])))/\
                 Decimal(self.make_number(a[0])))
     r=[]
     for j in range(steps+1):
         v=a[1]+j*a[0]
         rd=self.right_digits(v,self.exp)
         ld=self.left_digits(v)
         ad=self.all_digits(v,self.exp)
         f="%%%d.%d%s" % (ad[0],rd,ad[1])
	 #print v,rd,ld,ad,f
         r.append((self.make_number(v,self.exp),f,v))
     return list(r)

   def make_number(self,v,exp=0):
     """ Makes a number (string) from a value (float)
         containing the correct numbers of digits left/right from decimal point
     """
     rd=self.right_digits(v,exp)
     ld=self.left_digits(v)
     ad=self.all_digits(v,exp)
     f="%%%d.%d%s" % (ad[0],rd,ad[1])
     return (f % v).strip()
   
if __name__=='__main__':
# print scale_axis((12345678,12345699),7,30,'y')
# a=axis((2,222),7,120,'x')
# a=axis((-12,0.15),7,120,'x')
# a=axis((0.1,0.15),7,120,'x')
# a=axis((-120.,-18),7,120,'x')
# a=axis((1234567,1234568),7,120,'x')
# a=axis((12345678,12345699),7,120,'x') # nicht so gut
# a=axis((2,222),7,20,'y')
 a=axis(( 4.693037, 4.733559),7,27,'y')
 print a.scale,a.min,a.max,a.delta
 a=axis(( 20.277, 290.01),7,72,'x')
 print a.scale,a.min,a.max,a.delta

