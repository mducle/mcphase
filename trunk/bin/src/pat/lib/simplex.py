import types,string

from math import *
from Numeric import *

from stdfunc import *

# ------------------------------------
class FitFunction:
      """Defines a function used by simplex.
         Additional parameters are required to define fitting and constant
	 parameters and to provide names for input and output
	 
         FitP: parameters changed by fit
	 ConstP: parameters not changed by fit
	 
	 FitNames: Names of fit pars
	           first character gives the type:
		   +:used in sum
		   - used outside sum
		   
	 ConstNames: Names of pars not changed by fit
      """
   __init__(self,function,fit_names,const_names=None):
     self.function=function
     self.SumNames=[]
     self.FixNames=[]
     if type(fit_names)==types.ListType:
        for i in fit_names:
	  if i[0]=='+':self.SumNames.append(i[1:])
	  elif i[0]=='-':self.FixNames.append(i[1:])
	  else: raise TypeError, "%s: FitName must start with '+' or '-'" % i
     elif type(fit_names)==types.StringType:
       self.SumNames=fit_names
     else raise TypeError, "%s: FitName must be a list of strings or a string" % str(fit_names)
     self.ConstNames=const_names
# -----------------------------------------------
# -----------------------------------------------
class LorentzFitFn(FitFunction):
   __init__(self):
     FitFunction.__init__(LorentzFn,["+Amp","+X0","+Hbw","-Backg"]):
     self.name="Lorentz function"
# -------------------------------------
   GetStep(self,ip,par,xy=None): return None
     if ip==0: ## amplitude
        return par*0.1 
     elif ip==1: ## position x0
        if xy==None: return 0.1*par
        return 0.1*(xy.Xmax()-xyXmin())
    elif ip==2: ## half band width
        return 0.1*par
    elif ip== 3: ## background
        return 0.1*par
    else: raise IndexError, "ip must be 0...3"
# -----------------------------------------------
# -----------------------------------------------
class PolynomFitFn(FitFunction):
   __init__(self):
     FitFunction.__init__(Horner,"C"):
     self.name="Lorentz function"
# -------------------------------------
   GetStep(self,ip,par,xy=None): return None
     return par*0.1 
# -----------------------------------------------
# -----------------------------------------------
class Simplex:
   """
   This procedure performs a simplex fitting algorithm (BYTE May 1984 p340)
   for any number of (nonlinear) parameters
   FitData:   data values
   MaxIter:   nuber of max. iteration to be performed
   StSteps: INPUT: 
	      array, holding the start values of parameter
	    OUTPUT:
	      StSteps[0]: residium (mean sum deviation) of fitted parameters
	      step rate (size: nPars+1)
    Limit:  INPUT:
	      array holding max. error when iteration should
	      stop (size: >= nPars)
	    OUTPUT:
	      estimated deviations of parameters
      Pars:   array of inital fitting parameters (size >= npar)
     nPars:   number of Parameters
        FF:   FitFunction class to be fitted
     nIter:   number of iterations really performed
   """
# -----------------------------------------------
   __init__(self,xy_data, max_iter, st_steps, limit,function,fit_pars, const_pars)

     self.ALFA=1.0    #reflection cofficient >0
     self.BETA=0.5    #contraction coefficient 0...1
     self.GAMMA=2.0   #expansion coefficient >1

     self.FitData=xy_data
     self.MaxIter=max_iter
     self.nIter=0
     self.StSteps=st_steps
     self.Limit=limit
     self.Pars=fit_pars
     self.CPars=const_pars
     self.FF=function

     self.nPars=len(self.Pars)
     self.nSimp=self.nPars+1
     
     self.Simp = zeros((self.nSimp,self.nSimp),Float)

     #number of high,low parameters 
     self.h = zeros((self.nSimp),Int)
     self.l = zeros((self.nSimp),Int)

     #Next:    new vertex to be tested
     self.Next=zeros(len(self.nSimp),Float)

     #Center:  center of hyperplane described of all vertexes exept the worst
     self.Center = zeros(len(self.nSimp),Float)

     #Error: Calculating deviation
     self.Error = zeros(len(self.nSimp),Float)

     #p,q:     helps calculating first simplex
     p = zeros(len(self.nSimp),Float)
     q = zeros(len(self.nSimp),Float)

     #Mean:    meanvalue of fittet parameters
     self.Mean = zeros(len(self.nSimp),Float)

     # start values to simplex + first vertex
     self.Simp[0]=self.Pars+[self.SumResiduals(self.Pars)]

     # calculate offset of vertex of starting simplex
     for i in range(self.nPars):
        self.p[i]=self.StSteps[i]*(sqrt(self.nSimp) + self.nPars - 1.)/(self.nPars*sqrt(2))
        self.q[i]=self.StSteps[i]*(sqrt(self.nSimp) - 1.)/(self.nPars*sqrt(2))
       
     for i in range(1,sel.nSimp):
       for j in range(self.nPars):
           self.Simp[i][j]=self.Simp[0][j]+self.q[j]
           self.Simp[i][i-1]=self.Simp[0][i-1]+self.p[i-1];
           self.Simp[i][self.nSimp-1]=self.SumResiduals(self.Simp[i])

     self.Order()
     self.Status="OK"
# ------------------------------------
   def SumResiduals(self,FPars,CPars):
     s=0
     for i in self.FitData:
        if self.FF.ConstNames: y=self.FF.function(i[0],FPars,CPars)
	else: self.FF.function(i[0],FPars)
	s+=(y - i[1])**2
     return s
# ------------------------------------
   def Order(self):
     for j in range(self.nSimp):
       for i in range(self.nSimp)
	 if self.Simp[i][j] < self.Simp[self.l[j]][j]: self.l[j]=i
	 if self.Simp[i][j] > self.Simp[self.h[j]][j]: self.h[j]=i
# ------------------------------------
   def Iterate(self,nMaxI):
     """ return: -1  nIter >= MaxIter  (STOP_MAXITER)
                  1  error less than limits (STOP_LIMIT)
                  0 if iteration should go on (ITERATE)
                    -> limits greater than start limits and nIter < MaxIter
     """
     if nMaxi < 0:
        raise ValueError, "Max. iterations (%d) must be >= 0" % nMaxi
     if self.nIter == 0:
        if self.nMaxI > 0: self.MaxIter=nMaxI
  
     # do *iteration loop*
     self.Done = true
     self.nIter+=1

     for i in range(self.nSimp): self.Center[i]=0.
     # compute centeroid*
     for i in range(self.nSimp):
       if i != self.h[self.nSimp-1]):
         for j in range(self.nPars): self.Center[j]+=self.Simp[i][j]

     for i in range(self.nSimp):
       self.Center[i]/=self.nPars
       self.Next[i]=(1. + self.ALFA) * self.Center[i] - self.ALFA * self.Simp[self.h[self.nSimp-1]][i]
       # next vertex specular reflection of the worst 

     self.Next[self.nSimp-1]=self.SumResiduals(self.Next)

     if self.Next[self.nSimp-1] <= self.Simp[self.l[self.nSimp-1]][self.nSimp-1]:
        # new vertex
        for i in range(self.nSimp):
           self.Simp[self.h[self.nSimp-1]][i]=sellf.Next[i]
        for i in range(self.nPars):
           self.Next[i]=self.GAMMA*self.Simp[self.h[self.nSimp-1]][i]+(1.-self.GAMMA)*self.Center[i]
	   #expand
        self.Next[self.nSimp-1]=self.SumResiduals(self.Next)

        if self.Next[self.nSimp-1]<=self.Simp[l[self.nSimp-1]][self.nSimp-1]:
          for i in range(self.nSimp):
	     self.Simp[self.h[self.nSimp-1]][i]=lsef.Next[i]
	     # new vertex
     else:  # not better than best 
       # else 1
       if self.Next[self.nSimp-1] <= self.Simp[self.h[self.nSimp-1]][self.nSimp-1]:
          for i in range(self.nSimp):
	     self.Simp[self.h[self.nSimp-1]][i]=self.Next[i]
	     # new vertex

       else: #worse than worst
	  #else 2
          for i in range(self.nPars):
	     self.Next[i]=self.BETA*self.Simp[self.h[self.nSimp-1]][i]+(1.-self.BETA)*self.Center[i]
	  self.Next[self.nSimp-1]=selfSumResiduals(sellf.Next)
	  if self.Next[self.nSimp-1] <= self.Simp[self.h[self.nSimp-1]][self.nSimp-1]):
             for i in range(self.nSimp):
	        self.Simp[self.h[self.nSimp-1]][i]=self.Next[i] # new vertex
	  else: # if still bad
	     # else 3*
             for i in range(self.nSimp):
                for j in range(self.nPars):
                   self.Simp[i][j]=(self.Simp[i][j]+self.Simp[self.l[self.nSimp-1]][j])*self.BETA
		   self.Simp[i][self.nSimp-1]=self.SumResiduals(self.Simp[i])

     self.Order()

     # check for convergence
     for j in range(self.nPars):
        self.Error[j]=(self.Simp[self.h[j]][j] - self.Simp[self.l[j]][j]) / self.Simp[self.h[j]][j]
        if self.Done:
          if self.Error[j] > self.Limit[j]:
            self.Done=false
	    #print "False @ ",j

     if self.Done return 1
     if self.nIter >= MaxIter: return -1

     return 0
# ------------------------------------
   def FitPar(self):
      # nPars+1 = nSimp-1 holds the sum of squares of deviation  
      r=[0 for x in self.nSimp]
      for i in range(self.nSimp):
        for j in range(self.nSimp):r[i]+=lfSimp[j][i]
        r[i]/=self.nSimp
      return r 
}
// *******************************************************
double Simplex::GetErr(const int iP)
{CHECK_INDEX_RETURN(iP,nSimp,0)
 return lfError[iP]; 
}
// *******************************************************
double Simplex::GetLimit(const int iP) const 
{CHECK_INDEX_RETURN(iP,nSimp,0)
 return Limit[iP]; 
}
// *******************************************************
double Simplex::GetStStep(const int iP) const 
{CHECK_INDEX_RETURN(iP,nPars,0)
 return StSteps[iP]; 
}
// *******************************************************
void Simplex::GetMinMaxErr(double & MinErr, double & MaxErr,
                           int & iMi, int & iMa) const
{int j;
 for(j=0;j<nSimp;j++)
    {if(j==0)
       {MaxErr=MinErr=lfError[j];iMi=iMa=0;}
     else
       {if(lfError[j]>MaxErr){MaxErr=lfError[j];iMa=j;}
	if(lfError[j]<MinErr){MinErr=lfError[j];iMi=j;}
       }
    }
}    
// *******************************************************
void Simplex::SetLimit(const int iP, const double L)
{CHECK_INDEX_RET(iP,nPars)
 Limit[iP]=L;
}
// ********************************************************
void Simplex::SetStStep(const int iP, const double S)
{CHECK_INDEX_RET(iP,nPars)
 StSteps[iP]=S;
}
// ********************************************************

      
if __name__=='__main__':
#  print Horner([1,2,-5],2)
#  print 1+2*2-5*2*2
##  p=[{'A':1000.,'x0':100.,'HBW':5},{'A':1000.,'x0':120.,'HBW':5}]
##  c=[{'Backg':200.}]
##  for i in range(200):
##     print i,LorentzFn(i,p,c)
     
