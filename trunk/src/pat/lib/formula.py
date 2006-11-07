from types import *
import math,string,sys
from re import *

#$Log: formula.py,v $
#Revision 1.1  2004/08/04 11:13:51  herbie
#*** empty log message ***
#
CVS_ID="$Id: formula.py,v 1.1 2004/08/04 11:13:51 herbie Exp $"

class Formula:
  """Allowed formulas
     #: Colum number : 1 - n
     c# = 1.23 <+-*/^> c#	 1
     c# = c# <+-*/^> 1.234	 2
     c# = c# <+-*/^> c#	         3
     c# = FUN(c#)		 4

     FUN:  ABS, COS, EXP, LOG, SIN, TAN, ACOS, ASIN, ATAN, SQRT
     Analyzes the formula and set the variables:
     Src[0 or 1][{'type':-,'value':-},{'type':-,'value':-}]:
        where 'type' is one of 'fun', 'col' or 'val'
             'value' i an an int, float or function name
          Op: is + - * / or ^
     DestCol: is an int specifing the destination column
  """
  Funs={  'ABS': math.fabs,
          'COS': math.cos,
          'SIN': math.sin,
          'TAN': math.tan,
         'ACOS': math.acos,
         'ASIN': math.asin,
         'ATAN': math.atan,
         'SQRT': math.sqrt,
          'LOG': math.log,
        'LOG10': math.log10,
          'EXP': math.exp
       }
  def __init__(self,formula):
    """Analyzes the formula and set the variables:
       Src[0 or 1][{'type':-,'value':-},{'type':-,'value':-}]:
              where 'type' is one of 'fun', 'col' or 'val'
                   'value' i an an int, float or function name
             Op: is + - * / or ^
        DestCol: is an int specifing the destination column

        Input: formula: string
        Error: When an error during anlyzing is found
               IndexError or ValueError is raisen
    """
    if type(formula) != StringType:
       raise ValueError('Formula must be a string')
    self.formula=formula.upper().strip()
    self.DestCol=None
    self.Src=[{'type':None,'value':None},{'type':None,'value':None}]
    self.Op=None
    Error=1

    t=self.formula.split('=')
    if len(t) != 2:
      raise IndexError('Formula must have exact one "="')
    LeftF=t[0].strip()
    RightF=t[1].strip()
    self.DestCol=string.atoi(LeftF[1:])

    ReCol="C[1-9][0-9]*"
    ReOp="[\^\*\-\+\/]"
    ReB=" *"
    ReNo="-?(([0-9]+\.[0-9]*)|([0-9]+)|(\.[0-9]+))([eE][---+]?[0-9]+)?"
    ReFun="[A-Z][A-Z0-9]+ *\("

    Re="^"+ReCol+ReB+'='+ReB+ReCol+ReB+ReOp+ReB+ReCol+"$"
    # c# = c# <+-*/^> c#
    R=compile(Re)   
    if R.search(self.formula)!=None:
       O=compile(ReOp)
       t=O.split(RightF)
       self.Src[0]['value']=string.atoi(t[0].strip()[1:])
       self.Src[1]['value']=string.atoi(t[1].strip()[1:])
       self.Src[0]['type']='col'
       self.Src[1]['type']='col'
       self.Op=O.findall(RightF)[0]
       Error=0

    # c# = c# <+-*/^> 1.23456
    Re="^"+ReCol+ReB+'='+ReB+ReCol+ReB+ReOp+ReB+ReNo+"$"
    R=compile(Re)   
    if R.search(self.formula)!=None and Error != 0:
       O=compile(ReOp)
       self.Op=O.findall(RightF)[0]
       t=O.split(RightF)
       self.Src[0]['value']=string.atoi(t[0].strip()[1:])
       O=compile(ReCol+ReB+ReOp)
       t=O.split(RightF)
       self.Src[1]['value']=string.atof(t[1].strip())
       self.Src[0]['type']='col'
       self.Src[1]['type']='val'
       Error=0
 
   # c# = 1.23456 <+-*/^> c#
    Re="^"+ReCol+ReB+'='+ReB+ReNo+ReB+ReOp+ReB+ReCol+"$"
    R=compile(Re)   
    if R.search(self.formula)!=None and Error != 0:
       O=compile(ReB+ReOp+ReB+ReCol)
       t=O.split(RightF)
       self.Src[0]['value']=string.atof(t[0].strip())
       O=compile(ReOp)
       self.Op=O.findall(RightF)[-1]
       t=O.split(RightF)
       self.Src[1]['value']=string.atoi(t[-1].strip()[1:])
       self.Src[0]['type']='val'
       self.Src[1]['type']='col'
       Error=0

    # c# = FUN(c#)		 4
    Re="^"+ReCol+ReB+'='+ReB+ReFun+ReB+ReCol+ReB+"\)$"
    R=compile(Re)   
    if R.search(self.formula)!=None and Error != 0:
       O=compile('[\(\)]')
       t=O.split(RightF)
#       print t
       f=t[0].strip()
       if not (f in self.Funs.keys()):
          raise IndexError('%s: illegal function' % f)
       self.Src[0]['value']=self.Funs[f]
       self.Src[1]['value']=string.atoi(t[1].strip()[1:])
       self.Src[0]['type']='fun'
       self.Src[1]['type']='col'
       Error=0


    if Error:
         raise ValueError('Error analyzing formula: "%s"' % self.formula)    

  def __str__(self):
    return self.formula.__str__()+\
           ": Dcol=%d, Src:" % self.DestCol+str(self.Src)+\
           ", Op: %s" % self.Op

  def GetOp(self):
      """Returns the operation (+-*/^) after analyzing the formula
      """
      return self.Op
 
  def GetDest(self):
      """Returns the destination column number after analyzing the formula
      """
      return self.DestCol

  def GetSrc(self):
      """Returns the dictionary specifing both source parmeters
         after analyzing the formula:
         A list containing two dictionarys
         [{'type':-,'value':-},{'type':-,'value':-}]:
          where 'type' is one of 'col', 'val' or 'fun'
               'value'    is an    int, float or function name
        """
      return self.Src

if __name__=='__main__':
 F=Formula(sys.argv[1])
 print F
