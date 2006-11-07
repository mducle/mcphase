from types import *
from array import *
from formula import *
from stdfunc import *
import string,sys
#$Log: NumTab.py,v $
#Revision 1.1  2004/08/04 11:30:10  herbie
#*** empty log message ***
#
CVS_ID="$Id: NumTab.py,v 1.1 2004/08/04 11:30:10 herbie Exp $"

class NumTable:
   """Represents a numerical table where you can perform (mathematical) 
      operations (like a spreadsheet).
      It consists of a list of lines where the columns are seperated by
      blanks or tabs.
      Simple: A matrix of strings representing rows and columns of numbers.
      Only the columns involved in mathematical operations will be converted
      to floats then the operation will be performed and the result will be
      converted back to strings.
   """
# -----------------
   def __init__(self,line_list):
     """Makes a NumTable object from a list of lines containing
        numbers delimited by blank or tab.
	Each line (row) must have the same amount of numbers (columns)
	End of line can be \\n: UNIX or \\n\\r: DOS
	Parameters:
	  line_list: List of strings representing the data lines
     """
     self.LL=[]
 
     if type(line_list) != ListType or len(line_list)==0:       
       raise ValueError('Input must be a nonempty list')
     self.ColDesc=[]
     self.DeLim=None

     if line_list[0].find(DOS_END) != -1:
       self.DeLim=DOS_END
     else:
       if line_list[0][-1]==UNIX_END: self.DeLim=UNIX_END
     lc=0
     for i in line_list:
         if type(i) != StringType:
            raise ValueError (i,'Must be a string')

         l=i.strip().split()
         if lc==0:
           nco=len(l)
           for n in range(0,nco): self.LL.append([])
         nc=len(l)
         if nc != nco:
           raise ValueError('Number of columns changed %d: %d -> %d' %\
                     (lc+1,nc,nco))
         nco=nc
         for n in range(0,nc):
           try: t=string.atof(l[n])
           except: raise ValueError (l[n],'Must be a string representing a float')
           self.LL[n].append(l[n])
         lc+=1  
# -----------------
   def SetDelim(self,d):
     """ Changes the line delimiter.
         Parameters:
	 d: Type of line delimiter (DOS_END or UNIX_END)
	                           defined in stdfunc.py
     """
     if not (d==DOS_END or d==UNIX_END):
        raise ValueError("Unsupported delimiter")
     self.DeLim=d
# -----------------
   def SetColumnSpec(self):
     """Tries to find out the numerical type (int,double,exp) 
        and the size (precision) of a column.
	For internal use only.
     """
     for i in range(len(self.LL)):
       ct={'type':'int',  # 'int' or 'float' or 'exp'
    	   'size':-1,     # max number of digits (incl sign or .)
    	   'decd':0       # max number of digits after .
          }
       self.ColDesc.append(ct)
       for j in self.LL[i]:
         sz=len(j)
         dd=0
         tp='int'
         dp=j.find('.')
         ep=j.find('E')
         if ep == -1: ep=j.find('e')
         if dp != -1:
     	    tp='float'
     	    dd=len(j)-dp-1
         else:
           if ep!=-1: dp=ep+1
         if ep != -1:
     	   tp='gexp'
     	   dd=ep-dp-1
#         print i,j,sz
         if sz > self.ColDesc[i]['size']: self.ColDesc[i]['size']=sz
         if dd > self.ColDesc[i]['decd']: self.ColDesc[i]['decd']=dd
         if tp < self.ColDesc[i]['type']: self.ColDesc[i]['type']=tp
# -----------------
   def __str__(self): return self.LL.__str__()
# -----------------
   def __getitem__(self,iCol):
       """Returns the specified column
       """
       return self.LL[iCol]
# -----------------
   def __len__(self):
       """Returns the number of lines
       """
       return len(self.LL)
# -----------------
   def __setitem__(self,iCol,v):
       """Replaces the specified column with list v
       """
       self.TestColList(v,ListType,StringType)
       self.LL[iCol]=v
# -----------------
   def append(self,l):
       """Appends a column from a list of strings
       """
       self.TestColList(l,ListType,StringType)
       self.LL.append(l)
# -----------------
   def append_rows(self,nt):
       """Appends rows from an other NumTable
       """
       ol=len(self.LL)
       nl=len(nt.LL)
       if ol != nl:
         raise ValueError('Objects must have same numbers of columns (%d,%d)' % (ol,nl))
       for c in range(0,ol):
         for r in range(0,len(nt.LL[0])):
             self.LL[c].append(nt.LL[c][r])     
# -----------------
   def swap_col(self,a,b):
       """Swaps two columns
       """
       t=list(self.LL[a]); self.LL[a]=self.LL[b] ; self.LL[b]=t 
# -----------------
   def pop_row(self,i=-1):
       """Pops row i from the list of strings
       """
       for c in self.LL:
           c.pop(i)
# -----------------
   def Export(self,iCol,typ=ListType):
       """Exports a column to a list of floats
          Parameter:
            iCol: Column to export (0...nCols-1)
             typ:  ListType: a list of floats is exported
                  TupleType: a list of tuples (value,index) is exported
       """
       if not (typ == ListType or typ == TupleType):
          raise TypeError('typ must be list or tuple')
       l=[]
       ic=0
       for i in self.LL[iCol]:
         if typ == ListType: l.append(string.atof(i))
         else: l.append( (string.atof(i),ic+1) ); ic+=1
       return l
# -----------------
   def Import(self,iCol,l,format="%-18.11g"):
       """Imports a column from a list of floats
          if iCol=-1 list is appended
	  Parameters:
	    iCol: number of column to be imported (replaced)
	          -1: append column
	       l: list of numbers (float or int) to be imported
	  format: format specifiers for converting input list to string
       """
       self.TestColList(l)
       if iCol >= 0 and iCol < len(self.LL):
         for i in range(0,len(self.LL[iCol])):
           t=format % l[i]
           self.LL[iCol][i]=t.strip()
       else:
         tl=[]
         for i in l:
           t=format % i
           tl.append(t.strip())
         self.LL.append(tl)
# -----------------
   def Write(self,fp=sys.stdout,Align=0,Rows=(-1,-1,-1),Cols=(-1,-1),Vals=(-1,None,None),\
                  Sort=None,quiet=0):
       """ Row range: (start,end,step)
           Col range: (start,end)
                Sort: list with output order, 0 no output
           Val range: (colmn, min, max)
                      only only values between min and max are printed
               Align: columns are aligned by inserting blanks
       """
       ro=[0,len(self.LL[0]),1]
       if Rows[0] != -1: ro[0]=Rows[0]
       if Rows[1] != -1: ro[1]=Rows[1]+1
       if Rows[2] != -1: ro[2]=Rows[2]
#       print Rows,ro

       co=[0,len(self.LL)]
       if Cols[0] != -1: co[0]=Cols[0]
       if Cols[1] != -1: co[1]=Cols[1]+1

       vo=[0,min(self.Export(Vals[0])),max(self.Export(Vals[0]))]
       if Vals[0] != -1: vo[0]=Vals[0]
       if Vals[1] != None: vo[1]=Vals[1]
       if Vals[2] != None: vo[2]=Vals[2]

       if Sort != None and type(Sort) != ListType:
          raise TypeError('Sort must be a list')

#       print ro,Rows,co,Cols,Vals,vo
       if Align:
          self.SetColumnSpec()
       count=0
       for r in range(ro[0],ro[1],ro[2]):
          if not (Vals[0]==-1 and Vals[1]==None and Vals[2]==None): 
             v=string.atof(self.LL[vo[0]][r])
             if not (v >= vo[1] and v <= vo[2]): continue

          if Sort != None and Sort[r]==0: continue
       
          for c in range(co[0],co[1]):
            fill=''
#	    print "***",c
            if Sort == None: z=self.LL[c][r]  # unsorted
            else: z=self.LL[c][Sort[r]-1]
            fill=''
            if Align: fill=(self.ColDesc[c]['size']-len(z))*' '
            if not quiet: fp.write("%s%s  " % (fill,z))
	  count+=1
          if not quiet: fp.write("%s" % self.DeLim)
       return (count,range(co[0],co[1]))
# -----------------
   def TestColList(self,cl,typ1=None,typ2=None):
       """Tests the type combatibility of column list and input list
          and if input list is a list of strings representing numbers.
          Parameters:
	    cl: input list
	  typ1: type of input list
	  typ2: type of elements of input list 
	  Result:
             Raises an error if lists not compatible
	     length of input list != number of rows: ValueError
	     type of input list != typ1: ValueError
	     type of element != typ2: ValueError
	     element is a string NOT representing a number: ValueErrror
          For internal use only.
       """
       ll=len(cl)
       if ll != len(self.LL[0]):
         raise ValueError('Must have length %d but has %d' % (len(self.LL[0]),ll))
       if typ1 != None:
         if type(cl) != typ1:
           raise ValueError('Must be a %s' % typ1)
       if typ2 != None:
         for i in cl:
           if type(i) != typ2:
             raise ValueError (i,'Must be a %s' % typ2)
           if typ2==StringType:
             try: t=string.atof(i)
             except: raise ValueError (i,'Must be a string representing a float')
# -----------------
   def Calculate(self,form):
       """Performs a matematical calculation:
          #: Colum number : 1 - n
          c# = 1.23 <+-*/^> c#        1
          c# = c# <+-*/^> 1.234       2
          c# = c# <+-*/^> c#	      3
          c# = FUN(c#)  	      4
          FUN:  ABS, COS, EXP, LOG, SIN, TAN, ACOS, ASIN, ATAN, SQRT
	  Parameter:
	    form: a formula (string) from the list above (e.g.: "c2=c3*55.6")
	  Result:
	    The content of the column on left side of the formula will be
	    changed corresponding to the fromula.
       """
       OpD={'+':float.__add__,
            '-':float.__sub__,
            '*':float.__mul__,
            '/':float.__div__,
            '^':float.__pow__}

       F = Formula(form)
       Src=F.GetSrc()
       Op=F.GetOp()

       if Src[0]['type']=='fun':
         r=[]
         for i in self.Export(Src[1]['value']-1):
           r.append(Src[0]['value'](i))
         self.Import(F.GetDest()-1,r)
         return
      
       loop=None
       if Src[0]['type']=='col': loop=self.Export(Src[0]['value']-1)
       if Src[0]['type']=='val': loop=self.Export(Src[1]['value']-1)

       if loop == None:
          raise ValueError('Bad formula: "%s"' % form)

       r=[]
       c=0
       if Src[1]['type']=='col': right_l=self.Export(Src[1]['value']-1)
       if Src[1]['type']=='val': right=Src[1]['value']
       if Src[0]['type']=='val': left=Src[0]['value']

       for i in loop:
        if Src[1]['type']=='col': right=right_l[c]
        if Src[0]['type']=='col': left=i
        r.append(OpD[Op](left,right))
        c+=1

       self.Import(F.GetDest()-1,r)
# ==========================================
if __name__=='__main__':
 l1=[" 1 2 3 4 5 6\n",
     " 10 20 30 40 50 60  \n"]
 l2=["1 2 3 4 5 6\r\n",
     "10 20 30 40 50 60\r\n"]
 l3=[" 1 2.44 -33. 4 5 6\n",
     " 10.0  20.1    30   40.1E-3 50 60\n"]

 l4=[" 1     2.44    -33. 4       5  6000\n",
     " 10.0  20.1    30.1 40.1E-3 50 60\n",
     " 345.  1024.5  0.01 99.5    55 600\n"]

 N=NumTable(l4)
 print N
 N.Import(-1,array('d',[88,812434230.123,11]))
# print N
# l=N.Export(6,TupleType)
# print l
# l.sort()
# sl=[]
# for i in l: sl.append(i[1])
# print sl
# N.Write(Align=1,Vals=(0,1,10))
 N.Calculate("C1=100.+C2")
 N.Write(Align=1)
