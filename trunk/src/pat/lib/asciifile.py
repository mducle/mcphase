from types import *
#from pre import *
from re import *
from array import *
from stdfunc import *
from NumTab import *
from filepar import *
import string,sys,os

try: locale()['MOD_ftypes']
except: MOD_ftypes={"the":0,"sxs":0};
#print "asciifile", MOD_ftypes

#$Log: datafile.py,v $
#Revision 1.3  2006/01/04 14:41:36  herbie
#*** empty log message ***
#
#Revision 1.2  2005/12/19 09:13:05  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 08:57:36  herbie
#Initial revision
#
#CVS_ID="$Id: datafile.py,v 1.3 2006/01/04 14:41:36 herbie Exp herbie $"

FT_ASCII=10
FT_SPLINE=20
FT_SXSMulti=30
FT_THECAP=40
FILE_T={'10': ['Ascii',10],             # General Ascii file
        '20': ['Spline',20],            # Spline file; see spline.py
        '30': ['SXS Multi set',30],     # Soft x-ray spectroscopy; see sxsfile.py
        '40': ['THE Cap',40]            # Thermal expansion; thefile.py
       }

def CheckFileType(File):
   """ Checks if file is one of known file types
       and the delimiter.
       Returns a tuple containing file type and line type
   """
   fp=None
   if type(File) == StringType: fp=open(File,"r")
   if type(File) == FileType: fp=File
   if os.name!='posix' and type(File) == InstanceType: fp=File
   if fp == None: raise TypeError

   l=fp.readline()
   LineT=None
   FileT=None

   if l.find(DOS_END) != -1:
      LineT=LT_DOS
   else:
      if l[-1]==UNIX_END: LineT=LT_UNIX

   if l.find("SPLINETABLE")==0 : return (FT_SPLINE,LineT)
   if l.find("[SXSFileHeader]")==0 : return (FT_SXSMulti,LineT)
   if l.find("[FileHeader]")==0 : return (FT_THECAP,LineT)

   if type(File)==StringType: fp.close()
   if type(File) == FileType: fp.seek(0,0)
   if os.name!='posix' and type(File) == InstanceType: fp=File

   return (FT_ASCII,LineT)

class AsciiFile:
     """A data file that contains some text lines on top followed by
        a numeric table:
        Rows with the same number columns of numbers separated by blancs
	The following data object are availabe:
	
	   self.DeLim: One of DOS_END, UNIX_END (defined in stdfunc.py)
        self.LineType: One of LT_UNIX, LT_DOS (defined in stdfunc.py)
           self.FText: List of strings representing the header or text lines
            self.Nums: NubTable object representing the data columns
           self.ColId: List of strings representing the description of the data columns;
	               each list element the corresponding column
             self.iCx: Number of x column, if y=f(x) operation or plots performed
             self.iCy: Number of y column, if y=f(x) operation or plots performed
        self.FileName: filename
            self.Info: A string giving some (short) info text (for plotting)
     """
     def __init__(self,filename):
        """ Initialize variables with 'empty' values.
	    Vars will be set by Read-method
	    parameters:
	      filename: name of the data file read by Read-method later
	"""
        self.DeLim=None

        self.LineType=None 
        self.FileType=None
        self.nSets=0     # Number of sets in one file (if supported)

        self.FText=[]    # Text lines
        self.nTLines=0

        self.Nums=[]     # Data lines
        self.nNLines=0

        self.ColId=[]    # Name of colums

        self.iCx=None    # Number of x/y-Column
        self.iCy=None
        self.FileName=filename
        self.Info='No file readed'

     def Read(self,iSet=1,File=None):
        """ Read ascii file and set variables.
	    parameters:
	      iSet: Number of data set to be read.
	            Must be 1 for this class, can be changed in derived classes
	      File: None: self.FileName is used for reading
	            FileType: file pointer is used for reading. 
	"""
        # iSet must be 1 for pure ASCII files, can be used in derived classes
        if iSet!=1:
           raise ValueError("Unsupported # of data sets %d" % iSet)
        if File == None: fp=open(self.FileName,"r")
        if type(File) == FileType: fp=File; fp.seek(0,0)
        if os.name!='posix' and type(File) == InstanceType: fp=File
        buf=fp.readlines()  # Read all lines from file
        fp.close()

        cut=-self.GetDelim(buf)

        self.nTLines=self.GetTextLines(buf,self.FText,cut)

        self.Nums=NumTable(buf[self.nTLines:])
        self.nNLines=len(self.Nums[0])

        if self.nTLines > 0: self.Info=self.FText[0]
        else: self.Info='No info available'

        for i in range(0,len(self.Nums)):
           self.ColId.append('Column %d' % (i+1))

        self.FileType=FT_ASCII
        self.nSets=1
        self.iCx=0
        self.iCy=1

     def GetTextLines(self,lines,TList,cut):
        """Searches for text lines on top of the files.
	   parameters:
	     lines: lines of complete data files in wich will be searched
	     TList: text lines found will be appended to this list
	       cut: slice to be cut from end of line to remove line delimiters
	            -1: remove last char from end of line (UNIX)
		    -2: remove last two chars (DOS)
           return: Number of text lines found
	"""
#        ReNos="^[ ]*(-?(([0-9]+\.[0-9]*)|([0-9]+)|(\.[0-9]+))([eE][---+]?[0-9]+)?[ ]*)+"
        ReNos="^[ ]*(-?(([0-9]+\.[0-9]*)|([0-9]+)|(\.[0-9]+))([eE][---+]?[0-9]+)?\s+){2,}$"
        R=compile(ReNos)
        ic=0
        for i in lines:
           if R.search(i) == None: TList.append(i[:cut]); ic+=1
           else: break
        return ic

     def GetNumLines(self,lines,NList):
        """Searches for numeric lines containing the data columns
	     lines: lines of complete data files in wich will be searched
	     TList: numeric lines found will be appended to this list
           return: Number of numeric lines found
	"""
        ReNos="^[ ]*(-?(([0-9]+\.[0-9]*)|([0-9]+)|(\.[0-9]+))([eE][---+]?[0-9]+)?[ ]*)+"
        R=compile(ReNos)
        ic=0
        for i in lines:
           if R.search(i) != None: NList.append(i); ic+=1
           else: break
        return ic

     def GetDelim(self,line):
        """ Determine line delimiter
	    DOS: line ends <return><line feed>
	    UNIX: line ends wit <line-feed>
	    parameters: 
	     line: text line to be analyzed
	    return:
	     1: DOS
	     2: UNIX
	"""
        if line[0].find(DOS_END) != -1:
           self.DeLim=DOS_END
           self.LineType=LT_DOS
           return 2
        else:
           if line[0][-1]==UNIX_END:
             self.DeLim=UNIX_END
             self.LineType=LT_UNIX
             return 1

     def SetDelim(self,delim):
        """ Set line delimiter
	    parameters:
	        delim: line delimiter must be DOS_END or UNIX_END
	"""
        if not (delim == DOS_END or delim == UNIX_END):
           raise ValueError("Unsupported delimiter")
        self.DeLim=delim
        self.Nums.SetDelim(delim)

     def Info(self):
       """Assembles some info about data file
          return:
	     string containing info
       """
       r="File: %s (%s)%s" % (self.FileName, LINE_T[self.LineType],self.DeLim)
       r+="File Type: %s (%d)%s" % (FILE_T['%-d' % self.FileType][0],FILE_T['%-d' % self.FileType][1],self.DeLim)
       r+="%d data rows, %d columns%s" % (len(self.Nums[0]),len(self.Nums),self.DeLim)
       r+="%d text lines%s" % (len(self.FText),self.DeLim)
       for i in range(0,len(self.Nums)):
        a=self.Nums.Export(i,TupleType)
        r+="Column %d: %14.7g @ %d ... %14.7g @ %d%s" %\
           (i+1,min(a)[0],min(a)[1],max(a)[0],max(a)[1],self.DeLim)
       return r

     def SetXYCol(self,ix, iy):
       """Set column number to be considered as x and y.
          Necessary to plot data.
	  parameters:
	    ix: number of x-column: 0: first column
	    iy: number of y-column: 0: first column
       """
       if ix<0 or ix >= len(self.Nums) or\
          iy<0 or iy >= len(self.Nums):
           raise IndexError("Illegal column number %d or %d (0 ... %d)\n" %\
                            (ix,iy,len(self.Nums)))
       self.iCx=ix
       self.iCy=iy

     def ChgTextPar(self,PList):
       """Changes line in text if line looks like
            ParName=some_value to ParName=new_value
            Parameters:
            PList: list of string-tuples (ParName,new_value)
       """
       for i in range(0,len(self.FText)):
          for p in PList:
             if self.FText[i].find(p[0]+'=') == 0:
                self.FText[i]='%s=%s' % (p[0],p[1])

     def Write(self,fOut=sys.stdout,Sort=None,Comment=None):
        """ Write complete data file to destination
	    paramters:
	      fOut: Destination file-object
	      Sort: Sort order (see NumTab.py)
	"""
	c=''
	if Comment: c="#"
        for i in self.FText:
            fOut.write("%s%s%s" % (c,i,self.DeLim))
        return self.Nums.Write(fOut,Align=1,Sort=Sort)

     def InsTextPar(self,index,val):
       """Insert line in text after line were index is found
          or if index is an int after that index
       """
       if type(index) == StringType:
         found=None
         for i in range(0,len(self.FText)):
             if self.FText[i].find(index) != -1: found=i
         if found == None: raise IndexError("%s not fond in text lines" % val)
         self.FText.insert(found+1,val)
       elif type(index) == IntType:
         self.FText.insert(index+1,val)
       else: raise TypeError("Index must be string or int")

     def InfoText(self):
         return self.FileName

#     def Write(self,fOut=sys.stdout,Sort=None):
#        for i in self.FText:
#            fOut.write("%s%s" % (i,self.DeLim))
#        self.Nums.Write(fOut,Align=1,Sort=Sort)

if __name__=='__main__':
 F=AsciiFile(sys.argv[1])
 F.Read()
 print F.Nums
 print F.Info()
 F.Write()
