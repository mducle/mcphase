import string,types
import sys
from pprint import *

#$Log: filepar.py,v $
#Revision 1.2  2005/12/19 09:13:05  herbie
#*** empty log message ***
#
#Revision 1.1  2005/12/15 08:57:36  herbie
#Initial revision
#
CVS_ID="$Id: filepar.py,v 1.2 2005/12/19 09:13:05 herbie Exp $"

#NoFile = 'Input file not found'

class FilePar:
   """ Finds 'parameters' in a file or a list like:
       ---------
       conf_file = /etc/my_config
       Iterate = 1000
       # this is a comment
       ; this is also comment
       Retry = 5 # also a comment
       [Home]
       #This is a new topic
       dir=/home/herbie
       name = herbie
       -----------
       The file can be splitted in topics.
       Topic names are enclosed in []: e.g. [topic_name]
       Parameter names and topic names must NOT include blancs ' '.
       Comment starts with # or ;
       Topic lines must not contain a comment
       Parameter lines can contain a comment up to end of the line
       Empty lines are ignored
   """
   def __init__(self,source,ignore_comment=1,sep='='):
       """ __init__(self,filename,ignore_comment=1,sep='=')
          Opens the parameter file/list
          Parameters:
                  source : Name of the parameter file, file or list
           ignore_comment: 1: lines beginning with '#' are ignored
                           0: lines beginning with '#' are processed
                      sep: separator between name <sep> value
                                             host   =   localhost
               Error:
                 raises OSError if 'filename' can not be opened
       """
       self.content=[]
       self.ignore_comment=ignore_comment
       self.sep=sep
       self.buffer=None
       inp=None
       if type(source) == types.ListType: self.buffer=source
       if type(source) == types.StringType:
          inp=open(source,'r')
       if type(source) == types.FileType:
          inp=source
       if inp != None:
          self.buffer=inp.readlines() 
          inp.close()
       if self.buffer == None:
          raise TypeError('Source must be filename, file or list')
 
   def FindTopic(self,topic):
       """ FindTopic(self,topic)
           Write all parameters of a topic in the internal buffer
           ------------
           [topic]
           host = localhost
           retry = 5
           timeout = 10
           [next_topic]
           ....
          Parameters:
            topic: Name of the topic in that the parameters are searched
                   If '' the 'first parts of the file belonging to no topic
                   are processed
                   The [] in the topic name must be omitted
              Return:
                   Length of the buffer (no of matching lines found)
          For internal class use only!
          For searching use FindPar(...)
       """
       self.topic_content = []
       if topic == None or topic == '':
          for line in self.content:
             if line[0] == '[': break
             self.topic_content.append(line)
          return len(self.topic_content)

       n_topic='['+topic+']'
       found=0
       for line in self.content:
           if line == n_topic:
              found=1
              continue
           if found:
              if line[0] == '[': break
              self.topic_content.append(line)
       return len(self.topic_content)
       
      
   def FindPar(self,name,topic=None,until=None):
       """ FindPar(self,topic)
           Find all parameters in a topic
           ------------
           # Next 3 lines belong to topic = '' (empty)
           Figure=block
           height = 10
           depth = 5
           [my_site]
           # topic lines must not contain a comment
           host = localhost
           retry = 5
           timeout = 10
           [next_topic]
           ....
           ------------
          Parameters:
            topic: Name of the topic in that the parameters are searched
                   If '' the first parts of the file belonging to no topic
                   are processed
                   The [] in the topic name must be omitted
             name: Name of the parameter
            until: Searching stops when a line matches with this string
                   Can be used if only the header part of a file
                   contains parameters 
              Return:
             string: contents of the parameter, if found
               None: if no matching parameter can be found
               Error:
                    raises 'IndexError' if file contains no lines
       """
       if self.buffer == None or len(self.buffer) == 0 :
          raise IndexError('Source empty')

       for l in self.buffer:
           line=string.strip(l)
           if len(line) == 0: continue

           if self.ignore_comment :
              if line[0] == '#' or line[0] == ';': continue
           else: 
              if line[0] == '#':
                 t=line[1:]
                 line=string.strip(t)

           # #; can not be 1st char now
           # remove comment part of the line
           # a=5  #This is a
           f=string.find(line,'#')
           if f == -1:
              f=string.find(line,';')
           if f != -1:
              t=line[:f]
              line =string.strip(t)

           if until != None and line == until: break
           self.content.append(line)

       if not self.FindTopic(topic): return None
       for line in self.topic_content:
           tok=string.split(line,self.sep)
           if string.strip(tok[0]) == string.strip(name):
              return string.strip(tok[1])            

def ReadParDic(pd,filename):
    """ ReadParDic(pd,filename)
        Reads parameters from a file which are specified
        in the dictionary and performs a range check
        Parameters:
               pd: dictionary as described below
         filename: parameter file name
        The input dictionary must look like this:
        pd={  'topic1': {'name11': {'value':'', 'min': None,'max': None},
              'topic2': {'name21': {'value': -1,'min': 0,'max': None}, 
                           'Type': {'value': -1,'min': 0,'max': 1}, 
                           'Conn': {'value':  0,'min': 1,'max': 8}, 
                        'Current': {'value': -1,'min': 1,'max': 200} }
               }
        walks through the dictionary searches for 'names' in 'topic'
        and assigns the found value to 'value'.
        If min or max not None a value check is performed.
        return: a copy of input with values filled in
        errors: IndexError: Parameter file not readable
                LookupError: Parameter not found
                ValueError: Parameter out of given range
    """
    fp=FilePar(filename)
    r=pd
    for i in r.keys():
         for j in r[i].keys():
           try:
             p=fp.FindPar(j,i)
           except IndexError:
             raise IndexError
             return
           if p == None and r[i][j]['value'] != None:
             raise LookupError('Parameter [%s], %s not found' % (i,j))
           if p != None:  
            min=r[i][j]['min']
            max=r[i][j]['max']
            if min != None:
              if type(min)==types.IntType: p=string.atoi(p)
              if type(min)==types.FloatType: p=string.atof(p)
              if p < min : 
                raise ValueError('Parameter [%s], %s < min' % (i,j))
            if max != None:
              if p > max : 
                raise ValueError('Parameter [%s], %s > max' % (i,j))

            r[i][j]['value']=p
    return r

def ReadParDicNc(pd,filename):
    """ ReadParDic(pd,filename)
        Reads parameters from a file which are specified
        in the dictionary NO range check is performed
        Parameters:
               pd: dictionary as described below
         filename: parameter file name
        The input dictionary must look like this:
        pd={  'topic1': {'name11':  None},
              'topic2': {'name21': ''}, 
                           'Type': 5, 
                           'Conn': 1, 
                        'Current': 30} }
               }
        walks through the dictionary searches for 'names' in 'topic'
        and assigns the found value to 'value'.
        return: a copy of input with values filled in
        errors: IndexError: Parameter file not readable
                LookupError: Parameter not found
    """
    fp=FilePar(filename)
    r=pd
    for i in r.keys():
         for j in r[i].keys():
           try:
             p=fp.FindPar(j,i)
           except IndexError:
             raise IndexError
             return
           if p == None and r[i][j]['value'] != None:
             raise LookupError('Parameter [%s], %s not found' % (i,j))
           r[i][j]=p
    return r

def FindFromDic(source,dic,idx=0,topic=None,flag_ix=None):
    """ FindFomDic(filename,dic,idx=0,topic=None)
        Reads parameters from a source.
        The paramters are specified from the keys in the dictionary
        The the values in the dictionary get the paramerter values
        Parameters:
            source: parameter file name or list of strings
               dic: dictionary as described below
                ix: index of the value if values in the dictionary are lists
             topic: topic were parameters will be searched
           
        The input dictionary must look like this:
        dic={'name11': 1.},
             'name21': 0, 
             'name_a': '' 
            }
        or
        dic={'name11': [...,1.,...},
             'name21': [...,0,.... 
             'name_a': [...,'',...] 
            }
    """
    fip=FilePar(source)
    f=0
    for k,v in dic.iteritems():
      p=fip.FindPar(k,topic)
      if type(v)==types.ListType:
        if type(v[idx])==types.IntType: dic[k][idx]=string.atoi(p)
        if type(v[idx])==types.FloatType: dic[k][0]=string.atof(p)
        if type(v[idx])==types.StringType: dic[k][0]=p
        if flag_ix != None and p != None: f |= dic[k][flag_ix]
      else:
        if type(v)==types.IntType: dic[k]=string.atoi(p)
        if type(v)==types.FloatType: dic[k]=string.atof(p)
        if type(v)==types.StringType: dic[k]=p
    return f

if __name__=='__main__':
 print 'Usage:', sys.argv[0],'file topic name'

 f=FilePar(sys.argv[1])

 print 'File: %s Topic=%s Name=%s' % (sys.argv[1],sys.argv[2],sys.argv[3])

 try:
   r=f.FindPar(sys.argv[3],sys.argv[2])
 except InedxError,msg:
   print sys.argv[1],msg
   sys.exit()

 print 'Found *%s*' % r

 L=['[Topic]','value=6.78']
 fl=FilePar(L)
 print fl.FindPar('value','Topic')
