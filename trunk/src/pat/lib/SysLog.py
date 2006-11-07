import sys, syslog, string

#$Log: SysLog.py,v $
#Revision 1.4  2004/12/22 16:07:26  herbie
#*** empty log message ***
#
#Revision 1.3  2004/08/04 11:16:04  herbie
#log independenly to stderr or syslog added
#
#Revision 1.2  2004/05/26 13:53:13  herbie
#Log mask added
#
#Revision 1.1  2004/02/24 12:23:26  herbie
#Initial version
#

CVS_ID="$Id: SysLog.py,v 1.4 2004/12/22 16:07:26 herbie Exp $"

LOG_LEVEL=-1

def RepNpr(s):
    t=string.replace(s,'\n','\\n')
    u=string.replace(t,'\t','\\t')
    return string.replace(u,'\r','\\r')
   

def SLog(level, suffix, msg, StdErr=1):

    """syslog function: Logs messages to syslog and stderr
       SLog(level, suffix, msg, StdErr=0)
            level: (int) Message is logged if global variable LOG_LEVEL >= level
           suffix: (string) Preceeding text; identifyer appearing before message
              msg: (string) Text to be logged
           StdErr: 1: suffix and msg -> syslog
           	   2: suffix and msg -> stderr
           	   3: suffix and msg -> syslog & stderr
    """
    global LOG_LEVEL
    if LOG_LEVEL >= level:
       if StdErr & 1:
         syslog.openlog(suffix)
         syslog.syslog(msg)
         syslog.closelog()
       if StdErr & 1: sys.stderr.write(suffix+": "+msg+"\n")
       
class SysLog :

     """ Logs messages to syslog and stderr
         A 'global' log level can be defined.
         Only messages with message level >= 'global' log level or
          level & 'global' log level is 1 are logged
         Faciliy and priority options are set to default (USER) (see python syslog)
     """
     def __init__(self, level, suffix, stderr=1):
        """ __init__(self,level, suffix, stderr=1)
            Parameters:
              level:    (int) Defines the 'global' level
                     if positive interpretation as level:
                        Message is logged if 'global' variable LOG_LEVEL >= level
                     if negative interpretation as mask
                        Message is logged if 'global' variable LOG_LEVEL & abs(level) is 1 

             suffix: (string) Preceeding text; identifyer appearing before any message
             stderr: 1: suffix and msg -> syslog
                     2: suffix and msg -> stderr
                     3: suffix and msg -> syslog & stderr
        """
        self.LogLevel=level
        self.StdErr=stderr
        self.LogSuffix=suffix

     def NewSuffix(self, suffix):
         """ NewSuffix(self, suffix)
             Set a new log suffix
         """ 
         self.LogSuffix=suffix

     def SetLevel(self, level):
         """ SetLevel(self, level)
             Set a new 'global' log level
         """ 
         try: a=level+1
         except: return
         self.LogLevel=level

     def Log(self,level,msg):
         """ Log(self,level,msg)
             Logs suffix + msg to syslog
             Parameters:
               level: (int) Defines if message is logged
                      Depends on 'global' level
                     if 'global' level was positive interpretation as level:
                        Message is logged if global variable LOG_LEVEL >= level
                     if 'global' level was negative interpretation as mask
                        Message is logged if global variable LOG_LEVEL & abs(level) is 1 
                 msg: (string) Text to be logged
         """
         if self.LogLevel >= 0:
            if self.LogLevel >= level:
               if self.StdErr & 1:
                  syslog.openlog(self.LogSuffix)
                  syslog.syslog(msg)
                  syslog.closelog()
               if self.StdErr &2:
                  sys.stderr.write(self.LogSuffix+": "+msg+"\n")
         else:
            if (abs(self.LogLevel) & level) or (level == 0):
               if self.StdErr & 1:
                  syslog.openlog(self.LogSuffix)
                  syslog.syslog(msg)
                  syslog.closelog()
               if self.StdErr & 2:
                  sys.stderr.write(self.LogSuffix+": "+msg+"\n")
	    
if __name__=='__main__':
# S=SysLog(10,"Test Syslog-class")
# S.Log(1,"Level1");
# S.SetLevel(-2)
# S.Log(1,"Level-1");
# S.Log(3,"Level3");
# S.Log(2+4,"Level6");
# S.Log(4,"Level4");
# S.Log(0x6,"Level6x");
# S.SetLevel(0)
# S.Log(1,"Level1");
# S.SetLevel(-0x8)
# for i in range(0,16):S.Log(i,"Level%0x" % i) 
 S=SysLog(0,"Test Syslog-class")
# S.SetLevel(-1)
# S.Log(0,"Level0")
 S.Log(1,"Level1")
