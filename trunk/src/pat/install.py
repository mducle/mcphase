#!/usr/bin/env python

"""Usage:
   install.py
   
   Interactive installation of PAT tools.
   See README, INSTALL and LICENCE for details.
   
   Starting as normal user installs PAT for single user in
   her/his home directory.
   Starting as root installs PAT for all users in a system directory.
   
   Short installation instructions:
   
   1. Unpack your tar/zip-file (pat-v.v.vv-xxxx.tgz[zip])
      in a temporary directory.
   2. change to the created directory.
   3. Start ./install.py and answer the questions.
   4. a. Set your PATH environment variable to the pat-bin-path
      b. Set the environment variable PAT_LIB_PATH to the pat-lib-path
    5. If everthing is working you can remove the temorary directory.
   
   install.py -h   Print help message
   
"""
import sys,string,os,shutil,getpass
from getopt import *

#$Log: install.py,v $
#Revision 1.1  2005/12/15 11:59:23  herbie
#Initial revision
#
#$Id: install.py,v 1.1 2005/12/15 11:59:23 herbie Exp $

def calc_version(v_text):
   if v_text == None: return None
   try: l=v_text.split('-')
   except: return -1
   if l==[] or l[0]!='pat':
      return -2
   try: v=l[1].split('.')
   except: return -3
   s=0
   c=1000
   for i in v:
       try: s+=string.atoi(i)*c
       except: pass
       c/=10
   return s
# ---------------------------------
   
if len(sys.argv)<1: sys.exit(1)
if len(sys.argv)>1: sys.stderr.write(__doc__);sys.exit(1);
opt=getopt (sys.argv[1:], "h")
for o in opt[0]:
  if o[0]=='-h': print __doc__;sys.exit(1);

username=getpass.getuser()
euid=os.geteuid()
uid=os.getuid()

# Try to find an already installed version
try:
  inst_lib=os.environ['PAT_LIB_PATH']
except:
  inst_lib=None
  
installed_version=None  
if inst_lib:
   try: fp=open('%s/.version' % inst_lib ,"r")
   except: installed_version=None
   else:
      installed_version=fp.readline().strip()
      fp.close()
inst_version_value=calc_version(installed_version)

if username=="root" and euid==0 and uid==0:
   user="all users"
   install_dir="/user/local/pat"
else:
   user='user "%s"' % username
   try: install_dir=os.environ['HOME']+"/pat"
   except: install_dir=""

# get version
try:fp=open('bin/.version',"r")
except:
    raise IOError("No version avaliable in bin/.version")
version=fp.readline().strip()
fp.close()
current_version_value=calc_version(version)

if inst_version_value == None:
   prev="No previous version found"
else:
   if inst_lib[-1]=='/': t=inst_lib[:-1].rfind('/')
   else: t=inst_lib.rfind('/')
   if t!=-1:
       install_dir=inst_lib[:t]

   if current_version_value >= inst_version_value:
        prev="Previous version %s found in %s" % (installed_version,install_dir)
   if current_version_value < inst_version_value:
        prev="A newer version %s found in %s!!!\nIf you accept the path name the existing files will be replaced" % (installed_version,install_dir)

print "%s will be installed for %s" %(version,user)
print prev
print "Destination directory is %s" % install_dir
new_path=raw_input("Press enter to accept or give other path name:")
print new_path

nd=new_path.strip()
   
if nd!="": install_dir=nd
if install_dir[0]!="/":
   print "Please enter only absolute path names (beginning with /)"
   sys.exit(0)

print "%s will be installed to %s" %(version,install_dir)

if not os.path.exists(install_dir):
      print "%s does not exist, it will be created." % install_dir
      create=1
else:
      print os.listdir(install_dir)
      print "%s allready exists, its contents will be overridden." % install_dir
      create=0
      
g= raw_input("continue? [y/N]")
go=g.strip()
if not (go=="Y" or go=="y"):
   print "Installation canceled."
   sys.exit(1)
else:
   print "Installing ..."
   if create==1:
     os.makedirs(install_dir)
   for root, dirs, files in os.walk('.'):
     if root==".": 
       for s in dirs: 
	   cmd="cp -r %s %s" % (s,install_dir)
	   print cmd
	   os.system(cmd)
       for s in files:
           if s!="install.py": 
	      cmd="cp %s %s" % (s,install_dir)
	      print cmd
	      os.system(cmd)
   print "... finished"

   print "Dont forget to set the environment variable PAT_LIB_PATH"
   print "and add to your PATH %s:" % (install_dir+"/bin")
   print "If you are using bash, put the lines:"
   print 'export PAT_LIB_PATH="%s"' % (install_dir+"/lib")
   print "export PATH=$PATH:%s" % (install_dir+"/bin")

   if username=="root" and euid==0 and uid==0:
     print "in the system wide profile /etc/profile"
   
   else:
     print "in your .bashrc or .bash_profile"
