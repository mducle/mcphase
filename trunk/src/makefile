#*********************************************************************
#
# File:          Makefile
#
# Project:       McPhase
#
# Description:   Makefile for the program package - compiles and links all c++ code
#
# Author(s):     M. Rotter
#
# Original:	 August 2000
#
#**********************************************************************

# make looks in the following directories for files
vpath %.c  ./
vpath %.cpp ./
vpath %.h  ../include
vpath %.hpp ../include
vpath %.o  ./ ./Vector ../cf1ion_module

#---------------------------------------------------------------------
#  Operating system
#  (linux windows)
# LINUX = linux
 WINDOWS = windows
#
# STATIC = static
#
#---------------------------------------------------------------------
#               ... and his Flags
#
#                Includes
#

INCLDIRS =      -I. \
		-I../include \
		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include

include ./Makefile.Common


CFLAGS  = -fpermissive -c $(CCOPTS) $(INCLDIRS)

# PGPLLIBS =       -L$(PGPLOT_DIR)  -lcpgplot -lpgplot
ifdef LINUX
XLIBS= -ldl
#XLIBS=           -L/usr/X11/lib -lX11 -lXext -lc -lXpm -ldl
#EZLIB=           -ljpeg -ldl
endif
LFLAGS =          $(XLIBS) $(EZLIB) -lstdc++ -lm 


#---------------------------------------------------------------------
# RCS options
#
#               RCS lock options
RCSLOCK  =      co -l -r$(VERSION)
#               RCS check out options
RCSCO    =      co -r$(VERSION)
#               RCS check in options
RCSCI    =      ci -u -f -s"Rel" -r$(VERSION) -m"$(LOCKMSG)"


#---------------------------------------------------------------------
# All files needed for the library
#
#	File list for RCS control
#               all include files

#               source files


#--------------------------------------------------------------------
# What has to be made
#
#               Names of the include files in the 
#		installation directories
#

INCL = mcphas.h  htcalc.h testspin.hpp qvectors.hpp spincf.hpp par.hpp  jjjpar.hpp inipar.hpp mfcf.hpp vector.h martin.h
OBJS = jjjpar.o par.o physprop.o qvectors.o mcphas.o testspin.o spincf.o inipar.o mfcf.o martin.o ionpars.o

VOBS = ivector.o imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o cmatrix.o cmatmul.o \
	dludcomp.o dsvdcomp.o dortho.o cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
	chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o \
	matpack.o sortdbl.o sort2dbl.o sort3dbl.o

CFIELDOBJS= cfieldrout.o diahermx.o eingabe.o intensit.o komplex.o matrix.o stevens.o  \
	theta.o ortho.o mb11a.o minima.o spline.o va05a.o

LINKER= g++
#OBJS := $(patsubst %.c,%.o,$(wildcard ./src/*.c)) \
#	$(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp)) 
#---------------------------------------------------------------------
# Build mcphas
#

all :  martin.o cond  coq2jjj jjj2j mcdispit mcphasit spins addj singleion charges

# ****************************************************
# compile sources 

martin.o:  martin.c martin.h
	g++  $(CFLAGS)  $< 
	
paramete.o: paramete.cpp paramete.hpp martin.h vector.h 
	g++  $(CFLAGS)  $< 

par.o: par.cpp par.hpp jjjpar.hpp martin.h vector.h
	g++  $(CFLAGS)  $< 

jjjpar.o: jjjpar.cpp jjjpar.hpp martin.h vector.h myev.h
	g++  $(CFLAGS)  $< 

testspin.o: testspin.cpp testspin.hpp spincf.hpp martin.h
	g++  $(CFLAGS)  $< 

spincf.o: spincf.cpp spincf.hpp vector.h
	g++  $(CFLAGS)  $< 

chargedensity.o: chargedensity.cpp chargedensity.hpp vector.h
	g++  $(CFLAGS)  $< 

mfcf.o: mfcf.cpp mfcf.hpp vector.h
	g++  $(CFLAGS)  $< 

mdcf.o: mdcf.cpp mdcf.hpp vector.h
	g++  $(CFLAGS)  $< 

jq.o: jq.cpp jq.hpp vector.h
	g++  $(CFLAGS)  $< 

inipar.o: inipar.cpp inipar.hpp martin.h vector.h 
	g++  $(CFLAGS)  $< 

ionpars.o: ionpars.cpp ionpars.hpp
	g++ $(CFLAGS) $<

inimcdis.o: inimcdis.cpp inimcdis.hpp martin.h vector.h mfcf.hpp
	g++  $(CFLAGS)  $< 

physprop.o: physprop.cpp physprop.hpp spincf.hpp vector.h mfcf.hpp martin.h par.hpp
	g++  $(CFLAGS)  $< 

qvectors.o: qvectors.cpp qvectors.hpp spincf.hpp martin.h vector.h inipar.hpp myev.c
	g++  $(CFLAGS)  $<

mcphas.o: mcphas.c $(INCL) htcalc.c 
	g++ $(CFLAGS) $<

mcdisp.o: mcdisp.c  myev.c intcalc.c jq.hpp mdcf.hpp $(INCL) 
	g++ $(CFLAGS) $<

cond.o: cond.c paramete.hpp martin.h
	g++ $(CFLAGS) $<

coq2jjj.o: coq2jjj.c paramete.hpp martin.h
	g++ $(CFLAGS) $<

jjj2j.o: jjj2j.c par.hpp jjj2j.h jjjpar.hpp vector.h
	g++ $(CFLAGS) $<

spins.o: spins.c spincf.hpp  martin.h vector.h
	g++ $(CFLAGS) $<

singleion.o: singleion.c par.hpp  martin.h vector.h
	g++ $(CFLAGS) $<

charges.o: charges.c par.hpp  spincf.hpp martin.h vector.h
	g++ $(CFLAGS) $<

addj.o: addj.c par.hpp  martin.h vector.h
	g++ $(CFLAGS) $<

ifdef LINUX

inimenu.o: inimenu.c
	g++ $(CFLAGS) $<

endif

# ************************************************
# make  binaries (link .o files)

cond: cond.o paramete.o $(VOBS) martin.o
	$(LINKER)  -Wall  -o ../bin/cond $^ $(LFLAGS)

coq2jjj: coq2jjj.o paramete.o $(VOBS) martin.o
	$(LINKER)  -Wall  -o ../bin/coq2jjj $^ $(LFLAGS)

jjj2j: jjj2j.o par.o jjjpar.o  $(VOBS) martin.o $(CFIELDOBJS) ionpars.o
	$(LINKER)  -Wall  -o ../bin/jjj2j $^ $(LFLAGS)

mcdispit: mcdisp.o martin.o par.o jjjpar.o inimcdis.o jq.o mdcf.o mfcf.o ionpars.o  $(VOBS)  $(CFIELDOBJS)
	$(LINKER)  -Wall  -o ../bin/mcdispit   $^  $(LFLAGS)

mcphasit: $(OBJS) martin.o $(VOBS)  $(CFIELDOBJS) ionpars.o
	$(LINKER)  -Wall  -o ../bin/mcphasit  $^ $(LFLAGS)

spins: spincf.o spins.o martin.o  $(VOBS)
	$(LINKER)  -Wall  -o ../bin/spins $^  $(LFLAGS)

singleion:  singleion.o martin.o  jjjpar.o par.o $(VOBS) $(CFIELDOBJS) ionpars.o
	$(LINKER)  -Wall  -o ../bin/singleion $^  $(LFLAGS)

charges:  charges.o martin.o  jjjpar.o par.o $(VOBS) $(CFIELDOBJS) ionpars.o spincf.o chargedensity.o
	$(LINKER)  -Wall  -o ../bin/charges $^  $(LFLAGS)

addj: addj.o par.o martin.o jjjpar.o $(VOBS) $(CFIELDOBJS) ionpars.o
	$(LINKER)  -Wall  -o ../bin/addj $^  $(LFLAGS)


ifdef LINUX

inimenu: inimenu.o inipar.o martin.o $(VOBS)
	$(LINKER)  -Wall  -o ../bin/inimenu $^  $(LFLAGS)

endif

clean: 
	rm -vf ./*.o

