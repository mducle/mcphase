####### This section was automatically generated from
#######    /garbage/src/qt-1.31/Makefile for building on
#######    linux-gcc-shared on Tue Feb  3 16:22:57 CET 1998
#
# $Id: Makefile,v 1.3 1999/02/15 13:03:50 herbie Exp herbie $  
# $Log: Makefile,v $
# Revision 1.3  1999/02/15 13:03:50  herbie
# *** empty log message ***
# 
#
#
#

INCDIR = /usr/lib/qt/include
CFLAGS = -g -Wall -W -I/usr/X11R6/include
LIBCFLAGS = -fPIC
YACCCFLAGS = -Wno-unused -Wno-parentheses
LFLAGS = -L/usr/lib/qt/lib -lm 
CC = gcc
MOC = /usr/lib/qt/bin/moc
SHELL =	/bin/sh
AR = ar
ARFLAGS = rvs

####### Files
OLD_SOURCES=	OldCmd/*.* OldQts/*.*
#####
HEADERS=	WidMan.h 
SOURCES=        WidMan.cpp
OBJECTS=	$(SOURCES:.cpp=.o)
SRCMETA=	moc_WidMan.cpp moc_grafic.cpp moc_sfit.cpp moc_WidGen.cpp
OBJMETA=	$(SRCMETA:.cpp=.o)
TARGET=		XAusw Addspec Head Calc Fileop Deriv Integ Sort Select
TARGET+=	Remove Exchg Spline Look Zero Gap Merge Plot Dll Reduce
TARGET+=        Ggap Cgap Disp Fit Calib Polyfit Show Xport Xhelp Xmount
TARGET+=	Xfun Rgb Convtrd Convuth Ctime Theview Pfit 
#TARGET+=       Clpar fmac
#
CMD_LIB=        basicio.cpp browser.cpp cdata.cpp dfile.cpp 
CMD_LIB+=       filepar.cpp
CMD_LIB+=	formula.cpp simplex.cpp spline.cpp stdfunc.cpp thefunc.cpp 
CMD_LIB+=       strings.cpp matrix.cpp
CMD_HEAD=	$(CMD_LIB:.cpp=.h) stdinc.h
XA_LIB=		grafic.cpp sfit.cpp WidGen.cpp browser.cpp
XA_HEAD=        $(XA_LIB:.cpp=.h) rgb.h
# xtextf.h ausw.h moved to ORIG_HEAD
#######
#
# TEST PROGRAMS
#CDATA_OBJ      = stdfunc.o formula.o cdata.o
#FORMULA_OBJ    = stdfunc.o formula.o
#
CMD_SOURCES  =   addspec.cpp head.cpp calc.cpp fileop.cpp \
                 deriv.cpp integ.cpp sort.cpp select.cpp remove.cpp \
		 exchg.cpp mspline.cpp look.cpp zero.cpp gap.cpp reduce.cpp \
		 merge.cpp plot.cpp dll.cpp ggap.cpp cgap.cpp disp.cpp fit.cpp \
		 calib.cpp polyfit.cpp show.cpp xport.cpp xhelp.cpp \
		 xmount.cpp xfun.cpp rgb.cpp convuth.cpp convtrd.cpp u2th.cpp \
                 ctime.cpp theview.cpp \
                 clpar.cpp findmac.cpp pfit.cpp
MORE_HEADERS  =  formula.h cdata.h rgb.h

TEST_SOURCES  =   tfile.cpp tfile1.cpp filepar.new.cpp
TEST_HEAD     =  $(TEST_SOURCES:.cpp=.h)

ORIG_SOURCES = filepar.orig.cpp dfile.orig.cpp
ORIG_HEAD=     $(ORIG_SOURCES:.cpp=.h) xtextf.h ausw.h  

####### Implicit rules

vpath %.so   .:/usr/X11/lib:/usr/lib:/usr/local/lib:/lib
#vpath %.a    .:/usr/X11/lib:/usr/lib:/usr/local/lib:/lib:../XMeas
vpath %.a    .:/usr/X11/lib:/usr/lib:/usr/local/lib:/lib
 
.SUFFIXES: .cpp

.cpp.o:
	$(CC) -c $(CFLAGS) -I$(INCDIR) $<
moc_%.cpp: %.h
	$(MOC) -o $@ $<

####### Build rules

all: $(TARGET)

libs: libcommand.a libxausw.a

#$(TARGET): libs

$(TARGET):
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

libxausw.a libcommand.a libncommand:
	$(AR) $(ARFLAGS) $@ $?

####################################
depend:
	@makedepend -I$(INCDIR) $(SOURCES) $(MORE_SOURCES)\
	              $(CMD_LIB) $(CMD_HEAD) $(XA_LIB) $(XA_HEAD)\
		      $(CMD_SOURCES) $(MORE_HEADERS) 2> /dev/null

showfiles:
	@echo $(HEADERS) $(SOURCES) Makefile
####################################
clean:
	-rm -f *.o
	-rm -f $(SRCMETA) $(TARGET)
#ZIP & UNZIP SOURCES
zipall: 
	tar -cvzf ausw.tgz *.cpp *.h \
	                 *.hlp *.ini *.man *.list *.log README Makefile\
#zip_tu zip_home:
#	tar -cvzf $@.tgz $(CMD_SOURCES) $(CMD_LIB) $(CMD_HEAD) \
#			 $(XA_LIB) $(XA_HEAD) \
#			 $(SOURCES) $(HEADERS) \
#			 $(MORE_HEADERS) \
#	                 *.hlp *.ini *.man *.list \
#                         ./RCS
zip_tu:
	tar -cvzf ausw_tu.tgz $(CMD_SOURCES) $(CMD_LIB) $(CMD_HEAD) \
			  $(XA_LIB) $(XA_HEAD) \
			  $(SOURCES) $(HEADERS) \
                          $(TEST_SOURCES) $(TEST_HEAD)\
                          $(ORIG_SOURCES) $(ORIG_HEAD)\
			  $(MORE_HEADERS) \
	                  *.hlp *.ini *.man *.list *.log README
ziprcs_tu:
	tar -cvzf auswr_tu.tgz ./RCS
zip_home:
	tar -cvzf ausw_ho.tgz $(CMD_SOURCES) $(CMD_LIB) $(CMD_HEAD) \
			  $(XA_LIB) $(XA_HEAD) \
			  $(SOURCES) $(HEADERS) \
			  $(MORE_HEADERS) \
                          $(TEST_SOURCES) $(TEST_HEAD)\
                          $(ORIG_SOURCES) $(ORIG_HEAD)\
	                  *.hlp *.ini *.man *.list *.log README
ziprcs_home:
	tar -cvzf auswr_ho.tgz ./RCS
###########
unzip_tu:
	 tar -xvzpf ausw_ho.tgz
unzip_home:
	 tar -xvzpf ausw_tu.tgz 
unziprcs_tu:
	 tar -xvzpf auswr_ho.tgz
unziprcs_home:
	 tar -xvzpf auswr_tu.tgz 

############
#unzipall_tu:
#	 tar -xvzpf zipall_home.tgz
#unzipall_home:
#	 tar -xvzpf zipall_tu.tgz 
##########################
libcommand.a: basicio.o dfile.o thefunc.o
libcommand.a: stdfunc.o cdata.o formula.o spline.o 
libcommand.a: simplex.o strings.o matrix.o filepar.o
libxausw.a:   grafic.o moc_grafic.o 
libxausw.a:   sfit.o moc_sfit.o
libxausw.a:   WidGen.o moc_WidGen.o
libxausw.a:   browser.o moc_browser.o

####
Addspec: addspec.o -lcommand
Head:	head.o -lcommand
Calc:	calc.o -lcommand
Fileop:	fileop.o -lcommand
Deriv:	deriv.o -lcommand
Integ:	integ.o -lcommand
Sort:	sort.o -lcommand
Select:	select.o -lcommand
Remove:	remove.o -lcommand
Exchg:	exchg.o -lcommand
Spline:	mspline.o -lcommand
Look:	look.o -lcommand
Zero:	zero.o -lcommand
Gap:	gap.o -lcommand
Merge:	merge.o -lcommand
Dll:	dll.o -lcommand
Reduce:	reduce.o -lcommand
Ggap:	ggap.o  -lcommand
Cgap:	cgap.o  -lcommand
Calib:	calib.o -lcommand
Polyfit: polyfit.o -lcommand 
Convtrd: convtrd.o -lcommand
Convuth: convuth.o -lcommand
Ctime:   ctime.o -lcommand  
#Clpar:   clpar.o -lcommand
fmac:	findmac.o -lcommand -lpq -lcrypt
#U2TH:    u2th.o -lcommand -lmeas
tfile:   tfile.o -lcommand -lm
Pfit:   pfit.o -lcommand -lm
####
Plot:	plot.o libX11.so libqt.so -lxausw -lcommand
Disp:	disp.o -lcommand libX11.so libqt.so -lxausw
Fit:	fit.o -lxausw -lcommand libX11.so libqt.so 
Show:	show.o -lcommand libX11.so libqt.so -lxausw
Xport:	xport.o -lcommand libX11.so libqt.so -lxausw
Xhelp:	xhelp.o -lxausw -lcommand libX11.so libqt.so 
Xmount:	xmount.o -lxausw -lcommand libX11.so libqt.so 
Xfun:   xfun.o -lxausw -lcommand libX11.so libqt.so 
Theview: theview.o -lxausw -lcommand libX11.so libqt.so 
Rgb:	moc_rgb.o rgb.o libX11.so libqt.so
XAusw:	WidMan.o moc_WidMan.o libX11.so libqt.so -lxausw -lcommand
# DO NOT DELETE THIS LINE -- make depend depends on it.


WidGen.c.o: /usr/include/_G_config.h /usr/include/alloca.h
WidGen.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
WidGen.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
WidGen.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
WidGen.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
WidGen.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
WidGen.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
WidGen.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
WidGen.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
WidGen.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
WidGen.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
WidGen.c.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
WidGen.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
WidGen.c.o: /usr/include/limits.h /usr/include/linux/limits.h
WidGen.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
WidGen.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
WidGen.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
WidGen.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
WidGen.c.o: /usr/include/xlocale.h
WidGen.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
WidGen.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
WidGen.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
WidGen.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
WidGen.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
WidGen.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
WidGen.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
WidGen.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdatastream.h
WidGen.c.o: /usr/lib/qt/include/qdatetime.h /usr/lib/qt/include/qdialog.h
WidGen.c.o: /usr/lib/qt/include/qdir.h /usr/lib/qt/include/qevent.h
WidGen.c.o: /usr/lib/qt/include/qfile.h /usr/lib/qt/include/qfiledialog.h
WidGen.c.o: /usr/lib/qt/include/qfileinfo.h /usr/lib/qt/include/qfont.h
WidGen.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
WidGen.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
WidGen.c.o: /usr/lib/qt/include/qgdict.h /usr/lib/qt/include/qgeneric.h
WidGen.c.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
WidGen.c.o: /usr/lib/qt/include/qgmanager.h /usr/lib/qt/include/qintdict.h
WidGen.c.o: /usr/lib/qt/include/qiodevice.h /usr/lib/qt/include/qkeycode.h
WidGen.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlayout.h
WidGen.c.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
WidGen.c.o: /usr/lib/qt/include/qlistbox.h /usr/lib/qt/include/qmessagebox.h
WidGen.c.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
WidGen.c.o: /usr/lib/qt/include/qobjectdefs.h
WidGen.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
WidGen.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
WidGen.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
WidGen.c.o: /usr/lib/qt/include/qpointarray.h
WidGen.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
WidGen.c.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
WidGen.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
WidGen.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qstrlist.h
WidGen.c.o: /usr/lib/qt/include/qtableview.h /usr/lib/qt/include/qtimer.h
WidGen.c.o: /usr/lib/qt/include/qtooltip.h /usr/lib/qt/include/qwidget.h
WidGen.c.o: /usr/lib/qt/include/qwindowdefs.h /usr/lib/qt/include/qwmatrix.h
WidGen.c.o: WidGen.h cdata.h dfile.h formula.h grafic.h stdfunc.h stdinc.h
WidGen.c.o: strings.h thefunc.h
WidGen.o: /usr/include/_G_config.h /usr/include/alloca.h
WidGen.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
WidGen.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
WidGen.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
WidGen.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
WidGen.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
WidGen.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
WidGen.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
WidGen.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
WidGen.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
WidGen.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
WidGen.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
WidGen.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
WidGen.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
WidGen.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
WidGen.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
WidGen.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
WidGen.o: /usr/include/xlocale.h
WidGen.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
WidGen.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
WidGen.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
WidGen.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
WidGen.o: /usr/lib/qt/include/qarray.h /usr/lib/qt/include/qbrush.h
WidGen.o: /usr/lib/qt/include/qbutton.h /usr/lib/qt/include/qcollection.h
WidGen.o: /usr/lib/qt/include/qcolor.h /usr/lib/qt/include/qcursor.h
WidGen.o: /usr/lib/qt/include/qdatastream.h /usr/lib/qt/include/qdatetime.h
WidGen.o: /usr/lib/qt/include/qdialog.h /usr/lib/qt/include/qdir.h
WidGen.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfile.h
WidGen.o: /usr/lib/qt/include/qfiledialog.h /usr/lib/qt/include/qfileinfo.h
WidGen.o: /usr/lib/qt/include/qfont.h /usr/lib/qt/include/qfontinfo.h
WidGen.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
WidGen.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgeneric.h
WidGen.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
WidGen.o: /usr/lib/qt/include/qiodevice.h /usr/lib/qt/include/qlabel.h
WidGen.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
WidGen.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
WidGen.o: /usr/lib/qt/include/qobjectdefs.h /usr/lib/qt/include/qpaintdevice.h
WidGen.o: /usr/lib/qt/include/qpainter.h /usr/lib/qt/include/qpalette.h
WidGen.o: /usr/lib/qt/include/qpen.h /usr/lib/qt/include/qpoint.h
WidGen.o: /usr/lib/qt/include/qpointarray.h /usr/lib/qt/include/qpushbutton.h
WidGen.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
WidGen.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
WidGen.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qstrlist.h
WidGen.o: /usr/lib/qt/include/qtableview.h /usr/lib/qt/include/qwidget.h
WidGen.o: /usr/lib/qt/include/qwindowdefs.h /usr/lib/qt/include/qwmatrix.h
WidGen.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
WidMan.c.o: /usr/include/_G_config.h /usr/include/alloca.h
WidMan.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
WidMan.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
WidMan.c.o: /usr/include/bits/nan.h /usr/include/bits/pthreadtypes.h
WidMan.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
WidMan.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
WidMan.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
WidMan.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
WidMan.c.o: /usr/include/bits/wchar.h /usr/include/endian.h
WidMan.c.o: /usr/include/features.h /usr/include/gconv.h
WidMan.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/math.h
WidMan.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
WidMan.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
WidMan.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
WidMan.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
WidMan.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
WidMan.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
WidMan.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
WidMan.c.o: /usr/lib/qt/include/qbutton.h /usr/lib/qt/include/qbuttongroup.h
WidMan.c.o: /usr/lib/qt/include/qcheckbox.h /usr/lib/qt/include/qcollection.h
WidMan.c.o: /usr/lib/qt/include/qcolor.h /usr/lib/qt/include/qcursor.h
WidMan.c.o: /usr/lib/qt/include/qdatastream.h /usr/lib/qt/include/qdatetime.h
WidMan.c.o: /usr/lib/qt/include/qdialog.h /usr/lib/qt/include/qdir.h
WidMan.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfile.h
WidMan.c.o: /usr/lib/qt/include/qfiledialog.h /usr/lib/qt/include/qfileinfo.h
WidMan.c.o: /usr/lib/qt/include/qfont.h /usr/lib/qt/include/qfontinfo.h
WidMan.c.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
WidMan.c.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgdict.h
WidMan.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
WidMan.c.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qgmanager.h
WidMan.c.o: /usr/lib/qt/include/qgroupbox.h /usr/lib/qt/include/qintdict.h
WidMan.c.o: /usr/lib/qt/include/qiodevice.h /usr/lib/qt/include/qkeycode.h
WidMan.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlayout.h
WidMan.c.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
WidMan.c.o: /usr/lib/qt/include/qlistbox.h /usr/lib/qt/include/qmessagebox.h
WidMan.c.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
WidMan.c.o: /usr/lib/qt/include/qobjectdefs.h
WidMan.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpalette.h
WidMan.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
WidMan.c.o: /usr/lib/qt/include/qpushbutton.h
WidMan.c.o: /usr/lib/qt/include/qradiobutton.h /usr/lib/qt/include/qrect.h
WidMan.c.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
WidMan.c.o: /usr/lib/qt/include/qsocketnotifier.h
WidMan.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qstrlist.h
WidMan.c.o: /usr/lib/qt/include/qtableview.h /usr/lib/qt/include/qtimer.h
WidMan.c.o: /usr/lib/qt/include/qtooltip.h /usr/lib/qt/include/qwidget.h
WidMan.c.o: /usr/lib/qt/include/qwindowdefs.h WidMan.h ausw.h basicio.h
WidMan.c.o: filepar.h stdfunc.h stdinc.h strings.h
addspec.c.o: /usr/include/_G_config.h /usr/include/alloca.h
addspec.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
addspec.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
addspec.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
addspec.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
addspec.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
addspec.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
addspec.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
addspec.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
addspec.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
addspec.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
addspec.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
addspec.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
addspec.c.o: /usr/include/features.h /usr/include/gconv.h
addspec.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
addspec.c.o: /usr/include/libio.h /usr/include/limits.h
addspec.c.o: /usr/include/linux/limits.h /usr/include/math.h
addspec.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
addspec.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
addspec.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
addspec.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
addspec.c.o: /usr/include/xlocale.h
addspec.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
addspec.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
addspec.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
addspec.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
addspec.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
basicio.c.o: /usr/include/_G_config.h /usr/include/alloca.h
basicio.c.o: /usr/include/asm/ioctl.h /usr/include/asm/ioctls.h
basicio.c.o: /usr/include/asm/param.h /usr/include/asm/sigcontext.h
basicio.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
basicio.c.o: /usr/include/bits/environments.h /usr/include/bits/fcntl.h
basicio.c.o: /usr/include/bits/huge_val.h /usr/include/bits/ioctl-types.h
basicio.c.o: /usr/include/bits/ioctls.h /usr/include/bits/local_lim.h
basicio.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
basicio.c.o: /usr/include/bits/nan.h /usr/include/bits/posix1_lim.h
basicio.c.o: /usr/include/bits/posix2_lim.h /usr/include/bits/posix_opt.h
basicio.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/resource.h
basicio.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
basicio.c.o: /usr/include/bits/sigaction.h /usr/include/bits/sigcontext.h
basicio.c.o: /usr/include/bits/siginfo.h /usr/include/bits/signum.h
basicio.c.o: /usr/include/bits/sigset.h /usr/include/bits/sigstack.h
basicio.c.o: /usr/include/bits/sigthread.h /usr/include/bits/stat.h
basicio.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
basicio.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
basicio.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
basicio.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
basicio.c.o: /usr/include/endian.h /usr/include/fcntl.h
basicio.c.o: /usr/include/features.h /usr/include/gconv.h
basicio.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
basicio.c.o: /usr/include/libio.h /usr/include/limits.h
basicio.c.o: /usr/include/linux/limits.h /usr/include/linux/param.h
basicio.c.o: /usr/include/math.h /usr/include/signal.h /usr/include/stdio.h
basicio.c.o: /usr/include/stdlib.h /usr/include/sys/cdefs.h
basicio.c.o: /usr/include/sys/file.h /usr/include/sys/ioctl.h
basicio.c.o: /usr/include/sys/param.h /usr/include/sys/resource.h
basicio.c.o: /usr/include/sys/select.h /usr/include/sys/stat.h
basicio.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/ttydefaults.h
basicio.c.o: /usr/include/sys/types.h /usr/include/sys/ucontext.h
basicio.c.o: /usr/include/sys/wait.h /usr/include/sysexits.h
basicio.c.o: /usr/include/time.h /usr/include/ucontext.h /usr/include/unistd.h
basicio.c.o: /usr/include/wchar.h /usr/include/xlocale.h
basicio.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
basicio.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
basicio.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
basicio.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
basicio.c.o: basicio.h stdfunc.h stdinc.h strings.h
basicio.o: /usr/include/_G_config.h /usr/include/bits/pthreadtypes.h
basicio.o: /usr/include/bits/sched.h /usr/include/bits/stdio_lim.h
basicio.o: /usr/include/bits/types.h /usr/include/bits/wchar.h
basicio.o: /usr/include/features.h /usr/include/gconv.h
basicio.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/stdio.h
basicio.o: /usr/include/sys/cdefs.h /usr/include/wchar.h
basicio.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
basicio.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
browser.c.o: /usr/include/_G_config.h /usr/include/alloca.h
browser.c.o: /usr/include/asm/errno.h /usr/include/asm/ioctl.h
browser.c.o: /usr/include/asm/ioctls.h /usr/include/bits/confname.h
browser.c.o: /usr/include/bits/endian.h /usr/include/bits/environments.h
browser.c.o: /usr/include/bits/errno.h /usr/include/bits/fcntl.h
browser.c.o: /usr/include/bits/huge_val.h /usr/include/bits/ioctl-types.h
browser.c.o: /usr/include/bits/ioctls.h /usr/include/bits/mathcalls.h
browser.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
browser.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
browser.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
browser.c.o: /usr/include/bits/sigset.h /usr/include/bits/stat.h
browser.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
browser.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
browser.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
browser.c.o: /usr/include/bits/wordsize.h /usr/include/endian.h
browser.c.o: /usr/include/errno.h /usr/include/fcntl.h /usr/include/features.h
browser.c.o: /usr/include/gconv.h /usr/include/getopt.h
browser.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
browser.c.o: /usr/include/linux/errno.h /usr/include/math.h
browser.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
browser.c.o: /usr/include/sys/cdefs.h /usr/include/sys/file.h
browser.c.o: /usr/include/sys/ioctl.h /usr/include/sys/select.h
browser.c.o: /usr/include/sys/stat.h /usr/include/sys/sysmacros.h
browser.c.o: /usr/include/sys/time.h /usr/include/sys/ttydefaults.h
browser.c.o: /usr/include/sys/types.h /usr/include/time.h
browser.c.o: /usr/include/unistd.h /usr/include/wchar.h /usr/include/xlocale.h
browser.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
browser.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
browser.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
browser.c.o: /usr/lib/qt/include/qbutton.h /usr/lib/qt/include/qcollection.h
browser.c.o: /usr/lib/qt/include/qcolor.h /usr/lib/qt/include/qcursor.h
browser.c.o: /usr/lib/qt/include/qdialog.h /usr/lib/qt/include/qevent.h
browser.c.o: /usr/lib/qt/include/qfont.h /usr/lib/qt/include/qfontinfo.h
browser.c.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
browser.c.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgeneric.h
browser.c.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
browser.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlineedit.h
browser.c.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qlistbox.h
browser.c.o: /usr/lib/qt/include/qmessagebox.h /usr/lib/qt/include/qobject.h
browser.c.o: /usr/lib/qt/include/qobjectdefs.h
browser.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpalette.h
browser.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
browser.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
browser.c.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
browser.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
browser.c.o: /usr/lib/qt/include/qtimer.h /usr/lib/qt/include/qtooltip.h
browser.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
browser.c.o: basicio.h browser.h stdfunc.h stdinc.h strings.h
browser.o: /usr/include/_G_config.h /usr/include/alloca.h
browser.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
browser.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
browser.o: /usr/include/bits/nan.h /usr/include/bits/pthreadtypes.h
browser.o: /usr/include/bits/sched.h /usr/include/bits/select.h
browser.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
browser.o: /usr/include/bits/time.h /usr/include/bits/types.h
browser.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
browser.o: /usr/include/bits/wchar.h /usr/include/endian.h
browser.o: /usr/include/features.h /usr/include/gconv.h
browser.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/math.h
browser.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
browser.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
browser.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
browser.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
browser.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
browser.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
browser.o: /usr/lib/qt/include/qarray.h /usr/lib/qt/include/qbutton.h
browser.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
browser.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
browser.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
browser.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
browser.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
browser.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
browser.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qlineedit.h
browser.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qlistbox.h
browser.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
browser.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpalette.h
browser.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
browser.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
browser.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
browser.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
browser.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
browser.o: stdfunc.h stdinc.h strings.h
calc.c.o: /usr/include/_G_config.h /usr/include/alloca.h
calc.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
calc.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
calc.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
calc.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
calc.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
calc.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
calc.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
calc.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
calc.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
calc.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
calc.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
calc.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
calc.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
calc.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
calc.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
calc.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
calc.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
calc.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
calc.c.o: /usr/include/wchar.h /usr/include/xlocale.h
calc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
calc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
calc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
calc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
calc.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
calib.c.o: /usr/include/_G_config.h /usr/include/alloca.h
calib.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
calib.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
calib.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
calib.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
calib.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
calib.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
calib.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
calib.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
calib.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
calib.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
calib.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
calib.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
calib.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
calib.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
calib.c.o: /usr/include/linux/limits.h /usr/include/math.h
calib.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
calib.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
calib.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
calib.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
calib.c.o: /usr/include/xlocale.h
calib.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
calib.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
calib.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
calib.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
calib.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
cdata.c.o: /usr/include/_G_config.h /usr/include/alloca.h
cdata.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
cdata.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
cdata.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
cdata.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
cdata.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
cdata.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
cdata.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
cdata.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
cdata.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
cdata.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
cdata.c.o: /usr/include/ctype.h /usr/include/endian.h /usr/include/features.h
cdata.c.o: /usr/include/gconv.h /usr/include/gnu/stubs.h /usr/include/libio.h
cdata.c.o: /usr/include/limits.h /usr/include/linux/limits.h
cdata.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
cdata.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
cdata.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
cdata.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
cdata.c.o: /usr/include/xlocale.h
cdata.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
cdata.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
cdata.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
cdata.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
cdata.c.o: formula.h spline.h stdfunc.h stdinc.h strings.h
cdata.o: /usr/include/_G_config.h /usr/include/alloca.h
cdata.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
cdata.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
cdata.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
cdata.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
cdata.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
cdata.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
cdata.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
cdata.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
cdata.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
cdata.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
cdata.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
cdata.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
cdata.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
cdata.o: /usr/include/stdlib.h /usr/include/sys/cdefs.h
cdata.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
cdata.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
cdata.o: /usr/include/xlocale.h
cdata.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
cdata.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
cdata.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
cdata.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h formula.h
cdata.o: stdfunc.h stdinc.h strings.h
cgap.c.o: /usr/include/_G_config.h /usr/include/alloca.h
cgap.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
cgap.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
cgap.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
cgap.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
cgap.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
cgap.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
cgap.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
cgap.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
cgap.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
cgap.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
cgap.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
cgap.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
cgap.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
cgap.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
cgap.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
cgap.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
cgap.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
cgap.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
cgap.c.o: /usr/include/wchar.h /usr/include/xlocale.h
cgap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
cgap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
cgap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
cgap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
cgap.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
clpar.c.o: /usr/include/_G_config.h /usr/include/alloca.h
clpar.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
clpar.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
clpar.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
clpar.c.o: /usr/include/bits/nan.h /usr/include/bits/posix_opt.h
clpar.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
clpar.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
clpar.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
clpar.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
clpar.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
clpar.c.o: /usr/include/bits/wordsize.h /usr/include/endian.h
clpar.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
clpar.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/math.h
clpar.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
clpar.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
clpar.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
clpar.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
clpar.c.o: /usr/include/xlocale.h
clpar.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
clpar.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h stdfunc.h
clpar.c.o: stdinc.h strings.h
convtrd.c.o: /usr/include/_G_config.h /usr/include/alloca.h
convtrd.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
convtrd.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
convtrd.c.o: /usr/include/bits/nan.h /usr/include/bits/pthreadtypes.h
convtrd.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
convtrd.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
convtrd.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
convtrd.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
convtrd.c.o: /usr/include/bits/wchar.h /usr/include/endian.h
convtrd.c.o: /usr/include/features.h /usr/include/gconv.h
convtrd.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/math.h
convtrd.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
convtrd.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
convtrd.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
convtrd.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
convtrd.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
convtrd.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
convtrd.c.o: stdfunc.h stdinc.h strings.h
convuth.c.o: /usr/include/_G_config.h /usr/include/alloca.h
convuth.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
convuth.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
convuth.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
convuth.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
convuth.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
convuth.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
convuth.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
convuth.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
convuth.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
convuth.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
convuth.c.o: /usr/include/endian.h /usr/include/features.h
convuth.c.o: /usr/include/gconv.h /usr/include/gnu/stubs.h
convuth.c.o: /usr/include/libio.h /usr/include/limits.h
convuth.c.o: /usr/include/linux/limits.h /usr/include/math.h
convuth.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
convuth.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
convuth.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
convuth.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
convuth.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
convuth.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
convuth.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
convuth.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
convuth.c.o: cdata.h formula.h spline.h stdfunc.h stdinc.h strings.h
ctime.c.o: /usr/include/_G_config.h /usr/include/alloca.h
ctime.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
ctime.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
ctime.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
ctime.c.o: /usr/include/bits/nan.h /usr/include/bits/posix_opt.h
ctime.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
ctime.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
ctime.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
ctime.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
ctime.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
ctime.c.o: /usr/include/bits/wordsize.h /usr/include/endian.h
ctime.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
ctime.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/math.h
ctime.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
ctime.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
ctime.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
ctime.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
ctime.c.o: /usr/include/xlocale.h
ctime.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
ctime.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h stdfunc.h
ctime.c.o: stdinc.h strings.h
deriv.c.o: /usr/include/_G_config.h /usr/include/alloca.h
deriv.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
deriv.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
deriv.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
deriv.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
deriv.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
deriv.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
deriv.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
deriv.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
deriv.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
deriv.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
deriv.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
deriv.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
deriv.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
deriv.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
deriv.c.o: /usr/include/linux/limits.h /usr/include/math.h
deriv.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
deriv.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
deriv.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
deriv.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
deriv.c.o: /usr/include/xlocale.h
deriv.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
deriv.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
deriv.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
deriv.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
deriv.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
dfile.c.o: /usr/include/_G_config.h /usr/include/alloca.h
dfile.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
dfile.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
dfile.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
dfile.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
dfile.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
dfile.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
dfile.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
dfile.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
dfile.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
dfile.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
dfile.c.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
dfile.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
dfile.c.o: /usr/include/linux/limits.h /usr/include/math.h
dfile.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
dfile.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
dfile.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
dfile.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
dfile.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/float.h
dfile.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
dfile.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
dfile.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
dfile.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
dfile.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
dfile.o: /usr/include/_G_config.h /usr/include/alloca.h
dfile.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
dfile.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
dfile.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
dfile.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
dfile.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
dfile.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
dfile.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
dfile.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
dfile.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
dfile.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
dfile.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
dfile.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
dfile.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
dfile.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
dfile.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
dfile.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
dfile.o: /usr/include/xlocale.h
dfile.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
dfile.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
dfile.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
dfile.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
dfile.o: formula.h stdfunc.h stdinc.h strings.h thefunc.h
disp.c.o: /usr/include/_G_config.h /usr/include/alloca.h
disp.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
disp.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
disp.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
disp.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
disp.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
disp.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
disp.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
disp.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
disp.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
disp.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
disp.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
disp.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
disp.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
disp.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
disp.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
disp.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
disp.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
disp.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
disp.c.o: /usr/include/wchar.h /usr/include/xlocale.h
disp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
disp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
disp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
disp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
disp.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
disp.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
disp.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
disp.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
disp.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
disp.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
disp.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
disp.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
disp.c.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qlabel.h
disp.c.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
disp.c.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
disp.c.o: /usr/lib/qt/include/qobjectdefs.h /usr/lib/qt/include/qpaintdevice.h
disp.c.o: /usr/lib/qt/include/qpainter.h /usr/lib/qt/include/qpalette.h
disp.c.o: /usr/lib/qt/include/qpen.h /usr/lib/qt/include/qpixmap.h
disp.c.o: /usr/lib/qt/include/qpoint.h /usr/lib/qt/include/qpointarray.h
disp.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
disp.c.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
disp.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
disp.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
disp.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
disp.c.o: /usr/lib/qt/include/qwmatrix.h cdata.h dfile.h formula.h grafic.h
disp.c.o: stdfunc.h stdinc.h strings.h thefunc.h
dll.c.o: /usr/include/_G_config.h /usr/include/alloca.h
dll.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
dll.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
dll.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
dll.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
dll.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
dll.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
dll.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
dll.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
dll.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
dll.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
dll.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
dll.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
dll.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
dll.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
dll.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
dll.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
dll.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
dll.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
dll.c.o: /usr/include/wchar.h /usr/include/xlocale.h
dll.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
dll.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
dll.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
dll.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
dll.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
exchg.c.o: /usr/include/_G_config.h /usr/include/alloca.h
exchg.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
exchg.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
exchg.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
exchg.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
exchg.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
exchg.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
exchg.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
exchg.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
exchg.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
exchg.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
exchg.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
exchg.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
exchg.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
exchg.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
exchg.c.o: /usr/include/linux/limits.h /usr/include/math.h
exchg.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
exchg.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
exchg.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
exchg.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
exchg.c.o: /usr/include/xlocale.h
exchg.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
exchg.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
exchg.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
exchg.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
exchg.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
fileop.c.o: /usr/include/_G_config.h /usr/include/alloca.h
fileop.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
fileop.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
fileop.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
fileop.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
fileop.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
fileop.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
fileop.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
fileop.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
fileop.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
fileop.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
fileop.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
fileop.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
fileop.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
fileop.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
fileop.c.o: /usr/include/limits.h /usr/include/linux/limits.h
fileop.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
fileop.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
fileop.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
fileop.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
fileop.c.o: /usr/include/wchar.h /usr/include/xlocale.h
fileop.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
fileop.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
fileop.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
fileop.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
fileop.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
filepar.c.o: /usr/include/_G_config.h /usr/include/alloca.h
filepar.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
filepar.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
filepar.c.o: /usr/include/bits/nan.h /usr/include/bits/pthreadtypes.h
filepar.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
filepar.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
filepar.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
filepar.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
filepar.c.o: /usr/include/bits/wchar.h /usr/include/endian.h
filepar.c.o: /usr/include/features.h /usr/include/gconv.h
filepar.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/math.h
filepar.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
filepar.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
filepar.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
filepar.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
filepar.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
filepar.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
filepar.c.o: filepar.h stdfunc.h stdinc.h strings.h
filepar.o: /usr/include/_G_config.h /usr/include/bits/pthreadtypes.h
filepar.o: /usr/include/bits/sched.h /usr/include/bits/stdio_lim.h
filepar.o: /usr/include/bits/types.h /usr/include/bits/wchar.h
filepar.o: /usr/include/features.h /usr/include/gconv.h
filepar.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/stdio.h
filepar.o: /usr/include/sys/cdefs.h /usr/include/wchar.h
filepar.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
filepar.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
findmac.c.o: /usr/include/_G_config.h /usr/include/alloca.h
findmac.c.o: /usr/include/arpa/inet.h /usr/include/asm/ioctl.h
findmac.c.o: /usr/include/asm/ioctls.h /usr/include/asm/socket.h
findmac.c.o: /usr/include/asm/sockios.h /usr/include/bits/byteswap.h
findmac.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
findmac.c.o: /usr/include/bits/environments.h /usr/include/bits/fcntl.h
findmac.c.o: /usr/include/bits/in.h /usr/include/bits/ioctl-types.h
findmac.c.o: /usr/include/bits/ioctls.h /usr/include/bits/local_lim.h
findmac.c.o: /usr/include/bits/netdb.h /usr/include/bits/posix1_lim.h
findmac.c.o: /usr/include/bits/posix2_lim.h /usr/include/bits/posix_opt.h
findmac.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
findmac.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
findmac.c.o: /usr/include/bits/sockaddr.h /usr/include/bits/socket.h
findmac.c.o: /usr/include/bits/stat.h /usr/include/bits/stdio_lim.h
findmac.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
findmac.c.o: /usr/include/bits/uio.h /usr/include/bits/waitflags.h
findmac.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
findmac.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
findmac.c.o: /usr/include/ctype.h /usr/include/endian.h /usr/include/fcntl.h
findmac.c.o: /usr/include/features.h /usr/include/gconv.h
findmac.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
findmac.c.o: /usr/include/libio.h /usr/include/limits.h
findmac.c.o: /usr/include/linux/limits.h /usr/include/malloc.h
findmac.c.o: /usr/include/net/if.h /usr/include/netdb.h
findmac.c.o: /usr/include/netinet/in.h /usr/include/netinet/tcp.h
findmac.c.o: /usr/include/rpc/netdb.h /usr/include/stdint.h
findmac.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
findmac.c.o: /usr/include/sys/cdefs.h /usr/include/sys/ioctl.h
findmac.c.o: /usr/include/sys/select.h /usr/include/sys/socket.h
findmac.c.o: /usr/include/sys/stat.h /usr/include/sys/sysmacros.h
findmac.c.o: /usr/include/sys/time.h /usr/include/sys/ttydefaults.h
findmac.c.o: /usr/include/sys/types.h /usr/include/sys/uio.h
findmac.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
findmac.c.o: /usr/include/xlocale.h
findmac.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
findmac.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
findmac.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
findmac.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
findmac.c.o: stdinc.h strings.h
fit.c.o: /usr/include/_G_config.h /usr/include/alloca.h
fit.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
fit.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
fit.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
fit.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
fit.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
fit.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
fit.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
fit.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
fit.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
fit.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
fit.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
fit.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
fit.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
fit.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
fit.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
fit.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
fit.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
fit.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
fit.c.o: /usr/include/wchar.h /usr/include/xlocale.h
fit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
fit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
fit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
fit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
fit.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
fit.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
fit.c.o: /usr/lib/qt/include/qcheckbox.h /usr/lib/qt/include/qcollection.h
fit.c.o: /usr/lib/qt/include/qcolor.h /usr/lib/qt/include/qcursor.h
fit.c.o: /usr/lib/qt/include/qdialog.h /usr/lib/qt/include/qevent.h
fit.c.o: /usr/lib/qt/include/qfont.h /usr/lib/qt/include/qfontinfo.h
fit.c.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
fit.c.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgeneric.h
fit.c.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
fit.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlineedit.h
fit.c.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qlistbox.h
fit.c.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
fit.c.o: /usr/lib/qt/include/qobjectdefs.h /usr/lib/qt/include/qpaintdevice.h
fit.c.o: /usr/lib/qt/include/qpainter.h /usr/lib/qt/include/qpalette.h
fit.c.o: /usr/lib/qt/include/qpen.h /usr/lib/qt/include/qpixmap.h
fit.c.o: /usr/lib/qt/include/qpoint.h /usr/lib/qt/include/qpointarray.h
fit.c.o: /usr/lib/qt/include/qprogressbar.h /usr/lib/qt/include/qpushbutton.h
fit.c.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
fit.c.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
fit.c.o: /usr/lib/qt/include/qstack.h /usr/lib/qt/include/qstring.h
fit.c.o: /usr/lib/qt/include/qtableview.h /usr/lib/qt/include/qwidget.h
fit.c.o: /usr/lib/qt/include/qwindowdefs.h /usr/lib/qt/include/qwmatrix.h
fit.c.o: cdata.h dfile.h formula.h grafic.h sfit.h simplex.h stdfunc.h
fit.c.o: stdinc.h strings.h thefunc.h
formula.c.o: /usr/include/_G_config.h /usr/include/alloca.h
formula.c.o: /usr/include/asm/errno.h /usr/include/bits/endian.h
formula.c.o: /usr/include/bits/errno.h /usr/include/bits/huge_val.h
formula.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
formula.c.o: /usr/include/bits/nan.h /usr/include/bits/pthreadtypes.h
formula.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
formula.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
formula.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
formula.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
formula.c.o: /usr/include/bits/wchar.h /usr/include/ctype.h
formula.c.o: /usr/include/endian.h /usr/include/errno.h
formula.c.o: /usr/include/features.h /usr/include/gconv.h
formula.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
formula.c.o: /usr/include/linux/errno.h /usr/include/math.h
formula.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
formula.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
formula.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
formula.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
formula.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
formula.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
formula.c.o: formula.h stdfunc.h stdinc.h strings.h
formula.o: /usr/include/_G_config.h /usr/include/bits/pthreadtypes.h
formula.o: /usr/include/bits/sched.h /usr/include/bits/stdio_lim.h
formula.o: /usr/include/bits/types.h /usr/include/bits/wchar.h
formula.o: /usr/include/features.h /usr/include/gconv.h
formula.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/stdio.h
formula.o: /usr/include/sys/cdefs.h /usr/include/wchar.h
formula.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
formula.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h stdinc.h
gap.c.o: /usr/include/_G_config.h /usr/include/alloca.h
gap.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
gap.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
gap.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
gap.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
gap.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
gap.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
gap.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
gap.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
gap.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
gap.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
gap.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
gap.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
gap.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
gap.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
gap.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
gap.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
gap.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
gap.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
gap.c.o: /usr/include/wchar.h /usr/include/xlocale.h
gap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
gap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
gap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
gap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
gap.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
ggap.c.o: /usr/include/_G_config.h /usr/include/alloca.h
ggap.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
ggap.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
ggap.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
ggap.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
ggap.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
ggap.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
ggap.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
ggap.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
ggap.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
ggap.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
ggap.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
ggap.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
ggap.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
ggap.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
ggap.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
ggap.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
ggap.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
ggap.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
ggap.c.o: /usr/include/wchar.h /usr/include/xlocale.h
ggap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
ggap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
ggap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
ggap.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
ggap.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
grafic.c.o: /usr/include/_G_config.h /usr/include/alloca.h
grafic.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
grafic.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
grafic.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
grafic.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
grafic.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
grafic.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
grafic.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
grafic.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
grafic.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
grafic.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
grafic.c.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
grafic.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
grafic.c.o: /usr/include/limits.h /usr/include/linux/limits.h
grafic.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
grafic.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
grafic.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
grafic.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
grafic.c.o: /usr/include/xlocale.h
grafic.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
grafic.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
grafic.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
grafic.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
grafic.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
grafic.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
grafic.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
grafic.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
grafic.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
grafic.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
grafic.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
grafic.c.o: /usr/lib/qt/include/qgdict.h /usr/lib/qt/include/qgeneric.h
grafic.c.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
grafic.c.o: /usr/lib/qt/include/qgmanager.h /usr/lib/qt/include/qintdict.h
grafic.c.o: /usr/lib/qt/include/qkeycode.h /usr/lib/qt/include/qlabel.h
grafic.c.o: /usr/lib/qt/include/qlayout.h /usr/lib/qt/include/qlineedit.h
grafic.c.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qmultilinedit.h
grafic.c.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
grafic.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
grafic.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
grafic.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
grafic.c.o: /usr/lib/qt/include/qpointarray.h /usr/lib/qt/include/qprinter.h
grafic.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
grafic.c.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
grafic.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
grafic.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
grafic.c.o: /usr/lib/qt/include/qtimer.h /usr/lib/qt/include/qtooltip.h
grafic.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
grafic.c.o: /usr/lib/qt/include/qwmatrix.h cdata.h dfile.h formula.h grafic.h
grafic.c.o: stdfunc.h stdinc.h strings.h thefunc.h
grafic.o: /usr/include/_G_config.h /usr/include/alloca.h
grafic.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
grafic.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
grafic.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
grafic.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
grafic.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
grafic.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
grafic.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
grafic.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
grafic.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
grafic.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
grafic.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
grafic.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
grafic.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
grafic.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
grafic.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
grafic.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
grafic.o: /usr/include/xlocale.h
grafic.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
grafic.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
grafic.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
grafic.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
grafic.o: /usr/lib/qt/include/qarray.h /usr/lib/qt/include/qbrush.h
grafic.o: /usr/lib/qt/include/qbutton.h /usr/lib/qt/include/qcollection.h
grafic.o: /usr/lib/qt/include/qcolor.h /usr/lib/qt/include/qcursor.h
grafic.o: /usr/lib/qt/include/qdialog.h /usr/lib/qt/include/qevent.h
grafic.o: /usr/lib/qt/include/qfont.h /usr/lib/qt/include/qfontinfo.h
grafic.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
grafic.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgeneric.h
grafic.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
grafic.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlineedit.h
grafic.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qmultilinedit.h
grafic.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
grafic.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
grafic.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
grafic.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
grafic.o: /usr/lib/qt/include/qpointarray.h /usr/lib/qt/include/qpushbutton.h
grafic.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
grafic.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
grafic.o: /usr/lib/qt/include/qstack.h /usr/lib/qt/include/qstring.h
grafic.o: /usr/lib/qt/include/qtableview.h /usr/lib/qt/include/qwidget.h
grafic.o: /usr/lib/qt/include/qwindowdefs.h /usr/lib/qt/include/qwmatrix.h
grafic.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
head.c.o: /usr/include/_G_config.h /usr/include/alloca.h
head.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
head.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
head.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
head.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
head.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
head.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
head.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
head.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
head.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
head.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
head.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
head.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
head.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
head.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
head.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
head.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
head.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
head.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
head.c.o: /usr/include/wchar.h /usr/include/xlocale.h
head.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
head.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
head.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
head.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
head.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
integ.c.o: /usr/include/_G_config.h /usr/include/alloca.h
integ.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
integ.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
integ.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
integ.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
integ.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
integ.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
integ.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
integ.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
integ.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
integ.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
integ.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
integ.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
integ.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
integ.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
integ.c.o: /usr/include/linux/limits.h /usr/include/math.h
integ.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
integ.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
integ.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
integ.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
integ.c.o: /usr/include/xlocale.h
integ.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
integ.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
integ.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
integ.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
integ.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
look.c.o: /usr/include/_G_config.h /usr/include/alloca.h
look.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
look.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
look.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
look.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
look.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
look.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
look.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
look.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
look.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
look.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
look.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
look.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
look.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
look.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
look.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
look.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
look.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
look.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
look.c.o: /usr/include/wchar.h /usr/include/xlocale.h
look.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
look.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
look.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
look.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
look.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
matrix.c.o: /usr/include/_G_config.h /usr/include/alloca.h
matrix.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
matrix.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
matrix.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
matrix.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
matrix.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
matrix.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
matrix.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
matrix.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
matrix.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
matrix.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
matrix.c.o: /usr/include/ctype.h /usr/include/endian.h /usr/include/features.h
matrix.c.o: /usr/include/gconv.h /usr/include/gnu/stubs.h /usr/include/libio.h
matrix.c.o: /usr/include/limits.h /usr/include/linux/limits.h
matrix.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
matrix.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
matrix.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
matrix.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
matrix.c.o: /usr/include/xlocale.h
matrix.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
matrix.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
matrix.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
matrix.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
matrix.c.o: matrix.h stdfunc.h stdinc.h strings.h
matrix.o: /usr/include/_G_config.h /usr/include/alloca.h
matrix.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
matrix.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
matrix.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
matrix.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
matrix.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
matrix.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
matrix.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
matrix.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
matrix.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
matrix.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
matrix.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
matrix.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
matrix.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
matrix.o: /usr/include/stdlib.h /usr/include/sys/cdefs.h
matrix.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
matrix.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
matrix.o: /usr/include/xlocale.h
matrix.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
matrix.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
matrix.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
matrix.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
matrix.o: stdfunc.h stdinc.h strings.h
merge.c.o: /usr/include/_G_config.h /usr/include/alloca.h
merge.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
merge.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
merge.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
merge.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
merge.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
merge.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
merge.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
merge.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
merge.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
merge.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
merge.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
merge.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
merge.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
merge.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
merge.c.o: /usr/include/linux/limits.h /usr/include/math.h
merge.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
merge.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
merge.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
merge.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
merge.c.o: /usr/include/xlocale.h
merge.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
merge.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
merge.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
merge.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
merge.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
mspline.c.o: /usr/include/_G_config.h /usr/include/alloca.h
mspline.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
mspline.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
mspline.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
mspline.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
mspline.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
mspline.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
mspline.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
mspline.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
mspline.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
mspline.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
mspline.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
mspline.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
mspline.c.o: /usr/include/features.h /usr/include/gconv.h
mspline.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
mspline.c.o: /usr/include/libio.h /usr/include/limits.h
mspline.c.o: /usr/include/linux/limits.h /usr/include/math.h
mspline.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
mspline.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
mspline.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
mspline.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
mspline.c.o: /usr/include/xlocale.h
mspline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
mspline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
mspline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
mspline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
mspline.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
pfit.c.o: /usr/include/_G_config.h /usr/include/alloca.h
pfit.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
pfit.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
pfit.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
pfit.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
pfit.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
pfit.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
pfit.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
pfit.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
pfit.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
pfit.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
pfit.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
pfit.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
pfit.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
pfit.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
pfit.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
pfit.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
pfit.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
pfit.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
pfit.c.o: /usr/include/wchar.h /usr/include/xlocale.h
pfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
pfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
pfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
pfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
pfit.c.o: dfile.h formula.h simplex.h stdfunc.h stdinc.h strings.h thefunc.h
plot.c.o: /usr/include/_G_config.h /usr/include/alloca.h
plot.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
plot.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
plot.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
plot.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
plot.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
plot.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
plot.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
plot.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
plot.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
plot.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
plot.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
plot.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
plot.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
plot.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
plot.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
plot.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
plot.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
plot.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
plot.c.o: /usr/include/wchar.h /usr/include/xlocale.h
plot.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
plot.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
plot.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
plot.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
plot.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
plot.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
plot.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
plot.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
plot.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
plot.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
plot.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
plot.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
plot.c.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qlabel.h
plot.c.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
plot.c.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
plot.c.o: /usr/lib/qt/include/qobjectdefs.h /usr/lib/qt/include/qpaintdevice.h
plot.c.o: /usr/lib/qt/include/qpainter.h /usr/lib/qt/include/qpalette.h
plot.c.o: /usr/lib/qt/include/qpen.h /usr/lib/qt/include/qpixmap.h
plot.c.o: /usr/lib/qt/include/qpoint.h /usr/lib/qt/include/qpointarray.h
plot.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
plot.c.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
plot.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
plot.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
plot.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
plot.c.o: /usr/lib/qt/include/qwmatrix.h cdata.h dfile.h formula.h grafic.h
plot.c.o: stdfunc.h stdinc.h strings.h thefunc.h
polyfit.c.o: /usr/include/_G_config.h /usr/include/alloca.h
polyfit.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
polyfit.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
polyfit.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
polyfit.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
polyfit.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
polyfit.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
polyfit.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
polyfit.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
polyfit.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
polyfit.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
polyfit.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
polyfit.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
polyfit.c.o: /usr/include/features.h /usr/include/gconv.h
polyfit.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
polyfit.c.o: /usr/include/libio.h /usr/include/limits.h
polyfit.c.o: /usr/include/linux/limits.h /usr/include/math.h
polyfit.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
polyfit.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
polyfit.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
polyfit.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
polyfit.c.o: /usr/include/xlocale.h
polyfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
polyfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
polyfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
polyfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
polyfit.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
reduce.c.o: /usr/include/_G_config.h /usr/include/alloca.h
reduce.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
reduce.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
reduce.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
reduce.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
reduce.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
reduce.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
reduce.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
reduce.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
reduce.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
reduce.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
reduce.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
reduce.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
reduce.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
reduce.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
reduce.c.o: /usr/include/limits.h /usr/include/linux/limits.h
reduce.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
reduce.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
reduce.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
reduce.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
reduce.c.o: /usr/include/wchar.h /usr/include/xlocale.h
reduce.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
reduce.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
reduce.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
reduce.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
reduce.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
remove.c.o: /usr/include/_G_config.h /usr/include/alloca.h
remove.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
remove.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
remove.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
remove.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
remove.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
remove.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
remove.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
remove.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
remove.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
remove.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
remove.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
remove.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
remove.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
remove.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
remove.c.o: /usr/include/limits.h /usr/include/linux/limits.h
remove.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
remove.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
remove.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
remove.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
remove.c.o: /usr/include/wchar.h /usr/include/xlocale.h
remove.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
remove.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
remove.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
remove.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
remove.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
rgb.c.o: /usr/include/_G_config.h /usr/include/bits/pthreadtypes.h
rgb.c.o: /usr/include/bits/sched.h /usr/include/bits/stdio_lim.h
rgb.c.o: /usr/include/bits/types.h /usr/include/bits/wchar.h
rgb.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/gnu/stubs.h
rgb.c.o: /usr/include/libio.h /usr/include/stdio.h /usr/include/string.h
rgb.c.o: /usr/include/sys/cdefs.h /usr/include/wchar.h /usr/include/xlocale.h
rgb.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
rgb.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
rgb.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
rgb.c.o: /usr/lib/qt/include/qbitarray.h /usr/lib/qt/include/qbrush.h
rgb.c.o: /usr/lib/qt/include/qbutton.h /usr/lib/qt/include/qcolor.h
rgb.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdrawutil.h
rgb.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
rgb.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
rgb.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
rgb.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglobal.h
rgb.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlcdnumber.h
rgb.c.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
rgb.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
rgb.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
rgb.c.o: /usr/lib/qt/include/qpoint.h /usr/lib/qt/include/qpointarray.h
rgb.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrangecontrol.h
rgb.c.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
rgb.c.o: /usr/lib/qt/include/qscrollbar.h /usr/lib/qt/include/qshared.h
rgb.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstring.h
rgb.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
rgb.c.o: /usr/lib/qt/include/qwmatrix.h rgb.h
rgb.o: /usr/include/features.h /usr/include/gnu/stubs.h /usr/include/string.h
rgb.o: /usr/include/sys/cdefs.h /usr/include/xlocale.h
rgb.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
rgb.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
rgb.o: /usr/lib/qt/include/qbitarray.h /usr/lib/qt/include/qbrush.h
rgb.o: /usr/lib/qt/include/qbutton.h /usr/lib/qt/include/qcolor.h
rgb.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdrawutil.h
rgb.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
rgb.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
rgb.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
rgb.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglobal.h
rgb.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlcdnumber.h
rgb.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
rgb.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
rgb.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
rgb.o: /usr/lib/qt/include/qpoint.h /usr/lib/qt/include/qpointarray.h
rgb.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrangecontrol.h
rgb.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
rgb.o: /usr/lib/qt/include/qscrollbar.h /usr/lib/qt/include/qshared.h
rgb.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstring.h
rgb.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
rgb.o: /usr/lib/qt/include/qwmatrix.h
select.c.o: /usr/include/_G_config.h /usr/include/alloca.h
select.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
select.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
select.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
select.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
select.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
select.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
select.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
select.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
select.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
select.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
select.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
select.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
select.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
select.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
select.c.o: /usr/include/limits.h /usr/include/linux/limits.h
select.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
select.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
select.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
select.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
select.c.o: /usr/include/wchar.h /usr/include/xlocale.h
select.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
select.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
select.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
select.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
select.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
sfit.c.o: /usr/include/_G_config.h /usr/include/alloca.h
sfit.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
sfit.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
sfit.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
sfit.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
sfit.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
sfit.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
sfit.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
sfit.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
sfit.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
sfit.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
sfit.c.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
sfit.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
sfit.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
sfit.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
sfit.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
sfit.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
sfit.c.o: /usr/include/xlocale.h
sfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
sfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
sfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
sfit.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
sfit.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
sfit.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
sfit.c.o: /usr/lib/qt/include/qcheckbox.h /usr/lib/qt/include/qcollection.h
sfit.c.o: /usr/lib/qt/include/qcolor.h /usr/lib/qt/include/qcursor.h
sfit.c.o: /usr/lib/qt/include/qdialog.h /usr/lib/qt/include/qevent.h
sfit.c.o: /usr/lib/qt/include/qfont.h /usr/lib/qt/include/qfontinfo.h
sfit.c.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
sfit.c.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgeneric.h
sfit.c.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
sfit.c.o: /usr/lib/qt/include/qkeycode.h /usr/lib/qt/include/qlabel.h
sfit.c.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
sfit.c.o: /usr/lib/qt/include/qlistbox.h /usr/lib/qt/include/qmultilinedit.h
sfit.c.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
sfit.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
sfit.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
sfit.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
sfit.c.o: /usr/lib/qt/include/qpointarray.h /usr/lib/qt/include/qprinter.h
sfit.c.o: /usr/lib/qt/include/qprogressbar.h /usr/lib/qt/include/qpushbutton.h
sfit.c.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
sfit.c.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
sfit.c.o: /usr/lib/qt/include/qstack.h /usr/lib/qt/include/qstring.h
sfit.c.o: /usr/lib/qt/include/qtableview.h /usr/lib/qt/include/qtimer.h
sfit.c.o: /usr/lib/qt/include/qtooltip.h /usr/lib/qt/include/qwidget.h
sfit.c.o: /usr/lib/qt/include/qwindowdefs.h /usr/lib/qt/include/qwmatrix.h
sfit.c.o: cdata.h dfile.h formula.h grafic.h sfit.h simplex.h stdfunc.h
sfit.c.o: stdinc.h strings.h thefunc.h
sfit.o: /usr/include/_G_config.h /usr/include/alloca.h
sfit.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
sfit.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
sfit.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
sfit.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
sfit.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
sfit.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
sfit.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
sfit.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
sfit.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
sfit.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
sfit.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
sfit.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
sfit.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
sfit.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
sfit.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
sfit.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
sfit.o: /usr/include/xlocale.h
sfit.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
sfit.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
sfit.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
sfit.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
sfit.o: /usr/lib/qt/include/qarray.h /usr/lib/qt/include/qbrush.h
sfit.o: /usr/lib/qt/include/qbutton.h /usr/lib/qt/include/qcheckbox.h
sfit.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
sfit.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
sfit.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
sfit.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
sfit.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
sfit.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
sfit.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qlabel.h
sfit.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
sfit.o: /usr/lib/qt/include/qlistbox.h /usr/lib/qt/include/qmultilinedit.h
sfit.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
sfit.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
sfit.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
sfit.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
sfit.o: /usr/lib/qt/include/qpointarray.h /usr/lib/qt/include/qprogressbar.h
sfit.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
sfit.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
sfit.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
sfit.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
sfit.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
sfit.o: /usr/lib/qt/include/qwmatrix.h cdata.h dfile.h formula.h grafic.h
sfit.o: simplex.h stdfunc.h stdinc.h strings.h thefunc.h
show.c.o: /usr/include/_G_config.h /usr/include/alloca.h
show.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
show.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
show.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
show.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
show.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
show.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
show.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
show.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
show.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
show.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
show.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
show.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
show.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
show.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
show.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
show.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
show.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
show.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
show.c.o: /usr/include/wchar.h /usr/include/xlocale.h
show.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
show.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
show.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
show.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
show.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
show.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
show.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
show.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
show.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
show.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
show.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
show.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
show.c.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qlabel.h
show.c.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
show.c.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
show.c.o: /usr/lib/qt/include/qobjectdefs.h /usr/lib/qt/include/qpaintdevice.h
show.c.o: /usr/lib/qt/include/qpainter.h /usr/lib/qt/include/qpalette.h
show.c.o: /usr/lib/qt/include/qpen.h /usr/lib/qt/include/qpixmap.h
show.c.o: /usr/lib/qt/include/qpoint.h /usr/lib/qt/include/qpointarray.h
show.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
show.c.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
show.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
show.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
show.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
show.c.o: /usr/lib/qt/include/qwmatrix.h cdata.h dfile.h formula.h grafic.h
show.c.o: stdfunc.h stdinc.h strings.h thefunc.h
simplex.c.o: /usr/include/_G_config.h /usr/include/alloca.h
simplex.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
simplex.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
simplex.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
simplex.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
simplex.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
simplex.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
simplex.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
simplex.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
simplex.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
simplex.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
simplex.c.o: /usr/include/endian.h /usr/include/features.h
simplex.c.o: /usr/include/gconv.h /usr/include/gnu/stubs.h
simplex.c.o: /usr/include/libio.h /usr/include/limits.h
simplex.c.o: /usr/include/linux/limits.h /usr/include/math.h
simplex.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
simplex.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
simplex.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
simplex.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
simplex.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
simplex.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
simplex.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
simplex.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
simplex.c.o: cdata.h formula.h simplex.h stdfunc.h stdinc.h strings.h
simplex.o: /usr/include/_G_config.h /usr/include/alloca.h
simplex.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
simplex.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
simplex.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
simplex.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
simplex.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
simplex.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
simplex.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
simplex.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
simplex.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
simplex.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
simplex.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
simplex.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
simplex.o: /usr/include/linux/limits.h /usr/include/math.h
simplex.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sys/cdefs.h
simplex.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
simplex.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
simplex.o: /usr/include/xlocale.h
simplex.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
simplex.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
simplex.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
simplex.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
simplex.o: formula.h stdfunc.h stdinc.h strings.h
sort.c.o: /usr/include/_G_config.h /usr/include/alloca.h
sort.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
sort.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
sort.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
sort.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
sort.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
sort.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
sort.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
sort.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
sort.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
sort.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
sort.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
sort.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
sort.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
sort.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
sort.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
sort.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
sort.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
sort.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
sort.c.o: /usr/include/wchar.h /usr/include/xlocale.h
sort.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
sort.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
sort.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
sort.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
sort.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
spline.c.o: /usr/include/_G_config.h /usr/include/alloca.h
spline.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
spline.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
spline.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
spline.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
spline.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
spline.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
spline.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
spline.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
spline.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
spline.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
spline.c.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
spline.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
spline.c.o: /usr/include/limits.h /usr/include/linux/limits.h
spline.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
spline.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
spline.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
spline.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
spline.c.o: /usr/include/xlocale.h
spline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
spline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
spline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
spline.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
spline.c.o: cdata.h formula.h spline.h stdfunc.h stdinc.h strings.h
spline.o: /usr/include/_G_config.h /usr/include/alloca.h
spline.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
spline.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
spline.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
spline.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
spline.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
spline.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
spline.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
spline.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
spline.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
spline.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
spline.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
spline.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
spline.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
spline.o: /usr/include/stdlib.h /usr/include/sys/cdefs.h
spline.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
spline.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
spline.o: /usr/include/xlocale.h
spline.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
spline.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
spline.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
spline.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
spline.o: formula.h stdfunc.h stdinc.h strings.h
stdfunc.c.o: /usr/include/_G_config.h /usr/include/alloca.h
stdfunc.c.o: /usr/include/asm/errno.h /usr/include/bits/confname.h
stdfunc.c.o: /usr/include/bits/dirent.h /usr/include/bits/endian.h
stdfunc.c.o: /usr/include/bits/environments.h /usr/include/bits/errno.h
stdfunc.c.o: /usr/include/bits/huge_val.h /usr/include/bits/local_lim.h
stdfunc.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
stdfunc.c.o: /usr/include/bits/nan.h /usr/include/bits/posix1_lim.h
stdfunc.c.o: /usr/include/bits/posix2_lim.h /usr/include/bits/posix_opt.h
stdfunc.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
stdfunc.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
stdfunc.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
stdfunc.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
stdfunc.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
stdfunc.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
stdfunc.c.o: /usr/include/ctype.h /usr/include/dirent.h /usr/include/endian.h
stdfunc.c.o: /usr/include/errno.h /usr/include/features.h /usr/include/gconv.h
stdfunc.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
stdfunc.c.o: /usr/include/libio.h /usr/include/limits.h
stdfunc.c.o: /usr/include/linux/errno.h /usr/include/linux/limits.h
stdfunc.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
stdfunc.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
stdfunc.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
stdfunc.c.o: /usr/include/sys/types.h /usr/include/time.h
stdfunc.c.o: /usr/include/unistd.h /usr/include/wchar.h /usr/include/xlocale.h
stdfunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
stdfunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
stdfunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
stdfunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
stdfunc.c.o: cdata.h formula.h simplex.h stdfunc.h stdinc.h strings.h
stdfunc.o: /usr/include/_G_config.h /usr/include/alloca.h
stdfunc.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
stdfunc.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
stdfunc.o: /usr/include/bits/nan.h /usr/include/bits/pthreadtypes.h
stdfunc.o: /usr/include/bits/sched.h /usr/include/bits/select.h
stdfunc.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
stdfunc.o: /usr/include/bits/time.h /usr/include/bits/types.h
stdfunc.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
stdfunc.o: /usr/include/bits/wchar.h /usr/include/endian.h
stdfunc.o: /usr/include/features.h /usr/include/gconv.h
stdfunc.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/math.h
stdfunc.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sys/cdefs.h
stdfunc.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
stdfunc.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
stdfunc.o: /usr/include/xlocale.h
stdfunc.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
stdfunc.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h stdinc.h
stdfunc.o: strings.h
strings.c.o: /usr/include/_G_config.h /usr/include/alloca.h
strings.c.o: /usr/include/asm/errno.h /usr/include/bits/confname.h
strings.c.o: /usr/include/bits/dirent.h /usr/include/bits/endian.h
strings.c.o: /usr/include/bits/environments.h /usr/include/bits/errno.h
strings.c.o: /usr/include/bits/huge_val.h /usr/include/bits/local_lim.h
strings.c.o: /usr/include/bits/mathcalls.h /usr/include/bits/mathdef.h
strings.c.o: /usr/include/bits/nan.h /usr/include/bits/posix1_lim.h
strings.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
strings.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
strings.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
strings.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
strings.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
strings.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
strings.c.o: /usr/include/ctype.h /usr/include/dirent.h /usr/include/endian.h
strings.c.o: /usr/include/errno.h /usr/include/features.h /usr/include/gconv.h
strings.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
strings.c.o: /usr/include/libio.h /usr/include/linux/errno.h
strings.c.o: /usr/include/linux/limits.h /usr/include/math.h
strings.c.o: /usr/include/regex.h /usr/include/stdio.h /usr/include/stdlib.h
strings.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
strings.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
strings.c.o: /usr/include/sys/types.h /usr/include/time.h
strings.c.o: /usr/include/unistd.h /usr/include/wchar.h /usr/include/xlocale.h
strings.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
strings.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
strings.c.o: stdfunc.h stdinc.h strings.h
strings.o: stdinc.h
thefunc.c.o: /usr/include/_G_config.h /usr/include/alloca.h
thefunc.c.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
thefunc.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
thefunc.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
thefunc.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
thefunc.c.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
thefunc.c.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
thefunc.c.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
thefunc.c.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
thefunc.c.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
thefunc.c.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
thefunc.c.o: /usr/include/endian.h /usr/include/features.h
thefunc.c.o: /usr/include/gconv.h /usr/include/gnu/stubs.h
thefunc.c.o: /usr/include/libio.h /usr/include/limits.h
thefunc.c.o: /usr/include/linux/limits.h /usr/include/math.h
thefunc.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
thefunc.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
thefunc.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
thefunc.c.o: /usr/include/time.h /usr/include/wchar.h /usr/include/xlocale.h
thefunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
thefunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
thefunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
thefunc.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
thefunc.c.o: cdata.h dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
thefunc.o: /usr/include/_G_config.h /usr/include/alloca.h
thefunc.o: /usr/include/bits/endian.h /usr/include/bits/huge_val.h
thefunc.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
thefunc.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
thefunc.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
thefunc.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
thefunc.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
thefunc.o: /usr/include/bits/stdio_lim.h /usr/include/bits/time.h
thefunc.o: /usr/include/bits/types.h /usr/include/bits/waitflags.h
thefunc.o: /usr/include/bits/waitstatus.h /usr/include/bits/wchar.h
thefunc.o: /usr/include/bits/wordsize.h /usr/include/bits/xopen_lim.h
thefunc.o: /usr/include/endian.h /usr/include/features.h /usr/include/gconv.h
thefunc.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
thefunc.o: /usr/include/linux/limits.h /usr/include/math.h
thefunc.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sys/cdefs.h
thefunc.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
thefunc.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/wchar.h
thefunc.o: /usr/include/xlocale.h
thefunc.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
thefunc.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
thefunc.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
thefunc.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
thefunc.o: formula.h stdfunc.h stdinc.h strings.h
theview.c.o: /usr/include/_G_config.h /usr/include/alloca.h
theview.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
theview.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
theview.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
theview.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
theview.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
theview.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
theview.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
theview.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
theview.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
theview.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
theview.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
theview.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
theview.c.o: /usr/include/features.h /usr/include/gconv.h
theview.c.o: /usr/include/getopt.h /usr/include/gnu/stubs.h
theview.c.o: /usr/include/libio.h /usr/include/limits.h
theview.c.o: /usr/include/linux/limits.h /usr/include/math.h
theview.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
theview.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
theview.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
theview.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
theview.c.o: /usr/include/xlocale.h
theview.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
theview.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
theview.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
theview.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
theview.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
theview.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
theview.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
theview.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
theview.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
theview.c.o: /usr/lib/qt/include/qfontinfo.h
theview.c.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
theview.c.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgeneric.h
theview.c.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
theview.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlineedit.h
theview.c.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qmultilinedit.h
theview.c.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
theview.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
theview.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
theview.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
theview.c.o: /usr/lib/qt/include/qpointarray.h
theview.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
theview.c.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
theview.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
theview.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
theview.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
theview.c.o: /usr/lib/qt/include/qwmatrix.h cdata.h dfile.h formula.h grafic.h
theview.c.o: stdfunc.h stdinc.h strings.h thefunc.h
u2th.c.o: /usr/include/_G_config.h /usr/include/alloca.h
u2th.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
u2th.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
u2th.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
u2th.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
u2th.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
u2th.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
u2th.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
u2th.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
u2th.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
u2th.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
u2th.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
u2th.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
u2th.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
u2th.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
u2th.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
u2th.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
u2th.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
u2th.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
u2th.c.o: /usr/include/wchar.h /usr/include/xlocale.h
u2th.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
u2th.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
u2th.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
u2th.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
u2th.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
xfun.c.o: /usr/include/_G_config.h /usr/include/alloca.h
xfun.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
xfun.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
xfun.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
xfun.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
xfun.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
xfun.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
xfun.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
xfun.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
xfun.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
xfun.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
xfun.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
xfun.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
xfun.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
xfun.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
xfun.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
xfun.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
xfun.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
xfun.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
xfun.c.o: /usr/include/wchar.h /usr/include/xlocale.h
xfun.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
xfun.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
xfun.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
xfun.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
xfun.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
xfun.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
xfun.c.o: /usr/lib/qt/include/qcheckbox.h /usr/lib/qt/include/qcollection.h
xfun.c.o: /usr/lib/qt/include/qcolor.h /usr/lib/qt/include/qcursor.h
xfun.c.o: /usr/lib/qt/include/qdialog.h /usr/lib/qt/include/qevent.h
xfun.c.o: /usr/lib/qt/include/qfont.h /usr/lib/qt/include/qfontinfo.h
xfun.c.o: /usr/lib/qt/include/qfontmetrics.h /usr/lib/qt/include/qframe.h
xfun.c.o: /usr/lib/qt/include/qgarray.h /usr/lib/qt/include/qgeneric.h
xfun.c.o: /usr/lib/qt/include/qglist.h /usr/lib/qt/include/qglobal.h
xfun.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlineedit.h
xfun.c.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qlistbox.h
xfun.c.o: /usr/lib/qt/include/qmultilinedit.h /usr/lib/qt/include/qobject.h
xfun.c.o: /usr/lib/qt/include/qobjectdefs.h /usr/lib/qt/include/qpaintdevice.h
xfun.c.o: /usr/lib/qt/include/qpainter.h /usr/lib/qt/include/qpalette.h
xfun.c.o: /usr/lib/qt/include/qpen.h /usr/lib/qt/include/qpixmap.h
xfun.c.o: /usr/lib/qt/include/qpoint.h /usr/lib/qt/include/qpointarray.h
xfun.c.o: /usr/lib/qt/include/qprogressbar.h /usr/lib/qt/include/qpushbutton.h
xfun.c.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
xfun.c.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
xfun.c.o: /usr/lib/qt/include/qstack.h /usr/lib/qt/include/qstring.h
xfun.c.o: /usr/lib/qt/include/qtableview.h /usr/lib/qt/include/qwidget.h
xfun.c.o: /usr/lib/qt/include/qwindowdefs.h /usr/lib/qt/include/qwmatrix.h
xfun.c.o: cdata.h dfile.h formula.h grafic.h sfit.h simplex.h stdfunc.h
xfun.c.o: stdinc.h strings.h thefunc.h
xhelp.c.o: /usr/include/_G_config.h /usr/include/alloca.h
xhelp.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
xhelp.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
xhelp.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
xhelp.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
xhelp.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
xhelp.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
xhelp.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
xhelp.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
xhelp.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
xhelp.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
xhelp.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
xhelp.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
xhelp.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
xhelp.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
xhelp.c.o: /usr/include/linux/limits.h /usr/include/math.h
xhelp.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
xhelp.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
xhelp.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
xhelp.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
xhelp.c.o: /usr/include/xlocale.h
xhelp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
xhelp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
xhelp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
xhelp.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
xhelp.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
xhelp.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
xhelp.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
xhelp.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdatastream.h
xhelp.c.o: /usr/lib/qt/include/qdatetime.h /usr/lib/qt/include/qdialog.h
xhelp.c.o: /usr/lib/qt/include/qdir.h /usr/lib/qt/include/qevent.h
xhelp.c.o: /usr/lib/qt/include/qfile.h /usr/lib/qt/include/qfiledialog.h
xhelp.c.o: /usr/lib/qt/include/qfileinfo.h /usr/lib/qt/include/qfont.h
xhelp.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
xhelp.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
xhelp.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
xhelp.c.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qiodevice.h
xhelp.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlineedit.h
xhelp.c.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qmultilinedit.h
xhelp.c.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
xhelp.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
xhelp.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
xhelp.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
xhelp.c.o: /usr/lib/qt/include/qpointarray.h /usr/lib/qt/include/qpushbutton.h
xhelp.c.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
xhelp.c.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
xhelp.c.o: /usr/lib/qt/include/qstack.h /usr/lib/qt/include/qstring.h
xhelp.c.o: /usr/lib/qt/include/qstrlist.h /usr/lib/qt/include/qtableview.h
xhelp.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
xhelp.c.o: /usr/lib/qt/include/qwmatrix.h WidGen.h cdata.h dfile.h formula.h
xhelp.c.o: grafic.h stdfunc.h stdinc.h strings.h thefunc.h
xmount.c.o: /usr/include/_G_config.h /usr/include/alloca.h
xmount.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
xmount.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
xmount.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
xmount.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
xmount.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
xmount.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
xmount.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
xmount.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
xmount.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
xmount.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
xmount.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
xmount.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
xmount.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
xmount.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h
xmount.c.o: /usr/include/limits.h /usr/include/linux/limits.h
xmount.c.o: /usr/include/math.h /usr/include/stdio.h /usr/include/stdlib.h
xmount.c.o: /usr/include/string.h /usr/include/sys/cdefs.h
xmount.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
xmount.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
xmount.c.o: /usr/include/wchar.h /usr/include/xlocale.h
xmount.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
xmount.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
xmount.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
xmount.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
xmount.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
xmount.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
xmount.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
xmount.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdialog.h
xmount.c.o: /usr/lib/qt/include/qevent.h /usr/lib/qt/include/qfont.h
xmount.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
xmount.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
xmount.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
xmount.c.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qlabel.h
xmount.c.o: /usr/lib/qt/include/qlineedit.h /usr/lib/qt/include/qlist.h
xmount.c.o: /usr/lib/qt/include/qlistbox.h /usr/lib/qt/include/qmultilinedit.h
xmount.c.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
xmount.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
xmount.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
xmount.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
xmount.c.o: /usr/lib/qt/include/qpointarray.h
xmount.c.o: /usr/lib/qt/include/qpushbutton.h /usr/lib/qt/include/qrect.h
xmount.c.o: /usr/lib/qt/include/qregion.h /usr/lib/qt/include/qshared.h
xmount.c.o: /usr/lib/qt/include/qsize.h /usr/lib/qt/include/qstack.h
xmount.c.o: /usr/lib/qt/include/qstring.h /usr/lib/qt/include/qtableview.h
xmount.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
xmount.c.o: /usr/lib/qt/include/qwmatrix.h browser.h cdata.h dfile.h formula.h
xmount.c.o: grafic.h stdfunc.h stdinc.h strings.h thefunc.h
xport.c.o: /usr/include/_G_config.h /usr/include/alloca.h
xport.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
xport.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
xport.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
xport.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
xport.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
xport.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
xport.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
xport.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
xport.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
xport.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
xport.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
xport.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
xport.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
xport.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
xport.c.o: /usr/include/linux/limits.h /usr/include/math.h
xport.c.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/string.h
xport.c.o: /usr/include/sys/cdefs.h /usr/include/sys/select.h
xport.c.o: /usr/include/sys/sysmacros.h /usr/include/sys/types.h
xport.c.o: /usr/include/time.h /usr/include/unistd.h /usr/include/wchar.h
xport.c.o: /usr/include/xlocale.h
xport.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
xport.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
xport.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
xport.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h
xport.c.o: /usr/lib/qt/include/qapplication.h /usr/lib/qt/include/qarray.h
xport.c.o: /usr/lib/qt/include/qbrush.h /usr/lib/qt/include/qbutton.h
xport.c.o: /usr/lib/qt/include/qcollection.h /usr/lib/qt/include/qcolor.h
xport.c.o: /usr/lib/qt/include/qcursor.h /usr/lib/qt/include/qdatastream.h
xport.c.o: /usr/lib/qt/include/qdatetime.h /usr/lib/qt/include/qdialog.h
xport.c.o: /usr/lib/qt/include/qdir.h /usr/lib/qt/include/qevent.h
xport.c.o: /usr/lib/qt/include/qfile.h /usr/lib/qt/include/qfiledialog.h
xport.c.o: /usr/lib/qt/include/qfileinfo.h /usr/lib/qt/include/qfont.h
xport.c.o: /usr/lib/qt/include/qfontinfo.h /usr/lib/qt/include/qfontmetrics.h
xport.c.o: /usr/lib/qt/include/qframe.h /usr/lib/qt/include/qgarray.h
xport.c.o: /usr/lib/qt/include/qgeneric.h /usr/lib/qt/include/qglist.h
xport.c.o: /usr/lib/qt/include/qglobal.h /usr/lib/qt/include/qiodevice.h
xport.c.o: /usr/lib/qt/include/qlabel.h /usr/lib/qt/include/qlineedit.h
xport.c.o: /usr/lib/qt/include/qlist.h /usr/lib/qt/include/qmultilinedit.h
xport.c.o: /usr/lib/qt/include/qobject.h /usr/lib/qt/include/qobjectdefs.h
xport.c.o: /usr/lib/qt/include/qpaintdevice.h /usr/lib/qt/include/qpainter.h
xport.c.o: /usr/lib/qt/include/qpalette.h /usr/lib/qt/include/qpen.h
xport.c.o: /usr/lib/qt/include/qpixmap.h /usr/lib/qt/include/qpoint.h
xport.c.o: /usr/lib/qt/include/qpointarray.h /usr/lib/qt/include/qpushbutton.h
xport.c.o: /usr/lib/qt/include/qrect.h /usr/lib/qt/include/qregion.h
xport.c.o: /usr/lib/qt/include/qshared.h /usr/lib/qt/include/qsize.h
xport.c.o: /usr/lib/qt/include/qstack.h /usr/lib/qt/include/qstring.h
xport.c.o: /usr/lib/qt/include/qstrlist.h /usr/lib/qt/include/qtableview.h
xport.c.o: /usr/lib/qt/include/qwidget.h /usr/lib/qt/include/qwindowdefs.h
xport.c.o: /usr/lib/qt/include/qwmatrix.h WidGen.h cdata.h dfile.h formula.h
xport.c.o: grafic.h stdfunc.h stdinc.h strings.h thefunc.h
zero.c.o: /usr/include/_G_config.h /usr/include/alloca.h
zero.c.o: /usr/include/bits/confname.h /usr/include/bits/endian.h
zero.c.o: /usr/include/bits/environments.h /usr/include/bits/huge_val.h
zero.c.o: /usr/include/bits/local_lim.h /usr/include/bits/mathcalls.h
zero.c.o: /usr/include/bits/mathdef.h /usr/include/bits/nan.h
zero.c.o: /usr/include/bits/posix1_lim.h /usr/include/bits/posix2_lim.h
zero.c.o: /usr/include/bits/posix_opt.h /usr/include/bits/pthreadtypes.h
zero.c.o: /usr/include/bits/sched.h /usr/include/bits/select.h
zero.c.o: /usr/include/bits/sigset.h /usr/include/bits/stdio_lim.h
zero.c.o: /usr/include/bits/time.h /usr/include/bits/types.h
zero.c.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
zero.c.o: /usr/include/bits/wchar.h /usr/include/bits/wordsize.h
zero.c.o: /usr/include/bits/xopen_lim.h /usr/include/endian.h
zero.c.o: /usr/include/features.h /usr/include/gconv.h /usr/include/getopt.h
zero.c.o: /usr/include/gnu/stubs.h /usr/include/libio.h /usr/include/limits.h
zero.c.o: /usr/include/linux/limits.h /usr/include/math.h /usr/include/stdio.h
zero.c.o: /usr/include/stdlib.h /usr/include/string.h /usr/include/sys/cdefs.h
zero.c.o: /usr/include/sys/select.h /usr/include/sys/sysmacros.h
zero.c.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/unistd.h
zero.c.o: /usr/include/wchar.h /usr/include/xlocale.h
zero.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/limits.h
zero.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stdarg.h
zero.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/stddef.h
zero.c.o: /usr/lib/gcc-lib/i486-suse-linux/2.95.3/include/syslimits.h cdata.h
zero.c.o: dfile.h formula.h stdfunc.h stdinc.h strings.h thefunc.h
