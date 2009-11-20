#makefile for program cfield

OSTYPE=$(shell uname)
#  (linux windows)
ifeq ($(OSTYPE),Linux)
 LINUX = linux
else
 WINDOWS = windows
endif

vpath %.h  ../include
vpath %.hpp ../include
vpath %.o  ../src/vector

###########################################################################
# gcc version egcs-2.95.2 19991024 (release) on a Linux system
# Intel Pentium PC, use -m486 for 486 CPU or -mpentium for Pentium Systems
###########################################################################

# C compiler 
#CC= g++ -pipe
#COPTS= -ansi -O3 -fomit-frame-pointer -fforce-addr -m486 -x c 
#CC= g++ 
#COPTS= -ansi -O0 -x c 

# C++ compiler
#CCC= g++ -pipe
#CCC= g++ 
#CCOPTS= -O3 -fomit-frame-pointer -fforce-addr -mpentium -x c++ 
#CCOPTS= -O0 -Wall -DDEBUG
# ... define in makefile.common
include ../src/Makefile.common

#LIBL= /usr/lib/libfl.a
###########################################################################
### include directory: don't change !					  #
###########################################################################

INCLUDE= ../../include

###########################################################################
### general rules							  #
###########################################################################

.SUFFIXES: .cc .o

.c.o: 
	${CC} -c ${COPTS} -Wno-format -I ${INCLUDE} $<

.cc.o: 
	${CCC} -c ${CCOPTS} -Wno-format -I ${INCLUDE} $<

###########################################################################
###########################################################################



INCLDIRS =      -I. \
		-I../include \
		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include



all : cfieldrout.o cfield.so cfield module.o

VECOBJS= matpack.o sortdbl.o sort2dbl.o sort3dbl.o ivector.o \
		imatrix.o dvector.o dmatrix.o dmatmul.o cvector.o \
		cmatrix.o cmatmul.o dludcomp.o dsvdcomp.o dortho.o \
		cludcomp.o tred.o imtql.o imtql2.o rseigen.o cheigen.o \
		chgeigen.o matlapl.o moment.o matexph.o matsqth.o magicsquare.o

CFOBJS= cfield.o diahermx.o eingabe.o intensit.o komplex.o matrix.o mb11a.o minima.o \
	       ortho.o spline.o stevens.o theta.o va05a.o


${CFOBJS}: types.c

# used for internal module in mcphas and mcdisp
cfieldrout.o: cfieldrout.c $(VECOBJS) $(CFOBJS) 

#	$(CC) -o cfieldrout.o cfieldrout.c  $(INCLDIRS)


cfield.so: module.o cfield.o  $(VECOBJS) $(CFOBJS)
	$(CCC) $(CCOPTS) -Wall -o cfield.so $^ -shared $(LDFLAGS)

cfield: cfield.o  $(VECOBJS) $(CFOBJS)
	$(CCC) $(CCOPTS) -Wall -o ../cfield $^ $(LDFLAGS)

module.o: module.c 
	$(CCC) $(CCOPTS) -Wall -c module.c $(INCLDIRS) $(LDFLAGS)

clean:
	rm -vf ./*.o ./*.info ./*.parameter ./*.output

