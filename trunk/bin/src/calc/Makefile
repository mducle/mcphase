#makefile for program calc
#vpath %.h  ../include
#vpath %.hpp ../include
#vpath %.o  ../src/Vector

###########################################################################
# gcc version egcs-2.95.2 19991024 (release) on a Linux system
# Intel Pentium PC, use -m486 for 486 CPU or -mpentium for Pentium Systems
###########################################################################

# C compiler 
CC= g++ -pipe
COPTS= -ansi -O3 -fomit-frame-pointer -fforce-addr -m486 -x c 

# C++ compiler
CCC= g++ -pipe
CCOPTS= -O3 -fomit-frame-pointer -fforce-addr -mpentium -x c++ 

#LIBL= /usr/lib/libfl.a
###########################################################################
### include directory: don't change !					  #
###########################################################################

INCLUDE= ./

###########################################################################
### general rules							  #
###########################################################################

.SUFFIXES: .cc .o

.c.o: 
	${CC} -c ${COPTS}  -I ${INCLUDE} $<

.cc.o: 
	${CCC} -c ${CCOPTS}   -I ${INCLUDE} $<

###########################################################################
###########################################################################



INCLDIRS =      -I. \
		-I../include \
		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include


all : calc fileop


CFOBJS= cdata.o dfile.o formula.o simplex.o spline.o stdfunc.o strings.o thefunc.o 
PROGS = calc.o fileop.o



${CFOBJS}:


${PROGS}: 

# used for internal module in mcphas and mcdisp
#cfieldrout.o: cfieldrout.c

#	$(CC) -o cfieldrout.o cfieldrout.c  $(INCLDIRS)


#cfield.so: module.o cfield.o  $(VECOBJS) $(CFOBJS)
#	g++  -Wall -o cfield.so $^ -shared 

calc: calc.o  $(CFOBJS)
	g++  -Wall -o ../../bin/calc $^ 

fileop: fileop.o  $(CFOBJS)
	g++  -Wall -o ../../bin/fileop $^ 


clean:
	rm -vf ./*.o 

