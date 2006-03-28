# Makefile for all McPhase

include src/Makefile.common

cfdir = cf1ion_module
mcpdir = src
vecdir = src/vector

VECOBJ= matpack.o sortdbl.o sort2dbl.o sort3dbl.o\
        ivector.o imatrix.o \
        dvector.o dmatrix.o dmatmul.o \
        cvector.o cmatrix.o cmatmul.o \
        dludcomp.o dsvdcomp.o dortho.o \
        cludcomp.o \
        tred.o imtql.o imtql2.o rseigen.o cheigen.o chgeigen.o  \
        matlapl.o moment.o matexph.o matsqth.o \
        magicsquare.o

#export INCLDIRS =      -I. \
		-I../include \
		-I/usr/include/readline \
		-I/local/lib  \
		-I/usr/X11R6/include

#export INCLUDE = ./include

all: vector cfield mcphase

vector: 
	cd $(vecdir) && $(MAKE) 

cfield: vector
	cd $(cfdir) && $(MAKE)

mcphase: vector cfield
	cd $(mcpdir) && $(MAKE)

clean: 
	cd $(vecdir) && $(MAKE) clean
	cd $(cfdir) && $(MAKE) clean
	cd $(mcpdir) && $(MAKE) clean
	rm -vf bin/addj.exe bin/charges.exe bin/coq2jjj.exe \
		bin/mcdispit.exe bin/singleion.exe bin/cfield.exe \
		bin/cond.exe bin/jjj2j.exe bin/mcphasit.exe bin/spins.exe 
	rm -vf bin/addj bin/charges bin/coq2jjj \
		bin/mcdispit bin/singleion bin/cfield \
		bin/cond bin/jjj2j bin/mcphasit bin/spins 
