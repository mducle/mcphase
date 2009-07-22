# Figures operating system
ifeq (${OSTYPE},msys)
MSYS=1
endif
ifeq ($(shell uname),Linux)
 LINUX = LINUX
 LIBFORTRAN=-lgfortran
else
 MSYS=1
 LIBFORTRAN=-lg2c
endif

ifneq ($(shell gcc --version|head -n1|grep " 3."),)
 LIBFORTRAN=-lg2c
else
ifneq ($(shell gcc --version|head -n1|grep " 4."),)
 LIBFORTRAN=-lgfortran
endif
endif

# Determines if we have GOTOBLAS
ifneq ($(strip $(wildcard libs/GotoBLAS/libgoto.a)),)
gotoblas=1
endif

##################################################################################
#                              Declarations                                      #
##################################################################################
CC = g++
LINKER = g++
#------------------------------    MSYS    --------------------------------------#
ifdef MSYS
CFLAGS = -g -Wall --pedantic -D_WINDOWS
else
CFLAGS = -g -ansi -Wall --pedantic -fPIC
endif
#------------------------------    FAST    --------------------------------------#
ifdef fast
CFLAGS+= -O2 -march=native -pipe -fomit-frame-pointer
endif
#------------------------------   LINKER   --------------------------------------#
ifdef matlab
ifdef MSYS
LDFLAGS = -L. -lmwrefblas -lmwlapack -lmwarpack
else
LDFLAGS = -L/usr/local/MATLAB/bin/glnxa64 -Wl,-rpath,/usr/local/MATLAB/bin/glnxa64 -lut -lmwbinder -lmwblas -lmwlapack -Llibs -larpack $(LIBFORTRAN)
endif
else
ifdef gotoblas
LDFLAGS = -Llibs -Llibs/GotoBLAS -larpack -llapack -lgoto $(LIBFORTRAN) -lpthread
else
LDFLAGS = -Llibs -larpack -llapack -lblas $(LIBFORTRAN)
endif
endif
#------------------------------  32/64bit  --------------------------------------#
ifdef b32
LINKER += -m32

CFLAGS += -m32
LDFLAGS += -m32
endif
#------------------------------   DISTRO   --------------------------------------#
ifdef distrib
CFLAGS += -mtune=pentium3 -march=i686
STATIC = -static
endif
#------------------------------   LAPACK   --------------------------------------#
# Use lapack routines in maths.cpp by default
ifdef nomathslapack
MATHSCFLAGS = $(CFLAGS)
else
MATHSCFLAGS = $(CFLAGS) -DUSE_LAPACK
endif
##################################################################################

# List of object files to be created and their headers
OBJECT =       states               mmio njsyms cfp coulomb so_cf lovesey#maths_irbl
HEADER = maths states ic1ion icpars
MATHOB = maths               icpars                                       ic_diag

.PHONY: default all clean

default: all
all: ic1ion ic1ion.so

# Defines some variables
OBJ := $(foreach objs,$(OBJECT),$(objs).o)
HPP := $(foreach head,$(HEADER),$(head).hpp)
MOB := $(foreach mobj,$(MATHOB),$(mobj).o)

VECOBJ = ../src/vector/cheigen.o ../src/vector/chgeigen.o ../src/vector/cludcomp.o ../src/vector/cmatmul.o ../src/vector/cmatrix.o ../src/vector/cvector.o ../src/vector/dludcomp.o ../src/vector/dmatmul.o ../src/vector/dmatrix.o ../src/vector/dortho.o ../src/vector/dsvdcomp.o ../src/vector/dvector.o ../src/vector/imatrix.o ../src/vector/imtql.o ../src/vector/imtql2.o ../src/vector/ivector.o ../src/vector/magicsquare.o ../src/vector/matexph.o ../src/vector/matlapl.o ../src/vector/matpack.o ../src/vector/matsqth.o ../src/vector/moment.o ../src/vector/rseigen.o ../src/vector/sort2dbl.o ../src/vector/sort3dbl.o ../src/vector/sortdbl.o ../src/vector/tred.o

libs/libarpack.a:
	cd libs && make

ic1ion: libs/libarpack.a maths.o $(OBJ) ic_hmltn.o ic1ion.o ic_diag.o icpars.o ic1ion_module.o
ifdef integral
	$(LINKER) -o ic1ion maths.o $(OBJ) ic_hmltn.o ic1ion.o ic_diag.o icpars.o $(LDFLAGS) ic1ion_module.o $(VECOBJ) $(STATIC)
else
	$(LINKER) -o ic1ion maths.o $(OBJ) ic_hmltn.o ic1ion.o ic_diag.o icpars.o $(LDFLAGS) $(STATIC)
endif

ic1ion.o: ic1ion.cpp maths.hpp states.hpp ic1ion.hpp
ifdef integral
	$(CC) $(MATHSCFLAGS) -D_INTEGRAL -c $<
else
	$(CC) $(MATHSCFLAGS) -c $<
endif

ic_hmltn.o: ic_hmltn.cpp
ifdef arpackpp
	$(CC) $(CFLAGS) -DUSE_ARPACKPP -Ilibs/arpack++/include -c $<
else
ifdef arpack
	$(CC) $(CFLAGS) -DUSE_ARPACK -c $<
else
ifndef nomathslapack
	$(CC) $(CFLAGS) -DUSE_LAPACK -c $<
else
	@echo "Warning: Compiling without Lapack. Many routines will be slow and may not work!"
	$(CC) $(CFLAGS) -c $<
endif
endif
endif

ic1ion_module.o: ic1ion_module.cpp
ifdef integral
	$(CC) $(MATHSCFLAGS) -D_INTEGRAL -I../include -c $<
else
ifdef MSYS
	$(CC) $(MATHSCFLAGS) -DBUILD_DLL -c $<
else
	$(CC) $(MATHSCFLAGS) -I../include -c $<
endif
endif

ic1ion.so: libs/libarpack.a maths.o $(OBJ) ic_hmltn.o ic1ion.o ic_diag.o icpars.o ic1ion_module.o
ifdef MSYS
	$(LINKER) -o ic1ion.dll maths.o $(OBJ) ic_hmltn.o ic1ion.o ic_diag.o icpars.o ic1ion_module.o $(VECOBJ) -shared $(LDFLAGS) -Wl,--out-implib,libic1ion_module.a
else
	$(LINKER) -o ic1ion.so maths.o $(OBJ) ic_hmltn.o ic1ion.o ic_diag.o icpars.o ic1ion_module.o $(VECOBJ) -shared $(LDFLAGS)
endif

# Patern rules
# Syntax is
# targets ...: target-pattern: dep-patterns ...
#	commands
#	...

$(OBJ): %.o: %.cpp $(HPP)
	$(CC) $(CFLAGS) -c $<

$(MOB): %.o: %.cpp $(HPP)
	$(CC) $(MATHSCFLAGS) -c $<

# Implicit rules (suffix rules)
# NB. $@ = current target, $< = implied source, $* = filename without suffix.
# Syntax is
# .insuffix.outsuffix:
#	command
# For the implicit rules to work we must declare the suffixes used here
#.SUFFIXES: .tex .aux .toc .bbl .blg .log .dvi .ps .pdf .myt
#.tex.pdf:
#ifneq ($(shell if [ $*.aux -nt $*.bbl ]; then echo 1; fi),1)
#	$(LaTeX) $< $(endnot)
#endif
#	$(LaTeX) $< $(endnot)
#
#.aux.bbl:
#	$(BibTeX)  $* $(endnot)
#
#.tex.aux:
#	$(LaTeX) $< $(endnot)

clean:
	rm -vf *.o
	rm -vf *~
	rm -vf ic_hmltn

cleanall: clean
	cd libs && make clean
