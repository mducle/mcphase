CC = g++
LINKER = g++

# Figures operating system
ifeq (${OSTYPE},msys)
MSYS=1
endif
ifeq ($(shell uname),Linux)
 LINUX = LINUX
 LIBFORTRAN=-lgfortran
else
# Check for cygwin
ifneq (,$(findstring CYGWIN,$(shell uname)))
#ARCHFLAG+=-mno-cygwin
# To fix an error in spincf where min and max is used as variable names
#   and conflicts with definitions in <w32api/windef.h>
ARCHFLAG+=-DNOMINMAX
# Checks to see if we want to use the 64bit crosscompiler in Cygwin
ifdef cross64
CC=x86_64-w64-mingw32-g++ -m64 -pipe $(ARCHFLAG)
LINKER=x86_64-w64-mingw32-g++ -m64
endif
ifdef cross32
CC=i686-w64-mingw32-g++ -m32 -pipe $(ARCHFLAG)
LINKER=i686-w64-mingw32-g++ -m32
endif
endif
 MSYS=1
 LIBFORTRAN=-lg2c
endif

ifneq ($(shell $(CC) --version|head -n1|grep " 3."),)
 LIBFORTRAN=-lg2c
else
ifneq ($(shell $(CC) --version|head -n1|grep " 4."),)
 LIBFORTRAN=-lgfortran
endif
endif

# Determines if we have GOTOBLAS
ifneq ($(strip $(wildcard libs/GotoBLAS/libgoto.a)),)
ifndef cross64
gotoblas=1
endif
endif

##################################################################################
#                              Declarations                                      #
##################################################################################
#------------------------------    MSYS    --------------------------------------#
ifdef MSYS
CFLAGS = -g -Wall --pedantic -D_WINDOWS
else
CFLAGS = -g -ansi -Wall --pedantic -fPIC
endif
#------------------------------    FAST    --------------------------------------#
ifdef fast
ifneq ($(shell $(CC) --version|head -n1|grep " 3."),)
CFLAGS+= -O2 -pipe -fomit-frame-pointer
else
CFLAGS+= -O2 -march=native -pipe -fomit-frame-pointer
endif
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
LDFLAGS = -Llibs -Llibs/GotoBLAS -larpack -llapack -lgoto $(LIBFORTRAN)
ifeq (,$(findstring CYGWIN,$(shell uname)))
LDFLAGS += -lpthread
endif
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
##################################################################################

# List of object files to be created and their headers
OBJECT = maths states        icpars ic_hmltn ic_diag mmio njsyms cfp coulomb so_cf lovesey
HEADER = maths states ic1ion icpars

.PHONY: default all clean

default: all
all: ic1ion ic1ion.so

# Defines some variables
OBJ := $(foreach objs,$(OBJECT),$(objs).o)
HPP := $(foreach head,$(HEADER),$(head).hpp)

VECOBJ = ../src/vector/cheigen.o ../src/vector/chgeigen.o ../src/vector/cludcomp.o ../src/vector/cmatmul.o ../src/vector/cmatrix.o ../src/vector/cvector.o ../src/vector/dludcomp.o ../src/vector/dmatmul.o ../src/vector/dmatrix.o ../src/vector/dortho.o ../src/vector/dsvdcomp.o ../src/vector/dvector.o ../src/vector/imatrix.o ../src/vector/imtql.o ../src/vector/imtql2.o ../src/vector/ivector.o ../src/vector/magicsquare.o ../src/vector/matexph.o ../src/vector/matlapl.o ../src/vector/matpack.o ../src/vector/matsqth.o ../src/vector/moment.o ../src/vector/rseigen.o ../src/vector/sort2dbl.o ../src/vector/sort3dbl.o ../src/vector/sortdbl.o ../src/vector/tred.o

libs/libarpack.a:
	cd libs && make

ifdef integral
ic1ion: libs/libarpack.a $(OBJ) ic1ion.o ic1ion_module.o
	$(LINKER) -o ic1ion ic1ion.o $(OBJ) $(LDFLAGS) ic1ion_module.o $(VECOBJ) $(STATIC)
else
ic1ion: libs/libarpack.a $(OBJ) ic1ion.o
	$(LINKER) -o ic1ion ic1ion.o $(OBJ) $(LDFLAGS) $(STATIC)
endif

ic1ion.o: ic1ion.cpp maths.hpp states.hpp ic1ion.hpp
ifdef integral
	$(CC) $(MATHSCFLAGS) -D_INTEGRAL -c $<
else
	$(CC) $(MATHSCFLAGS) -c $<
endif

ic1ion_module.o: ic1ion_module.cpp
ifdef integral
	$(CC) $(MATHSCFLAGS) -D_INTEGRAL -I../include -c $<
else
ifdef MSYS
	$(CC) $(MATHSCFLAGS) -I../include -DBUILD_DLL -c $<
else
	$(CC) $(MATHSCFLAGS) -I../include -c $<
endif
endif

ic1ion.so: libs/libarpack.a $(OBJ) ic1ion_module.o ic1ion.o
ifdef MSYS
	$(LINKER) -o ic1ion.dll ic1ion_module.o ic1ion.o $(OBJ) $(VECOBJ) -shared $(LDFLAGS) -Wl,--out-implib,libic1ion_module.a
else
	$(LINKER) -o ic1ion.so  ic1ion_module.o ic1ion.o $(OBJ) $(VECOBJ) -shared $(LDFLAGS)
endif

# Patern rules
# Syntax is
# targets ...: target-pattern: dep-patterns ...
#	commands
#	...

$(OBJ): %.o: %.cpp $(HPP)
	$(CC) $(CFLAGS) -c $<

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
