# Makefile for all McPhase

include bin/src/Makefile.common

cfdir = bin/cf1ion_module
mcpdir = bin/src
vecdir = bin/src/vector

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
		bin/cond.exe bin/jjj2j.exe bin/mcphasit.exe bin/spins.exe \
                bin/chrgplt.exe bin/pointc.exe
	rm -vf bin/addj bin/charges bin/coq2jjj \
		bin/mcdispit bin/singleion bin/cfield \
		bin/cond bin/jjj2j bin/mcphasit bin/spins \
                bin/chrgplt bin/pointc \
                bin/cf1ion_module/cfield.so
